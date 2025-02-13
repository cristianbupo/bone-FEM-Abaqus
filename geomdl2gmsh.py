import os
import numpy as np
from geomdl import  operations, multi
import gmsh
import FreeCAD2gmsh as f2g


def createTransfiniteSurface(edges, meshSize, direction="Left", ignorePoint=0):
    adjacencies = createTransfiniteLine(edges, meshSize)
    nEdges = len(edges)
    transfinitePoints = []

    # Get the unique values from the adjacencies
    for i in range(nEdges):
        if i < nEdges - 1:
            transfinitePoints.extend(np.setdiff1d(adjacencies[i], adjacencies[i+1]))
        else:
            transfinitePoints.extend(np.setdiff1d(adjacencies[i], adjacencies[0]))

    if nEdges > 4:
        if isinstance(ignorePoint, int):
            transfinitePoints.pop(ignorePoint)
        elif isinstance(ignorePoint, list):
            for ip in sorted(ignorePoint, reverse=True):
                transfinitePoints.pop(ip)

    loop = gmsh.model.occ.addCurveLoop(edges)
    surface = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setTransfiniteSurface(surface, direction, transfinitePoints)


def createTransfiniteLine(edges, meshSize):
    adjacencies = []
    for i, edge in enumerate(edges):
        gmsh.model.mesh.setTransfiniteCurve(edge, meshSize[i])
        _, down = gmsh.model.getAdjacencies(1, edge)
        adjacencies.append(down)

    return adjacencies


def findT(curve, point):
    x = point[0]
    y = point[1]
    dom = curve.domain
    N = 200
    sample = np.linspace(dom[0], dom[1], N)
    points = curve.evaluate_list(sample)
    coordsX = np.array([j[0] for j in points])
    coordsY = np.array([j[1] for j in points])

    sqrdist = (coordsX-x)**2+(coordsY-y)**2
    k = np.argmin(sqrdist)

    t = sample[k]

    ders = curve.derivatives(t, order=1)
    X = ders[0][0]
    Y = ders[0][1]
    sqrdist = (X-x)**2+(Y-y)**2
    derivative = 2*(X-x)*ders[1][0]+2*(Y-y)*ders[1][1]

    maxiter = 100
    i = 0

    while np.sqrt(sqrdist) > 1e-8 and i < maxiter:
        i = i+1
        told = t
        t = told - sqrdist/derivative
        if t > 1 or t < 0:
            raise Exception("t out of bounds")
        ders = curve.derivatives(t, order=1)
        X = ders[0][0]
        Y = ders[0][1]
        sqrdist = (X-x)**2+(Y-y)**2
        derivative = 2*(X-x)*ders[1][0]+2*(Y-y)*ders[1][1]

    if i == maxiter:
        print("Max iterations reached")
        print("Iterations=", i)
    else:
        print("Converged in ", i, " iterations")
        print(f"Current distance: {np.sqrt(sqrdist):.2f}")

    return t


def borderPoints(curve):
    dom = curve.domain
    return curve.evaluate_list([dom[0], dom[1]])


def splitTwice(curve, points1):
    curves = multi.CurveContainer()
    curves.sample_size = 40

    T1 = findT(curve, min(points1[0], points1[1]))
    curve1, curve2 = operations.split_curve(curve, T1)
    T2 = findT(curve2, max(points1[0], points1[1]))
    curve3, curve4 = operations.split_curve(curve2, T2)

    curves.add(curve1)
    curves.add(curve3)
    curves.add(curve4)
    return curves


def differentPoints(curveA1, curveA2, curveB, curveC):
    pointsB = borderPoints(curveB)
    pointsA1 = borderPoints(curveA1)
    pointsA2 = borderPoints(curveA2)
    new_points = [point for point in pointsA2 + pointsA1 if point not in pointsB]
    curves = splitTwice(curveC, new_points)

    return curves


def mergeContainers(container1, container2, selectCurves=None):
    if selectCurves is not None:
        for i in selectCurves:
            container1.add(container2[i])
    else:
        for curve in container2:
            container1.add(curve)
    return container1


def pointsClose(p1, p2, rel_tol=1e-6):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2))) <= rel_tol


def pixelateBorder(entity):
    boundaryPoints = gmsh.model.getBoundary([(1, entity)], False, True)
    _, _, node_tags = gmsh.model.mesh.getElements(1, entity)
    nodes = node_tags[0][::-1]
    
    elementaryNodeTags = [boundaryPoints[1][1]]
    elementaryLineTags = []

    for i in range(1, len(nodes) - 1, 2):
        node = nodes[i]
        x, y, z = gmsh.model.mesh.getNode(node)[0]
        nodeTag = gmsh.model.occ.addPoint(x, y, z)
        elementaryNodeTags.append(nodeTag)

    elementaryNodeTags.append(boundaryPoints[0][1])

    for i, tag in enumerate(elementaryNodeTags[:-1]):
        lineTag = gmsh.model.occ.addLine(tag, elementaryNodeTags[i + 1])
        elementaryLineTags.append(lineTag)
    
    gmsh.model.occ.synchronize()

    for tag in elementaryLineTags:
        gmsh.model.addPhysicalGroup(1, [tag], -1, f"border_{tag}")

    return elementaryLineTags, elementaryNodeTags


def container2gmsh(bone, boneConfig, curves):

    numberElements = bone.mesh_vars.number_elements

    runFltk = boneConfig.runFltk

    pointsSet = getUniqueControlPoints(curves)
    pointTags = addPointsToModel(pointsSet)
    addCurvesToModel(curves, pointsSet, pointTags)
    gmsh.model.occ.synchronize()

    setGmshOptions()
    parameters = createTrasnfiniteParameters(numberElements)
    for param in parameters:
        createTransfiniteLine(param[0], param[1])

    gmsh.model.occ.remove(gmsh.model.occ.getEntities(0), True)
    gmsh.model.mesh.generate(1)

    createTransfiniteSurfaces(parameters)

    gmsh.model.mesh.generate(2)

    gmsh.model.mesh.renumberNodes()
    renumberElements()

    setSurfaceColors()
    addPhysicalGroups()

    if not boneConfig.skipWrite:
        writeMeshFiles(bone, boneConfig)

#    if runFltk:
#        gmsh.fltk.run()


def writeMeshFiles(bone, boneConfig):

    print("Writing .inp mesh")

    loadVars = bone.load_vars
    number_loads = loadVars.number_loads
    sigma = loadVars.load_extension
    rl = loadVars.number_loads

    inputPath = boneConfig.inputPath

    tags, coords, _ = gmsh.model.mesh.getNodes(-1, -1, False)
    allElements = gmsh.model.mesh.getElements()
    elemTypes, elemTags, elemNodeTags = allElements

    f2g.write_nodes(tags, coords, inputPath)
    f2g.write_connectivities(elemTypes, elemTags, elemNodeTags, inputPath)

    f2g.write_vtk(tags, coords, allElements, boneConfig.inputPath)

    physicalGroups_1D = gmsh.model.getPhysicalGroups(1)
    physicalGroups_2D = gmsh.model.getPhysicalGroups(2)

    f2g.writeBody(physicalGroups_2D, inputPath)
    f2g.writeBoundaries(physicalGroups_1D, inputPath)
    all2DElements = gmsh.model.mesh.getElements(2)

    if number_loads == 1:
        bone.load_vars.load_center = 0.0
        f2g.writeLoads(bone, boneConfig, [(1,3)], all2DElements, '')
    else:
        for i in range(0, number_loads):
            bone.load_vars.load_center = ((sigma - rl) * (2 * i - number_loads + 1)
                                        / (2 * number_loads - 2))

            f2g.writeLoads(bone, boneConfig, [(1,3)], all2DElements, str(i))

    f2g. writeParameters(bone, boneConfig, tags, all2DElements)


def getUniqueControlPoints(curves):
    pointsSet = []
    for curve in curves:
        for point in curve.ctrlpts:
            if not any(pointsClose(point, p) for p in pointsSet):
                pointsSet.append(tuple(point))
    return pointsSet


def addPointsToModel(pointsSet):
    pointTags = []
    for point in pointsSet:
        t = gmsh.model.occ.addPoint(point[0], point[1], point[2])
        pointTags.append(t)
    return pointTags


def addCurvesToModel(curves, pointsSet, pointTags):
    for curve in curves:
        ctrlPoints = curve.ctrlpts
        degree = curve.degree
        weights = curve.weights
        ctrlPointTags = []
        for ctrlPoint in ctrlPoints:
            for point, tag in zip(pointsSet, pointTags):
                if pointsClose(point, ctrlPoint):
                    ctrlPointTags.append(tag)
                    break
        gmsh.model.occ.addBSpline(ctrlPointTags, -1, degree, weights)


def setGmshOptions():
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
    gmsh.option.setNumber("Mesh.Algorithm", 8)


def createTrasnfiniteParameters(numberElements, pixelInfo=None):
    a2 = numberElements + 1
    a3 = 2 * numberElements + 1
    b = 2 * numberElements + 2

    parameters = [
        ([1, 2, 10, 5, 4, 3, 9], [a3, a3, b, a2, a3, a2, b], [1, 4, 5]),
        ([4, 13, 11, 12], [a3, a2, a3, a2]),
        ([3, 12, 14, 6], [a2, a2, a2, a2]),
        ([5, 8, 15, 13], [a2, a2, a2, a2]),
        ([11, 15, 7, 14], [a3, a2, a3, a2])
    ]

    if pixelInfo is not None:
        elements6, nodes6, elements7, nodes7, elements8, nodes8 = pixelInfo
        # Selection of nodes for border 6
        p16 = [3, 12, 14] + elements6
        p26 = [a2, a2, a2] + [1] * len(elements6)
        p36 = list(range(4, len(nodes6) + 2))

        # Selection of nodes for border 7
        p17 = [11, 15] + elements7 + [14]
        p27 = [a3, a2] + [1] * len(elements7) + [a2]
        p37 = list(range(3, len(nodes7) + 1))

        # Selection of nodes for border 8
        p18 = [5] + elements8 + [15, 13]
        p28 = [a2] + [1] * len(elements8) + [a2, a2]
        p38 = list(range(2, len(nodes8)))

        parameters = parameters[:-3]
        parameters.extend([
            (p16, p26, p36),
            (p17, p27, p37),
            (p18, p28, p38)
        ])

    return parameters


def createTransfiniteSurfaces(parameters):
    for param in parameters:
        if len(param) == 3:
            createTransfiniteSurface(param[0], param[1], ignorePoint=param[2])
        else:
            createTransfiniteSurface(param[0], param[1])


def renumberElements():
    allElements = gmsh.model.mesh.getElements()
    elemLineTags = allElements[1][0]
    elemQuadTags = allElements[1][1]
    elemPointTags = allElements[1][2]
    oldTags = np.concatenate((elemQuadTags, elemLineTags, elemPointTags))
    newTags = list(range(1, len(oldTags) + 1))
    gmsh.model.mesh.renumberElements(oldTags, newTags)


def setSurfaceColors():
    surfaces = gmsh.model.getEntities(2)
    for surface in surfaces:
        gmsh.model.setColor([surface], 64, 104, 177)
    gmsh.model.setColor([(2, 1)], 206, 166, 104)


def addPhysicalGroups():
    gmsh.model.addPhysicalGroup(2, [1], -1, "Hueso")
    gmsh.model.addPhysicalGroup(2, [2, 3, 4, 5], -1, "Cartilago")
    gmsh.model.addPhysicalGroup(1, [6, 7, 8], -1, "contorno2")
    gmsh.model.addPhysicalGroup(1, [1, 2], -1, "contorno1")
    gmsh.model.addPhysicalGroup(0, [2], -1, "pin")


def splitSymmetric(curve, t):
    curves = multi.CurveContainer()
    curve1, curve2 = operations.split_curve(curve, t)
    curve3, curve4 = operations.split_curve(curve2, (1-2*t)/(1-t))
    curves.add(curve1)
    curves.add(curve3)
    curves.add(curve4)
    return curves


def processSketchNurbs(local_sketch, boneConfig):
    renderRaw = boneConfig.renderRaw
    renderArea = boneConfig.renderArea
    renderMesh = boneConfig.renderMesh
    renderLength = boneConfig.renderLength

    curves0 = f2g.getBSplineGeom(local_sketch, 1000)
    curves1 = splitTwice(curves0[0], borderPoints(curves0[1]))
    curves2a = differentPoints(curves0[3], curves0[4], curves0[2], curves0[1])
    curves2b = differentPoints(curves0[5], curves0[6], curves0[2], curves1[1])

    curvesMesh0 = mergeContainers(curves2a, curves2b)
    curvesMesh0 = mergeContainers(curvesMesh0, curves1, selectCurves=[0, 2])
    curvesMesh0 = mergeContainers(curvesMesh0, curves0, selectCurves=[2, 3, 4, 5, 6, 7])
    curvesArea = mergeContainers(curves1, curves0, selectCurves=[1, 7])

    # curvesMesh0[13] is the curve to divide and delete
    curvesBottom = multi.CurveContainer()
    curveBottom1, curveBottom2 = operations.split_curve(curvesMesh0[13], 0.5)
    curvesBottom.add(curveBottom1)
    curvesBottom.add(curveBottom2)

    curvesMesh0 = mergeContainers(curvesBottom,
                                  curvesMesh0,
                                  selectCurves=[0, 1, 2, 3, 4, 5, 6,
                                                7, 8, 9, 10, 11, 12])

    render_flags = {
        renderRaw: curves0,
        renderMesh: curvesMesh0,
        renderArea: curvesArea,
        renderLength: splitSymmetric(curves0[0], 5/16)
    }

    for flag, curve in render_flags.items():
        if flag:
            curve.render()

    return curvesMesh0, curvesArea, render_flags[renderLength]
