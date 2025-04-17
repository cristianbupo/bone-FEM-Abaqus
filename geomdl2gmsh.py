import os
import numpy as np
from geomdl import  operations, multi
import gmsh
import FreeCAD2gmsh as f2g
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt


def createTransfiniteSurface(edges, meshSize, reverse=False, direction="Right", ignorePoint=0):
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

    if reverse:
        edges = [-edge for edge in edges] 

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

    return t


def borderPoints(curve):
    dom = curve.domain
    return curve.evaluate_list([dom[0], dom[1]])


def splitOnce(curve, point):
    curves = multi.CurveContainer()
    curves.sample_size = 40

    T = findT(curve, point)
    curve1, curve2 = operations.split_curve(curve, T)

    curves.add(curve1)
    curves.add(curve2)
    return curves


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


def extractCurves(container, selectCurves=None):
    curves = multi.CurveContainer()
    for i, curve in enumerate(container):
        if i in selectCurves:
            curves.add(curve)
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


def container2gmsh(bone, boneConfig, curvesMesh, curvesArea):

    # Main mesh

    numberElements = bone.mesh_vars.number_elements

    pointsSet = getUniqueControlPoints(curvesMesh)
    pointTags = addPointsToModel(pointsSet)
    addCurvesToModel(curvesMesh, pointsSet, pointTags)
    gmsh.model.occ.synchronize()

    setGmshOptions()
    parameters = createTrasnfiniteParameters(numberElements)
    for param in parameters:
        createTransfiniteLine(param[0], param[1])

    gmsh.model.occ.remove(gmsh.model.occ.getEntities(0), True)

    createTransfiniteSurfaces(parameters)
    
    addPhysicalGroups(boneConfig)

    gmsh.model.mesh.generate(2)

    gmsh.model.mesh.renumberNodes()
    renumberElements()
    setSurfaceColors()
    
    # Write mesh files

    inputPath = boneConfig.inputPath
    outputPath = boneConfig.outputPath

    tags, coords, _ = gmsh.model.mesh.getNodes(-1, -1, False)
    allElements = gmsh.model.mesh.getElements()
    elemTypes, elemTags, elemNodeTags = allElements

    f2g.write_nodes(tags, coords, inputPath, "nodos.inp", "restriccionesMultipunto.inp")
    f2g.write_connectivities(elemTypes, elemTags, elemNodeTags, inputPath, "conectividades.inp")
    f2g.write_vtk(tags, coords, allElements, inputPath, "malla.vtu")
    
    # physicalGroups_1D = gmsh.model.getPhysicalGroups(1)
    physicalGroups_2D = gmsh.model.getPhysicalGroups(2)
    physicalGroups = gmsh.model.getPhysicalGroups()

    f2g.writeBody(physicalGroups_2D, inputPath)
    lines = f2g.writeBoundaries(physicalGroups, inputPath)
    all2DElements = gmsh.model.mesh.getElements(2)

    listNElementLoads =[]

    # Create carga and resultado folder
    if not os.path.exists(os.path.join(outputPath, "carga")):
        os.makedirs(os.path.join(outputPath, "carga"))

    number_loads = bone.load_vars.number_loads
    number_steps = bone.time_vars.number_steps

    # Create carga and resultado PDV and VTM files
    srcPattern = os.path.join(outputPath, 'carga', 'carga*.vtp')
    destFile = os.path.join(outputPath, 'carga.vtm')
    createVTMbefore(srcPattern, destFile, number_loads, 1, shift=1)

    srcPattern = os.path.join(outputPath, 'carga', 'carga*.vtp')
    destFile = os.path.join(outputPath, 'carga.pvd')
    createPVDbefore(srcPattern, destFile, number_loads, 1, shift=1)

    srcPattern = os.path.join(outputPath, 'resultado', 'analisis*.vtu')
    destFile = os.path.join(outputPath, 'resultado.pvd')
    createPVDbefore(srcPattern, destFile, number_loads, 1, shift=1)
    
    loadCurve = curvesArea[1]

    h_vector, k_vector, r_vector = loadVectors(bone, loadCurve)

    listNElementLoads = []

    for i in range(number_loads):

        h = h_vector[i]
        k = k_vector[i]
        r = r_vector[i]
    
        nElementLoads = f2g.writeLoads(boneConfig, h, k, r, "contorno2", all2DElements, i)
        listNElementLoads.append(nElementLoads)


    writeParameters(bone, boneConfig, tags, all2DElements, lines, listNElementLoads)

    f2g.findNWriteMeshAdjacencies(all2DElements, inputPath)
    
    f2g.writeSteps(boneConfig, number_steps, number_loads)

    if boneConfig.runFltk:
        gmsh.fltk.run()

    if boneConfig.saveMsh:
        gmsh.write(os.path.join(inputPath,"malla.msh"))


def container2advanceMesh(bone, boneConfig, curvesAdvance):
    numberElements = bone.mesh_vars.number_elements

    a2, a3, _, _, _= meshLineElements(numberElements)

    horizontalElements = a2 + a3 + a2 - 2
    verticalElements = int(np.round(bone.geom_vars.total_advance // 0.1)) # each element is approximately 0.1 mm

    pointsSet = getUniqueControlPoints(curvesAdvance)
    pointTags = addPointsToModel(pointsSet)
    addCurvesToModel(curvesAdvance, pointsSet, pointTags)
    gmsh.model.occ.synchronize()

    setGmshOptions()

    parameters = [
        ([1, 2, 3, 4], [horizontalElements, verticalElements, horizontalElements, verticalElements]),
    ]

    for param in parameters:
        createTransfiniteLine(param[0], param[1])

    createTransfiniteSurfaces(parameters)

    gmsh.model.mesh.reverse([(2,7), (2,9)])

    gmsh.model.mesh.generate(2)
    

    if boneConfig.runFltk:
        gmsh.fltk.run()

    inputPath = boneConfig.inputPath
    # outputPath = boneConfig.outputPath

    tags, coords, _ = gmsh.model.mesh.getNodes(-1, -1, False)
    allElements = gmsh.model.mesh.getElements()
    elemTypes, elemTags, elemNodeTags = allElements

    f2g.write_nodes(tags, coords, inputPath, "nodosAdvance.inp")
    f2g.write_connectivities(elemTypes, elemTags, elemNodeTags, inputPath, "conectividadesAdvance.inp")
    f2g.write_vtk(tags, coords, allElements, inputPath, "mallaAdvance.vtu")
    if boneConfig.saveMsh:
        gmsh.write(os.path.join(inputPath,"mallaAvance.msh"))


def loadVectors(bone, loadCurve):

    number_loads = bone.load_vars.number_loads
    head_angle = bone.geom_vars.head_angle
    total_force = bone.load_vars.total_force
    # load_ext = bone.load_vars.load_extension
    # r = bone.load_vars.load_radius

    h_vector = np.zeros(number_loads)
    k_vector = np.zeros(number_loads)
    r_vector = np.ones(number_loads)

    maxLength, midLength, minLength = f2g.characteristicLengths(loadCurve)
    
    # if number_loads % 2 != 0: # Odd number of loads, ie. 5
    if True:                # Always active
        concave_length = minLength + 1 * (midLength - minLength) / 5
        convex_length = midLength + 3 * (maxLength - midLength) / 5


        if head_angle <= -15:
            load_ext = concave_length
        elif -15 < head_angle <= 0:
            load_ext = f2g.linear_interpolate(head_angle, -15, concave_length, 0, minLength)
        elif  0 < head_angle <= 15:
            load_ext = f2g.linear_interpolate(head_angle, 0, minLength, 15, convex_length)
        elif 15 < head_angle:
            load_ext = convex_length

        concave_radius = load_ext / 2
        convex_radius = load_ext / 4

        if head_angle <= -15:
            r = concave_radius
        elif -15 < head_angle <= 15:
            r = f2g.linear_interpolate(head_angle, -15, concave_radius, 15, convex_radius)
        elif 15 < head_angle:
            r = convex_radius

        k_vector = triangularPattern(number_loads, 0.5, 1.0)
        h_vector = np.linspace(-(load_ext-r)/2, (load_ext-r)/2, number_loads)
    else: # Even number of loads, ie. 6
        concaveExt = (minLength + 2*midLength)/3
        convexExt = (2*midLength + maxLength)/3

        concaveR = (maxLength - midLength)/3
        convexR = (maxLength - midLength)/3

        if head_angle <= -15:
            load_ext = concaveExt
            r = concaveR
        elif -15 < head_angle <= 15:
            load_ext = f2g.linear_interpolate(head_angle, -15, concaveExt, 15, convexExt)
            r = f2g.linear_interpolate(head_angle, -15, concaveR, 15, convexR)
        elif 15 < head_angle:
            load_ext = convexExt
            r = convexR

        k_vector = np.zeros(number_loads)
        k_vector[0:number_loads//2] = triangularPattern(number_loads//2, 0.5, 1.0)
        k_vector[number_loads//2:number_loads] = triangularPattern(number_loads//2, 0.5, 1.0)

        leftVal = load_ext/2-r
        rightVal = load_ext/2+r
        number = number_loads//2

        rightVal = load_ext/2+2*r

        h_vector[0:number_loads//2] = np.linspace(-rightVal, -leftVal, number)
        h_vector[number_loads//2:number_loads] = np.linspace(leftVal, rightVal, number)
    
    k_max = (3 * total_force)/(4 * r)
    k_vector = k_vector * k_max
    r_vector = r_vector * r

    return h_vector, k_vector, r_vector

def triangularPattern(size, min, max):
    vector = np.zeros(size)
    halfSize = size // 2
    for i in range(size):
        if i < halfSize:
            vector[i] = f2g.linear_interpolate(i, 0, min, halfSize, max)
        else:
            vector[i] = vector[size - i - 1]
    
    if size % 2 == 1:
        vector[halfSize] = max

    return vector


def writeOnFile(originFile, destinationFile, content):
    with open(originFile, 'r') as file, open(destinationFile, "w") as f:
        f_longBone_content = file.read()
        f.write(f_longBone_content.format(**content))


def writeParameters(bone, boneConfig, tags, all2DElements, lines, listNElementLoads):
    kOI = bone.oss_vars.kOI

    nLoads = bone.load_vars.number_loads
    nElems = len(all2DElements[1][0])
    numNode = len(tags)
    formatedNList = f"{' ,'.join(map(str, listNElementLoads))}"

    a2, a3, a4, b, c = meshLineElements(bone.mesh_vars.number_elements)

    originFile = os.path.join(boneConfig.masterPath, "parametros.txt")
    destinationFile = os.path.join(boneConfig.inputPath, "parametros.txt")
    
    # Get all physical groups
    physical_groups = gmsh.model.getPhysicalGroups()

    # Filter for 1D physical entities (dimension = 1)
    physical_1D_groups = [group for group in physical_groups if group[0] == 1]

    # Count the number of 1D physical entities
    num_1D_physical_entities = len(physical_1D_groups)

    content = {
        'stdWeight': bone.oss_vars.stdWeight,
        'nLoads': nLoads,
        'maxNElementLoads': max(listNElementLoads),
        'numNode': numNode,
        'nElems': nElems,
        'kOI': kOI,
        'nContornos': num_1D_physical_entities,
        'filasContorno1': (lines[0] + 5) // 6,
        'filasContorno2': (lines[1] + 5) // 6,
        'filasContorno3': (lines[2] + 5) // 6,
        'filasContorno4': (lines[3] + 5) // 6,
        'velocidad': bone.simulation_vars.growth_vel,
        'a2': a2-1,
        'a3': a3-1,
        'a4': a4-1,
        'b': b-1
    }

    with open(os.path.join(boneConfig.inputPath,'nFilasCargas.txt'), "w") as f:
        f.write(formatedNList)

    writeOnFile(originFile, destinationFile, content)


def createPVDbefore(srcPattern, destFile, numCombinations, numDigits, shift=0):
    # Ensure the destination directory exists
    destDir = os.path.dirname(destFile)
    os.makedirs(destDir, exist_ok=True)

    vtpFiles = []
    # Find all .vtp files matching the source pattern
    for i in range(numCombinations):
        vtpFile = srcPattern.replace('*', f"{i+shift:0{numDigits}d}")
        vtpFiles.append(vtpFile)

    # Create the root element
    vtkFile = ET.Element('VTKFile', type='Collection', version='1.0')

    collection = ET.SubElement(vtkFile, 'Collection')

    # Add DataSet entries
    for index, vtpFile in enumerate(vtpFiles):
        # Make the file path relative to the destination directory
        relativeInputPath = os.path.relpath(vtpFile, destDir)
        ET.SubElement(collection, 'DataSet', timestep=str(index), file=relativeInputPath)

    # Write the XML to the destination file
    tree = ET.ElementTree(vtkFile)
    tree.write(destFile, encoding='utf-8', xml_declaration=True)


def createVTMbefore(srcPattern, destFile, numCombinations, numDigits, shift = 0):
    # Ensure the destination directory exists
    destDir = os.path.dirname(destFile)
    os.makedirs(destDir, exist_ok=True)

    vtpFiles = []
    # Find all .vtp files matching the source pattern
    for i in range(numCombinations):
        vtpFile = srcPattern.replace('*', f"{i+shift:0{numDigits}d}")
        vtpFiles.append(vtpFile)

    # Create the root element
    vtkFile = ET.Element('VTKFile', type='vtkMultiBlockDataSet', version='1.0',
                         byte_order='LittleEndian', header_type='UInt32',
                         compressor='vtkZLibDataCompressor')

    vtkMultiBlockDataSet = ET.SubElement(vtkFile, 'vtkMultiBlockDataSet')

    # Add DataSet entries
    for index, vtpFile in enumerate(vtpFiles):
        # Make the file path relative to the destination directory
        relativeInputPath = os.path.relpath(vtpFile, destDir)
        ET.SubElement(vtkMultiBlockDataSet, 'DataSet', index=str(index + 1), file=relativeInputPath)

    # Write the XML to the destination file
    tree = ET.ElementTree(vtkFile)
    tree.write(destFile, encoding='utf-8', xml_declaration=True)


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


def meshLineElements(numberElements):
    a2 = numberElements + 1
    a3 = 2 * numberElements + 1
    a4 = 2 * numberElements + 1
    b = numberElements + 1
    c = 2 * numberElements + 1
    return a2, a3, a4, b, c

def createTrasnfiniteParameters(numberElements, pixelInfo=None):
    a2, a3, a4, b, c = meshLineElements(numberElements)

    parameters = [
        ([1, 2, 10, 5, 4, 3, 9], [a3, a3, b, a2, a3, a2, b], [1, 4, 5]),
        ([4, 13, 11, 12], [a3, a4, a3, a4]),
        ([3, 12, 14, 6], [a2, a4, a2, a4]),
        ([5, 8, 15, 13], [a2, a4, a2, a4]),
        ([11, 15, 7, 14], [a3, a2, a3, a2]),
        ([19, 26, 17, 8], [c, a4, c, a4]),
        ([18, 6, 16, 22], [c, a4, c, a4], True),
        ([21, 27, 19, 10], [c, b, c, b]),
        ([20, 9, 18, 24], [c, b, c, b], True),
        ([7, 17, 23, 16], [a3, c, a3, c]),
        ([25, 21, 2, 1, 20], [2*a3-1, c, a3, a3, c], [3])

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
            if isinstance(param[2], bool):
                createTransfiniteSurface(param[0], param[1], reverse=param[2])      
            else:
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

    physicalGroups = gmsh.model.getPhysicalGroups(2)
    colors = ([206, 166, 104], [64, 104, 177],[240, 100, 100])
    for i, physicalSurface in enumerate(physicalGroups):
        dim = physicalSurface[0]
        tag = physicalSurface[1]
        r = colors[i][0]
        g = colors[i][1]
        b = colors[i][2]
        surfaces = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)

        for surface in surfaces:
            gmsh.model.setColor([(2, surface)], r, g, b)


def addPhysicalGroups(boneConfig):
    gmsh.model.addPhysicalGroup(2, [1], -1, "Hueso")
    gmsh.model.addPhysicalGroup(2, range(2,6), -1, "Cartilago")

    r1 = range(6, 12)
    r2 = range(22,28)

    if boneConfig.capsule:
        gmsh.model.addPhysicalGroup(2, r1, -1, "Capsula")
        gmsh.model.addPhysicalGroup(1, r2, -1, "contorno0")
    else:
        gmsh.model.removeEntities([(2, i) for i in r1])
        gmsh.model.removeEntities([(1, i) for i in r2])
        gmsh.model.removeEntities([(1, i) for i in range(16,22)])
        gmsh.model.removeEntities([(0, i) for i in range(45,51)])

    gmsh.model.addPhysicalGroup(1, [1, 2], -1, "contorno1")
    gmsh.model.addPhysicalGroup(1, [6, 7, 8], -1, "contorno2")
    gmsh.model.addPhysicalGroup(1, [3, 4, 5], -1, "contorno3")
    gmsh.model.addPhysicalGroup(1, [7], -1, "contorno4")
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
    renderAdvance = boneConfig.renderAdvance

    curves0 = f2g.getBSplineGeom(local_sketch, 1000)
    curves1 = splitTwice(curves0[0], borderPoints(curves0[1]))
    curves2a = differentPoints(curves0[3], curves0[4], curves0[2], curves0[1])
    curves2b = differentPoints(curves0[5], curves0[6], curves0[2], curves1[1])

    curvesMesh0 = mergeContainers(curves2a, curves2b)
    curvesMesh0 = mergeContainers(curvesMesh0, curves1, selectCurves=[0, 2])
    curvesMesh0 = mergeContainers(curvesMesh0, curves0, selectCurves=[i for i in range(2, 21) if i != 8])
    curvesArea = mergeContainers(curves1, curves0, selectCurves=[1, 7])

    # curvesMesh0[13] is the curve to divide and delete
    curvesBottom = multi.CurveContainer()
    curveBottom1, curveBottom2 = operations.split_curve(curvesMesh0[13], 0.5)
    curvesBottom.add(curveBottom1)
    curvesBottom.add(curveBottom2)

    curvesMesh0 = mergeContainers(curvesBottom,
                                  curvesMesh0,
                                  selectCurves=[i for i in range(0, 26) if i != 13])

    curvesLength = splitSymmetric(curves0[0], 5/16)

    # Advance mesh

    curvesAdvance0 = multi.CurveContainer()
    curvesAdvance0.add(curves0[1])
    curvesAdvance0.add(curves0[8])

    curvesAdvance1 = splitTwice(curves0[0], borderPoints(curves0[8]))
    curveAdvance2 = splitOnce(curvesAdvance1[0], borderPoints(curves0[1])[0])[1]
    curveAdvance3 = splitOnce(curvesAdvance1[2], borderPoints(curves0[1])[1])[0]

    curvesAdvance1.add(curveAdvance2)
    curvesAdvance1.add(curveAdvance3)
    # curvesAdvance = mergeContainers(curvesAdvance0, curvesAdvance1)
    # curvesAdvance.add(curveAdvance2)

    curvesAdvance = multi.CurveContainer()
    curvesAdvance.add(curvesAdvance0[0])
    curvesAdvance.add(curvesAdvance1[4])
    curvesAdvance.add(curvesAdvance0[1])
    curvesAdvance.add(curvesAdvance1[3])


    render_flags = {
        'renderRaw': (renderRaw, curves0),
        'renderMesh': (renderMesh, curvesMesh0),
        'renderArea': (renderArea, curvesArea),
        'renderLength': (renderLength, curvesLength),
        'renderAdvance': (renderAdvance, curvesAdvance)
    }


    i = 1
    for key, (flag, container) in render_flags.items():
        if flag:
            drawContainer(container)
            plt.title(f"Container: {key}")
            plt.show()

            i += 1

    return curvesMesh0, curvesArea, curvesAdvance


def drawContainer(container, controlPol=False, knotEval=False):
    fig, ax = plt.subplots()

    for i, curve in enumerate(container):
        knotVector = curve.knotvector
        ctrlPoints = np.array(curve.ctrlpts)
        evaluated_points = np.array(curve.evalpts)
        
        # Plot the evaluated points
        plt.plot(evaluated_points[:, 0], evaluated_points[:, 1], label=f'curve {i}')
        
        # Plot the control points if controlPol is True
        if controlPol:
            plt.plot(ctrlPoints[:, 0], ctrlPoints[:, 1], 'o--', markersize=3, linewidth=0.5, label='Control Points', color='black')       
        # Evaluate and plot the spline at the knots
        if knotEval:
            knot_evaluations = np.array([curve.evaluate_single(knot) for knot in knotVector])
            plt.plot(knot_evaluations[:, 0], knot_evaluations[:, 1], 'x', label='Knot Evaluations', color='red')
        
    
    # Position the legend to the right of the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.0, frameon=False)
    
    # Set aspect ratio and turn off axis
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    return fig, ax