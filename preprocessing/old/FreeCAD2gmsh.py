import gmsh
import math
import numpy as np
from math import ceil
import vtk
import os
from geomdl import NURBS, multi

points = []


def calculateKnotVector(ctrlPoints, degree):
    NKnots = len(ctrlPoints) + degree + 1
    knots = [0.0] * NKnots
    for i in range(degree+1):
        knots[i] = 0.0
        knots[NKnots-1-i] = 1.0

    div = NKnots-2*degree
    for i in range(div-2):
        knots[i+degree+1] = (i+1)*1.0/(div-1)

    return knots


def getBSplineGeom(sketch, nPoints):
    # Get the geometry information of BSplines
    # global points
    Geometry = sketch.Geometry
    GeometryFacadeList = sketch.GeometryFacadeList

    curves = multi.CurveContainer()

    for i in range(sketch.GeometryCount):
        if not GeometryFacadeList[i].Construction:
            iGeometry = Geometry[i]
            iclass = iGeometry.__class__.__name__
            ctrlPoints = []
            weights = []

            if iclass == 'BSplineCurve':
                NbPoles = iGeometry.NbPoles
                degree = iGeometry.Degree

                j = 1
                k = i + 1
                while j <= NbPoles:
                    kGeometry = Geometry[k - NbPoles - 1]

                    ctrlPointCoords = kGeometry.Center
                    ctrlPointWeight = kGeometry.Radius

                    weights.append(ctrlPointWeight)

                    ctrlPoints.append(np.array([ctrlPointCoords.x,
                                                ctrlPointCoords.y,
                                                ctrlPointCoords.z]))

                    j += 1
                    k += 1

                knotVector = calculateKnotVector(ctrlPoints, degree)

            elif iclass == 'LineSegment':
                sp = iGeometry.StartPoint
                ep = iGeometry.EndPoint

                ctrlPoints = []
                ctrlPoints.append(np.array([sp.x, sp.y, sp.z]))
                ctrlPoints.append(np.array([ep.x, ep.y, ep.z]))
                degree = 1
                weights = [1, 1]
                knotVector = [1, 1]
                knotVector = calculateKnotVector(ctrlPoints, degree)

            ctrlPoints = np.array(ctrlPoints)

            curve = NURBS.Curve()
            curve.degree = degree
            curve.ctrlpts = ctrlPoints
            curve.weights = weights
            curve.delta = 1.0 / nPoints
            curve.knotvector = knotVector
            curves.add(curve)

    return curves

# MESH


def appendEntity(e, m):
    dim = e[0]
    tag = abs(e[1])
    bnd = gmsh.model.getBoundary([e])
    nod = gmsh.model.mesh.getNodes(dim, tag)
    ele = gmsh.model.mesh.getElements(dim, tag)
    m[(dim, tag)] = (bnd, nod, ele)

    if dim > 0:
        for b in bnd:
            appendEntity(b, m)
    return m

# Writing INP


def write_nodes(nodes, coords, inputPath):
    numNode = len(nodes)
    nodosPath = os.path.join(inputPath, "nodos.inp")
    with open(nodosPath, "w") as f:
        f.write("*NODE, NSET=N2\n")
        for n in range(numNode):
            f.write(f"{nodes[n]}, {coords[3 * n]}, {coords[3 * n + 1]}\n")


def write_connectivities(elemTypes, elemTags, elemNodeTags, inputPath):
    conectivitidadesVerPath = os.path.join(inputPath, "conectividades.inp")
    with open(conectivitidadesVerPath, "w") as h:
        h.write("*ELEMENT, TYPE=CPE4, ELSET=UEL\n")
        for i, _ in enumerate(elemTypes):
            for elem, item in enumerate(elemTags[i]):
                elemNum = int(item)
                if elemTypes[i] == 2:
                    h.write(f"{elemNum}, {elemNodeTags[i][3 * elem]}, "
                            f"{elemNodeTags[i][3 * elem + 1]}, "
                            f"{elemNodeTags[i][3 * elem + 2]}\n")
                elif elemTypes[i] == 3:
                    h.write(f"{elemNum}, {elemNodeTags[i][4 * elem]}, "
                            f"{elemNodeTags[i][4 * elem + 1]}, "
                            f"{elemNodeTags[i][4 * elem + 2]}, "
                            f"{elemNodeTags[i][4 * elem + 3]}\n")


def write_vtk(nodes, coords, allElements, inputPath):
    elemTypes = allElements[0]
    elemTags = allElements[1]
    elemNodeTags = allElements[2]
    vtk_points = vtk.vtkPoints()
    for i in range(len(nodes)):
        vtk_points.InsertNextPoint(coords[3 * i], coords[3 * i + 1], 0.0)

    unstructuredGrid = vtk.vtkUnstructuredGrid()
    unstructuredGrid.SetPoints(vtk_points)

    for i, item in enumerate(elemTypes):
        if item == 2:  # Triangle
            for elem_index, _ in enumerate(elemTags[i]):
                triangle = vtk.vtkTriangle()
                for j in range(3):  # Triangle has 3 vertices
                    start_index = elem_index * 3
                    node_id = int(elemNodeTags[i][start_index + j]) - 1
                    triangle.GetPointIds().SetId(j, node_id)
                unstructuredGrid.InsertNextCell(triangle.GetCellType(),
                                                triangle.GetPointIds())
        elif item == 3:  # Quadrangle
            for elem_index, _ in enumerate(elemTags[i]):
                quad = vtk.vtkQuad()
                for j in range(4):  # Quadrangle has 4 vertices
                    start_index = elem_index * 4
                    node_id = int(elemNodeTags[i][start_index + j]) - 1
                    quad.GetPointIds().SetId(j, node_id)
                unstructuredGrid.InsertNextCell(quad.GetCellType(),
                                                quad.GetPointIds())

    writer = vtk.vtkXMLUnstructuredGridWriter()
    vtuPath = os.path.join(inputPath, "malla.vtu")
    writer.SetFileName(vtuPath)
    writer.SetInputData(unstructuredGrid)
    writer.Write()


def writeLines(array, file):
    Ncols = 6
    N = len(array)
    N2 = Ncols*ceil(N/Ncols)
    for elem in range(N2):
        if elem < N:
            file.write(f"{array[elem]}, ")
        else:
            file.write("0, ")
        if (elem+1) % Ncols == 0:
            file.write("\n")


def writeBody(physicalGroups, inputPath):
    cuerpoPath = os.path.join(inputPath, "cuerpo.inp")
    gruposFisicosPath = os.path.join(inputPath, "gruposFisicos.txt")

    with open(cuerpoPath, "w") as f, open(gruposFisicosPath, "w") as g:
        g.write('Element Tag, Physical Group Tag\n')
        for dim, tag in physicalGroups:
            name = gmsh.model.getPhysicalName(dim, tag)
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
            f.write(f"*Nset, nset={name}\n")
            writeLines(nodeTags, f)

            f.write("*Elset, elset="+name+"\n")

            allElemTags = []
            for e in entities:
                _, elemTags, _ = gmsh.model.mesh.getElements(dim, e)
                for elem in range(0, len(elemTags[0])):
                    allElemTags.append(elemTags[0][elem])

            writeLines(allElemTags, f)

            if dim == 2:
                # Get the entities that belong to the physical group
                entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
                # For each entity
                for entity in entities:
                    # Get the elements that belong to the entity
                    _, element_tags, _ = gmsh.model.mesh.getElements(dim, entity)
                    # For each element
                    for element_tag in element_tags[0]:
                        # Write the physical group tag to the text file
                        g.write(f'{element_tag}, {tag}\n')


def writeBoundaries(physicalGroup, inputPath):
    contornoPath = os.path.join(inputPath, "contorno.inp")
    boundaryConditionsPath = os.path.join(inputPath, "condicionesContorno.inp")

    with open(contornoPath, "w") as f, open(boundaryConditionsPath, "w") as g:
        g.write("*Boundary\n")
        g.write("bottom, 2\n")
        for dim, tag in physicalGroup:
            name = gmsh.model.getPhysicalName(dim, tag)
            nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
            f.write("*Nset, nset="+name+"\n")
            writeLines(nodeTags, f)

            if name == "bottom":
                pinTag = int(nodeTags[(len(nodeTags)+1) // 2])
                g.write(f"{pinTag}, 1, 1\n")


def findLastPosition(A, B):
    # A and B are arrays, B is a 2-element array
    n = len(A)
    for i in range(1, n):
        if (A[i] == B[0] and A[i-1] == B[1]) or (A[i] == B[1] and A[i-1] == B[0]):
            return i
    if (A[0] == B[0] and A[n-1] == B[1]) or (A[0] == B[1] and A[n-1] == B[0]):
        return n
    return -1  # return -1 if no match found


def writeLoads(bone, boneConfig, physicalName, load_index):

    physicalGroups = gmsh.model.getPhysicalGroups()
    newEntities = []
    for dim, tag in physicalGroups:
        name = gmsh.model.getPhysicalName(dim, tag)
        if name == physicalName:
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            for entity in entities:
                newEntities.append((dim, entity))

    inputPath = boneConfig.inputPath

    h = - bone.load_vars.load_center  # Negative sign to match right as positive
    k = bone.load_vars.load_amplitude
    r = bone.load_vars.load_radius

    pNodCoords = []

    for i, entity in enumerate(newEntities):
        # Get contour elements and nodes of a single entity
        _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(1, entity)
        elemTags = elemTags[0]
        elemNodeTags = elemNodeTags[0]

        for nodeTag in elemNodeTags:
            coord, _, _, _ = gmsh.model.mesh.getNode(nodeTag)
            pNodCoords.append(coord)

    # Calculate midpoints and normals for the outer surface
    midpoints = vtk.vtkPoints()
    normals = vtk.vtkDoubleArray()
    load_magnitudes = vtk.vtkDoubleArray()
    normals.SetNumberOfComponents(3)  # Normal vectors are 3D
    normals.SetName("Load Orientation")
    load_magnitudes.SetName("Load Magnitude")

    n = len(newEntities)
    p = -r**2/(4*k)
    elemLenghts = np.zeros(n)
    distributionPath = os.path.join(inputPath, f"distribution{load_index}.txt")
    with open(distributionPath, "w") as g:
        for i in range(n):
            lenght = math.sqrt(
                (pNodCoords[2*i][0] - pNodCoords[2*i+1][0])**2 +
                (pNodCoords[2*i][1] - pNodCoords[2*i+1][1])**2
            )
            elemLenghts[i] = lenght

        for i in range(n):
            midpoint = [(pNodCoords[2*i][0] + pNodCoords[2*i+1][0])/2,
                        (pNodCoords[2*i][1] + pNodCoords[2*i+1][1])/2, 0]
            midpoints.InsertNextPoint(midpoint)

            normal = [(pNodCoords[2*i][1] - pNodCoords[2*i+1][1])/elemLenghts[i],
                      -(pNodCoords[2*i][0] - pNodCoords[2*i+1][0])/elemLenghts[i], 0]
            normals.InsertNextTuple(normal)

        L = np.sum(elemLenghts)

        x_im1 = -L/2
        for j in range(n):
            pressure_i = 0.0
            x_i = x_im1 + elemLenghts[j]
            if (h-r < x_i < h+r) or (h-r < x_im1 < h+r):
                a = max(x_im1, h-r)
                b = min(x_i, h+r)
                pressure_i = (((b-h)**3-(a-h)**3)/(12*p)+k*(b-a))/elemLenghts[j]
                g.write(f"{newEntities[j]}, {pressure_i}\n")
            load_magnitudes.InsertNextValue(pressure_i)
            x_im1 = x_i

        # Writing the vtp for the load distribution
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(midpoints)
        polydata.GetPointData().SetNormals(normals)
        polydata.GetPointData().AddArray(load_magnitudes)

        # Write vtkPolyData to a file
        writer = vtk.vtkXMLPolyDataWriter()
        file_path = os.path.join(inputPath, f"distribution{load_index}.vtp")
        writer.SetFileName(file_path)
        writer.SetInputData(polydata)
        writer.Write()


def writeParameters(bone, boneConfig, tags, all2DElements):
    kOI = bone.load_vars.kOI

    nElems = len(all2DElements[1][0])
    numNode = len(tags)
    parametrosPath = os.path.join(boneConfig.inputPath, "fortVars.for")

    with open(parametrosPath, "w") as f:
        f.write(f"      integer, parameter :: NUMNODE={numNode}, NELEMS={nElems}, "
                "dim=2, nnod=4\n"
                "      real*8 nodes(NUMNODE, dim)\n"
                f"      real*8, parameter :: kOI={kOI}\n"
                "      integer conectividades(NELEMS, nnod+1)\n"
                "      integer grupoFisico(NELEMS, 2)\n"
                "      common nodes\n"
                "      common conectividades\n"
                "      common grupoFisico\n")
