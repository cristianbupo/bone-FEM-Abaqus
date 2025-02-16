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
        f.write("*NODE,NSET=N2\n")
        for n in range(numNode):
            f.write(f"{nodes[n]}, {coords[3 * n]}, {coords[3 * n + 1]}\n")


def write_connectivities(elemTypes, elemTags, elemNodeTags, inputPath):
    conectivitidadesVerPath = os.path.join(inputPath, "conectividades.inp")
    with open(conectivitidadesVerPath, "w") as h:
        h.write("*ELEMENT,TYPE=U1,ELSET=UEL\n")
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


def writeBody(physicalGroups, writingPath):
    cuerpoPath = os.path.join(writingPath, "cuerpo.inp")
    gruposFisicosPath = os.path.join(writingPath, "gruposFisicos.txt")

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


def writeBoundaries(physicalGroups, inputPath):
    lines = []
    contornoPath = os.path.join(inputPath, "contorno.inp")
    boundaryConditionsPath = os.path.join(inputPath, "condicionesContorno.inp")

    with open(contornoPath, "w") as f, open(boundaryConditionsPath, "w") as g:
        g.write("*Boundary\n")
        g.write("contorno1, 11, 12, 0.0\n")
        for dim, tag in physicalGroups:
            name = gmsh.model.getPhysicalName(dim, tag)
            nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
            f.write("*NSET,NSET="+name+"\n")
            writeLines(nodeTags, f)

            lines.append(len(nodeTags))

            # if name == "contorno1":
                # pinTag = int(nodeTags[(len(nodeTags)+1) // 2])
                # g.write(f"{pinTag}, 1, 1\n")
    
    return lines


def findLastPosition(A, B):
    # A and B are arrays, B is a 2-element array
    n = len(A)
    for i in range(1, n):
        if (A[i] == B[0] and A[i-1] == B[1]) or (A[i] == B[1] and A[i-1] == B[0]):
            return i
    if (A[0] == B[0] and A[n-1] == B[1]) or (A[0] == B[1] and A[n-1] == B[0]):
        return n
    return -1  # return -1 if no match found


def writeLoads(bone, boneConfig, physicalName, all2DElements, load_index):

    physicalGroups = gmsh.model.getPhysicalGroups()
    physicalGroup = []
    for dim, tag in physicalGroups:
        name = gmsh.model.getPhysicalName(dim, tag)
        if name == physicalName:
            physicalGroup.append((dim, tag))

    inputPath = boneConfig.inputPath

    h = bone.load_vars.load_center
    k = bone.load_vars.load_amplitude
    r = bone.load_vars.load_radius

    elemTypes = all2DElements[0][0]
    allElemTags = all2DElements[1][0]
    allElemNodeTags = all2DElements[2][0]
    elemNodes = elemTypes + 1

    for dim, tag in physicalGroup:
        # Calculate load Vector
        # Get contour elements and nodes of physical group
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)

        pContourElements = []
        pElemTags = []
        pElemNodeTags = []
        pLoadFaces = []
        pNodCoords = []

        for i, entity in enumerate(entities):
            # Get contour elements and nodes of a single entity
            _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, entity)
            elemTags = elemTags[0]
            elemNodeTags = elemNodeTags[0]

            j = 0
            contourElements = np.zeros_like(elemTags)  # From body
            contourElementNodes = np.zeros(len(contourElements) * elemNodes)
            while j < len(elemTags):
                i = 0
                searchedVector = elemNodeTags[2*j:2*(j+1)]
                condition = all(
                    elem in allElemNodeTags[0:elemNodes] for elem in searchedVector
                )

                while i < len(allElemTags) and not condition:
                    i = i+1
                    slicedArr = allElemNodeTags[elemNodes*i:elemNodes*(i+1)]
                    condition = all(elem in slicedArr for elem in searchedVector)

                contourElements[j] = allElemTags[i]
                contourElementNodes[elemNodes*j:elemNodes*(j+1)] = (
                    allElemNodeTags[elemNodes*i:elemNodes*(i+1)])
                j = j + 1

            firstContourElement = elemNodeTags[0:2]
            firstBodyElement = contourElementNodes[0:elemNodes]

            loadFace = findLastPosition(firstBodyElement, firstContourElement)

            pContourElements.extend(contourElements)
            pElemTags.extend(elemTags)
            pElemNodeTags.extend(elemNodeTags)
            pLoadFaces.extend([loadFace] * len(elemTags))

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

        nElementLoads = 0
        npElementTags = len(pElemTags)
        p = -r**2/(4*k)
        elemLenghts = np.zeros(npElementTags)
        distributionPath = os.path.join(inputPath, f"carga{load_index}.inp")
        with open(distributionPath, "w") as g:
            for i in range(npElementTags):
                lenght = math.sqrt(
                    (pNodCoords[2*i][0] - pNodCoords[2*i+1][0])**2 +
                    (pNodCoords[2*i][1] - pNodCoords[2*i+1][1])**2
                )
                elemLenghts[i] = lenght

            for i in range(npElementTags):
                midpoint = [(pNodCoords[2*i][0] + pNodCoords[2*i+1][0])/2,
                            (pNodCoords[2*i][1] + pNodCoords[2*i+1][1])/2, 0]
                midpoints.InsertNextPoint(midpoint)

                normal = [(pNodCoords[2*i+1][1] - pNodCoords[2*i][1])/elemLenghts[i],
                          -(pNodCoords[2*i+1][0] - pNodCoords[2*i][0])/elemLenghts[i],
                          0]
                normals.InsertNextTuple(normal)

            g.write("*DLOAD\n")
            L = np.sum(elemLenghts)

            x_im1 = -L/2
            for j in range(npElementTags):
                pressure_i = 0.0
                x_i = x_im1 + elemLenghts[j]
                if (h-r < x_i < h+r) or (h-r < x_im1 < h+r):
                    a = max(x_im1, h-r)
                    b = min(x_i, h+r)
                    pressure_i = (((b-h)**3-(a-h)**3)/(12*p)+k*(b-a))/elemLenghts[j]
                    g.write(f"{pContourElements[j]}, U{pLoadFaces[j]}, {pressure_i}\n")
                    nElementLoads += 1
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

            return nElementLoads


def writeParameters(bone, boneConfig, tags, all2DElements, lines, listNElementLoads):
    kOI = bone.load_vars.kOI

    nLoads = bone.load_vars.number_loads
    nElems = len(all2DElements[1][0])
    numNode = len(tags)
    formatedNList = f"(/{' ,'.join(map(str, listNElementLoads))}/)"
    fillPath = os.path.join(boneConfig.masterPath, "conec.for")
    parametrosPath = os.path.join(boneConfig.inputPath, "conec.for")

    with open(fillPath, 'r') as file:
        f_longbone_content = file.read()

    modified_content = f_longbone_content.format(
        nLoads=nLoads,
        listNElementLoads=formatedNList,
        numNode=numNode,
        nElems=nElems,
        kOI=kOI,
        filasContorno1=(lines[0] + 5) // 6,
        filasContorno2=(lines[1] + 5) // 6
    )

    with open(parametrosPath, "w") as f:
        f.write(modified_content)