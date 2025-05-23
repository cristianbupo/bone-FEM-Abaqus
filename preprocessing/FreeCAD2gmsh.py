import gmsh
import math
import numpy as np
from math import ceil
import vtk
import os
import sketchUtils as su
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


def write_nodes(nodes, coords, inputPath, filename1, filename2=None):
    numNode = len(nodes)
    nodosPath = os.path.join(inputPath, filename1)
    
    with open(nodosPath, "w") as f:
        f.write("*NODE,NSET=N2\n")
        for n in range(numNode):
            f.write(f"{nodes[n]}, {coords[3 * n]}, {coords[3 * n + 1]}\n")

    if (filename2 != None):
        restriccionesPath = os.path.join(inputPath, filename2)
        with open(restriccionesPath, "w") as g:
            g.write("*MPC, user, mode=dof\n")
            for n in range(numNode):
                g.write(f"{nodes[n]}, {nodes[n]}, {nodes[n]}\n")


def write_connectivities(elemTypes, elemTags, elemNodeTags, inputPath, filename):
    conectivitidadesVerPath = os.path.join(inputPath, filename)
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


def write_vtk(nodes, coords, allElements, inputPath, filename):
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
    vtuPath = os.path.join(inputPath, filename)
    writer.SetFileName(vtuPath)
    writer.SetInputData(unstructuredGrid)
    writer.Write()


def writeLines(array, file):
    Ncols = 15
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

    with open(contornoPath, "w") as f:
        for dim, tag in physicalGroups:
            if dim < 3:
                name = gmsh.model.getPhysicalName(dim, tag)
                nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
                f.write("*NSET,NSET="+name+"\n")
                writeLines(nodeTags, f)

                lines.append(len(nodeTags))
    
    return lines


def characteristicLengths(loadCurve):
    knots = loadCurve.knotvector
    Nknots = len(knots)
    midKnotIndex = Nknots//2
    
    if knots[midKnotIndex] != 0.5:
        if knots[midKnotIndex + 1] == 0.5:
            knots[midKnotIndex] = 0.5
        elif knots[midKnotIndex - 1] == 0.5:
            knots[midKnotIndex] = 0.5


    loadKnots = knots[midKnotIndex-3 : midKnotIndex+4]

    maxLength = su.findLenght(loadCurve, loadKnots[0], loadKnots[6])
    midLength = su.findLenght(loadCurve, loadKnots[1], loadKnots[5])
    minLength = su.findLenght(loadCurve, loadKnots[2], loadKnots[4])

    return maxLength, midLength, minLength


def linear_interpolate(x, x0, y0, x1, y1):
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def findLastPosition(A, B):
    # A and B are arrays, B is a 2-element array
    n = len(A)
    for i in range(1, n):
        if (A[i] == B[0] and A[i-1] == B[1]) or (A[i] == B[1] and A[i-1] == B[0]):
            return i
    if (A[0] == B[0] and A[n-1] == B[1]) or (A[0] == B[1] and A[n-1] == B[0]):
        return n
    return -1  # return -1 if no match found


def findPhysicalGroup(physicalName):
    physicalGroups = gmsh.model.getPhysicalGroups()
    dim = None
    tag = None
    for d, t in physicalGroups:
        name = gmsh.model.getPhysicalName(d, t)
        if name == physicalName:
            dim = d
            tag = t
            break
    
    if dim is None or tag is None:
        raise ValueError(f"Physical group with name '{physicalName}' not found.")

    return dim, tag


def findConectivityInfo(dim, entity, allElemTags, allElemNodeTags, elemNodes):
    # Get contour elements and nodes of a single entity
    _, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, entity)
    elemTags = elemTags[0]
    nElementTags = len(elemTags)
    elemNodeTags = elemNodeTags[0]
    nodCoords = []
    elemLengths = np.zeros(nElementTags)

    for i, nodeTag in enumerate(elemNodeTags):
        coord, _, _, _ = gmsh.model.mesh.getNode(nodeTag)
        nodCoords.append(coord)

    j = 0
    contourElements = np.zeros_like(elemTags)  # From body
    contourElementNodes = np.zeros(len(contourElements) * elemNodes)
    for j in range(nElementTags):
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
        
        lenght = math.sqrt(
            (nodCoords[2*j][0] - nodCoords[2*j+1][0])**2 +
            (nodCoords[2*j][1] - nodCoords[2*j+1][1])**2
        )
        elemLengths[j] = lenght

    firstContourElement = elemNodeTags[0:2]
    firstBodyElement = contourElementNodes[0:elemNodes]

    loadFace = findLastPosition(firstBodyElement, firstContourElement)

    loadFaces = [loadFace] * len(elemTags)

    return contourElements, elemTags, elemNodeTags, loadFaces, nodCoords, elemLengths


def findConectivityInfoPhysical(physicalName, all2DElements):
    dim, tag = findPhysicalGroup(physicalName)
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)

    elemTypes = all2DElements[0][0]
    allElemTags = all2DElements[1][0]
    allElemNodeTags = all2DElements[2][0]
    elemNodes = elemTypes + 1

    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)

    pContourElements = []
    pElemTags = []
    pElemNodeTags = []
    pLoadFaces = []
    pNodCoords = []
    pElemLengths = []

    for entity in entities:
        (
            contourElements, elemTags, elemNodeTags, loadFaces, nodCoords, elemLengths
        ) = findConectivityInfo(dim, entity, allElemTags, allElemNodeTags, elemNodes)

        pContourElements.extend(contourElements)
        pElemTags.extend(elemTags)
        pElemNodeTags.extend(elemNodeTags)
        pLoadFaces.extend(loadFaces)
        pNodCoords.extend(nodCoords)
        pElemLengths.extend(elemLengths)

    return pContourElements, pElemTags, pElemNodeTags, pLoadFaces, pNodCoords, pElemLengths


def findMeshAdjacencies(all2DElements):
    allElemTags = all2DElements[1][0]
    allElemNodeTags = all2DElements[2][0]
    nElem = len(allElemTags)

    print(f"Number of 2D elements: {nElem}")

    # Reshape node_tags into quads
    quads = [allElemNodeTags[i:i+4] for i in range(0, len(allElemNodeTags), 4)]

    # Initialize adjacency array (nElem x 5)
    adjacencyArray = np.zeros((nElem, 5), dtype=int)
    adjacencyArray[:, 0] = np.arange(1, nElem + 1)  # First column: element IDs

    # Build edge: list of elements sharing that edge
    edge2elems = {}

    def orderedEdge(n1, n2):
        return tuple(sorted((n1, n2)))

    for eid, quad in enumerate(quads, start=1):  # Element IDs start at 1
        edges = [
            orderedEdge(quad[0], quad[1]),
            orderedEdge(quad[1], quad[2]),
            orderedEdge(quad[2], quad[3]),
            orderedEdge(quad[3], quad[0])
        ]
        for edge in edges:
            if edge not in edge2elems:
                edge2elems[edge] = []
            edge2elems[edge].append(eid)

    # Build element: neighbor per face
    for eid, quad in enumerate(quads, start=1):
        edges = [
            orderedEdge(quad[0], quad[1]),
            orderedEdge(quad[1], quad[2]),
            orderedEdge(quad[2], quad[3]),
            orderedEdge(quad[3], quad[0])
        ]
        neighbors = []
        for edge in edges:
            elems = edge2elems[edge]
            neighbor = [e for e in elems if e != eid]
            if neighbor:
                neighbors.append(neighbor[0])
            else:
                neighbors.append(0)  # Boundary face
        adjacencyArray[eid - 1, 1:] = neighbors  # Fill connectivity columns

    return adjacencyArray


def writeMeshAdjacencies(adjacencyArray, inputPath):
    # Write connectivity to file
    with open(os.path.join(inputPath, "adyacenciaElementos.inp"), "w") as f:
        f.write("elem, face1, face2, face3, face4\n")
        for row in adjacencyArray:
            f.write(", ".join(map(str, row)) + "\n")

    return adjacencyArray
    


def writeLoads(boneConfig, h, k, r, physicalName, all2DElements, loadIndex):
    inputPath = boneConfig.inputPath

    (
        pContourElements, pElemTags, _, pLoadFaces, pNodCoords, pElemLengths
    ) = findConectivityInfoPhysical(physicalName, all2DElements)

    loadElements, loadFaces, pressureDist, nElementLoads = calculatePressureDist(
        h, k, r, pContourElements, pLoadFaces, pElemLengths, pElemTags
    )

    writePressureDist(
        loadElements, loadFaces, pressureDist, nElementLoads, loadIndex, inputPath
    )

    writePressureDistVTP(
         boneConfig, pNodCoords, pContourElements, loadElements, pressureDist, pElemLengths, loadIndex
    )

    return nElementLoads


def calculatePressureDist(h, k, r, pContourElements, pLoadFaces, pElemLengths, pElemTags):
    L = np.sum(pElemLengths)
    loadElements = []
    loadFaces = []
    pressureDist =[]
    x_im1 = -L/2 # x_im1 is the previous x value
    p = -r**2/(4*k)
    
    nElementLoads = 0
    npElementTags = len(pElemTags)

    for i in range(npElementTags):
        x_i = x_im1 + pElemLengths[i]
        if (h-r < x_i < h+r) or (h-r < x_im1 < h+r):
            a = max(x_im1, h-r)
            b = min(x_i, h+r)            
            loadElements.append(pContourElements[i])
            loadFaces.append(pLoadFaces[i])
            pressureDist.append((((b-h)**3-(a-h)**3)/(12*p)+k*(b-a))/pElemLengths[i])
            nElementLoads += 1
        x_im1 = x_i
    return loadElements, loadFaces, pressureDist, nElementLoads


def writePressureDist(loadElements, loadFaces, pressureDist, nElementLoads, loadIndex, inputPath):
    with open(os.path.join(inputPath, "carga.inp"), "a") as g:
        g.write(f"*Step, name=LoadStep{loadIndex+1}\n")
        g.write("*DLOAD\n")
        for i in range(nElementLoads):
            g.write(f"{loadElements[i]}, {loadFaces[i]}, {pressureDist[i]}\n")
        g.write("*End Step\n")


def writePressureDistVTP(boneConfig, pNodCoords, pContourElements, loadElements, pressureDist, pElemLengths, loadIndex):
    outputPath = boneConfig.outputPath
    # Calculate midpoints and normals for the outer surface
    midpoints = vtk.vtkPoints()
    normals = vtk.vtkDoubleArray()
    load_magnitudes = vtk.vtkDoubleArray()
    normals.SetNumberOfComponents(3)  # Normal vectors are 3D
    normals.SetName("Load Orientation")
    load_magnitudes.SetName("Load Magnitude")

    j = 0
    for i, elementTag in enumerate(pContourElements):
        midpoint = [(pNodCoords[2*i][0] + pNodCoords[2*i+1][0])/2,
                    (pNodCoords[2*i][1] + pNodCoords[2*i+1][1])/2, 0]
        midpoints.InsertNextPoint(midpoint)

        normal = [(pNodCoords[2*i+1][1] - pNodCoords[2*i][1])/pElemLengths[i],
                    -(pNodCoords[2*i+1][0] - pNodCoords[2*i][0])/pElemLengths[i],
                    0]
        normals.InsertNextTuple(normal)


        if elementTag in loadElements:
            load_magnitudes.InsertNextValue(pressureDist[j])
            j += 1
        else:
            load_magnitudes.InsertNextValue(0.0)

    # Writing the vtp for the load distribution
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(midpoints)
    polydata.GetPointData().SetNormals(normals)
    polydata.GetPointData().AddArray(load_magnitudes)

    # Write vtkPolyData to a file
    writer = vtk.vtkXMLPolyDataWriter()
    file_path = os.path.join(outputPath, f"carga/carga{loadIndex+1}.vtp")
    writer.SetFileName(file_path)
    writer.SetInputData(polydata)
    writer.SetDataModeToAscii()  # Set the data mode to ASCII
    writer.Write()


def writeSteps(boneConfig, nSteps, nLoads):
    pasoPath = os.path.join(boneConfig.masterPath, "paso.inp")
    pasosPath = os.path.join(boneConfig.inputPath, "pasos.inp")

    with open(pasoPath, "r") as f:
        paso_content = f.read()

    with open(pasosPath, "w") as g:
        if (boneConfig.mode == "original"):
            for i in range(nSteps):
                formatted_content = paso_content.format(nStep=i,nLoads=nLoads,i=1)
                g.write(formatted_content)
        elif (boneConfig.mode == "new"):
            for i in range(2):
                if i == 0:
                    formatted_content = paso_content.format(nStep=i,nLoads=nLoads,i=i+1)
                else:
                    formatted_content = paso_content.format(nStep=i, nLoads=nSteps,i=i+1)
                g.write(formatted_content)
        else:
            raise ValueError("Invalid mode. Use 'original' or 'new'.")