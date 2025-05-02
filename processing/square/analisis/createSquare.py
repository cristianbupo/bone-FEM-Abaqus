import gmsh
import os
import math
import numpy as np
from math import ceil
from collections import defaultdict


def findLastPosition(A, B):
    # A and B are arrays, B is a 2-element array
    n = len(A)
    for i in range(1, n):
        if (A[i] == B[0] and A[i-1] == B[1]) or (A[i] == B[1] and A[i-1] == B[0]):
            return i
    if (A[0] == B[0] and A[n-1] == B[1]) or (A[0] == B[1] and A[n-1] == B[0]):
        return n
    return -1  # return -1 if no match found


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


def findNWriteMeshAdjacencies(all2DElements, inputPath):
    allElemTags = all2DElements[1][0]
    allElemNodeTags = all2DElements[2][0]
    nElem = len(allElemTags)

    print(f"Number of 2D elements: {nElem}")

    # Reshape node_tags into quads
    quads = [allElemNodeTags[i:i+4] for i in range(0, len(allElemNodeTags), 4)]

    # Initialize adjacency array (nElem x 5)
    adjacency_array = np.zeros((nElem, 5), dtype=int)
    adjacency_array[:, 0] = np.arange(1, nElem + 1)  # First column: element IDs

    # Build edge: list of elements sharing that edge
    edge_to_elems = {}

    def ordered_edge(n1, n2):
        return tuple(sorted((n1, n2)))

    for eid, quad in enumerate(quads, start=1):  # Element IDs start at 1
        edges = [
            ordered_edge(quad[0], quad[1]),
            ordered_edge(quad[1], quad[2]),
            ordered_edge(quad[2], quad[3]),
            ordered_edge(quad[3], quad[0])
        ]
        for edge in edges:
            if edge not in edge_to_elems:
                edge_to_elems[edge] = []
            edge_to_elems[edge].append(eid)

    # Build element: neighbor per face
    for eid, quad in enumerate(quads, start=1):
        edges = [
            ordered_edge(quad[0], quad[1]),
            ordered_edge(quad[1], quad[2]),
            ordered_edge(quad[2], quad[3]),
            ordered_edge(quad[3], quad[0])
        ]
        neighbors = []
        for edge in edges:
            elems = edge_to_elems[edge]
            neighbor = [e for e in elems if e != eid]
            if neighbor:
                neighbors.append(neighbor[0])
            else:
                neighbors.append(0)  # Boundary face
        adjacency_array[eid - 1, 1:] = neighbors  # Fill connectivity columns

    # Write connectivity to file
    with open(os.path.join(inputPath, "adyacenciaElementos.inp"), "w") as f:
        f.write("elem, face1, face2, face3, face4\n")
        for row in adjacency_array:
            f.write(", ".join(map(str, row)) + "\n")


script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Initialize Gmsh
gmsh.initialize()

# Create a new model
gmsh.model.add("square")

gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
gmsh.option.setNumber("Mesh.Algorithm", 8)

# Parameters
length = 1.0
height = 1.0
meshSize = 0.05  # Characteristic length for meshing

# Create the square points
p1 = gmsh.model.geo.addPoint(0, 0, 0, meshSize)
p2 = gmsh.model.geo.addPoint(length, 0, 0, meshSize)
p3 = gmsh.model.geo.addPoint(length, height, 0, meshSize)
p4 = gmsh.model.geo.addPoint(0, height, 0, meshSize)

# Create the square lines
l1 = gmsh.model.geo.addLine(p1, p4)
l2 = gmsh.model.geo.addLine(p4, p3)
l3 = gmsh.model.geo.addLine(p3, p2)
l4 = gmsh.model.geo.addLine(p2, p1)

# Create the curve loop and surface
cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
surface = gmsh.model.geo.addPlaneSurface([cl])

elementsBorder = 10
nodesBorder = elementsBorder + 1

# Set transfinite meshing for the square
gmsh.model.geo.mesh.setTransfiniteCurve(l1, nodesBorder)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, nodesBorder)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, nodesBorder)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, nodesBorder)
gmsh.model.geo.mesh.setTransfiniteSurface(surface)

# Add a physical group for the surface
gmsh.model.addPhysicalGroup(2, [surface], 2, "square")

# Physical groups
# gmsh.model.addPhysicalGroup(2, [surface], -1,"wedge")

gmsh.model.geo.synchronize()
# Mesh the model
gmsh.model.mesh.generate(2)

gmsh.model.mesh.reverse([(2,surface)])

nodes=gmsh.model.mesh.getNodes()

nNodes = len(nodes[0])

# Extract node tags and coordinates
node_tags = nodes[0]
node_coords = nodes[1]

# Combine node tags and coordinates into a list of tuples
node_data = [(node_tags[i], node_coords[3 * i], node_coords[3 * i + 1]) for i in range(len(node_tags))]

# Sort nodes sequentially by x-coordinate, then by y-coordinate
node_data_sorted = sorted(node_data, key=lambda x: (x[2]+x[1]/(2*elementsBorder)))

oldTags = [node[0] for node in node_data_sorted]
newTags = range(1,nNodes+1)

# Apply the renumbering using the mapping
gmsh.model.mesh.renumberNodes(oldTags, newTags)
gmsh.model.mesh.renumberElements()


nodes=gmsh.model.mesh.getNodes()

with open("nodos.inp", "w") as f, open("restriccionesMultipunto.inp", "w") as g:
    f.write("*NODE,NSET=N2\n")
    g.write("*MPC, user, mode=dof\n")
    for i, node in enumerate(node_data_sorted):
        f.write("{}, {}, {}\n".format(i+1,node[1],node[2]))
        g.write("{}, {}, {}\n".format(i+1,i+1,i+1))

elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2)
nElem=0
for j in elemTags:
    nElem=nElem+len(j)


with open("conectividades.inp", "w") as f, open("decoy.inp", "w") as g, open("conectividades_ver.inp","w") as h:
    f.write("*User element, nodes=3, type=U1, properties=2, coordinates=2, variables=24\n1, 2\n")
    f.write("*ELEMENT,TYPE=U1,ELSET=UEL\n")
    if elemTypes[0]==2:
        elType="CPS3"
    elif elemTypes[0]==3:
        elType="CPS4"
        
    g.write("*ELEMENT, TYPE={}, ELSET=decoy\n".format(elType))
    h.write("*ELEMENT, TYPE={}\n".format(elType))
    for i in range(0,len(elemTypes)):
        for elem in range(0,len(elemTags[i])):
            if elemTypes[i]==2:
                f.write("{}, {}, {}, {}\n".format(elemTags[i][elem],elemNodeTags[i][3*elem],elemNodeTags[i][3*elem+1],elemNodeTags[i][3*elem+2]))
                h.write("{}, {}, {}, {}\n".format(elemTags[i][elem],elemNodeTags[i][3*elem],elemNodeTags[i][3*elem+1],elemNodeTags[i][3*elem+2]))
                g.write("{}, {}, {}, {}\n".format(int(elemTags[i][elem])+nElem,elemNodeTags[i][3*elem],elemNodeTags[i][3*elem+1],elemNodeTags[i][3*elem+2]))
            elif elemTypes[i]==3:
                f.write("{}, {}, {}, {}, {}\n".format(elemTags[i][elem],elemNodeTags[i][4*elem],elemNodeTags[i][4*elem+1],elemNodeTags[i][4*elem+2],elemNodeTags[i][4*elem+3]))
                h.write("{}, {}, {}, {}, {}\n".format(elemTags[i][elem],elemNodeTags[i][4*elem],elemNodeTags[i][4*elem+1],elemNodeTags[i][4*elem+2],elemNodeTags[i][4*elem+3]))
                g.write("{}, {}, {}, {}, {}\n".format(int(elemTags[i][elem])+nElem,elemNodeTags[i][4*elem],elemNodeTags[i][4*elem+1],elemNodeTags[i][4*elem+2],elemNodeTags[i][4*elem+3]))

os.remove("conectividades_ver.inp")
os.remove("decoy.inp")


PGroups = gmsh.model.getPhysicalGroups(2)


with open("cuerpo.inp", "w") as f:
    for p in PGroups:
        name=gmsh.model.getPhysicalName(p[0],p[1]) 
        entities = gmsh.model.getEntitiesForPhysicalGroup(p[0],p[1])
        nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(p[0],p[1])[0]
        f.write("*Nset, nset="+name+"\n")
        N=len(nodeTags)
        N2=6*ceil(N/6)
        for elem in range(0,N2):
            if elem<N:
                f.write("{}, ".format(nodeTags[elem]))
            else:
                f.write("0, ")

            if (elem+1) % 6  == 0:
                f.write("\n")

        f.write("*Elset, elset="+name+"\n")
        
        allElemTags=[]
        for e in entities:
            elemTypes, elemTags, _ = gmsh.model.mesh.getElements(2,e)
            for elem in range (0,len(elemTags[0])):
                allElemTags.append(elemTags[0][elem])

        N=len(allElemTags)
        N2=6*ceil(N/6)
        for elem in range(0,N2):
            if elem<N:
                f.write("{}, ".format(allElemTags[elem]))
            else:
                f.write("0, ")
                
            if (elem+1) % 6  == 0:
                f.write("\n")

#BOUNDARIES

gmsh.model.addPhysicalGroup(1, [l1],-1, "contorno4")
gmsh.model.addPhysicalGroup(1, [l2],-1, "contorno3")
gmsh.model.addPhysicalGroup(1, [l3],-1, "contorno2")
gmsh.model.addPhysicalGroup(1, [l4],-1, "contorno1")

PGroups = gmsh.model.getPhysicalGroups(1)
with open("contorno.inp", "w") as f:
    for p in PGroups:
        name=gmsh.model.getPhysicalName(p[0],p[1])
        nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(p[0],p[1])[0]
        f.write("*NSET,NSET="+name+"\n")
        #entities = gmsh.model.getEntitiesForPhysicalGroup(2,p[1])
        N=len(nodeTags)
        N2=6*ceil(N/6)
        for elem in range(0,N2):
            if elem<N:
                f.write("{}, ".format(nodeTags[elem]))
            else:
                f.write("0, ")

            if (elem+1) % 6  == 0:
                f.write("\n")

gmsh.model.geo.synchronize()

all2DElements = gmsh.model.mesh.getElements(2)

elemTypes = all2DElements[0][0]
allElemTags = all2DElements[1][0]
allElemNodeTags = all2DElements[2][0]
elemNodes = elemTypes + 1


contourElements, _, _, loadFaces, _, _ = findConectivityInfo(1, 2, allElemTags, allElemNodeTags, elemNodes)

findNWriteMeshAdjacencies(all2DElements, ".")

with open("carga.inp", "w") as f:
    f.write("*Step, name=LoadStep1\n")
    f.write("*DLOAD\n")
    for i, elem in enumerate(contourElements):
        f.write(f"{elem}, {loadFaces[i]}, 1.0\n")


with open("master/parametros.txt") as f:
    text = f.read()

nElementsLoad = len(loadFaces)

default_values = defaultdict(lambda: 1,{ #missing values are set to 1
    "numNode": nNodes,
    "nElems": nElem,
    "nLoads": 1,
    "listNElementLoads": f"(/{nElementsLoad}/)",
    "maxNElementLoads": nElementsLoad,
})


formatted_text = text.format_map(default_values)

with open(os.path.join(".",'nFilasCargas.txt'), "w") as f:
    f.write(str(nElementsLoad))

# If you want to save the formatted text back to a file:
with open("parametros.txt", "w") as f:
    f.write(formatted_text)

with open("gruposFisicos.txt", "w") as f:
    f.write("Element Tag, Physical Group Tag\n")
    for i in range(nElem):
        f.write(f"{i+1},{1}\n")

gmsh.fltk.run()
# gmsh.write(u"mesh3.inp")
# Finalize Gmsh
gmsh.finalize()