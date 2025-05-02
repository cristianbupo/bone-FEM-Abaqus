import gmsh
import os
from math import ceil
import numpy as np

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

            if dim == 1:
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
    # boundaryConditionsPath = os.path.join(inputPath, "condicionesContorno.inp")

    with open(contornoPath, "w") as f:
        for dim, tag in physicalGroups:
            if dim < 3:
                name = gmsh.model.getPhysicalName(dim, tag)
                nodeTags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
                f.write("*NSET,NSET="+name+"\n")
                writeLines(nodeTags, f)

                lines.append(len(nodeTags))
    
    return lines

def writeOnFile(originFile, destinationFile, content):
    with open(originFile, 'r') as file, open(destinationFile, "w") as f:
        f_longBone_content = file.read()
        f.write(f_longBone_content.format(**content))


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

def findNWriteMeshAdjacencies1D(all1DElements, inputPath):
    allElemTags = all1DElements[1][0]  # Element IDs
    allElemNodeTags = all1DElements[2][0]  # Node connectivity for each element
    nElem = len(allElemTags)

    print(f"Number of 1D elements: {nElem}")

    # Initialize adjacency array (nElem x 3)
    # Columns: [Element ID, Neighbor 1, Neighbor 2]
    adjacency_array = np.zeros((nElem, 3), dtype=int)
    adjacency_array[:, 0] = np.arange(1, nElem + 1)  # First column: element IDs

    # Build node-to-element mapping
    node_to_elems = {}
    for eid, (n1, n2) in enumerate(zip(allElemNodeTags[::2], allElemNodeTags[1::2]), start=1):
        for node in (n1, n2):
            if node not in node_to_elems:
                node_to_elems[node] = []
            node_to_elems[node].append(eid)

    # Determine neighbors for each element
    for eid, (n1, n2) in enumerate(zip(allElemNodeTags[::2], allElemNodeTags[1::2]), start=1):
        neighbors = []
        for node in (n1, n2):
            neighbor_elems = [e for e in node_to_elems[node] if e != eid]
            neighbors.extend(neighbor_elems)

        # Remove duplicates and limit to 2 neighbors (1D elements can have at most 2 neighbors)
        neighbors = list(dict.fromkeys(neighbors))[:2]

        # Fill adjacency array
        adjacency_array[eid - 1, 1:1 + len(neighbors)] = neighbors

    # Write adjacency to file
    with open(os.path.join(inputPath, "adyacenciaElementos.inp"), "w") as f:
        f.write("elem, face1, face2, face3, face4\n")
        for row in adjacency_array:
            f.write(", ".join(map(str, row)) + ", 0, 0\n")