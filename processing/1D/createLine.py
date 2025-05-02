import gmsh
import os
import utils as ut
import numpy as np

gmsh.initialize()
gmsh.model.add("1D Line")

len1 = 2.0
len2 = 2.2
approxElementSize = 0.05

cartLenghts = [0.0,0.3,0.15,0.3,1.0]
cartLimits = np.cumsum(cartLenghts) + len1
print(cartLimits)

# Create 3 points
p1 = gmsh.model.geo.add_point(0, 0, 0)
p2 = gmsh.model.geo.add_point(len1, 0, 0)
p3 = gmsh.model.geo.add_point(len1 + len2, 0, 0)

# Create 2 lines
l1 = gmsh.model.geo.add_line(p1, p2)
l2 = gmsh.model.geo.add_line(p2, p3)

n1 = round(len1/approxElementSize)+1
n2 = round(len2/approxElementSize)+1

gmsh.model.geo.synchronize()

# Create 2 transfinite lines
gmsh.model.geo.mesh.setTransfiniteCurve(l1, n1)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, n2)

gmsh.model.geo.synchronize()

# Mesh 1D
gmsh.model.mesh.generate(1)

gmsh.model.addPhysicalGroup(1, [l1], -1, "Hueso")
gmsh.model.addPhysicalGroup(1, [l2], -1, "Cartilago")
gmsh.model.addPhysicalGroup(0, [p1], -1, "contorno1")
gmsh.model.addPhysicalGroup(0, [p2], -1, "contorno2")
gmsh.model.addPhysicalGroup(0, [p3], -1, "contorno3")

os.chdir(os.path.join(os.path.dirname(__file__),'analisis1D'))

# nodes
nodes = gmsh.model.mesh.getNodes()
elements = gmsh.model.mesh.getElements()

node_connectivity = {}
for i, elem in enumerate(elements[1][0]):  # Loop through element tags
    node1 = elements[2][0][2 * i]  # First node of the element
    node2 = elements[2][0][2 * i + 1]  # Second node of the element
    if node1 not in node_connectivity:
        node_connectivity[node1] = []
    if node2 not in node_connectivity:
        node_connectivity[node2] = []
    node_connectivity[node1].append(node2)
    node_connectivity[node2].append(node1)

# Renumber nodes based on connectivity
visited_nodes = set()
node_id_mapping = {}
new_node_id = 1

def traverse_node(node):
    global new_node_id
    if node not in visited_nodes:
        visited_nodes.add(node)
        node_id_mapping[node] = new_node_id
        new_node_id += 1
        for neighbor in node_connectivity[node]:
            traverse_node(neighbor)

# Start traversal from the first node
traverse_node(list(node_connectivity.keys())[0])

# Apply the new node IDs using gmsh.model.mesh.renumberNodes
old_node_tags = list(node_id_mapping.keys())
new_node_tags = list(node_id_mapping.values())
gmsh.model.mesh.renumberNodes(old_node_tags, new_node_tags)

# Update element connectivity with new node IDs
updated_elements = []
for i, elem in enumerate(elements[1][0]):  # Loop through element tags
    node1 = elements[2][0][2 * i]  # First node of the element
    node2 = elements[2][0][2 * i + 1]  # Second node of the element
    updated_elements.append((elem, node_id_mapping[node1], node_id_mapping[node2]))

# Renumber elements based on connectivity
element_id_mapping = {old_id: new_id + 1 for new_id, (old_id, _, _) in enumerate(updated_elements)}

# Apply the new element IDs using gmsh.model.mesh.renumberElements
old_element_tags = [elem[0] for elem in updated_elements]
new_element_tags = [element_id_mapping[elem[0]] for elem in updated_elements]
gmsh.model.mesh.renumberElements(old_element_tags, new_element_tags)

# gmsh.fltk.run()
# elemType = elem[0][0]
# elemTags = elem[1][0]
# elemCon = elem[2][0]
# nodeTags = nod[0]
# nodeCoords = nod[1]

# nodes
nodes = gmsh.model.mesh.getNodes()
elements = gmsh.model.mesh.getElements()

# ...existing code...

# Extract node tags and coordinates
node_data = [(nodes[0][i], nodes[1][3 * i], nodes[1][3 * i + 1], nodes[1][3 * i + 2]) for i in range(len(nodes[0]))]

# Sort nodes by their tags
sorted_nodes = sorted(node_data, key=lambda x: x[0])  # Sort by node tag

# Overwrite the nodes variable with sorted data
nodes = (
    [node[0] for node in sorted_nodes],  # Sorted node tags
    [coord for node in sorted_nodes for coord in node[1:]],  # Flattened sorted coordinates
    nodes[2]  # Keep the third element (parametric coordinates) unchanged
)

# Write sorted nodes to a file (optional)
# with open("sorted_nodes.inp", "w") as f:
#     f.write("*NODE,NSET=N2\n")
#     for node in sorted_nodes:
#         f.write(f"{node[0]}, {node[1]}, {node[2]}, {node[3]}\n")


gmsh.model.geo.synchronize()

with open("nodos.inp", "w") as f, open("restriccionesMultipunto.inp", "w") as g:
    f.write("*NODE,NSET=N2\n")
    g.write("*MPC, user, mode=dof\n")
    for i in range(len(nodes[0])):
        f.write(f"{nodes[0][i]}, {nodes[1][3*i]}\n")
        g.write(f"{nodes[0][i]}, {nodes[0][i]}, {nodes[0][i]}\n")


with open("conectividades.inp", "w") as f:
    f.write("*ELEMENT,TYPE=U1,ELSET=UEL\n")
    for i in range(len(elements[1][0])):
        f.write(f"{elements[1][0][i]}, {elements[2][0][2*i]}, {elements[2][0][2*i+1]}\n")

# ut.writeBody(gmsh.model.getPhysicalGroups(1), ".")
lines = ut.writeBoundaries(gmsh.model.getPhysicalGroups(0), ".")
# gmsh.fltk.run()

# Filter for 1D physical entities (dimension = 1)
physical_0D_groups = gmsh.model.getPhysicalGroups(0)

# Count the number of 1D physical entities
num_0D_physical_entities = len(physical_0D_groups)

print("Number of 0D physical entities:", num_0D_physical_entities)


content = {
    'stdWeight': 0.0,
    'nLoads': 1,
    'maxNElementLoads': 1,
    'numNode': len(nodes[0]),
    'nElems': len(elements[1][0]),
    'kOI': 0.0,
    'nContornos': num_0D_physical_entities,
    'filasContorno1': (lines[0] + 5) // 6,
    'filasContorno2': (lines[1] + 5) // 6,
    'filasContorno3': (lines[2] + 5) // 6,
    'velocidad': 0.0,
    'a2': 1,
    'a3': 1,
    'a4': 1,
    'b': 1
}

# with open('nFilasCargas.txt', "w") as f:
#     f.write("0")

print("current directory:", os.getcwd())
originFile = os.path.join("master", "parametros.txt")
destinationFile = "parametros.txt"

ut.writeOnFile(originFile, destinationFile, content)

# Classify nodes and elements based on cartLimits
node_classifications = []
element_classifications = []

# Classify nodes
for i, coord in enumerate(nodes[1][::3]):  # Extract x-coordinates of nodes
    for j, limit in enumerate(cartLimits):
        if coord < limit:
            node_classifications.append((nodes[0][i], j))  # Node ID and classification
            break
    else:
        node_classifications.append((nodes[0][i], len(cartLimits)))  # Last category

# Classify elements
for i, elem in enumerate(elements[1][0]):  # Loop through element tags
    node1 = elements[2][0][2 * i]  # First node of the element
    node2 = elements[2][0][2 * i + 1]  # Second node of the element
    x1 = nodes[1][int(3 * (node1 - 1))]  # x-coordinate of the first node
    x2 = nodes[1][int(3 * (node2 - 1))]  # x-coordinate of the second node
    center = (x1 + x2) / 2  # Center of the element

    for j, limit in enumerate(cartLimits):
        if center < limit:
            element_classifications.append((elem, j))  # Element ID and classification
            break
    else:
        element_classifications.append((elem, len(cartLimits)))  # Last category

# Export node classifications
with open("gruposFisicosN.txt", "w") as f:
    f.write("Node Tag, Physical Group Tag\n")
    for node_id, classification in node_classifications:
        f.write(f"{node_id}, {classification}\n")

# Export element classifications
with open("gruposFisicos.txt", "w") as f:
    f.write("Element Tag, Physical Group Tag\n")
    for elem_id, classification in element_classifications:
        f.write(f"{int(elem_id)}, {classification}\n")

# ut.findNWriteMeshAdjacencies1D(gmsh.model.mesh.getElements(1), ".")

gmsh.finalize()