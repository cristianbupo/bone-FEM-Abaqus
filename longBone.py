import FreeCADPath
import json
import FreeCAD as App
import FreeCADGui as Gui
import FreeCAD2gmsh as f2g
import geomdl2gmsh as g2g
import shutil
import os
import gmsh
from sketchUtils import setConstraintValue
import numpy as np
import argparse
from collections import defaultdict, deque

# Get the existing system PATH
env = os.environ.copy()

# Add the required Abaqus and Intel Fortran paths
additional_paths = [
    r"C:\\SIMULIA\\Commands",  # Abaqus Command Path
    r"C:\\Program Files (x86)\\IntelSWTools\\compilers_and_libraries_2017.2.187\\windows\\bin",  # Intel Fortran Compiler Path
    r"C:\\Program Files (x86)\\Microsoft Visual Studio 12.0\\VC\\bin",  # Visual Studio
]

# Append them to the PATH
env["PATH"] = ";".join(additional_paths) + ";" + env["PATH"]

class myConfigObject:
    def __init__(self, d=None):
        if d is not None:
            for key, value in d.items():
                if isinstance(value, dict):
                    setattr(self, key, myConfigObject(value))
                else:
                    setattr(self, key, value)

    def __iter__(self):
        for key, value in self.__dict__.items():
            yield key, value


class myObject:
    def __init__(self, d=None):
        if d is not None:
            for key, value in d.items():
                if isinstance(value, dict) and 'default' in value:
                    default_value = value['default']
                    setattr(self, key, default_value)
                elif isinstance(value, dict):
                    setattr(self, key, myObject(value))
                else:
                    setattr(self, key, value)

    def __iter__(self):
        for key, value in self.__dict__.items():
            yield key, value


def getBoneData(configFile):
    with open(configFile, 'r') as file:
        data = json.load(file)

    bone = myObject(data['bone'])
    boneLimits = myConfigObject(data['bone'])
    boneConfig = myObject(data['boneConfig'])

    base_name = os.path.basename(configFile)
    no_extension = os.path.splitext(base_name)[0]

    boneConfig.inputPath = os.path.join(boneConfig.inputPath, no_extension)
    boneConfig.outputPath = os.path.join(boneConfig.outputPath, no_extension)
    
    return bone, boneLimits, boneConfig


def copyFile(origin, destination):
    # Create the destination directory if it does not exist
    destination_dir = os.path.dirname(destination)
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir, exist_ok=True)

    # Copy the file if the origin exists
    if os.path.exists(origin):
        shutil.copy(origin, destination)


def clear_folder(folder):

    # Clear the folder by deleting its contents
    # If the folder does not exist, create it

    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)
    else:
        for item in os.listdir(folder):
            item_path = os.path.join(folder, item)
            try:
                if os.path.isfile(item_path) or os.path.islink(item_path):
                    os.unlink(item_path)  # Removes files and symlinks
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)  # Removes directories and their contents
            except Exception as e:
                print(f'Failed to delete {item_path}. Reason: {e}')


def modifySketch(defaultGeometry=False):

    # Modify sketch
    # Create and write mesh

    Gui.setupWithoutGUI()  # Initialize FreeCAD without GUI
#    if boneConfig.capsule:
    fileName = os.path.join('CADs', 'longBone.FCStd')
#    else:
#        fileName = os.path.join('CADs', 'longBoneCapsule.FCStd')

    doc = App.openDocument(fileName)
    sketch = doc.getObject('Sketch')

    print(defaultGeometry)
    if defaultGeometry:
        for key, value in bone.geom_vars:
            current_value = sketch.getDatum(key).Value
            setattr(bone.geom_vars, key, current_value)

            print(f"Value for {key} was set to {current_value}")
    else:
        for key, value in bone.geom_vars:
            min_value = getattr(boneLimits.geom_vars, key).min_value
            max_value = getattr(boneLimits.geom_vars, key).max_value
            setConstraintValue(sketch, key, value, min_value, max_value)
            current_value = sketch.getDatum(key)

            if current_value != value:
                print(f"Value for {key} was not set correctly. Expected {value}, "
                    f"but got {current_value}.")
            else:
                print(f"Value for {key} was set correctly: {value}")

    # Save CAD file
    doc.saveAs(boneConfig.inputPath + '/longBone.FCStd')

    curvesMesh, curvesArea, curvesAdvance = g2g.processSketchNurbs(sketch, boneConfig)

    gmsh.initialize()
    elem1, nod1 = g2g.container2gmsh(bone, boneConfig, curvesMesh, curvesArea)
    gmsh.write(boneConfig.inputPath + '/longBone.msh')

    elem2, nod2 = g2g.container2advanceMesh(bone, boneConfig, curvesAdvance)
    gmsh.write(boneConfig.inputPath + '/advanceMesh.msh')

    intersectMeshes(bone, boneConfig, elem1, nod1, elem2, nod2)

    gmsh.finalize()

    return curvesArea


def intersectMeshes(bone, boneConfig, elem1, nod1, elem2, nod2):

    # n_e = bone.mesh_vars.number_elements
    # a2, a3, _, _, _ = g2g.meshLineElements(n_e)
    
    # elemType = elem[0][0]
    # elemTags = elem[1][0]
    # elemCon = elem[2][0]
    # nodeTags = nod[0]
    # nodeCoords = nod[1]
    
    # Create separate views for node and element data
    node_view_tag = gmsh.view.add("Node_Data")
    gmsh.view.addModelData(node_view_tag, 0, "Secondary model", "NodeData", nod2[0], [[tag] for tag in nod2[0]])

    element_view_tag = gmsh.view.add("Element_Data")
    gmsh.view.addModelData(element_view_tag, 0, "Secondary model", "ElementData", elem2[1][0], [[tag] for tag in elem2[1][0]])

    print("Started to calculate centroids")
    centroids = elementsCentroids(elem1, nod1)
    print("Finished calculating centroids")

    centroidsx = [[centroid[0]] for centroid in centroids]
    centroidsy = [[centroid[1]] for centroid in centroids]

    centroids_x_view_tag = gmsh.view.add("Centroids X")
    gmsh.view.addModelData(centroids_x_view_tag, 0, "Main model", "ElementData", elem1[1][0], centroidsx)

    centroids_y_view_tag = gmsh.view.add("Centroids Y")
    gmsh.view.addModelData(centroids_y_view_tag, 0, "Main model", "ElementData", elem1[1][0], centroidsy)

    output_file = boneConfig.inputPath + "/Combined_Data.pos"
    gmsh.view.write(node_view_tag, output_file)
    gmsh.view.write(element_view_tag, output_file, append=True)
    gmsh.view.write(centroids_x_view_tag, output_file, append=True)
    gmsh.view.write(centroids_y_view_tag, output_file, append=True)

    
    node_matrix, elem_matrix = organizeMeshByConnectivity(elem2)
    
    # Save node matrix
    np.savetxt(boneConfig.inputPath + "/nodes.csv", node_matrix, fmt='%d', delimiter=",", comments='')

    # Save element matrix
    np.savetxt(boneConfig.inputPath + "/elements.csv", elem_matrix, fmt='%d', delimiter=",", comments='')

    midIndex = int(node_matrix.shape[1]-1)//2

    # Verify the middle line is correct-----
    middleLine, middleCoords = getMiddleLine(elem2, nod2)
    middleLine.reverse()
    middleCoords.reverse()

    if  middleLine != node_matrix[:,midIndex].tolist():
        print("Middle line is not correct")
    else:
        print("Middle line is correct")

    referenceLine = node_matrix[:,midIndex].tolist()
    referenceCoords = [nod2[1][3*(node-1)+1] for node in referenceLine]

    maxLength = referenceCoords[0] - referenceCoords[-1]

    plate_1 = bone.geom_vars.plate_1
    plate_2 = bone.geom_vars.plate_2
    cart_thick = bone.geom_vars.cart_thick
    max_length = cart_thick - plate_2

    zonesLimits = [0.0, 0.25, 0.35, 0.6, max_length]+referenceCoords[-1]
    
    closest_indexes = []
    for limit in zonesLimits:
        # Find the index of the closest value in referenceCoords
        closest_index = np.argmin(np.abs(np.array(referenceCoords) - limit))
        closest_indexes.append(closest_index)

    print("Closest indexes:", closest_indexes)
    
    extracted_rows = node_matrix[closest_indexes, :]
    
    np.savetxt(boneConfig.inputPath + "/extracted_rows.csv", extracted_rows, fmt='%d', delimiter=",", comments='')

    # Create a new view for extracted rows
    extracted_rows_view_tag = gmsh.view.add("Extracted_Rows")

    # Prepare the data for the extracted rows
    extracted_rows_data = [[1 if node in extracted_rows else 0] for node in node_matrix.flatten()]

    # Add the extracted rows data to the view
    gmsh.view.addModelData(
        extracted_rows_view_tag,
        0,
        "Secondary model",
        "NodeData",
        node_matrix.flatten(),
        extracted_rows_data
    )

    # Write the extracted rows data to the output file
    gmsh.view.write(extracted_rows_view_tag, output_file, append=True)

    # centroids contains the centroids of the elements in elem1
    # centroids = elementsCentroids(elem1, nod1)

    # centroidsx = [[centroid[0]] for centroid in centroids]
    # centroidsy = [[centroid[1]] for centroid in centroids]

    print("Started to check if the elements and nodes are below the line")

    # Initialize lists for marking elements and nodes
    is_below_list_elements = [[] for _ in range(len(elem1[1][0]))]  # For elements
    is_below_list_nodes = [[] for _ in range(len(nod1[0]))]  # For nodes

    # Iterate over each row in extracted_rows
    for row_index in range(len(extracted_rows)):
        coordinatesRow = get_coordinates_from_extracted_row(row_index, extracted_rows, nod2)

        # Mark elements below the current row
        for i in range(len(elem1[1][0])):
            if is_below_list_elements[i] == []:  # Mark only unmarked elements
                centroid = centroids[i]
                x = centroid[0]
                y = centroid[1]

                # Check if the element's centroid is below the line
                is_below = is_point_below_line(x, y, coordinatesRow)

                if is_below:
                    is_below_list_elements[i] = row_index  # Mark with the current row index

                # Assign the element's nodes to the nodes list
                element_nodes = elem1[2][0][4 * i:4 * (i + 1)]  # Get the 4 nodes of the element
                for node in element_nodes:
                    # if is_below_list_nodes[int(node) - 1] == []:  # Mark only unmarked nodes
                    is_below_list_nodes[int(node) - 1] = row_index 
    
    # Mark remaining unmarked elements and nodes with the next number
    next_number = len(extracted_rows)
    for i in range(len(is_below_list_elements)):
        if is_below_list_elements[i] == []:  # If the element is still unmarked
            is_below_list_elements[i] = next_number
            element_nodes = elem1[2][0][4 * i:4 * (i + 1)]
            for node in element_nodes:
                is_below_list_nodes[int(node) - 1] = next_number

    # Create a Gmsh view for the marked elements
    is_below_view_tag_elements = gmsh.view.add("Is_Below_Line_Elements")
    gmsh.view.addModelData(
        is_below_view_tag_elements,
        0,
        "Main model",
        "ElementData",
        elem1[1][0],
        [[value] for value in is_below_list_elements]
    )

    # Create a Gmsh view for the marked nodes
    is_below_view_tag_nodes = gmsh.view.add("Is_Below_Line_Nodes")
    gmsh.view.addModelData(
        is_below_view_tag_nodes,
        0,
        "Main model",
        "NodeData",
        nod1[0],
        [[value] for value in is_below_list_nodes]
    )

    # Write the marked elements and nodes data to the output file
    gmsh.view.write(is_below_view_tag_elements, output_file, append=True)
    gmsh.view.write(is_below_view_tag_nodes, output_file, append=True)

    # Save the node markings to a file
    tipoCartilagoNodesPath = os.path.join(boneConfig.inputPath, 'gruposFisicosN.txt')
    with open(tipoCartilagoNodesPath, "w") as g:
        g.write('Node Tag, Physical Group Tag\n')
        for i in range(len(nod1[0])):
            g.write(f'{nod1[0][i]}, {is_below_list_nodes[i]}\n')

    # Save the element markings to a file
    tipoCartilagoElementsPath = os.path.join(boneConfig.inputPath, 'gruposFisicos.txt')
    with open(tipoCartilagoElementsPath, "w") as g:
        g.write('Element Tag, Physical Group Tag\n')
        for i in range(len(elem1[1][0])):
            g.write(f'{elem1[1][0][i]}, {is_below_list_elements[i]}\n')

def get_coordinates_from_extracted_row(row_index, extracted_rows, nod):
    """
    Extract the coordinates of nodes from a specific row of the extracted_rows matrix.

    Args:
        row_index (int): The index of the row in the extracted_rows matrix.
        extracted_rows (np.ndarray): A 2D array where each row contains node IDs.
        nod (list): A list where nod[1] contains the flattened node coordinates.

    Returns:
        list: A list of tuples containing the coordinates (x, y, z) of the nodes in the specified row.
    """
    if row_index < 0 or row_index >= extracted_rows.shape[0]:
        raise IndexError("Row index is out of bounds for the extracted_rows matrix.")

    # Get the node IDs from the specified row
    node_ids = extracted_rows[row_index]

    # Retrieve the coordinates for each node ID
    coordinates = []
    for node_id in node_ids:
        # Node IDs are 1-based, so adjust for 0-based indexing in the nod array
        x = nod[1][3 * (node_id - 1)]
        y = nod[1][3 * (node_id - 1) + 1]
        z = nod[1][3 * (node_id - 1) + 2]
        coordinates.append((x, y, z))

    return coordinates


def is_point_below_line(x, y, coordinatesRow):
    """
    Detect if a point (x, y) lies below the line created by a row of coordinatesRow.
    If x is outside the range of the line, assume an infinite horizontal line from the nearest point.

    Args:
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        coordinatesRow (np.ndarray): 2D array where each row represents a point (x, y).

    Returns:
        bool: True if the point lies below the line, False otherwise.
    """
    # Ensure coordinatesRow has at least two points to form a line
    if len(coordinatesRow) < 2:
        raise ValueError("coordinatesRow must contain at least two points to form a line.")

    # Sort coordinatesRow by x-coordinate to ensure proper line formation
    coordinatesRow = sorted(coordinatesRow, key=lambda point: point[0])
    
    # Iterate through consecutive points in coordinatesRow to find the segment containing x
    for i in range(len(coordinatesRow) - 1):
        x1, y1, _ = coordinatesRow[i]
        x2, y2, _ = coordinatesRow[i + 1]

        # Check if x lies between x1 and x2
        if x1 <= x <= x2 or x2 <= x <= x1:
            # Calculate the y-value on the line at x using linear interpolation
            y_on_line = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
            return y < y_on_line  # True if the point is below the line

    # If x is outside the range, assume an infinite horizontal line from the nearest point
    if x < coordinatesRow[0][0]:  # Left of the first point
        return y < coordinatesRow[0][1]
    elif x > coordinatesRow[-1][0]:  # Right of the last point
        return y < coordinatesRow[-1][1]

    return False


def findClosestIndexes(referenceCoords, zonesLimits):
    closest_indexes = []
    for limit in zonesLimits:
        # Find the index of the closest value in referenceCoords
        closest_index = np.argmin(np.abs(np.array(referenceCoords) - limit))
        closest_indexes.append(closest_index)
    return closest_indexes


def getMiddleLine(elem, nod):
    elemType = elem[0][0]
    # elemTags = elem[1][0]
    # elemCon = elem[2][0]
    nodeTags = nod[0]
    nodeCoords = nod[1]

    middleLine = []
    middleLineCoords = []
    for i in range(len(nodeTags)):
        if abs(nodeCoords[3*i]) < 1e-6:
            middleLine.append(nodeTags[i])
            middleLineCoords.append(nodeCoords[3*i+1])
    pass
    pass

    # Sort nodeTags and nodeCoords according to nodeCoords
    sorted_pairs = sorted(zip(middleLine, middleLineCoords), key=lambda x: x[1])  # Sort by nodeCoords
    middleLine2, middleLineCoords2 = zip(*sorted_pairs)  # Unzip the sorted pairs back into separate lists

    # Convert back to lists if needed
    middleLine = list(middleLine2)
    middleLineCoords = list(middleLineCoords2)
    
    return middleLine, middleLineCoords

# ----------------------------------------------------------------------
# elem  : gmsh "Elements" block, shaped (something, 1, 4) :
#         elem[1][0]  – element tags
#         elem[2][0]  – connectivity flattened (4 nodes per quad)
# ----------------------------------------------------------------------
def organizeMeshByConnectivity(elem):
    # ------------------------------------------------------------------ #
    # 1.  helpers                                                        #
    # ------------------------------------------------------------------ #
    elemTags = elem[1][0]                # (n_elem,)
    elemCon  = elem[2][0]                # (4*n_elem,)

    # element → [n1,n2,n3,n4]
    elem2nodes = {tag: elemCon[i*4:(i+1)*4] for i, tag in enumerate(elemTags)}

    # node → [e1,e2,…]   (needed to find neighbours)
    node2elems = defaultdict(list)
    for e, nds in elem2nodes.items():
        for n in nds:
            node2elems[n].append(e)

    # neighbour list : two quads share **exactly one edge** (two nodes)
    neigh = defaultdict(list)
    for e, nds in elem2nodes.items():
        for n in nds:
            for m in node2elems[n]:
                if m != e and len(set(nds) & set(elem2nodes[m])) == 2:
                    neigh[e].append(m)

    # ------------------------------------------------------------------ #
    # 2.  breadth-first traversal – build sparse (i,j) ↦ elemTag map     #
    # ------------------------------------------------------------------ #
    elem_matrix_dict = {}
    q        = deque([(min(elemTags), 0, 0)])   # start in lower-left quad
    visited  = set()

    while q:
        e, i, j = q.popleft()
        if e in visited:
            continue
        visited.add(e)
        elem_matrix_dict[(i, j)] = e

        n1, n2, n3, n4 = elem2nodes[e]

        for m in neigh[e]:
            if m in visited:
                continue
            shared = frozenset(elem2nodes[e]) & frozenset(elem2nodes[m])

            if shared == frozenset({n1, n2}):      # bottom edge
                q.append((m, i+1, j))
            elif shared == frozenset({n2, n3}):    # right edge
                q.append((m, i,   j+1))
            elif shared == frozenset({n3, n4}):    # top edge
                q.append((m, i-1, j))
            elif shared == frozenset({n4, n1}):    # left edge
                q.append((m, i,   j-1))

    # dense element matrix ---------------------------------------------
    rows = [ij[0] for ij in elem_matrix_dict]
    cols = [ij[1] for ij in elem_matrix_dict]
    rmin, rmax = min(rows), max(rows)
    cmin, cmax = min(cols), max(cols)

    elem_matrix = np.full((rmax-rmin+1, cmax-cmin+1), -1, dtype=int)
    for (i, j), tag in elem_matrix_dict.items():
        elem_matrix[i - rmin, j - cmin] = tag

    # ------------------------------------------------------------------ #
    # 3.  node matrix – first pass (may miss the very top boundary)      #
    # ------------------------------------------------------------------ #
    node_matrix = np.full((elem_matrix.shape[0]+1,
                           elem_matrix.shape[1]+1), -1, dtype=int)

    for i in range(elem_matrix.shape[0]):
        for j in range(elem_matrix.shape[1]):
            e = elem_matrix[i, j]
            if e == -1:
                continue
            n1, n2, n3, n4 = elem2nodes[e]
            node_matrix[i,   j]   = n1
            node_matrix[i,   j+1] = n2
            node_matrix[i+1, j+1] = n3
            node_matrix[i+1, j]   = n4

    # ------------------------------------------------------------------ #
    # 4.  insert missing TOP boundary row                                #
    # ------------------------------------------------------------------ #
    boundary_nodes = {n for n, els in node2elems.items() if len(els) == 1}
    have_top_row  = node_matrix[0, 0] in boundary_nodes

    if not have_top_row:
        top_nodes = []                            # will hold  (n_cols + 1) nodes
        for j, e in enumerate(elem_matrix[0]):    # first element-row
            n1, n2, n3, n4 = elem2nodes[e]
            if j == 0:
                top_nodes.append(n4)              # first corner
            top_nodes.append(n3)                  # right corner of this quad

        node_matrix = np.vstack([np.array(top_nodes, dtype=int),
                                 node_matrix])    # add at the top

    # ------------------------------------------------------------------ #
    # 5.  drop duplicate rows (fixes “swapped last two rows”)            #
    # ------------------------------------------------------------------ #
    unique = []
    seen   = set()
    for r in node_matrix:
        t = tuple(r)
        if t not in seen:
            unique.append(r)
            seen.add(t)
    node_matrix = np.vstack(unique)

    return node_matrix, elem_matrix


def elementsCentroids(elem, nod):
    """
    Calculate the centroids of the elements.
    For quadrilateral elements, the centroid is calculated using the intersection
    of the diagonals formed by the centroids of the triangles derived from the quad.
    """
    elemType = elem[0][0]
    elemTags = elem[1][0]
    elemCon = elem[2][0]
    nodeCoords = nod[1]

    centroids = []

    for i in range(len(elemTags)):
        start = (elemType + 1) * i
        end = start + elemType + 1
        nodes = elemCon[start:end]
        coords = [nodeCoords[int(3 * (node - 1)):int(3 * (node - 1) + 3)] for node in nodes]

        # Step 1: Divide the quad into 4 triangles
        triangles = [
            [coords[0], coords[1], coords[2]],  # Triangle 1
            [coords[0], coords[3], coords[2]],  # Triangle 2
            [coords[1], coords[2], coords[3]],  # Triangle 3
            [coords[0], coords[1], coords[3]],  # Triangle 4
        ]

        # Step 2: Calculate centroids of the triangles
        triangle_centroids = []
        for tri in triangles:
            cx = sum([point[0] for point in tri]) / 3
            cy = sum([point[1] for point in tri]) / 3
            triangle_centroids.append((cx, cy))

        # Step 3: Calculate the intersection of the diagonals
        centroid = calculateIntersection(
            triangle_centroids[0], triangle_centroids[1],
            triangle_centroids[2], triangle_centroids[3]
        )

        centroids.append(centroid)

    return centroids


def calculateIntersection(p1, p2, p3, p4):
    """
    Calculate the intersection point of two lines formed by points (p1, p2) and (p3, p4).
    Returns the intersection point as (x, y) or None if the lines are parallel or coincident.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    # Line intersection formula
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    ix = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom
    iy = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom

    return (ix, iy)


def copyAnalysisFiles():

    # Copy analysis files
    
    inputPath = boneConfig.inputPath
    # outputPath = boneConfig.outputPath

    # Define the source files and their destination filenames
    jobName = f"analisis{bone.simulation_vars.case_string}"
    files_to_copy = {
        #"master/run.bat": jobName + ".bat",
        "master/analisis.inp": jobName + ".inp",
        "master/debug.bat": jobName + "Debug.bat",
        "master/run.ps1":  "run.ps1",
        "master/propiedades.csv": "propiedades.csv",
        # "master/abaqus_v6.env": "abaqus_v6.env"
        "master/condicionesContorno.inp": "condicionesContorno.inp"
    }

    # Copy each file to the inputPath with the new name
    for src, dest in files_to_copy.items():
        src_path = os.path.join(os.getcwd(), src)
        dest_path = os.path.join(inputPath, dest)
        dest_dir = os.path.dirname(dest_path)

        # Ensure the destination directory exists
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir, exist_ok=True)
                
        shutil.copy(src_path, dest_path)

    formatFile(
        os.path.join('master','condicionesContorno.inp'),
        os.path.join(boneConfig.inputPath, 'condicionesContorno.inp'),
        '**'
    )

    formatFile(
        os.path.join('master','run.bat'),
        os.path.join(boneConfig.inputPath, jobName + ".bat"),
        jobName
    )

    with open (os.path.join(boneConfig.inputPath,'runCommands.txt'), 'a') as file:
        file.write(f"{inputPath},{jobName}\n")


def formatFile(inputFile, outputFile, fillText):
    with open(inputFile, 'r') as f:
        text = f.read()

    text = text.format(fillText)

    with open(outputFile, 'w') as g:
        g.write(text)



def setupSimulations():

    # Set up simulations parameters

    N = bone.simulation_vars.number_analyzes

    OI_threshold_vec = np.linspace(
        boneLimits.oss_vars.OI_threshold.min_value,
        boneLimits.oss_vars.OI_threshold.max_value,
        N
    )

    kOI_threshold_vec = np.linspace(
        boneLimits.oss_vars.kOI.min_value,
        boneLimits.oss_vars.kOI.max_value,
        N
    )

    if N > 1:
        for i in range(1, N+1):
            bone.simulation_vars.case_id = i
            bone.simulation_vars.case_string = str(i).zfill(3)
            bone.oss_vars.kOI = kOI_threshold_vec[i-1]
            copyAnalysisFiles()

    else: 
        copyAnalysisFiles()


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run bone FEM Abaqus simulations.")
    parser.add_argument(
        '--configs',
        nargs='+',
        required=True,
        help="List of configuration file paths (e.g., 'diffusionConvexSteps.json diffusionConcaveSteps.json')"
    )
    parser.add_argument(
        '--defaultGeometry',
        action='store_true',
        help='If set, it uses the saved geometry'
    )

    # parser.add_argument(
    #     '--capsule',
    #     action='store_true',
    #     help="If set, the model will be enclosed in a capsule"
    # )

    args = parser.parse_args()

    # Iterate over the provided configuration files
    for config_file in args.configs:
        print(f"Processing configuration: {config_file}")
        bone, boneLimits, boneConfig = getBoneData(os.path.join('loadCases',config_file))
        clear_folder(boneConfig.inputPath)

        # Create or clear the runCommands.txt file
        with open(os.path.join(boneConfig.inputPath, 'runCommands.txt'), 'w') as file:
            pass

        # Run the simulation steps
        modifySketch(defaultGeometry=args.defaultGeometry)  # Modify sketch, create and write mesh
        setupSimulations()  # Set up simulation parameters

    