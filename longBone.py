import FreeCADPath
import json
import FreeCAD as App
import FreeCADGui as Gui
import macros.FreeCAD.sliderWindow as sw
import geomdl2gmsh as g2g
import shutil
import os
import gmsh
from datetime import datetime
from sketchUtils import setConstraintValue
import subprocess

# Get the existing system PATH
env = os.environ.copy()

# Add the required Abaqus and Intel Fortran paths
additional_paths = [
    r"C:\SIMULIA\Commands",  # Abaqus Command Path
    r"C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.2.187\windows\bin",  # Intel Fortran Compiler Path
    r"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin",  # Visual Studio
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


def getBoneData():
    with open('combined_vars.json', 'r') as file:
        data = json.load(file)

    with open('configuration.json', 'r') as file:
        config_data = json.load(file)

    bone = myObject(data)
    boneLimits = myConfigObject(data)
    boneConfig = myObject(config_data)

    return bone, boneLimits, boneConfig


def modifySketch(bone, boneLimits, boneConfig):
    Gui.setupWithoutGUI()  # Initialize FreeCAD without GUI
    fileName = 'CADs/longBone.FCStd'
    doc = App.openDocument(fileName)
    sketch = doc.getObject('Sketch')

    gmsh.initialize()  # Initialize gmsh

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

    curvesMesh, curvesArea, _ = g2g.processSketchNurbs(sketch, boneConfig)

    g2g.container2gmsh(bone, boneConfig, curvesMesh, curvesArea)

    gmsh.finalize()  # Finalize gmsh

    return curvesArea


def myFunction(bone, boneLimits, boneConfig, sketch):

    gmsh.initialize()  # Initialize gmsh

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

    if boneConfig.deleteBackUp:
        backup_dir = os.path.join("currentBackUp")
        shutil.rmtree(backup_dir, ignore_errors=True)
        os.makedirs(backup_dir, exist_ok=True)

        currentLoadHistory_dir = os.path.join("currentLoadHistory")
        shutil.rmtree(currentLoadHistory_dir, ignore_errors=True)
        os.makedirs(currentLoadHistory_dir, exist_ok=True)

    curvesMesh, _, _ = g2g.processSketchNurbs(
        sketch,
        renderRaw=boneConfig.renderRaw,
        renderMesh=boneConfig.renderMesh,
        renderArea=boneConfig.renderArea,
        renderLength=boneConfig.renderLength
    )

    if boneConfig.pickleIn:
        # backup fiñes
        backup_dir = "currentBackUp"
        destination_dir = "current"

        os.makedirs(backup_dir, exist_ok=True)
        os.makedirs(destination_dir, exist_ok=True)

        for file_name in os.listdir(backup_dir):
            full_file_name = os.path.join(backup_dir, file_name)
            if os.path.isfile(full_file_name):
                copyFile(full_file_name, destination_dir)

        print(f"All files have been copied from {backup_dir} to {destination_dir}")

    g2g.container2gmsh(
        bone.load_vars.load_center,
        bone.load_vars.load_amplitude,
        bone.load_vars.load_radius,
        bone.load_vars.kOI,
        curvesMesh,
        inputPath='current',
        numberElements=bone.mesh_vars.number_elements,
        runFltk=boneConfig.runFltk,
        skipWrite=boneConfig.skipWrite,
        writeVTK=boneConfig.writeVTK
    )

    if boneConfig.pickleOut:
        # backup files
        source_dir = "current"
        backup_dir = os.path.join("currentBackUp")

        os.makedirs(source_dir, exist_ok=True)
        os.makedirs(backup_dir, exist_ok=True)

        for file_name in os.listdir(source_dir):
            if os.path.isfile(os.path.join(source_dir, file_name)):
                copyFile(os.path.join(source_dir, file_name), backup_dir)

        print(f"Files have been copied to {backup_dir}")

    gmsh.finalize()  # Finalize gmsh

    base_path = "D:\\bone-FEM-results"
    id_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    dir_name = os.path.join(base_path, f"test_{id_time}")
    
    saveResults(dir_name, 0)


def oldRunAbaqus(bone, boneConfig):
    # Construct the command string
    boneG = bone.geom_vars
    boneL = bone.load_vars
    command = (
        f'current\\analisis.bat "-ne {bone.mesh_vars.number_elements}'
        f'-rx {boneG.radius_x} -ry {boneG.radius_y} '
        f'-ha {boneG.head_angle} -bw {boneG.bone_width} -hh {boneG.head_angle} '
        f'-bh {boneG.bone_height} -ca {boneG.curve_angle} '
        f'{"--showLoad" if boneConfig.showLoad else ""} '
        # f'{"--renderRaw" if boneConfig.renderRaw else ""} '
        # f'{"--renderMesh" if boneConfig.renderMesh else ""} '
        # f'{"--renderArea" if boneConfig.renderArea else ""} '
        # f'{"--renderLength" if boneConfig.renderLength else ""} '
        f'{"--runFltk" if boneConfig.runFltk else ""} '
        f'{"--skipWrite" if boneConfig.skipWrite else ""} '
        f'{"--writeVTK" if boneConfig.writeVTK else ""} '
        f'{"--runMeshing" if boneConfig.runMeshing else ""} '
        f'{"--pickleIn" if boneConfig.pickleIn else ""} '
        f'{"--pickleOut" if boneConfig.pickleOut else ""} '
        f'{"--deleteBackUp" if boneConfig.deleteOutput else ""} '
        f'-hl {boneL.load_center} -kl {boneL.load_amplitude} -rl {boneL.load_radius} '
        f'-koi {boneL.kOI}" {str(boneConfig.runAbq).lower()}'
    )

    # Remove extra spaces
    command = ' '.join(command.split())

    # Run the command
    # skipcq: PYL-W1510, BAN-B602
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Print the output and error (if any)
    print("Output:\n", result.stdout)
    print("Error:\n", result.stderr)

    return command

def saveResults(dir_name, index):
    copyFile("current\\analisis.vtu", f"{dir_name}\\analisis{index}.vtu")
    copyFile("current\\malla.vtu", f"{dir_name}\\malla{index}.vtu")
    copyFile("current\\carga.vtp", f"{dir_name}\\carga{index}.vtp")


def copyFile(origin, destination):
    # Create the destination directory if it does not exist
    destination_dir = os.path.dirname(destination)
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir, exist_ok=True)

    # Copy the file if the origin exists
    if os.path.exists(origin):
        shutil.copy(origin, destination)


def clear_folder(folder):
    if not os.path.exists(folder):
        print(f'Folder {folder} ist.')
        return

    for item in os.listdir(folder):
        item_path = os.path.join(folder, item)
        try:
            if os.path.isfile(item_path) or os.path.islink(item_path):
                os.unlink(item_path)  # Removes files and symlinks
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)  # Removes directories and their contents
        except Exception as e:
            print(f'Failed to delete {item_path}. Reason: {e}')

def runAnalysis(boneConfig, boneLimits, bone):
    
    # Clear the input folder
    clear_folder(boneConfig.inputPath)

    inputPath = boneConfig.inputPath
    # outputPath = boneConfig.outputPath

    # clear_folder(inputPath)

    # Define the source files and their destination filenames
    files_to_copy = {
#        "master/abaqus_v6.env": "abaqus_v6.env",
        "master/analisis.inp": "analisis.inp",
        "master/debug.bat": "debug.bat",
        "master/general2DElastic.for": "user.for",
        "master/run.bat": "run.bat",
        "master/propiedades.txt": "propiedades.txt"
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

    modifySketch(bone, boneLimits, boneConfig)
    process = subprocess.run("run.bat", shell=True, capture_output=True,
                            text=True, cwd=boneConfig.inputPath, env=env)

    print(process.stdout)
    print(process.stderr)


def singleAnalysis():
    bone, boneLimits, boneConfig = getBoneData()

    setattr(boneConfig, 'deleteOutput', True)
    setattr(boneConfig, 'runAbq', True)
    # setattr(boneConfig, 'runFltk', True)
    setattr(bone.mesh_vars, 'number_elements', 30)
    setattr(bone.geom_vars, 'bone_width', 2.2)
    setattr(bone.geom_vars, 'head_height', 3.0)

    # Output concave single load case
    folder = r'C:\Users\crist\OneDrive\Documents\results\singleAnalysis'
    setattr(boneConfig, 'inputPath', folder)
    setattr(boneConfig, 'outputPath', folder)

    setattr(bone.geom_vars, 'radius_x', 1.4)
    setattr(bone.geom_vars, 'radius_y', 1.5)
    setattr(bone.geom_vars, 'head_angle', 25.0)
    setattr(bone.geom_vars, 'cart_thick', 2.5)
    setattr(bone.geom_vars, 'curve_angle', 0.0)
    setattr(bone.geom_vars, 'head_height', 1)

    setattr(bone.load_vars, 'number_loads', 5)
    setattr(bone.time_vars, 'number_steps', 5)
    runAnalysis(boneConfig, boneLimits, bone)



def multipleAnalysis():
    bone, boneLimits, boneConfig = getBoneData()
    # setattr(boneConfig, 'runFltk', True)
    # setattr(boneConfig, 'writeVTK', True)

    setattr(boneConfig, 'deleteOutput', True)
    setattr(boneConfig, 'runAbq', True)
    setattr(bone.mesh_vars, 'number_elements', 15)
    setattr(bone.geom_vars, 'bone_width', 2.2)
    setattr(bone.geom_vars, 'head_height', 3.0)

    # Output concave single load case
    folder = r'C:\Users\crist\OneDrive\Documents\results\concaveCase\singleLoad'
    setattr(boneConfig, 'inputPath', folder)
    setattr(boneConfig, 'outputPath', folder)
    setattr(bone.geom_vars, 'radius_x', 1.65)
    setattr(bone.geom_vars, 'radius_y', 1.2)
    setattr(bone.geom_vars, 'head_angle', -15.0)
    setattr(bone.geom_vars, 'cart_thick', 2.0)
    setattr(bone.geom_vars, 'curve_angle', 9.0)
    setattr(bone.load_vars, 'number_loads', 5)
    runAnalysis(boneConfig, boneLimits, bone)

    # Output concave double load case
    folder = r'C:\Users\crist\OneDrive\Documents\results\concaveCase\doubleLoad'
    setattr(boneConfig, 'inputPath', folder)
    setattr(boneConfig, 'outputPath', folder)
    setattr(bone.load_vars, 'number_loads', 6)
    runAnalysis(boneConfig, boneLimits, bone)

    # Output convex single load case
    folder = r'C:\Users\crist\OneDrive\Documents\results\convexCase\singleLoad'
    setattr(boneConfig, 'inputPath', folder)
    setattr(boneConfig, 'outputPath', folder)

    setattr(bone.geom_vars, 'radius_x', 1.4)
    setattr(bone.geom_vars, 'radius_y', 1.5)
    setattr(bone.geom_vars, 'head_angle', 25.0)
    setattr(bone.geom_vars, 'cart_thick', 2.5)
    setattr(bone.geom_vars, 'curve_angle', 0.0)
    setattr(bone.load_vars, 'number_loads', 5)
    runAnalysis(boneConfig, boneLimits, bone)

    # Output convex double load case
    folder = r'C:\Users\crist\OneDrive\Documents\results\convexCase\doubleLoad'
    setattr(boneConfig, 'inputPath', folder)
    setattr(boneConfig, 'outputPath', folder)
    setattr(bone.load_vars, 'number_loads', 6)
    runAnalysis(boneConfig, boneLimits, bone)

if __name__ == '__main__':
    singleAnalysis()
