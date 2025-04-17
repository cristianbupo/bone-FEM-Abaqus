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


def modifySketch():

    # Modify sketch
    # Create and write mesh

    Gui.setupWithoutGUI()  # Initialize FreeCAD without GUI
#    if boneConfig.capsule:
    fileName = os.path.join('CADs', 'longBone.FCStd')
#    else:
#        fileName = os.path.join('CADs', 'longBoneCapsule.FCStd')

    doc = App.openDocument(fileName)
    sketch = doc.getObject('Sketch')

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
    g2g.container2gmsh(bone, boneConfig, curvesMesh, curvesArea)
    # gmsh.write(boneConfig.inputPath + '/longBone.msh')
    gmsh.finalize()

    gmsh.initialize()
    # g2g.container2advanceMesh(bone, boneConfig, curvesAdvance)
    # gmsh.write(boneConfig.inputPath + '/advanceMesh.msh')
    gmsh.finalize()

    return curvesArea


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
        "master/propiedades.txt": "propiedades.txt",
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
        modifySketch()  # Modify sketch, create and write mesh
        setupSimulations()  # Set up simulation parameters

    