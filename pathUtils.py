import sys
import os
import json
# from dotenv import load_dotenv


def set_freecad_paths():
    condaPath = "/home/cristian/miniconda3/envs/fcenv/"
    freecadPath = "/usr/lib/freecad/"

    defPath = condaPath
    # load_dotenv(".env")
    FREECAD_BIN = defPath + "bin"
    FREECAD_EXT = defPath + "Ext"
    FREECAD_LIB = defPath + "lib"
    FREECAD_MOD = defPath + "Mod"
    sys.path.append(FREECAD_BIN)
    sys.path.append(FREECAD_EXT)
    sys.path.append(FREECAD_LIB)
    sys.path.append(FREECAD_MOD)

def load_paths(config_file):
    with open(config_file, 'r') as f:
        return json.load(f)

paths = load_paths('config/paths.json')

workspaceFolder = paths['workspaceFolder']
resultsFolder = paths['resultsFolder']
salomePath = paths['salomePath']
paraviewPath = paths['paraviewPath']
medReaderPath = paths['medReaderPath']
temporaryResultsFolder = paths['temporaryResultsFolder']

inputFolder = os.path.join(temporaryResultsFolder , 'input')
outputFolder = os.path.join(temporaryResultsFolder , 'output')
logFolder = os.path.join(temporaryResultsFolder , 'log')

os.makedirs(resultsFolder, exist_ok=True)
os.makedirs(temporaryResultsFolder, exist_ok=True)
os.makedirs(inputFolder, exist_ok=True)
os.makedirs(outputFolder, exist_ok=True)
os.makedirs(logFolder, exist_ok=True)

def getBoneData(jsonConfigurationFile):
    with open(jsonConfigurationFile, 'r') as file:
        data = json.load(file)

    bone = data['bone']
    boneConfig = data['boneConfig']

    return bone, boneConfig