import os
import shutil
import subprocess
import sys
import time
import longBone as lb
from parameterAnalysis import args


def main(args):
    print("Deleting old results")
    for filename in os.listdir("current"):
        file_path = os.path.join("current", filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

    print("Meshing and applying loads")
    lb.longBoneFunction(args, sketch)

    calculate(run_abaqus)


def calculate(run_abaqus):
    if run_abaqus.lower() == "true":
        print("Running Abaqus")
        shutil.copy("analisisMaster.inp", "current/analisis.inp")
        shutil.copy("userMaster.for", "current/user.for")
        shutil.copy("posMaster.pvsm", "current/pos.pvsm")

        os.chdir("current")
        subprocess.run(["../runAbaqus.bat"])
        os.chdir("..")

        print("Waiting for Abaqus to finish")
        while os.path.exists("current/analisis.lck"):
            time.sleep(1)
        print("Lock file removed")
        print("Abaqus finished")

def show_help(params):
    subprocess.run([sys.executable, "longBone.py", params])

if __name__ == "__main__":
    main()