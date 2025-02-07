import os
import subprocess


def runAbaqusMinimal():
    print("Running Abaqus")
    os.chdir("current")
    cmd = r'abaqus job=analisis.inp user=user.for ask_delete=OFF cpus=1'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Print the output and error (if any)
    print("Output:\n", result.stdout)
    print("Error:\n", result.stderr)
    os.chdir("..")


if __name__ == '__main__':
    runAbaqusMinimal()