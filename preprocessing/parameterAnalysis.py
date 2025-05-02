import os
import subprocess
import shutil
import time
import averageIO as avio
from datetime import datetime


class args:
    def __init__(self, ne=12, rx=2.0, ry=2.5, ha=20.0, bw=2.4, hh=4.5, bh=-1.0, ca=15,
                 hl=0.0, kl=6.85, rl=0.81, koi=0.0, runAbq=False,
                 showLoad=False, renderRaw=False, renderMesh=False, renderArea=False,
                 renderLength=False, runFltk=False, skipWrite=False, writeVTK=False,
                 runMeshing=False, pickleIn=False, pickleOut=False,
                 deleteBackUp=False):

        # Mesh parameters
        self.ne = ne

        # Geometry parameters
        self.radius_x = rx
        self.radius_y = ry
        self.head_angle = ha
        self.bone_width = bw
        self.head_height = hh
        self.bone_height = bh
        self.curve_angle = ca

        # Load parameters
        self.load_center = hl
        self.load_amplitude = kl
        self.load_radius = rl
        self.kOI = koi

        # Analysis parameters
        self.runAbq = runAbq
        self.showLoad = showLoad
        self.renderRaw = renderRaw
        self.renderMesh = renderMesh
        self.renderArea = renderArea
        self.renderLength = renderLength
        self.runFltk = runFltk
        self.skipWrite = skipWrite
        self.writeVTK = writeVTK
        self.runMeshing = runMeshing
        self.pickleIn = pickleIn

        # Output parameters
        self.pickleOut = pickleOut
        self.deleteBackUp = deleteBackUp

    def default_run(self):
        # Construct the command string
        command = (
            f'analisis.bat "-ne {self.ne} -rx {self.radius_x} -ry {self.radius_y} '
            f'-ha {self.head_angle} -bw {self.bone_width} -hh {self.head_height} '
            f'-bh {self.bone_height} -ca {self.curve_angle} '
            f'{"--showLoad" if self.showLoad else ""} '
            f'{"--renderRaw" if self.renderRaw else ""} '
            f'{"--renderMesh" if self.renderMesh else ""} '
            f'{"--renderArea" if self.renderArea else ""} '
            f'{"--renderLength" if self.renderLength else ""} '
            f'{"--runFltk" if self.runFltk else ""} '
            f'{"--skipWrite" if self.skipWrite else ""} '
            f'{"--writeVTK" if self.writeVTK else ""} '
            f'{"--runMeshing" if self.runMeshing else ""} '
            f'{"--pickleIn" if self.pickleIn else ""} '
            f'{"--pickleOut" if self.pickleOut else ""} '
            f'{"--deleteBackUp" if self.deleteBackUp else ""} '
            f'-hl {self.load_center} -kl {self.load_amplitude} -rl {self.load_radius} '
            f'-koi {self.kOI}" {str(self.runAbq).lower()}'
        )

        return run_command(command)


def exportResults(originFolder, destinationFolder, fileID):
    analysis_dir = os.path.join(originFolder, "analisis.vtu")
    mesh_dir = os.path.join(originFolder, "malla.vtu")
    dist_dir = os.path.join(originFolder, "carga.vtp")

    if os.path.exists(analysis_dir):
        shutil.copy(analysis_dir, os.path.join(destinationFolder,
                                               f"analisis{fileID}.vtu"))

    if os.path.exists(mesh_dir):
        shutil.copy(mesh_dir, os.path.join(destinationFolder,
                                           f"malla{fileID}.vtu"))

    if os.path.exists(dist_dir):
        shutil.copy(dist_dir, os.path.join(destinationFolder,
                                           f"carga{fileID}.vtu"))


def run_command(command):
    # Remove extra spaces
    command = ' '.join(command.split())

    # Run the command
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Print the output and error (if any)
    print("Output:\n", result.stdout)
    print("Error:\n", result.stderr)

    return command


def parametricAnalysisKoi():
    # Define base path
    base_path = "D:\\bone-FEM-results"
    # Get the current date and time to use in the ID string
    id_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    #  Define the directory
    dir_name = os.path.join(base_path, f"test_{id_time}")

    # Check if the directory exists
    if not os.path.exists(dir_name):
        # Create the directory
        os.makedirs(dir_name)

    with open(f"log/manyResults_{id_time}.txt", 'w') as log_file:
        # Get the current timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_file.write(f"{timestamp} - Starting the loop\n\n")

        # Record the start time
        start_time = time.time()
        previous_time = start_time

        N = 5
        # Loop over values of i
        for i in range(0, N):
            # Format the koi value as a decimal with two digits
            koi_val = f"{1.0*i/(N-1):.2f}"

            string = f"\n case{i}:\n"

            params = args()
            params.kOI = koi_val
            params.runAbq = True
            command = params.run()
            string += f"Command: {command}\n"

            analysis_dir = "current\\analisis.vtu"
            if os.path.exists(analysis_dir):
                shutil.copy(analysis_dir, f"{dir_name}\\analisis{i}.vtu")
                string += ("Analysis completed successfully\n")
            else:
                string += (f"Error: current\\analisis.vtu"
                           f" does not exist for case{i}\n")

            if os.path.exists("current\\malla.vtu"):
                shutil.copy("current\\malla.vtu", f"{dir_name}\\malla{i}.vtu")

            distribution_dir = "current\\carga.vtp"
            if os.path.exists(distribution_dir):
                shutil.copy(distribution_dir, f"{dir_name}\\carga{i}.vtp")
            else:
                string += (f"Error: current\\carga.vtp"
                           f" does not exist for case{i}\n")

            # Calculate the elapsed time
            elapsed_time = time.time() - start_time
            step_time = time.time() - previous_time
            previous_time = time.time()

            # Get the current timestamp again
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            string += (f"{timestamp} (Elapsed time: {elapsed_time:.2f} "
                       f"seconds, Step time:{step_time:.2f})\n")

            log_file.write(string)
            print(string)

    # Calculate the average of the IO values
    avio.calculate_average(dir_name)


def parametricAnalysisVariable():
    # Define base path
    base_path = "D:\\bone-FEM-results"
    # Get the current date and time to use in the ID string
    id_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Define the directory
    dir_name = os.path.join(base_path, f"test_{id_time}")

    # Check if the directory exists
    if not os.path.exists(dir_name):
        # Create the directory
        os.makedirs(dir_name)

    with open(f"log/manyResults_{id_time}.txt", 'w') as log_file:
        # Get the current timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_file.write(f"{timestamp} - Starting the loop\n\n")

        # Record the start time
        start_time = time.time()
        previous_time = start_time

        # Define the number of steps, cases increase by N^3
        N = 5
        # Loop over values of radius_x, radius_y, and head_max (hm_val)
        index = 0  # Initialize the index
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    # Format the radius_x, radius_y, and head_max values
                    ha_ran = [-20.0, 20.0]
                    rx_ran = [1.2, 2.5]
                    ry_ran = [1.2, 2.5]

                    ha_val = ha_ran[0] + (ha_ran[1] - ha_ran[0]) * i / (N - 1)
                    rx_val = ry_ran[0] + (ry_ran[1] - ry_ran[0]) * j / (N - 1)
                    ry_val = rx_ran[0] + (rx_ran[1] - rx_ran[0]) * k / (N - 1)

                    string = f"\n case{index} - i, j, k = {i}, {j}, {k}:\n"

                    params = args(rx=rx_val, ry=ry_val, ha=ha_val)
                    params.writeVTK = True
                    params.runAbq = True
                    command = params.run()
                    string += f"Command: {command}\n"

                    analysis_dir = "current\\analisis.vtu"
                    if os.path.exists(analysis_dir):
                        shutil.copy(analysis_dir, f"{dir_name}\\analisis{index}.vtu")
                        string += ("Analysis completed successfully\n")
                    else:
                        string += (f"Error: current\\analisis.vtu"
                                   f" does not exist for case{index}\n")

                    if os.path.exists("current\\malla.vtu"):
                        shutil.copy("current\\malla.vtu",
                                    f"{dir_name}\\malla{index}.vtu")

                    distribution_dir = "current\\carga.vtp"
                    if os.path.exists(distribution_dir):
                        shutil.copy(distribution_dir,
                                    f"{dir_name}\\carga{index}.vtp")
                    else:
                        string += (f"Error: current\\carga.vtp"
                                   f" does not exist for case{index}\n")

                    # Calculate the elapsed time
                    elapsed_time = time.time() - start_time
                    step_time = time.time() - previous_time
                    previous_time = time.time()

                    # Get the current timestamp again
                    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                    string += (f"{timestamp} (Elapsed time: {elapsed_time:.2f} "
                               f"seconds, Step time:{step_time:.2f})\n")

                    log_file.write(string)
                    print(string)
                    index += 1

    # Calculate the average of the IO values
    avio.calculate_average(dir_name)


def notMain(sigma, rl, N):

    destinationFolder = "currentLoadHistory"
    shutil.rmtree(destinationFolder, ignore_errors=True)
    os.makedirs(destinationFolder, exist_ok=True)

    case = args()
    index = 0

    for i in range(0, 5):
        case.load_center = (sigma - rl)*(2 * i - N + 1)/(2*N - 2)
        # Calculate the 5 load cases with different load_center (hl) values
        if i == 0:
            # Run the first analysis
            case.pickleOut = True  # Store the data
            case.deleteBackUp = True  # Delete the backup directory
            case.writeVTK = True  # Write the VTK mesh
            case.runMeshing = True  # Run the meshing
            case.pickleIn = True
            case.runAbq = True
            case.default_run()  # Run the loading
        else:
            # Run the other analyses with the same mesh
            case.pickleOut = False  # Do not store the data
            case.deleteBackUp = False  # Do not delete the backup directory
            case.writeVTK = True  # Do not write the VTK mesh
            case.runMeshing = False  # Do not run the meshing
            case.pickleIn = True  # Run the loading
            case.runAbq = True
            case.default_run()  # Run the loading

        exportResults("current", destinationFolder, index)
        index += 1


def notMain2(sigma, rl, N):

    destinationFolder = "currentLoadHistory"
    shutil.rmtree(destinationFolder, ignore_errors=True)
    os.makedirs(destinationFolder, exist_ok=True)

    case = args()

    for i in range(0, 5):
        case.load_center = (sigma - rl)*(2 * i - N + 1)/(2*N - 2)
        # Calculate the 5 load cases with different load_center (hl) values
        if i == 0:
            # Run the first analysis
            case.pickleOut = True  # Store the data
            case.deleteBackUp = True  # Delete the backup directory
            case.writeVTK = True  # Write the VTK mesh
            case.runMeshing = True  # Run the meshing
            case.pickleIn = True
            case.runAbq = True
            case.default_run()  # Run the loading
        else:
            # Run the other analyses with the same mesh
            case.pickleOut = False  # Do not store the data
            case.deleteBackUp = False  # Do not delete the backup directory
            case.writeVTK = True  # Do not write the VTK mesh
            case.runMeshing = False  # Do not run the meshing
            case.pickleIn = True  # Run the loading
            case.runAbq = True
            case.default_run()  # Run the loading

def main():
    print("lol xd")

if __name__ == "__main__":
    main()