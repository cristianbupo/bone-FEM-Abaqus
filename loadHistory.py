import argparse
import subprocess
import shutil
import os

# Create the argument parser
parser = argparse.ArgumentParser(description="Parse sigma, h, k, and r values.")

# Add arguments
parser.add_argument('--sigma', type=float, required=True, help='Value for sigma')
# parser.add_argument('--h', type=float, required=True, help='Value for h')
parser.add_argument('--k', type=float, required=True, help='Value for k')
parser.add_argument('--r', type=float, required=True, help='Value for r')

# Parse the arguments
args = parser.parse_args()

# Print the values
print(f"sigma: {args.sigma}")
# print(f"h: {args.h}")
print(f"k: {args.k}")
print(f"r: {args.r}")

sigma = args.sigma
kl = args.k
rl = args.r

cmd = ('analisis.bat "-ne 12 -rx 2.0 -ry 2.0 -hm 1.8 -bw 2.4 -hh 4.5 -bh -1.0 -ca 15 '
       '--writeVTK -hl {} -kl {} -rl {} -koi 0.25" true')

# For convex case
L = 5.967507723310594  # mm
N = 10

for i in range(0, N):
    hl = (sigma - rl)*(2 * i - N + 1)/(2*N - 2)
    print(cmd.format(hl, kl, rl))
    subprocess.run(cmd.format(hl, kl, rl))

    analysis_dir = "current\\analisis.vtu"
    if os.path.exists(analysis_dir):
        shutil.copy(analysis_dir, f"currentLoadHistory\\result{i}.vtu")

    mesh_dir = "current\\malla.vtu"
    if os.path.exists(mesh_dir):
        shutil.copy(mesh_dir, f"currentLoadHistory\\malla{i}.vtu")

    dist_dir = "current\\carga.vtp"
    if os.path.exists(dist_dir):
        shutil.copy(dist_dir, f"currentLoadHistory\\carga{i}.vtp")

    dist_dir_2 = "current\\carga2.inp"
    if os.path.exists(dist_dir_2):
        shutil.copy(dist_dir_2, f"currentLoadHistory\\carga{i}.inp")