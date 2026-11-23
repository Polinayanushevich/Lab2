#!/bin/sh
#SBATCH --partition=cascade
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --threads-per-core=1
#SBATCH --time=00:40:00
#SBATCH --job-name=lab4_bench
#SBATCH --output=lab2_bench_%j.out
#SBATCH --error=lab2_bench_%j.err

echo "Host: $(hostname -s)"
module purge
module load oneapi/2022.3.1/ALL

for NP in 1 2 4 8 16; do
    echo "===== np=$NP ====="
    mpirun -np $NP -bind-to core ./yacobi3d
    echo
done
