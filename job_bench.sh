#!/bin/sh
#SBATCH --partition=cascade
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --threads-per-core=1
#SBATCH --time=00:20:00
#SBATCH --job-name=lab4_bench
#SBATCH --output=lab4_%j.out
#SBATCH --error=lab4_%j.err

echo $(hostname -s)

module purge
module load oneapi/2022.3.1/ALL

echo "=== 1 ядро ==="
mpirun -np 1 ./yacobi3d
echo "=== 2 ядра ==="
mpirun -np 2 ./yacobi3d
echo "=== 4 ядра ==="
mpirun -np 4 ./yacobi3d
echo "=== 8 ядер ==="
mpirun -np 8 ./yacobi3d
echo "=== 16 ядер ==="
mpirun -np 16 ./yacobi3d
