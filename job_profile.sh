#!/bin/sh
#SBATCH --partition=cascade
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --threads-per-core=1
#SBATCH --time=00:10:00
#SBATCH --job-name=lab2_profile
#SBATCH --output=lab2_profile_%j.out
#SBATCH --error=lab2_profile_%j.err

echo $(hostname -s)

module purge
module load oneapi/2022.3.1/ALL

export PATH=$HOME/mpe2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/mpe2/lib:$LD_LIBRARY_PATH

mpirun -np 16 ./yacobi3d_prof
