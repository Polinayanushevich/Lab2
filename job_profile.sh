#!/bin/sh
#SBATCH --partition=cascade
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --threads-per-core=1
#SBATCH --time=00:20:00
#SBATCH --job-name=lab4_profile
#SBATCH --output=lab4_profile_%j.out
#SBATCH --error=lab4_profile_%j.err

echo "Host: $(hostname -s)"
module purge
module load oneapi/2022.3.1/ALL

# MPE (если ставили локально)
export PATH=$HOME/mpe2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/mpe2/lib:$LD_LIBRARY_PATH

# Профилировка: пишем лог
mpirun -np 16 ./yacobi3d_prof -mpilog -logfile=lab4.slog2
echo "Лог: lab4.slog2 — открыть: jumpshot lab4.slog2"
