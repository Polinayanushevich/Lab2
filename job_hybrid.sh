#!/bin/sh
#SBATCH --partition=cascade
#SBATCH --nodes=4                 # <--- 4 узла
#SBATCH --ntasks=64               # Всего 64 MPI-процесса (максимум)
#SBATCH --cpus-per-task=16        # Даем по 16 ядер на узел
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00
#SBATCH --job-name=hybrid_bench
#SBATCH --output=hybrid_%j.out
#SBATCH --error=hybrid_%j.err

module purge
module load oneapi/2022.3.1/ALL

# Функция для запуска
run_test() {
    local P=$1
    local T=$2
    echo "===== Запуск: Процессов (MPI) = $P, Потоков (OpenMP) = $T ====="
    export OMP_NUM_THREADS=$T
    mpirun -np $P -bind-to core ./yacobi3d_omp
    echo ""
}

run_test 64 1
run_test 32 2
run_test 16 4
run_test 8  8
run_test 4  16
