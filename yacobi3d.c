#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define X0  -1.0
#define Y0  -1.0
#define Z0  -1.0
#define DX   2.0
#define DY   2.0
#define DZ   2.0

#define A_PARAM   1e5
#define EPSILON   1e-8
#define MAX_ITER  50000

#define NX 350
#define NY 350
#define NZ 350

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double hx = DX / NX;
    double hy = DY / NY;
    double hz = DZ / NZ;

    double cx = 1.0 / (hx * hx);
    double cy = 1.0 / (hy * hy);
    double cz = 1.0 / (hz * hz);
    double cc = 2.0 * (cx + cy + cz) + A_PARAM;

    int base = (NZ - 1) / size;
    int rem  = (NZ - 1) % size;
    int local_nz, k_start;

    if (rank < rem) {
        local_nz = base + 1;
        k_start  = rank * local_nz + 1;
    } else {
        local_nz = base;
        k_start  = rem * (base + 1) + (rank - rem) * base + 1;
    }

    int nx1 = NX + 1;
    int ny1 = NY + 1;
    int nz_loc = local_nz + 2;
    int total = nx1 * ny1 * nz_loc;

    double *phi     = (double*)calloc(total, sizeof(double));
    double *phi_new = (double*)calloc(total, sizeof(double));

    for (int kl = 0; kl < nz_loc; kl++) {
        int kg;
        if (kl == 0)            kg = k_start - 1;
        else if (kl == nz_loc-1) kg = k_start + local_nz;
        else                     kg = k_start + kl - 1;

        double z = Z0 + kg * hz;
        for (int j = 0; j < ny1; j++) {
            double y = Y0 + j * hy;
            for (int i = 0; i < nx1; i++) {
                double x = X0 + i * hx;
                int idx = kl * nx1 * ny1 + j * nx1 + i;

                if (kg == 0 || kg == NZ || j == 0 || j == NY || i == 0 || i == NX) {
                    phi[idx] = x*x + y*y + z*z;
                }
                phi_new[idx] = phi[idx];
            }
        }
    }

    int rank_down = (rank > 0)         ? rank - 1 : MPI_PROC_NULL;
    int rank_up   = (rank < size - 1)  ? rank + 1 : MPI_PROC_NULL;

    int plane_size = nx1 * ny1;
    double *buf_sd = (double*)malloc(plane_size * sizeof(double));
    double *buf_rd = (double*)malloc(plane_size * sizeof(double));
    double *buf_su = (double*)malloc(plane_size * sizeof(double));
    double *buf_ru = (double*)malloc(plane_size * sizeof(double));

    MPI_Request reqs[4];
    MPI_Status  stats[4];

    double t_start = MPI_Wtime();
    int iter = 0;
    double max_diff = 1.0;

    while (max_diff > EPSILON && iter < MAX_ITER) {
        memcpy(buf_sd, phi + 1        * nx1 * ny1, plane_size * sizeof(double));
        memcpy(buf_su, phi + local_nz * nx1 * ny1, plane_size * sizeof(double));

        MPI_Isend(buf_sd, plane_size, MPI_DOUBLE, rank_down, 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Irecv(buf_rd, plane_size, MPI_DOUBLE, rank_down, 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Isend(buf_su, plane_size, MPI_DOUBLE, rank_up,   1, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(buf_ru, plane_size, MPI_DOUBLE, rank_up,   0, MPI_COMM_WORLD, &reqs[3]);

        double local_max = 0.0;
        for (int kl = 2; kl <= local_nz - 1; kl++) {
            int kg = k_start + kl - 1;
            double z = Z0 + kg * hz;
            for (int j = 1; j < ny1 - 1; j++) {
                double y = Y0 + j * hy;
                for (int i = 1; i < nx1 - 1; i++) {
                    double x = X0 + i * hx;
                    int idx = kl * nx1 * ny1 + j * nx1 + i;
                    double rho = 6.0 - A_PARAM * phi[idx];
                    phi_new[idx] = ( cx * (phi[idx-1] + phi[idx+1])
                                   + cy * (phi[idx-nx1] + phi[idx+nx1])
                                   + cz * (phi[idx-nx1*ny1] + phi[idx+nx1*ny1])
                                   - rho ) / cc;
                    double d = fabs(phi_new[idx] - phi[idx]);
                    if (d > local_max) local_max = d;
                }
            }
        }

        MPI_Waitall(4, reqs, stats);

        memcpy(phi + 0            * nx1 * ny1, buf_rd, plane_size * sizeof(double));
        memcpy(phi + (local_nz+1) * nx1 * ny1, buf_ru, plane_size * sizeof(double));

        int border_ks[2] = {1, local_nz};
        for (int b = 0; b < 2; b++) {
            int kl = border_ks[b];
            if (kl < 1 || kl > local_nz) continue;
            int kg = k_start + kl - 1;
            double z = Z0 + kg * hz;
            for (int j = 1; j < ny1 - 1; j++) {
                double y = Y0 + j * hy;
                for (int i = 1; i < nx1 - 1; i++) {
                    double x = X0 + i * hx;
                    int idx = kl * nx1 * ny1 + j * nx1 + i;
                    double rho = 6.0 - A_PARAM * phi[idx];
                    phi_new[idx] = ( cx * (phi[idx-1] + phi[idx+1])
                                   + cy * (phi[idx-nx1] + phi[idx+nx1])
                                   + cz * (phi[idx-nx1*ny1] + phi[idx+nx1*ny1])
                                   - rho ) / cc;
                    double d = fabs(phi_new[idx] - phi[idx]);
                    if (d > local_max) local_max = d;
                }
            }
        }

        MPI_Allreduce(&local_max, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        memcpy(phi, phi_new, total * sizeof(double));
        iter++;
        if (rank == 0 && iter % 1000 == 0) {
            printf("Iter %5d | max_diff = %e\n", iter, max_diff);
        }
    }

    double t_end = MPI_Wtime();
    if (rank == 0) {
        if (iter >= MAX_ITER) {
            printf("\nПревышено максимальное число итераций (%d)\n", MAX_ITER);
        } else {
            printf("\nСходимость достигнута за %d итераций\n", iter);
        }
        printf("Время счёта: %.3f сек\n", t_end - t_start);
    }

    double max_err = 0.0;
    for (int kl = 1; kl <= local_nz; kl++) {
        int kg = k_start + kl - 1;
        double z = Z0 + kg * hz;
        for (int j = 0; j < ny1; j++) {
            double y = Y0 + j * hy;
            for (int i = 0; i < nx1; i++) {
                double x = X0 + i * hx;
                double exact = x*x + y*y + z*z;
                int idx = kl * nx1 * ny1 + j * nx1 + i;
                double err = fabs(phi[idx] - exact);
                if (err > max_err) max_err = err;
            }
        }
    }
    double global_err;
    MPI_Allreduce(&max_err, &global_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0) printf("Максимальная ошибка: %e\n", global_err);

    free(phi); free(phi_new);
    free(buf_sd); free(buf_rd); free(buf_su); free(buf_ru);

    MPI_Finalize();
    return 0;
}
