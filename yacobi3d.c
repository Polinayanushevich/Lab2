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

#define A_PARAM   1e-5
#define EPSILON   1e-6
#define MAX_ITER  500000

#define NX 128
#define NY 128
#define NZ 128

int main(int argc, char** argv) {
    int rank, size;
    double hx, hy, hz;
    double cx, cy, cz, cc;
    int base, rem, local_nz, k_start;
    int nx1, ny1, nz_loc, total;
    double *phi, *phi_new;
    int kl, kg, j, i, idx, b;
    double z, y, x, rho, d, exact;
    int rank_down, rank_up;
    int plane;
    double *buf_sd, *buf_su, *buf_rd, *buf_ru;
    MPI_Request reqs[4];
    MPI_Status  stats[4];
    double t_start, t_end;
    int iter;
    double max_diff, local_max;
    int border_ks[4];
    double max_err, global_err;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hx = DX / (NX - 1);
    hy = DY / (NY - 1);
    hz = DZ / (NZ - 1);

    cx = 1.0 / (hx * hx);
    cy = 1.0 / (hy * hy);
    cz = 1.0 / (hz * hz);
    cc = 2.0 * (cx + cy + cz) + A_PARAM;

    base = (NZ - 1) / size;
    rem  = (NZ - 1) % size;

    if (rank < rem) {
        local_nz = base + 1;
        k_start  = rank * local_nz + 1;
    } else {
        local_nz = base;
        k_start  = rem * (base + 1) + (rank - rem) * base + 1;
    }

    nx1 = NX + 1;
    ny1 = NY + 1;
    nz_loc = local_nz + 2;
    total = nx1 * ny1 * nz_loc;

    phi     = (double*)calloc(total, sizeof(double));
    phi_new = (double*)calloc(total, sizeof(double));

    for (kl = 0; kl < nz_loc; kl++) {
        if (kl == 0)             kg = k_start - 1;
        else if (kl == nz_loc-1) kg = k_start + local_nz;
        else                     kg = k_start + kl - 1;

        z = Z0 + kg * hz;
        for (j = 0; j < ny1; j++) {
            y = Y0 + j * hy;
            for (i = 0; i < nx1; i++) {
                x = X0 + i * hx;
                idx = kl * nx1 * ny1 + j * nx1 + i;
                if (kg == 0 || kg == NZ || j == 0 || j == NY || i == 0 || i == NX)
                    phi[idx] = x*x + y*y + z*z;
                else
                    phi[idx] = 0.0;
                phi_new[idx] = phi[idx];
            }
        }
    }

    rank_down = (rank > 0)        ? rank - 1 : MPI_PROC_NULL;
    rank_up   = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

    plane = nx1 * ny1;
    buf_sd = (double*)malloc(plane * sizeof(double));
    buf_su = (double*)malloc(plane * sizeof(double));
    buf_rd = (double*)malloc(plane * sizeof(double));
    buf_ru = (double*)malloc(plane * sizeof(double));

    t_start = MPI_Wtime();
    iter = 0;
    max_diff = 1.0;

    while (max_diff > EPSILON && iter < MAX_ITER) {
        memcpy(buf_sd, phi + 1          * nx1 * ny1, plane * sizeof(double));
        memcpy(buf_su, phi + local_nz   * nx1 * ny1, plane * sizeof(double));

        MPI_Isend(buf_sd, plane, MPI_DOUBLE, rank_up,   0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(buf_su, plane, MPI_DOUBLE, rank_down, 1, MPI_COMM_WORLD, &reqs[1]);
        MPI_Irecv(buf_rd, plane, MPI_DOUBLE, rank_down, 1, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(buf_ru, plane, MPI_DOUBLE, rank_up,   0, MPI_COMM_WORLD, &reqs[3]);

        local_max = 0.0;

        /* Фаза 1: глубоко внутренние слои, не зависят от гало */
        for (kl = 3; kl <= local_nz - 2; kl++) {
            kg = k_start + kl - 1;
            z = Z0 + kg * hz;
            for (j = 1; j < ny1 - 1; j++) {
                y = Y0 + j * hy;
                for (i = 1; i < nx1 - 1; i++) {
                    x = X0 + i * hx;
                    idx = kl * nx1 * ny1 + j * nx1 + i;
                    rho = 6.0 - A_PARAM * (x*x + y*y + z*z);
                    phi_new[idx] = ( cx * (phi[idx-1] + phi[idx+1])
                                   + cy * (phi[idx-nx1] + phi[idx+nx1])
                                   + cz * (phi[idx-nx1*ny1] + phi[idx+nx1*ny1])
                                   - rho ) / cc;
                    d = fabs(phi_new[idx] - phi[idx]);
                    if (d > local_max) local_max = d;
                }
            }
        }

        MPI_Waitall(4, reqs, stats);

        memcpy(phi + 0                * nx1 * ny1, buf_rd, plane * sizeof(double));
        memcpy(phi + (local_nz + 1)   * nx1 * ny1, buf_ru, plane * sizeof(double));

        /* Фаза 2: приграничные слои (нужны свежие гало) */
        border_ks[0] = 1;
        border_ks[1] = 2;
        border_ks[2] = local_nz - 1;
        border_ks[3] = local_nz;
        for (b = 0; b < 4; b++) {
            kl = border_ks[b];
            if (kl < 1 || kl > local_nz) continue;
            kg = k_start + kl - 1;
            z = Z0 + kg * hz;
            for (j = 1; j < ny1 - 1; j++) {
                y = Y0 + j * hy;
                for (i = 1; i < nx1 - 1; i++) {
                    x = X0 + i * hx;
                    idx = kl * nx1 * ny1 + j * nx1 + i;
                    rho = 6.0 - A_PARAM * (x*x + y*y + z*z);
                    phi_new[idx] = ( cx * (phi[idx-1] + phi[idx+1])
                                   + cy * (phi[idx-nx1] + phi[idx+nx1])
                                   + cz * (phi[idx-nx1*ny1] + phi[idx+nx1*ny1])
                                   - rho ) / cc;
                    d = fabs(phi_new[idx] - phi[idx]);
                    if (d > local_max) local_max = d;
                }
            }
        }

        MPI_Allreduce(&local_max, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        memcpy(phi, phi_new, total * sizeof(double));
        iter++;
        if (rank == 0 && iter % 1000 == 0)
            printf("Iter %5d | max_diff = %e\n", iter, max_diff);
    }

    t_end = MPI_Wtime();
    if (rank == 0) {
        if (iter >= MAX_ITER) printf("\nПревышено MAX_ITER (%d)\n", MAX_ITER);
        else                  printf("\nСходимость за %d итераций\n", iter);
        printf("Время счёта: %.3f сек\n", t_end - t_start);
    }

    max_err = 0.0;
    for (kl = 1; kl <= local_nz; kl++) {
        kg = k_start + kl - 1;
        z = Z0 + kg * hz;
        for (j = 0; j < ny1; j++) {
            y = Y0 + j * hy;
            for (i = 0; i < nx1; i++) {
                x = X0 + i * hx;
                exact = x*x + y*y + z*z;
                idx = kl * nx1 * ny1 + j * nx1 + i;
                d = fabs(phi[idx] - exact);
                if (d > max_err) max_err = d;
            }
        }
    }
    MPI_Allreduce(&max_err, &global_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0) printf("Максимальная ошибка: %e\n", global_err);

    free(phi); free(phi_new);
    free(buf_sd); free(buf_rd); free(buf_su); free(buf_ru);

    MPI_Finalize();
    return 0;
}
