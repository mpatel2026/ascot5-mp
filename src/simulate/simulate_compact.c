/**
 * @file simulate_compact.c
 * @brief Compaction of the active marker arrays on GPU
 *
 * The GPU simulation loops keep all markers in one array and run every
 * kernel over the first n_mrk elements. When enough markers have finished,
 * the running ones are moved to the front (keeping their order) and n_mrk is
 * reduced, so that the kernels no longer spend time on finished markers.
 * Finished markers stay in the back of the arrays untouched; at the end n_mrk
 * is restored and the end states are written back by their queue index as
 * before.
 *
 * Every per-marker array has to move with its marker. The random numbers do
 * not: they are still generated for the full original array, and each marker
 * takes the values of its original slot, so the results are bit-identical to
 * a run without compaction.
 */
#include <stdlib.h>
#include "simulate_compact.h"

SIMULATE_COMPACT_DEFINE_PERMUTE(permute_real, real)
SIMULATE_COMPACT_DEFINE_PERMUTE(permute_int, int)
SIMULATE_COMPACT_DEFINE_PERMUTE(permute_integer, integer)
SIMULATE_COMPACT_DEFINE_PERMUTE(permute_err, a5err)

/**
 * @brief Allocate compaction buffers for n markers
 *
 * @param c compaction data
 * @param n number of markers (size of the marker arrays)
 */
void simulate_compact_init(simulate_compact_data* c, int n) {
    c->n     = n;
    c->perm  = (int*) malloc(n * sizeof(int));
    c->slot  = (int*) malloc(n * sizeof(int));
    c->tmp_r = (real*) malloc(n * sizeof(real));
    c->tmp_i = (int*) malloc(n * sizeof(int));
    c->tmp_l = (integer*) malloc(n * sizeof(integer));
    c->tmp_e = (a5err*) malloc(n * sizeof(a5err));
    for(int i = 0; i < n; i++) {
        c->perm[i] = i;
        c->slot[i] = i;
    }
    GPU_MAP_TO_DEVICE(c->perm[0:n], c->slot[0:n])
    GPU_MAP_CREATE(c->tmp_r[0:n], c->tmp_i[0:n], c->tmp_l[0:n], c->tmp_e[0:n])
}

/**
 * @brief Free compaction buffers
 *
 * @param c compaction data
 */
void simulate_compact_free(simulate_compact_data* c) {
    int n = c->n;
    GPU_MAP_DELETE_DEVICE(c->perm[0:n], c->slot[0:n], c->tmp_r[0:n],
                          c->tmp_i[0:n], c->tmp_l[0:n], c->tmp_e[0:n])
    free(c->perm);
    free(c->slot);
    free(c->tmp_r);
    free(c->tmp_i);
    free(c->tmp_l);
    free(c->tmp_e);
}

/**
 * @brief Build the permutation that moves running markers to the front
 *
 * Copies the running flags of the active markers to the host, builds the
 * permutation (running markers first, then finished ones, both in their
 * current order), copies it to the device, and permutes the slot array. The
 * caller must have waited for the GPU queue and then permute every
 * per-marker array with the same permutation.
 *
 * @param c compaction data
 * @param running running flags of the markers (device array)
 * @param n_active number of active markers
 *
 * @return number of running markers, i.e. the new number of active markers
 */
int simulate_compact_plan(simulate_compact_data* c, integer* running,
                          int n_active) {
    GPU_UPDATE_FROM_DEVICE(running[0:n_active])
    int k = 0;
    for(int j = 0; j < n_active; j++) {
        if(running[j] > 0) {
            c->perm[k++] = j;
        }
    }
    int n_running = k;
    for(int j = 0; j < n_active; j++) {
        if(running[j] <= 0) {
            c->perm[k++] = j;
        }
    }
    GPU_UPDATE_TO_DEVICE(c->perm[0:n_active])
    permute_int(c->slot, c->tmp_i, c->perm, n_active);
    return n_running;
}

/**
 * @brief Permute a real per-marker array with the current permutation
 *
 * @param c compaction data
 * @param a array to permute (device array)
 * @param n_active number of active markers before compaction
 */
void simulate_compact_real(simulate_compact_data* c, real* a, int n_active) {
    permute_real(a, c->tmp_r, c->perm, n_active);
}

/**
 * @brief Permute all fields of a guiding center marker array
 *
 * @param c compaction data
 * @param p marker array
 * @param n_active number of active markers before compaction
 */
void simulate_compact_gc(simulate_compact_data* c, particle_simd_gc* p,
                         int n_active) {
    real* reals[] = {
        p->r, p->phi, p->z, p->ppar, p->mu, p->zeta, p->mass, p->charge,
        p->time, p->B_r, p->B_phi, p->B_z, p->B_r_dr, p->B_phi_dr, p->B_z_dr,
        p->B_r_dphi, p->B_phi_dphi, p->B_z_dphi, p->B_r_dz, p->B_phi_dz,
        p->B_z_dz, p->weight, p->cputime, p->rho, p->theta, p->mileage};
    for(size_t k = 0; k < sizeof(reals) / sizeof(reals[0]); k++) {
        permute_real(reals[k], c->tmp_r, c->perm, n_active);
    }
    permute_int(p->bounces, c->tmp_i, c->perm, n_active);
    integer* integers[] = {p->id, p->endcond, p->walltile, p->running,
                           p->index};
    for(size_t k = 0; k < sizeof(integers) / sizeof(integers[0]); k++) {
        permute_integer(integers[k], c->tmp_l, c->perm, n_active);
    }
    permute_err(p->err, c->tmp_e, c->perm, n_active);
}

/**
 * @brief Permute all fields of a particle marker array
 *
 * @param c compaction data
 * @param p marker array
 * @param n_active number of active markers before compaction
 */
void simulate_compact_fo(simulate_compact_data* c, particle_simd_fo* p,
                         int n_active) {
    real* reals[] = {
        p->r, p->phi, p->z, p->p_r, p->p_phi, p->p_z, p->mass, p->charge,
        p->time, p->B_r, p->B_phi, p->B_z, p->B_r_dr, p->B_phi_dr, p->B_z_dr,
        p->B_r_dphi, p->B_phi_dphi, p->B_z_dphi, p->B_r_dz, p->B_phi_dz,
        p->B_z_dz, p->weight, p->cputime, p->rho, p->theta, p->mileage};
    for(size_t k = 0; k < sizeof(reals) / sizeof(reals[0]); k++) {
        permute_real(reals[k], c->tmp_r, c->perm, n_active);
    }
    int* ints[] = {p->znum, p->anum, p->bounces};
    for(size_t k = 0; k < sizeof(ints) / sizeof(ints[0]); k++) {
        permute_int(ints[k], c->tmp_i, c->perm, n_active);
    }
    integer* integers[] = {p->id, p->endcond, p->walltile, p->running,
                           p->index};
    for(size_t k = 0; k < sizeof(integers) / sizeof(integers[0]); k++) {
        permute_integer(integers[k], c->tmp_l, c->perm, n_active);
    }
    permute_err(p->err, c->tmp_e, c->perm, n_active);
}

/**
 * @brief Pick each active marker's random numbers, interleaved layout
 *
 * rnd_full holds k numbers per original slot, slot s at [k*s, k*s+k). Sets
 * rnd[k*j + l] = rnd_full[k*slot[j] + l] for the active markers j.
 *
 * @param c compaction data
 * @param rnd random numbers for the active markers (device array)
 * @param rnd_full random numbers for all slots (device array)
 * @param n_active number of active markers
 * @param k random numbers per marker
 */
void simulate_compact_rnd_interleaved(simulate_compact_data* c, real* rnd,
                                      const real* rnd_full, int n_active,
                                      int k) {
    int n = c->n;
    int* slot = c->slot;
    GPU_PARALLEL_LOOP_ALL_LEVELS_PRESENT(rnd[0:k*n], rnd_full[0:k*n],
                                         slot[0:n])
    for(int j = 0; j < n_active; j++) {
        for(int l = 0; l < k; l++) {
            rnd[k*j + l] = rnd_full[k*slot[j] + l];
        }
    }
}

/**
 * @brief Pick each active marker's random numbers, blocked layout
 *
 * rnd_full holds k blocks of n numbers, slot s at [l*n + s]. Sets
 * rnd[l*n_active + j] = rnd_full[l*n + slot[j]] for the active markers j.
 *
 * @param c compaction data
 * @param rnd random numbers for the active markers (device array)
 * @param rnd_full random numbers for all slots (device array)
 * @param n_active number of active markers
 * @param k random numbers per marker
 */
void simulate_compact_rnd_blocked(simulate_compact_data* c, real* rnd,
                                  const real* rnd_full, int n_active, int k) {
    int n = c->n;
    int* slot = c->slot;
    GPU_PARALLEL_LOOP_ALL_LEVELS_PRESENT(rnd[0:k*n], rnd_full[0:k*n],
                                         slot[0:n])
    for(int j = 0; j < n_active; j++) {
        for(int l = 0; l < k; l++) {
            rnd[l*n_active + j] = rnd_full[l*n + slot[j]];
        }
    }
}
