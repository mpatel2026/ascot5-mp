/**
 * @file simulate_compact.h
 * @brief Compaction of the active marker arrays on GPU
 */
#ifndef SIMULATE_COMPACT_H
#define SIMULATE_COMPACT_H

#include "../ascot5.h"
#include "../error.h"
#include "../offload.h"
#include "../particle.h"

/**
 * @brief Whether the GPU simulation loops compact their marker arrays
 *
 * On GPU all markers are loaded at once and every kernel runs over the whole
 * array, so finished markers keep costing time: a warp of 32 threads costs as
 * much as a full one as long as one of its markers is running. Compaction
 * moves the running markers to the front of the arrays and shrinks the
 * number of markers the kernels run over. Results are unchanged.
 */
#ifndef A5_GPU_COMPACT
#define A5_GPU_COMPACT 1
#endif

/**
 * @brief Compact when fewer than this fraction of the active markers run
 */
#ifndef A5_GPU_COMPACT_FRAC
#define A5_GPU_COMPACT_FRAC 0.75
#endif

/**
 * @brief Do not compact below this many active markers
 */
#ifndef A5_GPU_COMPACT_MIN
#define A5_GPU_COMPACT_MIN 256
#endif

/**
 * @brief Work buffers for compaction, one element per marker
 */
typedef struct {
    int n;          /**< Total number of markers (array size)            */
    int* perm;      /**< Permutation: new position j takes old perm[j]   */
    int* slot;      /**< Original array position of the marker at each
                         position; random numbers are taken by slot      */
    real* tmp_r;    /**< Temporary for real arrays                       */
    int* tmp_i;     /**< Temporary for int arrays                        */
    integer* tmp_l; /**< Temporary for integer arrays                    */
    a5err* tmp_e;   /**< Temporary for error flags                       */
} simulate_compact_data;

/**
 * @brief Define a function that permutes the first n elements of an array
 *
 * Defines name(a, tmp, perm, n) which sets a[j] = a[perm[j]] for j < n, using
 * tmp as scratch space. All arrays must be mapped to the device.
 */
#define SIMULATE_COMPACT_DEFINE_PERMUTE(name, T)                               \
static void name(T* a, T* tmp, const int* perm, int n) {                      \
    GPU_PARALLEL_LOOP_ALL_LEVELS_PRESENT(a[0:n], tmp[0:n], perm[0:n])          \
    for(int j = 0; j < n; j++) {                                               \
        tmp[j] = a[perm[j]];                                                   \
    }                                                                          \
    GPU_PARALLEL_LOOP_ALL_LEVELS_PRESENT(a[0:n], tmp[0:n])                     \
    for(int j = 0; j < n; j++) {                                               \
        a[j] = tmp[j];                                                         \
    }                                                                          \
}

void simulate_compact_init(simulate_compact_data* c, int n);
void simulate_compact_free(simulate_compact_data* c);
int simulate_compact_plan(simulate_compact_data* c, integer* running,
                          int n_active);
void simulate_compact_real(simulate_compact_data* c, real* a, int n_active);
void simulate_compact_gc(simulate_compact_data* c, particle_simd_gc* p,
                         int n_active);
void simulate_compact_fo(simulate_compact_data* c, particle_simd_fo* p,
                         int n_active);
void simulate_compact_rnd_interleaved(simulate_compact_data* c, real* rnd,
                                      const real* rnd_full, int n_active,
                                      int k);
void simulate_compact_rnd_blocked(simulate_compact_data* c, real* rnd,
                                  const real* rnd_full, int n_active, int k);

#endif
