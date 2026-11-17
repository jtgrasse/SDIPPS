//------------------------------------------------------------------------------
// SPEX_SDIPPS/Include/SDIPPS.h: SDIPPS column replacement for SPEX LU
//------------------------------------------------------------------------------

// SDIPPS (Sparse Diagonal Integer-Preserving Push-and-Swap) replaces one
// column of A in an updatable SPEX LU factorization
//
//      P A Q = L D^{-1} U,     D = diag (rho_{-1} rho_0, ..., rho_{n-2} rho_{n-1})
//
// with rho_{-1} = 1.  The exiting column is pushed from its position xi
// toward the end, one frame at a time, by diagonal pushes (swap rows and
// columns k, k+1) and column pushes (swap columns k, k+1).  The frame being
// pushed is carried in two seed vectors.  When the push reaches the end (or
// the landing position of an early exit), the entering column is put in by
// one left-looking REF forward solve.  The update is done in place.
//
// The factorization F keeps SPEX's updatable format, so SPEX_update_solve and
// SPEX_update_tsolve can be used on it afterwards, and SPEX_update_lu_colrep
// and SDIPPS_LU_col_rep can be mixed.  A vector of L or U may carry a pending
// scale (SPEX_vector->scale): its true entries are x[p] * scale.  SDIPPS
// reads and writes that convention.

#ifndef SDIPPS_H
#define SDIPPS_H

#include "SPEX.h"

//------------------------------------------------------------------------------
// options
//------------------------------------------------------------------------------

// Which push to take at frame k.  With sL, sU the seed entries at position
// k+1 (the supports) and pi = rho~_{k-1} rho_{k+1} + sL sU:
//
//  - if sL * sU == 0, only a diagonal push is possible (then pi != 0);
//  - if pi == 0, only a column push is possible (then sL * sU != 0);
//  - otherwise both are possible, and push_choice decides.  A diagonal push
//    keeps both seeds going; a column push restarts the L seed.

typedef enum
{
    SDIPPS_PREFER_DIAGONAL = 0,     // diagonal push whenever pi != 0
    SDIPPS_PREFER_COLUMN = 1        // column push whenever sL * sU != 0
}
SDIPPS_push_choice ;

// How a backtracking step that is a pure scaling is done.  That happens in a
// diagonal push when the other half's support is 0: the new frame is then the
// old one times lambda = rho~_{k-1} / rho_k.  Along a run of diagonal pushes
// with sL * sU = 0, lambda stays the same, so EAGER and LAZY reduce it
// (SDIPPS_Reduce) once, at the first push of the run.  The same choice is
// used for the frames past an early exit, which scale by one ratio.

typedef enum
{
    SDIPPS_SCALING_NONE = 0,    // the regular backtrack: every entry times
                                // rho~_{k-1}, divided by rho_k
    SDIPPS_SCALING_EAGER = 1,   // every entry divided by q, then times p,
                                // with p / q the reduced ratio; no pending
                                // scales are left
    SDIPPS_SCALING_LAZY = 2     // p / q is multiplied into the vector's
                                // pending scale; no entry is touched
}
SDIPPS_scaling ;

typedef struct
{
    SDIPPS_push_choice push_choice ;
    SDIPPS_scaling scaling ;

    // Stop pushing as soon as the exiting row of U is empty past the
    // diagonal.  The entering column then lands at that position, and the
    // later frames only change by one scale factor.
    bool early_exit ;

    // For benchmarking the seed history: do every RwSOP in full, on every
    // live row of the seed, at every push.  The seed history (deferring the
    // rows whose RwSOP is a pure scaling until they are read) is then not
    // used.  The result is the same either way.
    bool eager_rwsop ;
}
SDIPPS_options ;

// The default options: prefer diagonal pushes, the regular backtrack, no
// early exit, the seed history.  Passing NULL for the options to
// SDIPPS_LU_col_rep gives these.
#define SDIPPS_DEFAULT_OPTIONS \
    { SDIPPS_PREFER_DIAGONAL, SDIPPS_SCALING_NONE, false, false }

//------------------------------------------------------------------------------
// statistics (optional output)
//------------------------------------------------------------------------------

// Counts for one call of SDIPPS_LU_col_rep.  Gathering them costs O(1) per
// operation counted, whether or not they are returned.

typedef struct
{
    // where the exiting column starts and the entering column lands
    int64_t start ;         // xi: starting position of the exiting column
    int64_t landing ;       // position of the entering column (n-1 unless
                            // the update exited early)

    // the pushes, and why each one was taken
    int64_t n_diagonal ;    // number of diagonal pushes
    int64_t n_column ;      // number of column pushes
    int64_t n_diagonal_forced ;  // sL sU = 0: only a diagonal push possible
    int64_t n_diagonal_chosen ;  // both possible; diagonal preferred
    int64_t n_column_forced ;    // pi = 0: only a column push possible
    int64_t n_column_chosen ;    // both possible; column preferred

    // backtracking steps that are a pure scaling (in a diagonal push, the
    // other half's support is 0), by half; eager and lazy scaling act here
    int64_t n_pure_L ;
    int64_t n_pure_U ;

    // SDIPPS_Reduce (eager and lazy scaling only): one call per run of
    // pure-scaling pushes, and one after an early exit
    int64_t n_reduce ;          // calls
    int64_t n_reduced ;         // calls whose ratio had a common factor
    int64_t reduce_bits_saved ; // total over the calls of bits(num) +
                                // bits(den) - bits(p) - bits(q)

    // the seeds mu^L, mu^U: the largest number of nonzeros and the largest
    // total bits (of the stored values, at live positions) they reach
    int64_t max_nnz_mu_L, max_nnz_mu_U ;
    int64_t max_bits_mu_L, max_bits_mu_U ;

    // seed row updates, both halves: by the full RwSOP formula (s_p T[i,k+1]
    // != 0), and by a pure scaling (rho over rho, done at once with
    // eager_rwsop, or several frames together by the history update)
    int64_t n_rwsop_rows ;
    int64_t n_scaled_rows ;
}
SDIPPS_stats ;

//------------------------------------------------------------------------------
// SDIPPS_LU_col_rep: replace one column of A in F
//------------------------------------------------------------------------------

// Same calling convention as SPEX_update_lu_colrep.  If F is not updatable it
// is converted first (SPEX_factorization_convert).  On output F is the
// updatable factorization of the new A.  vk is not modified; call
// SPEX_update_matrix_colrep (A, vk, col, option) to update A itself.
//
// Returns SPEX_SINGULAR if the new A is singular.  As with
// SPEX_update_lu_colrep, F is undefined if this function fails for any
// reason, including SPEX_SINGULAR; refactorize in that case.

SPEX_info SDIPPS_LU_col_rep
(
    SPEX_factorization F,       // updatable LU factorization of A; modified
    const SPEX_matrix vk,       // n-by-1 SPEX_DYNAMIC_CSC SPEX_MPZ matrix:
                                // the entering column, in A's row order,
                                // with vk->scale = A->scale and
                                // vk->v[0]->scale = 1
    int64_t col,                // column of A replaced by vk, 0 <= col < n
    const SDIPPS_options *sdipps_option,    // NULL: SDIPPS_DEFAULT_OPTIONS
    SDIPPS_stats *stats,        // NULL: not returned
    const SPEX_options option
) ;

//------------------------------------------------------------------------------
// SDIPPS_Reduce: reduce a pivot-over-pivot ratio
//------------------------------------------------------------------------------

// ratio = num / den in lowest terms.  This is where a different reduction
// method can be tried.  The result must be in lowest terms (canonical):
//  - eager scaling divides each entry by den(ratio) first, which is exact
//    only because num and den have no common factor (x * num / den is an
//    integer, so den divides x);
//  - lazy scaling makes it part of a SPEX_vector's pending scale, an mpq_t,
//    and GMP requires every mpq_t to be canonical.

SPEX_info SDIPPS_Reduce
(
    mpq_t ratio,                // output: num / den, in lowest terms
    const mpz_t num,
    const mpz_t den             // nonzero
) ;

#endif
