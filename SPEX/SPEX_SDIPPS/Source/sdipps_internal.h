//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_internal.h: internal definitions for SDIPPS
//------------------------------------------------------------------------------

// Not for user applications; use SDIPPS.h.
//
// Notation (0-based positions, as in SPEX).  rho_j is the OLD pivot of frame
// j (F->rhos, not changed until the update finishes) and rho~_j the new one;
// rho_{-1} = rho~_{-1} = 1.  Frame j is column j of L together with row j of
// U.  A "half" is either L (columns of L, indexed by row) or U (rows of U,
// indexed by column); every push routine serves both halves.
//
// Slots.  While the exiting column is pushed from position k to k+1:
//      slots < k   hold new frames (final),
//      slot  k     is free (its old frame has already been read),
//      slots > k   hold old frames.
// New frame k is written into slot k.  The vectors are SPEX's own
// (F->L->v[j], F->U->v[j]); nothing is copied out of F.
//
// Seeds.  mu^L and mu^U hold the frame being pushed (the exiting frame), as
// true values, by position.  h[i] is the last frame whose RwSOP mu[i] has been
// through.  RwSOP updates only the rows it must; the others are a pure
// pivot-over-pivot scaling and are brought up to date (sdipps_catch_up) when
// they are next read.  At frame k only positions > k are live.
//
// Lazy signs.  A column push negates every later frame.  Instead of doing
// that, theta[j] = (-1)^(number of column pushes at frames <= j) is recorded;
// stored old values are never negated.  At frame k every old value read
// carries the sign theta[k-1], which cancels in all but a few formulas.
//
// Pending scales.  As in SPEX, the true entries of a SPEX_vector v are
// x[p] * v->scale, x[p] / den(v->scale) is an integer, and x[0] * v->scale
// is the pivot.  A vector's scale is applied (sdipps_apply_scale) before its
// entries are read.  New vectors are written with scale 1, except when lazy
// scaling (SDIPPS_SCALING_LAZY) moves an old vector and multiplies its scale.

#ifndef SDIPPS_INTERNAL_H
#define SDIPPS_INTERNAL_H

#include "spex_util_internal.h"
#include "../Include/SDIPPS.h"

//------------------------------------------------------------------------------
// workspace
//------------------------------------------------------------------------------

typedef struct
{
    mpz_t *mu ;         // size n: the seed, by position
    int64_t *h ;        // size n: history of each mu[i]
    int64_t *list ;     // positions where mu has been nonzero since the last
                        // load, in no particular order (a superset of the
                        // nonzero pattern of mu)
    bool *listed ;      // size n: listed[i] is true if i is in list
    int64_t nlist ;     // length of list
    int64_t nlive ;     // number of nonzero mu[i] at live positions
    int64_t bits ;      // total bits of the mu[i] at live positions
    int64_t max_nlive ; // the largest nlive and bits reached (statistics)
    int64_t max_bits ;
}
sdipps_seed ;

typedef struct
{
    SPEX_vector *v ;    // F->L->v or F->U->v: slot j holds frame j
    int64_t *perm ;     // F->P_perm or F->Q_perm: position -> original index
    int64_t *perm_inv ; // F->Pinv_perm or F->Qinv_perm: original -> position
    sdipps_seed seed ;
}
sdipps_half ;

typedef struct
{
    int64_t n ;
    int64_t col ;               // the exiting (and entering) column of A
    int64_t xi ;                // starting position of the exiting column
    sdipps_half L, U ;
    mpz_t *rho ;                // old pivots, F->rhos->x.mpz
    mpz_t *rho_new ;            // size n: rho~_j for xi <= j <= landing
    mpz_t one ;                 // rho~_{-1} = 1
    signed char *theta ;        // theta[-1..n-1], lazy column-push signs
    signed char *theta_mem ;    // the allocation behind theta

    // the current frame
    mpz_t sL, sU ;              // supports: mu^L[k+1], mu^U[k+1]
    mpz_t pi ;                  // rho~_{k-1} rho_{k+1} + theta_{k-1} sL sU
    mpz_t f ;                   // L~(k+1,k) of a column push
    mpz_t tmp ;

    // eager and lazy scaling: lambda = rho~_{k-1} / rho_k, reduced, kept
    // along a run of diagonal pushes in which it does not change (see
    // sdipps_backtrack.c)
    bool lambda_known ;
    mpq_t lambda ;

    // EnterColumn: the forward solve x and its history eta, by position
    mpz_t *x ;
    int64_t *eta ;

    // marking rows when two patterns are merged
    int64_t *mark ;
    int64_t stamp ;

    SDIPPS_options opt ;
    SDIPPS_stats stats ;        // counts gathered during the update
    SPEX_options option ;
}
sdipps_work ;

// rho~_j: 1 for j = -1, the old pivot for j < xi (those frames do not change),
// and the new pivot otherwise.
static inline mpz_ptr sdipps_rho_new (sdipps_work *W, int64_t j)
{
    if (j < 0) return (W->one) ;
    if (j < W->xi) return (W->rho [j]) ;
    return (W->rho_new [j]) ;
}

// true if x == 0.  mpz_sgn is a GMP macro that cannot fail, so it is used
// directly rather than through SPEX's wrapper.
static inline bool sdipps_is_zero (const mpz_t x)
{
    return (mpz_sgn (x) == 0) ;
}

// number of bits of |x| (0 for x = 0), for the statistics; mpz_sizeinbase
// cannot fail and does not allocate
static inline int64_t sdipps_bits (const mpz_t x)
{
    return ((mpz_sgn (x) == 0) ? 0 : (int64_t) mpz_sizeinbase (x, 2)) ;
}

//------------------------------------------------------------------------------
// workspace (sdipps_workspace.c)
//------------------------------------------------------------------------------

SPEX_info sdipps_work_alloc
(
    sdipps_work *W,             // workspace to allocate (all fields zeroed)
    SPEX_factorization F,
    int64_t col,
    const SDIPPS_options *sdipps_option,
    const SPEX_options option
) ;

void sdipps_work_free
(
    sdipps_work *W
) ;

//------------------------------------------------------------------------------
// the pushes
//------------------------------------------------------------------------------

// exchange rows and columns k, k+1
SPEX_info sdipps_diagonal_push
(
    sdipps_work *W,
    int64_t k
) ;

// exchange columns k, k+1
SPEX_info sdipps_column_push
(
    sdipps_work *W,
    int64_t k
) ;

// new frame k of half H from old frame k+1 and the seed:
//  pivot rho~_k, the extra entry (if any), then for rows i >= k+2
//  (rho~_{k-1} T[i,k+1] + theta_{k-1} s_o mu[i]) / rho_k
SPEX_info sdipps_backtrack
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t k,
    const mpz_t s_o,            // the other half's support
    int64_t extra_index,        // original index of the extra entry
    const mpz_t extra_value     // NULL: no extra entry
) ;

// move the seed of half H from frame k-1 to frame k:
//  mu[i] = eps (rho_{k+1} mu[i] - s_p T[i,k+1]) / rho_k   for i >= k+2
SPEX_info sdipps_rwsop
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t k,
    const mpz_t s_p,            // this half's support
    int eps                     // +1, or -1 for the U half of a column push
) ;

//------------------------------------------------------------------------------
// the seeds (sdipps_seed.c)
//------------------------------------------------------------------------------

// mu <- sign * (off-diagonal entries of v), all at history `frame`
SPEX_info sdipps_seed_load
(
    sdipps_half *H,
    const SPEX_vector v,        // an old frame of H, with scale 1
    int sign,
    int64_t frame
) ;

// position j is no longer live
SPEX_info sdipps_seed_retire
(
    sdipps_seed *seed,
    int64_t j
) ;

// record that mu[i] may have changed between zero and nonzero
void sdipps_seed_track
(
    sdipps_seed *seed,
    int64_t i,
    bool was_nonzero
) ;

// record the seed's current nnz and bits if they are its largest so far
void sdipps_seed_note_max
(
    sdipps_seed *seed
) ;

// history update: bring mu[i] up to frame m
//  mu[i] = theta_m theta_h rho_{m+1} mu[i] / rho_{h+1},  h = h[i]
SPEX_info sdipps_catch_up
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t i,
    int64_t m
) ;

// SDIPPS_Reduce, counted in W->stats
SPEX_info sdipps_reduce
(
    sdipps_work *W,
    mpq_t ratio,
    const mpz_t num,
    const mpz_t den
) ;

//------------------------------------------------------------------------------
// the swap
//------------------------------------------------------------------------------

// one SLIP LU iteration at position `landing` with the entering column a
SPEX_info sdipps_enter_column
(
    sdipps_work *W,
    const SPEX_vector a,
    int64_t landing
) ;

// after an early exit: frames past `landing` scale by rho~_landing/rho_landing
SPEX_info sdipps_rescale_untouched
(
    sdipps_work *W,
    int64_t landing
) ;

//------------------------------------------------------------------------------
// SPEX_vector utilities (sdipps_vector.c)
//------------------------------------------------------------------------------

// apply v's pending scale to its entries, so that x[p] are true values
SPEX_info sdipps_apply_scale
(
    SPEX_vector v,
    const mpz_t pivot           // the true pivot of v's frame
) ;

// set x[0] so that x[0] * scale is the given pivot
SPEX_info sdipps_set_pivot
(
    SPEX_vector v,
    const mpz_t pivot
) ;

// v <- (pivot_index : pivot), scale 1, with room for `need` entries
SPEX_info sdipps_start_vector
(
    SPEX_vector v,
    int64_t pivot_index,
    const mpz_t pivot,
    int64_t need,
    const SPEX_options option
) ;

// make room for `need` entries
SPEX_info sdipps_reserve
(
    SPEX_vector v,
    int64_t need,
    const SPEX_options option
) ;

// append (index, value) to v if value != 0; room must have been reserved
SPEX_info sdipps_append
(
    SPEX_vector v,
    int64_t index,
    const mpz_t value
) ;

// remove the off-diagonal entry of v at original index `index`, if present
SPEX_info sdipps_remove_entry
(
    SPEX_vector v,
    int64_t index
) ;

#endif
