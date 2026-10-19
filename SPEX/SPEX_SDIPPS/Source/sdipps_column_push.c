//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_column_push: exchange columns k, k+1
//------------------------------------------------------------------------------

// Column push (APCPU, [Escobedo & Moreno-Centeno 2017, Fig. 4]): swap
// positions k and k+1 of Q only.  The exiting column moves to k+1; the
// exiting frame's row stays at k.  Requires sU != 0.
//
//  - new pivot:       rho~_k = sU
//  - new L column k:  Backtrack with the other support sU, and the extra
//                     entry L~(k+1,k) = pi / rho_k (zero when pi = 0)
//  - L seed:          restarts as minus old column k+1 of L (new frame k+1
//                     of L is old frame k+1 negated)
//  - new U row k:     the U seed itself (the exiting row, at frame k-1)
//  - U seed:          RwSOP with eps = -1
//  - frames k+2..n-1 are negated: recorded lazily, theta_k = -theta_{k-1}

#include "sdipps_internal.h"

SPEX_info sdipps_column_push
(
    sdipps_work *W,
    int64_t k
)
{
    SPEX_info info ;
    W->stats.n_column++ ;
    ASSERT (!sdipps_is_zero (W->sU)) ;

    sdipps_half *L = &(W->L), *U = &(W->U) ;
    int64_t n = W->n ;
    int theta = W->theta [k-1] ;

    //--------------------------------------------------------------------------
    // permutation, sign and pivot
    //--------------------------------------------------------------------------

    int64_t *Q = U->perm, *Qinv = U->perm_inv ;
    int64_t t = Q [k] ; Q [k] = Q [k+1] ; Q [k+1] = t ;
    Qinv [Q [k]] = k ; Qinv [Q [k+1]] = k+1 ;

    W->theta [k] = -theta ;
    SPEX_MPZ_SET (W->rho_new [k], W->sU) ;

    //--------------------------------------------------------------------------
    // lower half
    //--------------------------------------------------------------------------

    SPEX_CHECK (sdipps_apply_scale (L->v [k+1], W->rho [k+1])) ;

    // L~(k+1,k) = pi / rho_k, then backtrack into slot k
    SPEX_MPZ_DIVEXACT (W->f, W->pi, W->rho [k]) ;
    SPEX_CHECK (sdipps_backtrack (W, L, k, W->sU, L->perm [k+1], W->f)) ;

    // restart the L seed from old column k+1 of L: mu^L = -theta_{k-1} L(:,k+1)
    if (k < n-2)
    {
        SPEX_CHECK (sdipps_seed_load (L, L->v [k+1], -theta, k)) ;
    }

    //--------------------------------------------------------------------------
    // upper half
    //--------------------------------------------------------------------------

    // RwSOP below reads old row k+1 of U
    SPEX_CHECK (sdipps_apply_scale (U->v [k+1], W->rho [k+1])) ;

    // new row k of U is the seed: pivot sU, then mu^U[i] for i >= k+2
    sdipps_seed *seed = &(U->seed) ;
    int64_t need = 1 ;
    for (int64_t s = 0 ; s < seed->nlist ; s++)
    {
        int64_t i = seed->list [s] ;
        if (i >= k+2 && !sdipps_is_zero (seed->mu [i])) { need++ ; }
    }
    SPEX_vector V = U->v [k] ;                  // the free slot
    SPEX_CHECK (sdipps_start_vector (V, Q [k], W->sU, need, W->option)) ;
    for (int64_t s = 0 ; s < seed->nlist ; s++)
    {
        int64_t i = seed->list [s] ;
        if (i < k+2 || sdipps_is_zero (seed->mu [i])) { continue ; }
        SPEX_CHECK (sdipps_catch_up (W, U, i, k-1)) ;
        SPEX_CHECK (sdipps_append (V, Q [i], seed->mu [i])) ;
    }

    // move the U seed to frame k (eps = -1 carries the column swap)
    SPEX_CHECK (sdipps_rwsop (W, U, k, W->sU, -1)) ;

    // the backtracking ratio changes, which ends any run of pure scalings
    W->lambda_known = false ;
    return (SPEX_OK) ;
}
