//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_diagonal_push: exchange rows and columns k, k+1
//------------------------------------------------------------------------------

// Diagonal push (APDPU, [Escobedo & Moreno-Centeno 2017, Fig. 5]): swap
// positions k and k+1 of both P and Q.  The exiting frame moves to k+1.
//
//  - new pivot:    rho~_k = pi / rho_k
//  - new frame k:  Backtrack of each half (old frame k+1 and the seed)
//  - the seeds:    RwSOP of each half, eps = +1
//  - frames k+2..n-1 do not change, and no sign changes: theta_k = theta_{k-1}
//
// For the L half the other half's support is sU, and the new frame gets the
// extra entry L~(k+1,k) = sU (the exiting row is now at position k+1).  For
// the U half the other support is sL; its entry U~(k,k+1) belongs to the
// exiting column and is not stored.

#include "sdipps_internal.h"

// Backtrack and RwSOP for one half.
static SPEX_info sdipps_diagonal_push_half
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t k,
    const mpz_t s_p,            // this half's support
    const mpz_t s_o,            // the other half's support
    int64_t extra_index         // original index of the extra entry, or -1
)
{
    SPEX_info info ;
    SPEX_vector T = H->v [k+1] ;                // old frame k+1
    bool pure_scaling = sdipps_is_zero (s_o) ;

    // T's entries are read unless it is only rescaled lazily and RwSOP does
    // not need it either (s_p = 0); in that case its pending scale is kept.
    bool keep_scale = pure_scaling &&
        (W->opt.scaling == SDIPPS_SCALING_LAZY) &&
        sdipps_is_zero (s_p) ;
    if (!keep_scale)
    {
        SPEX_CHECK (sdipps_apply_scale (T, W->rho [k+1])) ;
    }

    // the extra entry of the L half is L~(k+1,k) = s_o
    mpz_srcptr extra = (extra_index >= 0) ? s_o : NULL ;

    if (pure_scaling)
    {
        if (H == &(W->L)) W->stats.n_pure_L++ ; else W->stats.n_pure_U++ ;
        // Backtrack rescales T in place and moves it, so RwSOP reads T first
        SPEX_CHECK (sdipps_rwsop (W, H, k, s_p, +1)) ;
        SPEX_CHECK (sdipps_backtrack (W, H, k, s_o, extra_index, extra)) ;
    }
    else
    {
        // Backtrack reads the seed before RwSOP moves it to frame k
        SPEX_CHECK (sdipps_backtrack (W, H, k, s_o, extra_index, extra)) ;
        SPEX_CHECK (sdipps_rwsop (W, H, k, s_p, +1)) ;
    }
    return (SPEX_OK) ;
}

SPEX_info sdipps_diagonal_push
(
    sdipps_work *W,
    int64_t k
)
{
    SPEX_info info ;
    W->stats.n_diagonal++ ;

    //--------------------------------------------------------------------------
    // permutations, sign and pivot
    //--------------------------------------------------------------------------

    int64_t *P = W->L.perm, *Pinv = W->L.perm_inv ;
    int64_t *Q = W->U.perm, *Qinv = W->U.perm_inv ;
    int64_t t ;
    t = P [k] ; P [k] = P [k+1] ; P [k+1] = t ;
    Pinv [P [k]] = k ; Pinv [P [k+1]] = k+1 ;
    t = Q [k] ; Q [k] = Q [k+1] ; Q [k+1] = t ;
    Qinv [Q [k]] = k ; Qinv [Q [k+1]] = k+1 ;

    W->theta [k] = W->theta [k-1] ;

    // rho~_k = pi / rho_k
    SPEX_MPZ_DIVEXACT (W->rho_new [k], W->pi, W->rho [k]) ;

    //--------------------------------------------------------------------------
    // the two halves
    //--------------------------------------------------------------------------

    // L: this support sL, other support sU, extra entry at the exiting row
    SPEX_CHECK (sdipps_diagonal_push_half (W, &(W->L), k, W->sL, W->sU,
        P [k+1])) ;

    // U: this support sU, other support sL, no extra entry
    SPEX_CHECK (sdipps_diagonal_push_half (W, &(W->U), k, W->sU, W->sL, -1)) ;

    // A push with sL sU != 0 changes the backtracking ratio, which ends a
    // run of pure-scaling pushes.
    if (!sdipps_is_zero (W->sL) && !sdipps_is_zero (W->sU))
    {
        W->lambda_known = false ;
    }
    return (SPEX_OK) ;
}
