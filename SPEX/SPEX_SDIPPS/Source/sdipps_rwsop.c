//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_rwsop: move the seed of one half to frame k
//------------------------------------------------------------------------------

// RwSOP (switch of originating pivot, [Escobedo & Moreno-Centeno 2017,
// Thm 8]) moves the exiting frame past old frame k+1.  For each half it is
// one extra IPGE iteration on the seed:
//
//      mu[i] <- eps (rho_{k+1} mu[i] - s_p T[i,k+1]) / rho_k,   i >= k+2
//
// with T[:,k+1] the old frame k+1 of this half and s_p this half's support.
// eps = -1 for the U half of a column push (the column swap negates).
//
// Only the rows with s_p T[i,k+1] != 0 need real work.  Every other row is a
// pure scaling and is left for sdipps_catch_up, except row k+2: it is next
// frame's support, which must be current.
//
// theta[k] must already be set, and old frame k+1 must have scale 1 if
// s_p != 0.

#include "sdipps_internal.h"

SPEX_info sdipps_rwsop
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t k,
    const mpz_t s_p,            // this half's support
    int eps                     // +1, or -1 for the U half of a column push
)
{
    SPEX_info info ;
    if (k+2 >= W->n) { return (SPEX_OK) ; }     // no rows past the support

    sdipps_seed *seed = &(H->seed) ;
    mpz_t *mu = seed->mu ;

    //--------------------------------------------------------------------------
    // rows where T[i,k+1] != 0 (none if s_p = 0)
    //--------------------------------------------------------------------------

    if (!sdipps_is_zero (s_p))
    {
        SPEX_vector T = H->v [k+1] ;
        for (int64_t p = 1 ; p < T->nz ; p++)
        {
            int64_t i = H->perm_inv [T->i [p]] ;
            SPEX_CHECK (sdipps_catch_up (W, H, i, k-1)) ;
            bool was_nonzero = !sdipps_is_zero (mu [i]) ;
            seed->bits -= sdipps_bits (mu [i]) ;
            SPEX_MPZ_MUL (mu [i], mu [i], W->rho [k+1]) ;
            SPEX_MPZ_SUBMUL (mu [i], s_p, T->x [p]) ;
            SPEX_MPZ_DIVEXACT (mu [i], mu [i], W->rho [k]) ;
            if (eps < 0)
            {
                SPEX_MPZ_NEG (mu [i], mu [i]) ;
            }
            seed->h [i] = k ;
            seed->bits += sdipps_bits (mu [i]) ;
            sdipps_seed_track (seed, i, was_nonzero) ;
        }
    }

    //--------------------------------------------------------------------------
    // row k+2 is the next support: bring it to frame k
    //--------------------------------------------------------------------------

    // If row k+2 was updated above, it is already at frame k.  Otherwise its
    // update at this frame is a pure scaling, and the catch-up through
    // frame k does it (with theta[k] theta[k-1] = eps).
    SPEX_CHECK (sdipps_catch_up (W, H, k+2, k)) ;
    return (SPEX_OK) ;
}
