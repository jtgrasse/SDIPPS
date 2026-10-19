//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_seed: the seed vectors and their history
//------------------------------------------------------------------------------

// A seed mu carries one half of the frame being pushed (the exiting frame),
// by position.  Each push moves it one frame further by RwSOP, which is one
// extra IPGE iteration:
//
//      mu[i] <- eps (rho_{k+1} mu[i] - s_p T[i,k+1]) / rho_k
//
// When s_p T[i,k+1] = 0 this is a pure scaling, so it is not done right away.
// h[i] records the last frame mu[i] has been through, and sdipps_catch_up
// applies all the skipped scalings at once when mu[i] is next read: the
// rho's telescope, and the signs eps collect into theta_m theta_h.

#include "sdipps_internal.h"

//------------------------------------------------------------------------------
// sdipps_seed_load: mu <- sign * (off-diagonal entries of v), at history frame
//------------------------------------------------------------------------------

// The pivot of v is not loaded: its position is no longer live.

SPEX_info sdipps_seed_load
(
    sdipps_half *H,
    const SPEX_vector v,        // an old frame of H, with scale 1
    int sign,
    int64_t frame
)
{
    SPEX_info info ;
    sdipps_seed *seed = &(H->seed) ;

    // clear the previous seed
    for (int64_t t = 0 ; t < seed->nlist ; t++)
    {
        int64_t i = seed->list [t] ;
        SPEX_MPZ_SET_UI (seed->mu [i], 0) ;
        seed->listed [i] = false ;
    }
    seed->nlist = 0 ;
    seed->nlive = 0 ;
    seed->bits = 0 ;

    // load v
    for (int64_t p = 1 ; p < v->nz ; p++)
    {
        if (sdipps_is_zero (v->x [p])) { continue ; }    // explicit zero
        int64_t i = H->perm_inv [v->i [p]] ;
        SPEX_MPZ_SET (seed->mu [i], v->x [p]) ;
        if (sign < 0)
        {
            SPEX_MPZ_NEG (seed->mu [i], seed->mu [i]) ;
        }
        seed->h [i] = frame ;
        seed->listed [i] = true ;
        seed->list [seed->nlist++] = i ;
        seed->nlive++ ;
        seed->bits += sdipps_bits (seed->mu [i]) ;
    }
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_seed_retire: position j is no longer live
//------------------------------------------------------------------------------

SPEX_info sdipps_seed_retire
(
    sdipps_seed *seed,
    int64_t j
)
{
    SPEX_info info ;
    if (!sdipps_is_zero (seed->mu [j]))
    {
        seed->bits -= sdipps_bits (seed->mu [j]) ;
        SPEX_MPZ_SET_UI (seed->mu [j], 0) ;
        seed->nlive-- ;
    }
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_seed_note_max: keep the largest nnz and bits reached (statistics)
//------------------------------------------------------------------------------

void sdipps_seed_note_max
(
    sdipps_seed *seed
)
{
    if (seed->nlive > seed->max_nlive) seed->max_nlive = seed->nlive ;
    if (seed->bits > seed->max_bits) seed->max_bits = seed->bits ;
}

//------------------------------------------------------------------------------
// sdipps_seed_track: mu[i] may have changed between zero and nonzero
//------------------------------------------------------------------------------

void sdipps_seed_track
(
    sdipps_seed *seed,
    int64_t i,
    bool was_nonzero
)
{
    bool now_nonzero = !sdipps_is_zero (seed->mu [i]) ;
    if (now_nonzero && !was_nonzero)
    {
        seed->nlive++ ;
        if (!seed->listed [i])
        {
            seed->listed [i] = true ;
            seed->list [seed->nlist++] = i ;
        }
    }
    else if (was_nonzero && !now_nonzero)
    {
        seed->nlive-- ;
    }
}

//------------------------------------------------------------------------------
// sdipps_catch_up: history update, bring mu[i] up to frame m
//------------------------------------------------------------------------------

//      mu[i] <- theta_m theta_h rho_{m+1} mu[i] / rho_{h+1},  h = h[i]
//
// These are the pure scalings of frames h+1..m at once (pivot cancellation:
// rho_{h+2} ... rho_{m+1} over rho_{h+1} ... rho_m).  The signs: an RwSOP
// carries eps = -1 exactly at a column push, where theta changes sign.

SPEX_info sdipps_catch_up
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t i,
    int64_t m
)
{
    SPEX_info info ;
    sdipps_seed *seed = &(H->seed) ;
    if (sdipps_is_zero (seed->mu [i])) { return (SPEX_OK) ; }
    int64_t h = seed->h [i] ;
    if (h >= m) { return (SPEX_OK) ; }      // already current

    seed->bits -= sdipps_bits (seed->mu [i]) ;
    SPEX_MPZ_MUL (seed->mu [i], seed->mu [i], W->rho [m+1]) ;
    SPEX_MPZ_DIVEXACT (seed->mu [i], seed->mu [i], W->rho [h+1]) ;
    if (W->theta [m] != W->theta [h])
    {
        SPEX_MPZ_NEG (seed->mu [i], seed->mu [i]) ;
    }
    seed->bits += sdipps_bits (seed->mu [i]) ;
    seed->h [i] = m ;
    return (SPEX_OK) ;
}
