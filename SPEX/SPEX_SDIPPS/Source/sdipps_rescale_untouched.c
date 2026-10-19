//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_rescale_untouched: frames past an early exit
//------------------------------------------------------------------------------

// After an early exit at `landing`, the exiting row of U is empty past the
// diagonal, so every later frame is the old one times
//
//      r = rho~_landing / rho_landing      [Chen 2022, Thm 3.1.4]
//
// (the lazy signs cancel: both the old frames and rho_landing carry
// theta_{landing-1}).  As in backtracking, the scaling option decides how:
//      NONE   each entry times rho~_landing, divided by rho_landing;
//      EAGER  r reduced to p / q by SDIPPS_Reduce, each entry divided by q,
//             then times p (exact, since gcd (p, q) = 1);
//      LAZY   r reduced to p / q and multiplied into each vector's pending
//             scale; only the pivots are computed.
//
// The pivots of these frames are rewritten directly in F->rhos.

#include "sdipps_internal.h"

// x = x * num / den, exactly: multiplying first, or (divide_first) dividing
// first, which is exact only when num / den is in lowest terms
static SPEX_info sdipps_scale_entry (mpz_t x, const mpz_t num, const mpz_t den,
    bool divide_first)
{
    SPEX_info info ;
    if (divide_first)
    {
        SPEX_MPZ_DIVEXACT (x, x, den) ;
        SPEX_MPZ_MUL (x, x, num) ;
    }
    else
    {
        SPEX_MPZ_MUL (x, x, num) ;
        SPEX_MPZ_DIVEXACT (x, x, den) ;
    }
    return (SPEX_OK) ;
}

#undef  SPEX_FREE_ALL
#define SPEX_FREE_ALL                       \
{                                           \
    SPEX_mpq_clear (r) ;                    \
}

SPEX_info sdipps_rescale_untouched
(
    sdipps_work *W,
    int64_t landing
)
{
    SPEX_info info ;
    mpq_t r ;
    SPEX_mpq_set_null (r) ;
    SPEX_MPQ_INIT (r) ;

    mpz_ptr rho_new = W->rho_new [landing] ;    // rho~_landing
    mpz_ptr rho_old = W->rho [landing] ;        // rho_landing

    if (W->opt.scaling == SDIPPS_SCALING_LAZY)
    {
        SPEX_CHECK (sdipps_reduce (W, r, rho_new, rho_old)) ;
        for (int64_t j = landing+1 ; j < W->n ; j++)
        {
            SPEX_vector Lj = W->L.v [j], Uj = W->U.v [j] ;
            // rho_j = rho_j * r
            SPEX_MPZ_DIVEXACT (W->rho [j], W->rho [j], SPEX_MPQ_DEN (r)) ;
            SPEX_MPZ_MUL (W->rho [j], W->rho [j], SPEX_MPQ_NUM (r)) ;
            SPEX_MPQ_MUL (Lj->scale, Lj->scale, r) ;
            SPEX_MPQ_MUL (Uj->scale, Uj->scale, r) ;
            SPEX_CHECK (sdipps_set_pivot (Lj, W->rho [j])) ;
            SPEX_CHECK (sdipps_set_pivot (Uj, W->rho [j])) ;
        }
    }
    else
    {
        // every entry (pivots in rhos included) times num / den: for NONE,
        // num / den = rho~_landing / rho_landing, multiplying first; for
        // EAGER, the reduced p / q, dividing first (exact: gcd (p, q) = 1)
        bool eager = (W->opt.scaling == SDIPPS_SCALING_EAGER) ;
        mpz_ptr num = rho_new, den = rho_old ;
        if (eager)
        {
            SPEX_CHECK (sdipps_reduce (W, r, rho_new, rho_old)) ;
            num = SPEX_MPQ_NUM (r) ;
            den = SPEX_MPQ_DEN (r) ;
        }
        for (int64_t j = landing+1 ; j < W->n ; j++)
        {
            SPEX_vector Lj = W->L.v [j], Uj = W->U.v [j] ;
            SPEX_CHECK (sdipps_apply_scale (Lj, W->rho [j])) ;
            SPEX_CHECK (sdipps_apply_scale (Uj, W->rho [j])) ;
            SPEX_CHECK (sdipps_scale_entry (W->rho [j], num, den, eager)) ;
            for (int64_t p = 1 ; p < Lj->nz ; p++)
            {
                SPEX_CHECK (sdipps_scale_entry (Lj->x [p], num, den, eager)) ;
            }
            for (int64_t p = 1 ; p < Uj->nz ; p++)
            {
                SPEX_CHECK (sdipps_scale_entry (Uj->x [p], num, den, eager)) ;
            }
            SPEX_MPZ_SET (Lj->x [0], W->rho [j]) ;
            SPEX_MPZ_SET (Uj->x [0], W->rho [j]) ;
        }
    }

    SPEX_FREE_ALL ;
    return (SPEX_OK) ;
}
