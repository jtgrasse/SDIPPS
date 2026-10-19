//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_enter_column: put the entering column in
//------------------------------------------------------------------------------

// The exiting column has been pushed to position `landing`.  Frames
// 0..landing-1 do not involve it, so the entering column a replaces it with
// one iteration of left-looking SLIP LU [Lourenco et al. 2019, Alg. 1]: a REF
// forward solve of P a through columns 0..landing-1 of the new L.  Then
//
//      x[landing]       is the new pivot rho~_landing (0 iff A is singular),
//      x[landing+1:n-1] is column `landing` of the new L,
//      x[0:landing-1]   is column `landing` of the new U above the diagonal,
//
// and row `landing` of the new U is the pivot alone (trivially when landing
// is n-1; by the early-exit condition otherwise).
//
// eta[i] is the history of x[i]: the last pivot rho~_{eta[i]} it was
// multiplied by.  Bring(i, m) catches x[i] up to iteration m.

#include "sdipps_internal.h"

// x[i] <- rho~_m x[i] / rho~_{eta[i]}, eta[i] <- m
static SPEX_info sdipps_bring (sdipps_work *W, int64_t i, int64_t m)
{
    SPEX_info info ;
    if (W->eta [i] < m)
    {
        if (!sdipps_is_zero (W->x [i]))
        {
            SPEX_MPZ_MUL (W->x [i], W->x [i], sdipps_rho_new (W, m)) ;
            SPEX_MPZ_DIVEXACT (W->x [i], W->x [i],
                sdipps_rho_new (W, W->eta [i])) ;
        }
        W->eta [i] = m ;
    }
    return (SPEX_OK) ;
}

SPEX_info sdipps_enter_column
(
    sdipps_work *W,
    const SPEX_vector a,                // the entering column, A's row order
    int64_t landing
)
{
    SPEX_info info ;
    int64_t n = W->n ;
    int64_t *P = W->L.perm, *Pinv = W->L.perm_inv ;
    mpz_t *x = W->x ;

    //--------------------------------------------------------------------------
    // x = P a, history -1 (rho~_{-1} = 1)
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; i < n ; i++)
    {
        SPEX_MPZ_SET_UI (x [i], 0) ;
        W->eta [i] = -1 ;
    }
    for (int64_t p = 0 ; p < a->nz ; p++)
    {
        SPEX_MPZ_SET (x [Pinv [a->i [p]]], a->x [p]) ;
    }

    //--------------------------------------------------------------------------
    // forward solve through columns 0..landing-1 of the new L
    //--------------------------------------------------------------------------

    for (int64_t j = 0 ; j < landing ; j++)
    {
        if (sdipps_is_zero (x [j])) { continue ; }
        SPEX_CHECK (sdipps_bring (W, j, j-1)) ;     // x[j] is final

        SPEX_vector Lj = W->L.v [j] ;
        SPEX_CHECK (sdipps_apply_scale (Lj, sdipps_rho_new (W, j))) ;
        for (int64_t p = 1 ; p < Lj->nz ; p++)
        {
            if (sdipps_is_zero (Lj->x [p])) { continue ; }  // explicit zero
            int64_t i = Pinv [Lj->i [p]] ;
            SPEX_CHECK (sdipps_bring (W, i, j-1)) ;
            // x[i] = (rho~_j x[i] - L(i,j) x[j]) / rho~_{j-1}
            SPEX_MPZ_MUL (x [i], x [i], sdipps_rho_new (W, j)) ;
            SPEX_MPZ_SUBMUL (x [i], Lj->x [p], x [j]) ;
            SPEX_MPZ_DIVEXACT (x [i], x [i], sdipps_rho_new (W, j-1)) ;
            W->eta [i] = j ;
        }
    }

    // bring the rest to iteration landing-1, and count the new L column
    int64_t count = 0 ;
    for (int64_t i = landing ; i < n ; i++)
    {
        SPEX_CHECK (sdipps_bring (W, i, landing-1)) ;
        if (i > landing && !sdipps_is_zero (x [i])) { count++ ; }
    }

    if (sdipps_is_zero (x [landing]))
    {
        return (SPEX_SINGULAR) ;
    }

    //--------------------------------------------------------------------------
    // frame `landing`
    //--------------------------------------------------------------------------

    SPEX_MPZ_SET (W->rho_new [landing], x [landing]) ;

    SPEX_vector V = W->L.v [landing] ;
    SPEX_CHECK (sdipps_start_vector (V, P [landing], x [landing], 1 + count,
        W->option)) ;
    for (int64_t i = landing+1 ; i < n ; i++)
    {
        SPEX_CHECK (sdipps_append (V, P [i], x [i])) ;
    }

    ASSERT (W->U.perm [landing] == W->col) ;
    SPEX_CHECK (sdipps_start_vector (W->U.v [landing], W->col, x [landing], 1,
        W->option)) ;

    //--------------------------------------------------------------------------
    // column `landing` of U above the diagonal: append to rows 0..landing-1
    //--------------------------------------------------------------------------

    for (int64_t j = 0 ; j < landing ; j++)
    {
        if (sdipps_is_zero (x [j])) { continue ; }
        SPEX_vector Uj = W->U.v [j] ;
        SPEX_CHECK (sdipps_apply_scale (Uj, sdipps_rho_new (W, j))) ;
        SPEX_CHECK (sdipps_reserve (Uj, Uj->nz + 1, W->option)) ;
        SPEX_CHECK (sdipps_append (Uj, W->col, x [j])) ;
    }
    return (SPEX_OK) ;
}
