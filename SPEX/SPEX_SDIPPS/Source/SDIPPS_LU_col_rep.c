//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/SDIPPS_LU_col_rep: SDIPPS column replacement
//------------------------------------------------------------------------------

// Purpose: replace column `col` of A by vk in the updatable SPEX LU
// factorization F, in place.  Outline:
//
//  1. Frames 0..xi-1 do not change, except that the exiting column leaves
//     their rows of U.  The exiting frame (slot xi) becomes the two seeds.
//  2. For k = xi, ..., n-2: push the exiting column from position k to k+1
//     by a diagonal or a column push.  With early_exit, stop at the first k
//     where the exiting row of U (mu^U) is empty past the diagonal.
//  3. The entering column lands at the final position by one REF forward
//     solve; with an early exit, the frames past it are rescaled.
//
// See SDIPPS.h for the interface and sdipps_internal.h for the conventions.

#define SPEX_FREE_ALL                       \
{                                           \
    sdipps_work_free (&W) ;                 \
}

#include "sdipps_internal.h"

// Choose the push at frame k, and count why it was chosen.  A diagonal push
// is possible iff pi != 0, and a column push iff sU != 0.  If sL * sU == 0
// then pi = rho~_{k-1} rho_{k+1} is nonzero, so a diagonal push is always
// possible there.
static char sdipps_choose_push (sdipps_work *W)
{
    bool zero_support = sdipps_is_zero (W->sL) || sdipps_is_zero (W->sU) ;
    if (zero_support)
    {
        W->stats.n_diagonal_forced++ ;  // only the diagonal push is possible
        return ('D') ;
    }
    if (sdipps_is_zero (W->pi))
    {
        W->stats.n_column_forced++ ;    // only the column push is possible
        return ('C') ;
    }
    // both are possible
    if (W->opt.push_choice == SDIPPS_PREFER_COLUMN)
    {
        W->stats.n_column_chosen++ ;
        return ('C') ;
    }
    W->stats.n_diagonal_chosen++ ;
    return ('D') ;
}

SPEX_info SDIPPS_LU_col_rep
(
    SPEX_factorization F,       // updatable LU factorization of A; modified
    const SPEX_matrix vk,       // n-by-1 dynamic CSC MPZ, the entering column
    int64_t col,                // column of A replaced by vk, 0 <= col < n
    const SDIPPS_options *sdipps_option,    // NULL: SDIPPS_DEFAULT_OPTIONS
    SDIPPS_stats *stats,        // NULL: not returned
    const SPEX_options option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    sdipps_work W ;
    memset (&W, 0, sizeof (sdipps_work)) ;   // so SPEX_FREE_ALL is safe

    if (!spex_initialized ( )) { return (SPEX_PANIC) ; }

    SPEX_REQUIRE (vk, SPEX_DYNAMIC_CSC, SPEX_MPZ) ;
    if (F == NULL || F->kind != SPEX_LU_FACTORIZATION || F->L == NULL)
    {
        return (SPEX_INCORRECT_INPUT) ;
    }
    int64_t n = F->L->n ;
    if (col < 0 || col >= n || vk->m != n || vk->n != 1)
    {
        return (SPEX_INCORRECT_INPUT) ;
    }
    SPEX_vector a = vk->v [0] ;
    int r ;
    SPEX_MPQ_CMP_UI (&r, a->scale, 1, 1) ;
    if (r != 0)
    {
        return (SPEX_INCORRECT_INPUT) ;     // vk->v[0]->scale must be 1
    }
    if (a->nz == 0)
    {
        return (SPEX_SINGULAR) ;            // a zero column; F is unchanged
    }

    // make sure F is updatable (L and U as SPEX_DYNAMIC_CSC, Qinv_perm built)
    SPEX_CHECK (SPEX_factorization_convert (F, true, option)) ;

    //--------------------------------------------------------------------------
    // allocate the workspace
    //--------------------------------------------------------------------------

    SPEX_CHECK (sdipps_work_alloc (&W, F, col, sdipps_option, option)) ;
    int64_t xi = W.xi ;

    //--------------------------------------------------------------------------
    // frames 0..xi-1 are final, except for the exiting column's entries
    //--------------------------------------------------------------------------

    // Remove the exiting column from rows 0..xi-1 of U.  The entering
    // column's entries are put in by sdipps_enter_column.
    for (int64_t j = 0 ; j < xi ; j++)
    {
        SPEX_CHECK (sdipps_remove_entry (W.U.v [j], col)) ;
    }

    // The exiting frame is the starting seed, at history xi-1.
    SPEX_CHECK (sdipps_apply_scale (W.L.v [xi], W.rho [xi])) ;
    SPEX_CHECK (sdipps_apply_scale (W.U.v [xi], W.rho [xi])) ;
    SPEX_CHECK (sdipps_seed_load (&W.L, W.L.v [xi], +1, xi-1)) ;
    SPEX_CHECK (sdipps_seed_load (&W.U, W.U.v [xi], +1, xi-1)) ;
    sdipps_seed_note_max (&W.L.seed) ;
    sdipps_seed_note_max (&W.U.seed) ;

    //--------------------------------------------------------------------------
    // push the exiting column from xi toward n-1
    //--------------------------------------------------------------------------

    int64_t landing = n-1 ;
    for (int64_t k = xi ; k < n-1 ; k++)
    {
        // position k (last frame's support) is no longer live
        if (k > xi)
        {
            SPEX_CHECK (sdipps_seed_retire (&W.L.seed, k)) ;
            SPEX_CHECK (sdipps_seed_retire (&W.U.seed, k)) ;
        }

        // early exit: the exiting row of U is empty past the diagonal
        if (W.opt.early_exit && W.U.seed.nlive == 0)
        {
            landing = k ;
            break ;
        }

        // the supports (current to frame k-1, see sdipps_rwsop) and
        // pi = rho~_{k-1} rho_{k+1} + theta_{k-1} sL sU
        SPEX_MPZ_SET (W.sL, W.L.seed.mu [k+1]) ;
        SPEX_MPZ_SET (W.sU, W.U.seed.mu [k+1]) ;
        SPEX_MPZ_MUL (W.pi, sdipps_rho_new (&W, k-1), W.rho [k+1]) ;
        if (W.theta [k-1] > 0)
        {
            SPEX_MPZ_ADDMUL (W.pi, W.sL, W.sU) ;
        }
        else
        {
            SPEX_MPZ_SUBMUL (W.pi, W.sL, W.sU) ;
        }

        if (sdipps_choose_push (&W) == 'D')
        {
            SPEX_CHECK (sdipps_diagonal_push (&W, k)) ;
        }
        else
        {
            SPEX_CHECK (sdipps_column_push (&W, k)) ;
        }
        sdipps_seed_note_max (&W.L.seed) ;
        sdipps_seed_note_max (&W.U.seed) ;
    }

    //--------------------------------------------------------------------------
    // swap: the entering column lands at position `landing`
    //--------------------------------------------------------------------------

    ASSERT (F->Q_perm [landing] == col) ;
    SPEX_CHECK (sdipps_enter_column (&W, a, landing)) ;
    if (landing < n-1)
    {
        SPEX_CHECK (sdipps_rescale_untouched (&W, landing)) ;
    }

    // the new pivots replace the old ones at positions xi..landing
    for (int64_t j = xi ; j <= landing ; j++)
    {
        SPEX_MPZ_SWAP (W.rho [j], W.rho_new [j]) ;
    }

    if (stats != NULL)
    {
        W.stats.start = xi ;
        W.stats.landing = landing ;
        W.stats.max_nnz_mu_L = W.L.seed.max_nlive ;
        W.stats.max_nnz_mu_U = W.U.seed.max_nlive ;
        W.stats.max_bits_mu_L = W.L.seed.max_bits ;
        W.stats.max_bits_mu_U = W.U.seed.max_bits ;
        (*stats) = W.stats ;
    }

    SPEX_FREE_ALL ;
    return (SPEX_OK) ;
}
