//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_workspace: allocate and free the SDIPPS workspace
//------------------------------------------------------------------------------

// The workspace points into F for the factors, permutations and old pivots
// (the update is in place), and owns the seeds, the new pivots, the lazy sign
// array and the forward-solve vector.  W must be zeroed on input, so that
// sdipps_work_free is safe at any point.  If sdipps_work_alloc fails, the
// caller frees W.

#include "sdipps_internal.h"

// allocate one seed of size n
static SPEX_info sdipps_seed_alloc (sdipps_seed *seed, int64_t n)
{
    seed->mu     = spex_create_mpz_array (n) ;
    seed->h      = (int64_t *) SPEX_malloc (n * sizeof (int64_t)) ;
    seed->list   = (int64_t *) SPEX_malloc (n * sizeof (int64_t)) ;
    seed->listed = (bool *) SPEX_calloc (n, sizeof (bool)) ;
    seed->nlist  = 0 ;
    seed->nlive  = 0 ;
    if (!seed->mu || !seed->h || !seed->list || !seed->listed)
    {
        return (SPEX_OUT_OF_MEMORY) ;
    }
    return (SPEX_OK) ;
}

static void sdipps_seed_free (sdipps_seed *seed, int64_t n)
{
    spex_free_mpz_array (&(seed->mu), n) ;
    SPEX_FREE (seed->h) ;
    SPEX_FREE (seed->list) ;
    SPEX_FREE (seed->listed) ;
}

SPEX_info sdipps_work_alloc
(
    sdipps_work *W,             // workspace to allocate (all fields zeroed)
    SPEX_factorization F,
    int64_t col,
    const SDIPPS_options *sdipps_option,
    const SPEX_options option
)
{
    SPEX_info info ;
    int64_t n = F->L->n ;

    //--------------------------------------------------------------------------
    // options and the factorization
    //--------------------------------------------------------------------------

    SDIPPS_options defaults = SDIPPS_DEFAULT_OPTIONS ;
    W->opt = (sdipps_option != NULL) ? (*sdipps_option) : defaults ;
    W->option = option ;
    W->n = n ;
    W->col = col ;
    W->xi = F->Qinv_perm [col] ;

    W->L.v = F->L->v ;  W->L.perm = F->P_perm ;  W->L.perm_inv = F->Pinv_perm ;
    W->U.v = F->U->v ;  W->U.perm = F->Q_perm ;  W->U.perm_inv = F->Qinv_perm ;
    W->rho = F->rhos->x.mpz ;

    //--------------------------------------------------------------------------
    // arrays
    //--------------------------------------------------------------------------

    SPEX_CHECK (sdipps_seed_alloc (&(W->L.seed), n)) ;
    SPEX_CHECK (sdipps_seed_alloc (&(W->U.seed), n)) ;
    W->rho_new   = spex_create_mpz_array (n) ;
    W->x         = spex_create_mpz_array (n) ;
    W->eta       = (int64_t *) SPEX_malloc (n * sizeof (int64_t)) ;
    W->mark      = (int64_t *) SPEX_calloc (n, sizeof (int64_t)) ;
    W->theta_mem = (signed char *) SPEX_malloc ((n+1) * sizeof (signed char)) ;
    if (!W->rho_new || !W->x || !W->eta || !W->mark || !W->theta_mem)
    {
        return (SPEX_OUT_OF_MEMORY) ;
    }

    // theta[-1..n-1] = +1: no column pushes yet
    W->theta = W->theta_mem + 1 ;
    for (int64_t j = -1 ; j < n ; j++)
    {
        W->theta [j] = 1 ;
    }
    W->stamp = 0 ;

    //--------------------------------------------------------------------------
    // scalars
    //--------------------------------------------------------------------------

    SPEX_MPZ_INIT (W->one) ;
    SPEX_MPZ_INIT (W->sL) ;
    SPEX_MPZ_INIT (W->sU) ;
    SPEX_MPZ_INIT (W->pi) ;
    SPEX_MPZ_INIT (W->f) ;
    SPEX_MPZ_INIT (W->tmp) ;
    SPEX_MPQ_INIT (W->lambda) ;
    SPEX_MPZ_SET_UI (W->one, 1) ;
    W->lambda_known = false ;

    return (SPEX_OK) ;
}

void sdipps_work_free
(
    sdipps_work *W
)
{
    int64_t n = W->n ;
    sdipps_seed_free (&(W->L.seed), n) ;
    sdipps_seed_free (&(W->U.seed), n) ;
    spex_free_mpz_array (&(W->rho_new), n) ;
    spex_free_mpz_array (&(W->x), n) ;
    SPEX_FREE (W->eta) ;
    SPEX_FREE (W->mark) ;
    SPEX_FREE (W->theta_mem) ;
    W->theta = NULL ;
    SPEX_mpz_clear (W->one) ;
    SPEX_mpz_clear (W->sL) ;
    SPEX_mpz_clear (W->sU) ;
    SPEX_mpz_clear (W->pi) ;
    SPEX_mpz_clear (W->f) ;
    SPEX_mpz_clear (W->tmp) ;
    SPEX_mpq_clear (W->lambda) ;
}
