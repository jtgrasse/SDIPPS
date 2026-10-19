//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_backtrack: new frame k of one half
//------------------------------------------------------------------------------

// Backtracking ([Escobedo & Moreno-Centeno 2017, Thm 7]) builds new frame k
// of half H from old frame k+1 (T) and the seed (mu, the exiting frame):
//
//      pivot rho~_k,  the extra entry (if any),  and for rows i >= k+2
//      (rho~_{k-1} T[i,k+1] + theta_{k-1} s_o mu[i]) / rho_k
//
// where s_o is the OTHER half's support.  Two cases:
//
//  - s_o != 0 (seed merge): the new frame is written into the free slot k,
//    over the rows of T and of the seed.  T must have scale 1.
//
//  - s_o == 0 (pure scaling): the new frame is T times
//    lambda = rho~_{k-1} / rho_k.  T itself is moved to slot k (the free
//    vector goes to slot k+1).  How the entries are scaled is the scaling
//    option:
//      NONE   each entry times rho~_{k-1}, divided by rho_k;
//      EAGER  each entry divided by q, then times p, with p / q = lambda
//             in lowest terms (q divides every entry: x p / q is an
//             integer and gcd (p, q) = 1) [Chen 2022, Sec. 5.4];
//      LAZY   p / q is multiplied into T's pending scale; no entry is
//             touched.
//    For NONE and EAGER, T must have scale 1.  Since T is changed, the
//    caller runs this half's RwSOP (which reads T) first.
//
// lambda is constant along a run of diagonal pushes with sL sU = 0: there
// rho~_k = lambda rho_{k+1}, so rho~_k / rho_{k+1} = lambda again.  For EAGER
// and LAZY it is therefore reduced (SDIPPS_Reduce) once per run and kept in
// W->lambda; the push routines forget it when a run ends.

#include "sdipps_internal.h"

SPEX_info sdipps_backtrack
(
    sdipps_work *W,
    sdipps_half *H,
    int64_t k,
    const mpz_t s_o,            // the other half's support
    int64_t extra_index,        // original index of the extra entry
    const mpz_t extra_value     // NULL: no extra entry
)
{
    SPEX_info info ;
    SPEX_vector T = H->v [k+1] ;                // old frame k+1
    mpz_ptr pivot = W->rho_new [k] ;            // rho~_k
    mpz_ptr rho_new_prev = sdipps_rho_new (W, k-1) ;

    //--------------------------------------------------------------------------
    // pure scaling: new frame k = lambda * old frame k+1
    //--------------------------------------------------------------------------

    if (sdipps_is_zero (s_o))
    {
        // (the extra entry, if any, is s_o = 0 here)
        SDIPPS_scaling scaling = W->opt.scaling ;
        if (scaling != SDIPPS_SCALING_NONE && !W->lambda_known)
        {
            // the first pure-scaling push of a run: reduce lambda
            SPEX_CHECK (sdipps_reduce (W, W->lambda, rho_new_prev,
                W->rho [k])) ;
            W->lambda_known = true ;
        }

        if (scaling == SDIPPS_SCALING_LAZY)
        {
            SPEX_MPQ_MUL (T->scale, T->scale, W->lambda) ;
            SPEX_CHECK (sdipps_set_pivot (T, pivot)) ;
        }
        else if (scaling == SDIPPS_SCALING_EAGER)
        {
            for (int64_t p = 1 ; p < T->nz ; p++)
            {
                SPEX_MPZ_DIVEXACT (T->x [p], T->x [p],
                    SPEX_MPQ_DEN (W->lambda)) ;
                SPEX_MPZ_MUL (T->x [p], T->x [p], SPEX_MPQ_NUM (W->lambda)) ;
            }
            SPEX_MPZ_SET (T->x [0], pivot) ;
        }
        else
        {
            for (int64_t p = 1 ; p < T->nz ; p++)
            {
                SPEX_MPZ_MUL (T->x [p], T->x [p], rho_new_prev) ;
                SPEX_MPZ_DIVEXACT (T->x [p], T->x [p], W->rho [k]) ;
            }
            SPEX_MPZ_SET (T->x [0], pivot) ;
        }

        // T's pivot index is H->perm[k] (positions k, k+1 were swapped)
        ASSERT (T->i [0] == H->perm [k]) ;
        H->v [k+1] = H->v [k] ;
        H->v [k] = T ;
        return (SPEX_OK) ;
    }

    //--------------------------------------------------------------------------
    // seed merge: rows of T and rows of the seed, written into slot k
    //--------------------------------------------------------------------------

    sdipps_seed *seed = &(H->seed) ;
    int theta = W->theta [k-1] ;
    int64_t stamp = ++(W->stamp) ;

    // mark the rows of T, and count the entries of the new frame
    int64_t need = 2 ;                          // the pivot and the extra
    for (int64_t p = 1 ; p < T->nz ; p++)
    {
        W->mark [H->perm_inv [T->i [p]]] = stamp ;
        need++ ;
    }
    for (int64_t t = 0 ; t < seed->nlist ; t++)
    {
        int64_t i = seed->list [t] ;
        if (i >= k+2 && W->mark [i] != stamp && !sdipps_is_zero (seed->mu [i]))
        {
            need++ ;
        }
    }

    SPEX_vector V = H->v [k] ;                  // the free slot
    SPEX_CHECK (sdipps_start_vector (V, H->perm [k], pivot, need, W->option)) ;
    if (extra_value != NULL)
    {
        SPEX_CHECK (sdipps_append (V, extra_index, extra_value)) ;
    }

    // rows of T: (rho~_{k-1} T[i,k+1] + theta_{k-1} s_o mu[i]) / rho_k
    for (int64_t p = 1 ; p < T->nz ; p++)
    {
        int64_t i = H->perm_inv [T->i [p]] ;
        SPEX_CHECK (sdipps_catch_up (W, H, i, k-1)) ;
        SPEX_MPZ_MUL (W->tmp, rho_new_prev, T->x [p]) ;
        if (theta > 0)
        {
            SPEX_MPZ_ADDMUL (W->tmp, s_o, seed->mu [i]) ;
        }
        else
        {
            SPEX_MPZ_SUBMUL (W->tmp, s_o, seed->mu [i]) ;
        }
        SPEX_MPZ_DIVEXACT (W->tmp, W->tmp, W->rho [k]) ;
        SPEX_CHECK (sdipps_append (V, T->i [p], W->tmp)) ;
    }

    // rows of the seed only (T[i,k+1] = 0): theta_{k-1} s_o mu[i] / rho_k
    for (int64_t t = 0 ; t < seed->nlist ; t++)
    {
        int64_t i = seed->list [t] ;
        if (i < k+2 || W->mark [i] == stamp || sdipps_is_zero (seed->mu [i]))
        {
            continue ;
        }
        SPEX_CHECK (sdipps_catch_up (W, H, i, k-1)) ;
        SPEX_MPZ_MUL (W->tmp, s_o, seed->mu [i]) ;
        if (theta < 0)
        {
            SPEX_MPZ_NEG (W->tmp, W->tmp) ;
        }
        SPEX_MPZ_DIVEXACT (W->tmp, W->tmp, W->rho [k]) ;
        SPEX_CHECK (sdipps_append (V, H->perm [i], W->tmp)) ;
    }
    return (SPEX_OK) ;
}
