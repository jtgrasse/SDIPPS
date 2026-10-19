//------------------------------------------------------------------------------
// SPEX_SDIPPS/Test/sdipps_lp_refactor.c: refactorization baseline
//------------------------------------------------------------------------------

/* The baseline for sdipps_lp_test: it walks through the same last n bases of
 * a simplex run, but instead of updating, it factorizes the basis matrix
 * from scratch with SPEX LU every y updates, and solves with it.  Comparing
 * its rows with sdipps_lp_test's rows for the same updates shows how far the
 * updated factors drift from fresh ones (size, number growth, solve time).
 *
 * Usage:
 *
 *      sdipps_lp_refactor pivot_log [-n N] [-y Y] [-o metrics.csv] [-v] [-f]
 *
 *          -n N    the last N updates (default 50)
 *          -y Y    refactorize every Y updates (default 5)
 *          -o F    write the metrics to F (default sdipps_lp_refactor.csv)
 *          -v      check every solve exactly (default: only the last one)
 *          -f      forward solve only: no transposed solve, except one (not
 *                  timed) for the last basis, to check it
 *
 * With T updates in the log, the basis after update T-N is factorized, and
 * then the basis after every Y-th update from there, and after update T.
 * Each factorization is followed by the same solves as in sdipps_lp_test:
 * B x = b with SPEX_update_solve and B^T y = c_B with SPEX_update_tsolve.
 * The solves of the last basis (or, with -v, of every basis) are checked
 * exactly, and y against QSopt_ex's dual solution if the log has it.
 *
 * Output CSV columns, one row per factorization:
 *
 *      update          the basis is the one after this update (1-based)
 *      factor_time     seconds for SPEX_lu_analyze + SPEX_lu_factorize
 *      solve_time      seconds for SPEX_update_solve,  B x = b
 *      tsolve_time     seconds for SPEX_update_tsolve, B^T y = c_B (empty
 *                      with -f)
 *      nnz_L, nnz_U    nonzero entries of L and U, pivots included (SPEX
 *                      LU also stores entries that cancelled to zero; those
 *                      are not counted)
 *      bits_LU         total bits of all integers the factorization stores
 *      bits_A          total bits of the entries of the basis matrix
 *      max_bits_L, max_bits_U, max_bits_rho, max_bits_scale   largest
 *                      off-diagonal entry of L and of U, largest pivot, and
 *                      largest pending scale (always 0 here), in bits
 */

#include "sdipps_pivot_log.h"

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                          \
{                                               \
    if (out != NULL) fclose (out) ;             \
    free (basis) ;                              \
    free (slot_of) ;                            \
    SPEX_matrix_free (&b, option) ;             \
    SPEX_matrix_free (&x, option) ;             \
    SPEX_matrix_free (&y, option) ;             \
    SPEX_factorization_free (&F, option) ;      \
    lp_free (&L) ;                              \
    SPEX_FREE (option) ;                        \
    SPEX_finalize ( ) ;                         \
}

int main (int argc, char *argv [])
{
    pivot_log L ;
    memset (&L, 0, sizeof (pivot_log)) ;
    FILE *out = NULL ;
    int64_t *basis = NULL, *slot_of = NULL ;
    SPEX_matrix b = NULL, x = NULL, y = NULL ;
    SPEX_factorization F = NULL ;
    SPEX_options option = NULL ;

    //--------------------------------------------------------------------------
    // command line
    //--------------------------------------------------------------------------

    int64_t nrun = 50, every = 5 ;
    const char *out_name = "sdipps_lp_refactor.csv" ;
    bool verify_all = false, forward_only = false ;
    if (argc < 2)
    {
        printf ("usage: %s pivot_log [-n N] [-y Y] [-o metrics.csv] [-v] "
            "[-f]\n", argv [0]) ;
        return (1) ;
    }
    for (int a = 2 ; a < argc ; a++)
    {
        if (strcmp (argv [a], "-n") == 0 && a+1 < argc)
        {
            nrun = atoll (argv [++a]) ;
        }
        else if (strcmp (argv [a], "-y") == 0 && a+1 < argc)
        {
            every = atoll (argv [++a]) ;
        }
        else if (strcmp (argv [a], "-o") == 0 && a+1 < argc)
        {
            out_name = argv [++a] ;
        }
        else if (strcmp (argv [a], "-v") == 0)
        {
            verify_all = true ;
        }
        else if (strcmp (argv [a], "-f") == 0)
        {
            forward_only = true ;
        }
        else
        {
            printf ("unknown option %s\n", argv [a]) ;
            return (1) ;
        }
    }
    if (every < 1)
    {
        printf ("-y must be at least 1\n") ;
        return (1) ;
    }

    //--------------------------------------------------------------------------
    // read the log (before SPEX_initialize; see sdipps_pivot_log.h)
    //--------------------------------------------------------------------------

    double t0 = SUITESPARSE_TIME ;
    if (!lp_read (argv [1], &L))
    {
        lp_free (&L) ;
        return (1) ;
    }
    double t_read = SUITESPARSE_TIME - t0 ;

    SPEX_TRY (SPEX_initialize ( )) ;
    SPEX_TRY (SPEX_create_default_options (&option)) ;

    int64_t T = L.nupdates ;
    if (nrun > T) nrun = T ;
    if (nrun < 0) nrun = 0 ;
    int64_t first = T - nrun ;

    printf ("SDIPPS LP refactorization baseline: %s\n", argv [1]) ;
    printf ("A: %" PRId64 " x %" PRId64 ", nnz %" PRId64 "; %" PRId64
        " updates in the log (read in %.2f s)\n", L.m, L.N, L.nz, T, t_read) ;
    printf ("factorizing the bases after updates %" PRId64 ", %" PRId64
        ", ... (every %" PRId64 "), and %" PRId64 "\n", first, first+every,
        every, T) ;

    //--------------------------------------------------------------------------
    // walk through the bases, factorizing every `every` updates
    //--------------------------------------------------------------------------

    if (!lp_basis_start (&L, &basis, &slot_of))
    {
        FREE_WORKSPACE ;
        return (1) ;
    }
    SPEX_TRY (lp_rhs_b (&b, &L, option)) ;

    out = fopen (out_name, "w") ;
    if (out == NULL)
    {
        perror (out_name) ;
        FREE_WORKSPACE ;
        return (1) ;
    }
    fprintf (out, "update,factor_time,solve_time,tsolve_time,nnz_L,nnz_U,"
        "bits_LU,bits_A,max_bits_L,max_bits_U,max_bits_rho,max_bits_scale\n");

    int nfail = 0 ;
    for (int64_t k = 0 ; k <= T ; k++)
    {
        // here the basis is the one after update k (k = 0: the start)
        bool factorize = (k >= first) &&
            ((k - first) % every == 0 || k == T) ;
        if (factorize)
        {
            double t_factor, t_solve, t_tsolve ;
            lp_metrics M ;
            SPEX_factorization_free (&F, option) ;
            SPEX_TRY (lp_factorize (&F, &t_factor, &L, basis, option)) ;
            lp_factor_metrics (F, &M) ;
            int64_t bits_A = lp_basis_bits (&L, basis) ;

            bool do_tsolve = !forward_only || k == T ;
            SPEX_TRY (lp_solves (&x, &y, &t_solve, &t_tsolve, do_tsolve, F,
                b, &L, basis, option)) ;
            printf ("after update %" PRId64 ": factorized in %.3f s, solve "
                "%.3f s", k, t_factor, t_solve) ;
            if (!forward_only) printf (", transposed solve %.3f s", t_tsolve) ;
            printf ("\n") ;
            if (verify_all || k == T)
            {
                bool ok ;
                SPEX_TRY (lp_check (&ok, x, y, b, &L, basis, k == T,
                    option)) ;
                if (!ok) nfail++ ;
            }
            SPEX_matrix_free (&x, option) ;
            SPEX_matrix_free (&y, option) ;

            char tsolve_field [32] = "" ;      // empty with -f
            if (!forward_only) snprintf (tsolve_field, 32, "%.6f", t_tsolve) ;
            fprintf (out, "%" PRId64 ",%.6f,%.6f,%s", k, t_factor, t_solve,
                tsolve_field) ;
            fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64,
                M.nnz_L, M.nnz_U, M.bits, bits_A) ;
            fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64
                "\n", M.max_bits_L, M.max_bits_U, M.max_bits_rho,
                M.max_bits_scale) ;
            fflush (out) ;
        }

        // move to the basis after update k+1
        if (k < T)
        {
            if (lp_update_slot (&L, slot_of, k) < 0)
            {
                FREE_WORKSPACE ;
                return (1) ;
            }
            lp_basis_apply (&L, basis, slot_of, k) ;
        }
    }

    printf ("metrics written to %s\n", out_name) ;
    printf ("%s\n", (nfail == 0) ? "all checks passed" : "TEST FAILED") ;
    FREE_WORKSPACE ;
    return ((nfail == 0) ? 0 : 1) ;
}
