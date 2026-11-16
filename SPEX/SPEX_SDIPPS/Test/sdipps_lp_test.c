//------------------------------------------------------------------------------
// SPEX_SDIPPS/Test/sdipps_lp_test.c: SDIPPS on the bases of a simplex run
//------------------------------------------------------------------------------

/* A stress test of the update: replays the last n basis changes of a simplex
 * run with SDIPPS_LU_col_rep (or SPEX's update), solving after each one, and
 * writes metrics to a CSV file.
 *
 * Usage:
 *
 *      sdipps_lp_test pivot_log [-n N] [-o metrics.csv] [-s] [-c] [-g | -l]
 *                     [-e] [-r] [-v] [-f]
 *
 *          -n N    replay the last N updates (default 50)
 *          -o F    write the metrics to F (default sdipps_lp_metrics.csv)
 *          -s      use SPEX's update, SPEX_update_lu_colrep, instead of
 *                  SDIPPS_LU_col_rep (for comparison; -c -g -l -e -r are
 *                  ignored)
 *          -c      prefer column pushes (default: prefer diagonal pushes)
 *          -g      eager backtracking scaling (SDIPPS_SCALING_EAGER)
 *          -l      lazy backtracking scaling (SDIPPS_SCALING_LAZY)
 *          -e      early exit
 *          -r      eager RwSOP, without the seed history (benchmark)
 *          -v      check every solve exactly (default: only the last one)
 *          -f      forward solve only: no transposed solve, except one (not
 *                  timed, after the metrics) after the last update, to
 *                  check it
 *
 * See sdipps_pivot_log.h for the log format and the integral basis matrices
 * B D.  Steps:
 *
 *  1. With T updates in the log, the first T-n are replayed as bookkeeping.
 *  2. That basis matrix is factorized with SPEX LU.
 *  3. Each of the last n updates is done with the update method, replacing
 *     the leaving column's slot.  The metrics of the new factorization are
 *     taken, then B x = b and B^T y = c_B are solved with SPEX_update_solve
 *     and SPEX_update_tsolve.  (The metrics come before the solves because
 *     the solves apply pending scales in place.)
 *  4. The solves after the last update (or, with -v, every update) are
 *     checked exactly; at the last update y must also equal QSopt_ex's dual
 *     solution, if the log has it.
 *
 * Output CSV columns, one row per update:
 *
 *      update          the update's number in the log (1-based)
 *      update_time     seconds for the update
 *      solve_time      seconds for SPEX_update_solve,  B x = b
 *      tsolve_time     seconds for SPEX_update_tsolve, B^T y = c_B (its right
 *                      hand side is dense, so this solve is not fast); empty
 *                      with -f
 *      nnz_L, nnz_U    nonzero entries of L and U, pivots included (SPEX
 *                      LU also stores entries that cancelled to zero; those
 *                      are not counted)
 *      bits_LU         total bits of all integers the factorization stores
 *                      (entries of L and U, pending scales, pivots), right
 *                      after the update
 *      bits_LU_solved  the same after the solves of this row, which apply
 *                      pending scales in place (SPEX_update_solve: those of
 *                      the L columns it uses; SPEX_update_tsolve: U rows)
 *      bits_A          total bits of the entries of the basis matrix
 *      max_bits_L, max_bits_U   largest off-diagonal entry of L / U (as
 *                      stored), in bits, after the update
 *      max_bits_rho    largest pivot, in bits
 *      max_bits_scale  largest pending scale p/q of a vector: bits (p) +
 *                      bits (q) (0 if none)
 *      pushes_to_end   position of the exiting column before the update, as
 *                      the number of pushes to the last position: n-1-xi
 *      exit_nnz_L, exit_nnz_U   nonzeros of the exiting frame before the
 *                      update: its column of L and row of U, pivot included
 *      enter_nnz       nonzeros of the entering column
 *
 * SDIPPS only (from SDIPPS_stats; empty with -s):
 *
 *      pushes_skipped  pushes saved by an early exit: n-1-landing
 *      n_diag, n_col   number of diagonal and column pushes
 *      diag_forced     diagonal pushes with sL sU = 0 (the only choice)
 *      diag_chosen     diagonal pushes where both were possible
 *      col_forced      column pushes with pi = 0 (the only choice)
 *      col_chosen      column pushes where both were possible
 *      pure_L, pure_U  backtracking steps that were a pure scaling, by half
 *      reduce_calls    calls of SDIPPS_Reduce (eager and lazy scaling only)
 *      reduce_effective  calls where the ratio had a common factor
 *      reduce_bits_saved total bits removed by the reductions
 *      max_nnz_muL, max_nnz_muU, max_bits_muL, max_bits_muU   the largest
 *                      number of nonzeros and total bits each seed reached
 *      rwsop_rows      seed rows updated by the full RwSOP formula
 *      scaled_rows     seed rows updated by a pure scaling (with -r, one per
 *                      row per frame; otherwise by the history update)
 */

#include "sdipps_pivot_log.h"
#include "SDIPPS.h"

// The SDIPPS_stats columns of a row (empty fields if s is NULL)
static void write_sdipps_stats (FILE *out, const SDIPPS_stats *s, int64_t m)
{
    if (s == NULL)
    {
        for (int t = 0 ; t < 18 ; t++) fprintf (out, ",") ;
        return ;
    }
    fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64, (m-1) - s->landing,
        s->n_diagonal, s->n_column) ;
    fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64,
        s->n_diagonal_forced, s->n_diagonal_chosen, s->n_column_forced,
        s->n_column_chosen) ;
    fprintf (out, ",%" PRId64 ",%" PRId64, s->n_pure_L, s->n_pure_U) ;
    fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64, s->n_reduce,
        s->n_reduced, s->reduce_bits_saved) ;
    fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64,
        s->max_nnz_mu_L, s->max_nnz_mu_U, s->max_bits_mu_L,
        s->max_bits_mu_U) ;
    fprintf (out, ",%" PRId64 ",%" PRId64, s->n_rwsop_rows, s->n_scaled_rows) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                          \
{                                               \
    if (out != NULL) fclose (out) ;             \
    free (basis) ;                              \
    free (slot_of) ;                            \
    SPEX_matrix_free (&b, option) ;             \
    SPEX_matrix_free (&x, option) ;             \
    SPEX_matrix_free (&y, option) ;             \
    SPEX_matrix_free (&vk, option) ;            \
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
    SPEX_matrix b = NULL, x = NULL, y = NULL, vk = NULL ;
    SPEX_factorization F = NULL ;
    SPEX_options option = NULL ;

    //--------------------------------------------------------------------------
    // command line
    //--------------------------------------------------------------------------

    SDIPPS_options opt = SDIPPS_DEFAULT_OPTIONS ;
    int64_t nrun = 10 ;
    const char *out_name = "sdipps_lp_metrics.csv" ;
    bool verify_all = false, use_spex = false, forward_only = false ;
    if (argc < 2)
    {
        printf ("usage: %s pivot_log [-n N] [-o metrics.csv] [-s] [-c]"
            " [-g | -l] [-e] [-r] [-v] [-f]\n", argv [0]) ;
        return (1) ;
    }
    for (int a = 2 ; a < argc ; a++)
    {
        if (strcmp (argv [a], "-n") == 0 && a+1 < argc)
        {
            nrun = atoll (argv [++a]) ;
        }
        else if (strcmp (argv [a], "-o") == 0 && a+1 < argc)
        {
            out_name = argv [++a] ;
        }
        else if (strcmp (argv [a], "-s") == 0)
        {
            use_spex = true ;
        }
        else if (strcmp (argv [a], "-c") == 0)
        {
            opt.push_choice = SDIPPS_PREFER_COLUMN ;
        }
        else if (strcmp (argv [a], "-g") == 0)
        {
            opt.scaling = SDIPPS_SCALING_EAGER ;
        }
        else if (strcmp (argv [a], "-l") == 0)
        {
            opt.scaling = SDIPPS_SCALING_LAZY ;
        }
        else if (strcmp (argv [a], "-e") == 0)
        {
            opt.early_exit = true ;
        }
        else if (strcmp (argv [a], "-r") == 0)
        {
            opt.eager_rwsop = true ;
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

    int64_t m = L.m, T = L.nupdates ;
    if (nrun > T) nrun = T ;
    if (nrun < 0) nrun = 0 ;
    int64_t first = T - nrun ;      // updates first..T-1 are run (0-based)

    printf ("SDIPPS LP test: %s\n", argv [1]) ;
    printf ("A: %" PRId64 " x %" PRId64 ", nnz %" PRId64 "; %" PRId64
        " updates in the log (read in %.2f s)\n", m, L.N, L.nz, T, t_read) ;
    if (use_spex)
    {
        printf ("update method: SPEX_update_lu_colrep\n") ;
    }
    else
    {
        const char *scaling = (opt.scaling == SDIPPS_SCALING_EAGER) ? "eager" :
            (opt.scaling == SDIPPS_SCALING_LAZY) ? "lazy" : "regular" ;
        printf ("update method: SDIPPS_LU_col_rep (prefer %s pushes, %s "
            "scaling, early exit %s, seed history %s)\n",
            (opt.push_choice == SDIPPS_PREFER_COLUMN) ? "column" : "diagonal",
            scaling, opt.early_exit ? "on" : "off",
            opt.eager_rwsop ? "off" : "on") ;
    }

    //--------------------------------------------------------------------------
    // replay updates 0..first-1 as bookkeeping, then factorize
    //--------------------------------------------------------------------------

    if (!lp_basis_start (&L, &basis, &slot_of))
    {
        FREE_WORKSPACE ;
        return (1) ;
    }
    for (int64_t k = 0 ; k < first ; k++)
    {
        if (lp_update_slot (&L, slot_of, k) < 0)
        {
            FREE_WORKSPACE ;
            return (1) ;
        }
        lp_basis_apply (&L, basis, slot_of, k) ;
    }

    double t_factor ;
    SPEX_TRY (lp_factorize (&F, &t_factor, &L, basis, option)) ;
    SPEX_TRY (lp_rhs_b (&b, &L, option)) ;
    lp_metrics M ;
    lp_factor_metrics (F, &M) ;
    printf ("starting basis (after update %" PRId64 "): SPEX LU in %.3f s, "
        "nnz(L) %" PRId64 ", nnz(U) %" PRId64 ", %" PRId64 " bits\n",
        first, t_factor, M.nnz_L, M.nnz_U, M.bits) ;
    printf ("running updates %" PRId64 "..%" PRId64 "\n", first+1, T) ;

    out = fopen (out_name, "w") ;
    if (out == NULL)
    {
        perror (out_name) ;
        FREE_WORKSPACE ;
        return (1) ;
    }
    fprintf (out, "update,update_time,solve_time,tsolve_time,"
        "nnz_L,nnz_U,bits_LU,bits_LU_solved,bits_A,"
        "max_bits_L,max_bits_U,max_bits_rho,max_bits_scale,"
        "pushes_to_end,exit_nnz_L,exit_nnz_U,enter_nnz,"
        "pushes_skipped,n_diag,n_col,diag_forced,diag_chosen,col_forced,"
        "col_chosen,pure_L,pure_U,reduce_calls,reduce_effective,"
        "reduce_bits_saved,max_nnz_muL,max_nnz_muU,max_bits_muL,"
        "max_bits_muU,rwsop_rows,scaled_rows\n") ;

    //--------------------------------------------------------------------------
    // the last nrun updates
    //--------------------------------------------------------------------------

    int nfail = 0 ;
    double total_update = 0, total_solve = 0, total_tsolve = 0 ;
    for (int64_t k = first ; k < T ; k++)
    {
        int64_t slot = lp_update_slot (&L, slot_of, k) ;
        if (slot < 0)
        {
            FREE_WORKSPACE ;
            return (1) ;
        }
        // position of the exiting column in the factorization, before the
        // update, as the number of pushes to the last position
        int64_t xi = F->Qinv_perm [slot] ;
        int64_t pushes_to_end = (m-1) - xi ;
        int64_t exit_nnz_L = lp_vector_nnz (F->L->v [xi]) ;
        int64_t exit_nnz_U = lp_vector_nnz (F->U->v [xi]) ;
        SPEX_TRY (lp_column (&vk, &L, L.entering [k], option)) ;
        int64_t enter_nnz = lp_vector_nnz (vk->v [0]) ;

        // the update: the entering column replaces column `slot`
        SPEX_info status ;
        SDIPPS_stats stats ;
        t0 = SUITESPARSE_TIME ;
        if (use_spex)
        {
            status = SPEX_update_lu_colrep (F, vk, slot, option) ;
        }
        else
        {
            status = SDIPPS_LU_col_rep (F, vk, slot, &opt, &stats, option) ;
        }
        double t_update = SUITESPARSE_TIME - t0 ;
        if (status != SPEX_OK)
        {
            printf ("update %" PRId64 ": the update failed (%d)\n", k+1,
                (int) status) ;
            FREE_WORKSPACE ;
            return (1) ;
        }
        SPEX_matrix_free (&vk, option) ;
        lp_basis_apply (&L, basis, slot_of, k) ;

        // metrics of the updated factorization (before the solves)
        lp_factor_metrics (F, &M) ;
        int64_t bits_A = lp_basis_bits (&L, basis) ;

        // B x = b and B^T y = c_B (with -f, B x = b only)
        double t_solve, t_tsolve ;
        SPEX_TRY (lp_solves (&x, &y, &t_solve, &t_tsolve, !forward_only, F, b,
            &L, basis, option)) ;

        // bits after the solves, which apply pending scales in place
        lp_metrics M_solved ;
        lp_factor_metrics (F, &M_solved) ;

        // with -f, solve B^T y = c_B after the last update anyway (not timed,
        // after the metrics), to check it against QSopt_ex's y
        if (forward_only && k == T-1)
        {
            SPEX_matrix cB = NULL ;
            SPEX_TRY (lp_rhs_cB (&cB, &L, basis, option)) ;
            SPEX_info tsolve_status = SPEX_update_tsolve (&y, F, cB, option) ;
            SPEX_matrix_free (&cB, option) ;
            SPEX_TRY (tsolve_status) ;
        }
        if (verify_all || k == T-1)
        {
            bool ok ;
            printf ("update %" PRId64 ":\n", k+1) ;
            SPEX_TRY (lp_check (&ok, x, y, b, &L, basis, k == T-1, option)) ;
            if (!ok) nfail++ ;
        }
        SPEX_matrix_free (&x, option) ;
        SPEX_matrix_free (&y, option) ;

        char tsolve_field [32] = "" ;      // empty with -f
        if (!forward_only) snprintf (tsolve_field, 32, "%.6f", t_tsolve) ;
        fprintf (out, "%" PRId64 ",%.6f,%.6f,%s", k+1, t_update, t_solve,
            tsolve_field) ;
        fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64
            ",%" PRId64, M.nnz_L, M.nnz_U, M.bits, M_solved.bits, bits_A) ;
        fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64,
            M.max_bits_L, M.max_bits_U, M.max_bits_rho, M.max_bits_scale) ;
        fprintf (out, ",%" PRId64 ",%" PRId64 ",%" PRId64 ",%" PRId64,
            pushes_to_end, exit_nnz_L, exit_nnz_U, enter_nnz) ;
        write_sdipps_stats (out, use_spex ? NULL : &stats, m) ;
        fprintf (out, "\n") ;
        fflush (out) ;
        total_update += t_update ;
        total_solve += t_solve ;
        if (!forward_only) total_tsolve += t_tsolve ;
    }

    printf ("%" PRId64 " updates: total time %.3f s updating, %.3f s in solves,"
        " %.3f s in transposed solves\n", nrun, total_update, total_solve,
        total_tsolve) ;
    printf ("metrics written to %s\n", out_name) ;
    printf ("%s\n", (nfail == 0) ? "all checks passed" : "TEST FAILED") ;
    FREE_WORKSPACE ;
    return ((nfail == 0) ? 0 : 1) ;
}
