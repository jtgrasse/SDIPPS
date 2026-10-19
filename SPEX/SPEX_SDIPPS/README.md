# SPEX_SDIPPS: SDIPPS column replacement for SPEX LU

`SDIPPS_LU_col_rep` replaces one column of A in an updatable SPEX LU
factorization `P A Q = L D^{-1} U`, in place. The exiting column is pushed
to the end one frame at a time by diagonal and column pushes. The frame being
pushed is carried in two seed vectors (`mu^L`, `mu^U`). The entering column
is then put in by one left-looking REF forward solve on the last column.

It has the same calling convention as `SPEX_update_lu_colrep`, plus options
and optional statistics (`Include/SDIPPS.h`):

```c
SDIPPS_options opt = SDIPPS_DEFAULT_OPTIONS ;   // or NULL for the defaults
SDIPPS_stats stats ;                             // or NULL
info = SDIPPS_LU_col_rep (F, vk, col, &opt, &stats, option) ;
```

## Options

| Option | Values | Effect |
|---|---|---|
| `push_choice` | `SDIPPS_PREFER_DIAGONAL` (default), `SDIPPS_PREFER_COLUMN` | Which push to take when both are possible (`sL sU != 0` and `pi != 0`). If `sL sU == 0`, only a diagonal push is possible; if `pi == 0`, only a column push. |
| `scaling` | `SDIPPS_SCALING_NONE` (default), `SDIPPS_SCALING_EAGER`, `SDIPPS_SCALING_LAZY` | How a backtracking step that is a pure scaling is done. That happens in a diagonal push when the other half's support is 0; the new frame is then the old one times `lambda = rho~_{k-1}/rho_k`. **NONE:** the regular backtrack, with every entry multiplied by `rho~_{k-1}` and divided by `rho_k`. **EAGER:** `lambda` is reduced to `p/q` by `SDIPPS_Reduce` once per run of such pushes, and every entry is divided by `q` and then multiplied by `p`. Dividing first is exact because `SDIPPS_Reduce` returns lowest terms (Chen 2022, §5.4), and it keeps the intermediate numbers smaller. No pending scales are left; the only saving comes from the reduction. **LAZY:** the reduced `p/q` goes into the vector's pending scale and no entry is touched. The frames past an early exit are scaled the same way. |
| `early_exit` | `false` (default), `true` | Stop pushing once the exiting row of U is empty past the diagonal. The later frames are only rescaled. |

These parts are always on:
- **Seed history:** RwSOP updates only the rows it must, and the rest catch up later by the pivot-cancellation history update.
- **Lazy sign tracking:** the `theta` array.

## Layout

| File | Algorithm block |
|---|---|
| `Include/SDIPPS.h` | public interface |
| `Source/SDIPPS_LU_col_rep.c` | main loop, choice of push |
| `Source/SDIPPS_Reduce.c` | reduction of a pivot-over-pivot ratio (GMP GCD for now) |
| `Source/sdipps_diagonal_push.c` | diagonal push |
| `Source/sdipps_column_push.c` | column push |
| `Source/sdipps_backtrack.c` | backtracking (new frame k), lazy scaling |
| `Source/sdipps_rwsop.c` | RwSOP (moves a seed one frame) |
| `Source/sdipps_seed.c` | seed load and the history update (catch-up) |
| `Source/sdipps_enter_column.c` | the swap: one REF forward solve |
| `Source/sdipps_rescale_untouched.c` | frames past an early exit |
| `Source/sdipps_vector.c` | `SPEX_vector` helpers (pending scale, append, ...) |
| `Source/sdipps_workspace.c` | workspace |
| `Source/sdipps_internal.h` | conventions and internal prototypes |

`sdipps_internal.h` explains the conventions in full: 0-based positions, old and new pivots, slots, seeds, lazy signs and pending scales.

## Building into SPEX

No SPEX file is edited. SPEX's `CMakeLists.txt` already compiles every
`SPEX*/Source/*.c` into libspex, so this folder is built with SPEX.
`make local` in `SPEX/` (or `make library` here) puts `SDIPPS_LU_col_rep` in
`SPEX/build/libspex.a` and `libspex.so`. CMake must re-run to pick up a
newly added source file; `make local` does that. Programs include
`SPEX_SDIPPS/Include/SDIPPS.h`.

The output factorization uses SPEX's own updatable format, including pending
scales (`true entry = x[p] * scale`). So `SPEX_update_solve` and
`SPEX_update_tsolve` work on it unchanged, and SDIPPS and
`SPEX_update_lu_colrep` can be used on the same factorization in any order.
Because those scales are GMP `mpq_t` values, which must be canonical,
`SDIPPS_Reduce` must return the ratio in lowest terms; its comments say where
another method plugs in.

## Testing

```
make            # build libspex (with SDIPPS) and the test, run all tests
make run        # run the tests only
```

Usage of `Test/sdipps_test`:

```
Test/sdipps_test matrix_file update_file [-c] [-g | -l] [-e] [-n N]
    -c  prefer column pushes    -g  eager scaling    -l  lazy scaling
    -e  early exit
    -n N  frame check only if n <= N (default 200)
```

The program reads a matrix (SPEX triplet format) and factors it with SPEX LU.
It then applies the column replacements in `update_file` in sequence:

```
number_of_updates
col nz          (1-based column, number of nonzeros)
row value       (nz lines, 1-based row)
...
```

After every update it checks three things:
- **Solve:** `SPEX_update_solve` solves `A x = b` exactly.
- **Transposed solve:** `SPEX_update_tsolve` solves `A^T x = b` exactly.
- **Frames:** every pivot and entry of F equals a from-scratch REF LU of `P A Q` without pivoting.

An update that makes A singular must return `SPEX_SINGULAR`. The test
confirms this, then refactorizes and continues.

`make run` runs every option combination on three data sets in `Test/Data`:
- the 4×4 `ExampleMats/example`;
- a random 12×12 matrix with ±1 entries, where column pushes are common;
- `ExampleMats/10teams` (177×177).

Each has a sequence of 5–12 updates, including one singular update.

## Stress tests on a simplex run (pivot log)

These tests replay the last N basis changes of a simplex run recorded by QSopt_ex. They stress the update on real LP bases; testing updates between exact checks is a separate test.

```
make lp            # SDIPPS on the last N updates          -> <name>_sdipps.csv
make lp_spex       # the same with SPEX_update_lu_colrep    -> <name>_spex.csv
make lp_refactor   # SPEX LU from scratch every Y updates   -> <name>_refactor.csv
make lp LOG=Test/Data/cycle_bases.txt N=100     # another log, more updates
```

The defaults are `LOG=Test/Data/co9_bases.txt`, `N=50` and `Y=5`. `<name>` is the log's file name without `_bases.txt`, for example `co9`.

The programs behind these targets:
```
Test/sdipps_lp_test     pivot_log [-n N] [-o out.csv] [-s] [-c] [-l] [-e] [-v]
Test/sdipps_lp_refactor pivot_log [-n N] [-y Y] [-o out.csv] [-v]
```

- **`-s`** runs SPEX's own update (`SPEX_update_lu_colrep`) in place of SDIPPS, with everything else unchanged.
- **`-c`, `-g`, `-l`, `-e`** are the SDIPPS options: prefer column pushes, eager scaling, lazy scaling, early exit.
- **`-v`** checks every solve exactly; by default only the last one is checked.

**The pivot log** (`Test/sdipps_pivot_log.c/.h`, shared by both programs) holds:
- A, b and c as exact rationals, with the logical columns included;
- the starting basis;
- one `U leaving entering` line per basis change;
- the final exact solution (x, y, d, z), of which only y is used.

A is made integral column by column: each column is multiplied by the LCM of its denominators. Every basis matrix is therefore `B D` with D diagonal, and every entering column arrives in a consistent integer scale. An entering column takes the leaving column's slot (its column of the basis matrix).

**`sdipps_lp_test`** proceeds as follows:
1. Replay the first T−N updates as bookkeeping only.
2. Factorize that basis with SPEX LU.
3. For each of the last N updates: run the update, record the factorization's metrics, then solve `B x = b` (`SPEX_update_solve`) and `Bᵀ y = c_B` (`SPEX_update_tsolve`). The metrics come before the solves, because the solves apply pending scales in place.

**`sdipps_lp_refactor`** is the baseline. It walks through the same bases, but factorizes from scratch with SPEX LU at the starting basis, every Y updates, and the last basis, and runs the same solves. Comparing its rows with the update test's rows for the same update numbers shows how far the updated factors drift from fresh ones.

**Both programs check** the solves at the last basis exactly against the basis matrix built from scratch. With `-v` they check every solve. At the last basis, y must also equal QSopt_ex's exact dual solution from the log.

**CSV columns:**

| Column | Meaning |
|---|---|
| `update` | the update's number in the log; the factorized basis is the one after it |
| `update_time` / `factor_time` | seconds for the update / for `SPEX_lu_analyze` + `SPEX_lu_factorize` |
| `solve_time` | seconds for `SPEX_update_solve`, `B x = b` |
| `tsolve_time` | seconds for `SPEX_update_tsolve`, `Bᵀ y = c_B` (dense right-hand side, so not fast) |
| `nnz_L`, `nnz_U` | nonzero entries of L and U, pivots included. SPEX LU also stores entries that cancelled to zero; those are not counted. |
| `bits_LU` | total bits of every integer the factorization stores: entries of L and U, pending scales, pivots. Taken right after the update (or factorization). |
| `bits_LU_solved` | (update test only) the same, after that row's solves. The solves apply pending scales in place: `SPEX_update_solve` those of the L columns it uses, `SPEX_update_tsolve` those of the U rows. With lazy scaling, the difference from `bits_LU` is the size that was still pending. |
| `bits_A` | total bits of the entries of the basis matrix `B D` |
| `pushes_to_end` | (update test only) position of the exiting column before the update, as the number of pushes to the last position, `n−1−xi` |
| `max_bits_L`, `max_bits_U` | largest off-diagonal entry of L / U as stored, in bits |
| `max_bits_rho` | largest pivot, in bits |
| `max_bits_scale` | largest pending scale p/q of a vector, as bits(p) + bits(q); 0 if there is none |
| `exit_nnz_L`, `exit_nnz_U` | (update test only) nonzeros of the exiting frame before the update: its column of L and row of U, pivot included |
| `enter_nnz` | (update test only) nonzeros of the entering column |

**SDIPPS-only columns** of the update test, from `SDIPPS_stats`; they are empty for SPEX's update (`-s`):

| Column | Meaning |
|---|---|
| `pushes_skipped` | pushes saved by an early exit, `n−1−landing` |
| `n_diag`, `n_col` | number of diagonal and column pushes |
| `diag_forced` / `diag_chosen` | diagonal pushes with `sL sU = 0` (the only choice) / where both pushes were possible |
| `col_forced` / `col_chosen` | column pushes with `pi = 0` (the only choice) / where both were possible |
| `pure_L`, `pure_U` | backtracking steps that were a pure scaling, by half: where eager and lazy scaling act |
| `reduce_calls` | calls of `SDIPPS_Reduce`: one per run of pure-scaling pushes plus one after an early exit; eager and lazy scaling only |
| `reduce_effective` | calls whose ratio had a common factor |
| `reduce_bits_saved` | total bits removed by those reductions: bits(num) + bits(den) − bits(p) − bits(q) |
| `max_nnz_muL`, `max_nnz_muU` | the largest number of nonzeros each seed reached during the update |
| `max_bits_muL`, `max_bits_muU` | the largest total bits each seed reached, as stored: a deferred entry is counted at its last history |
