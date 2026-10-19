#!/bin/bash
#-------------------------------------------------------------------------------
# SPEX/SPEX_SDIPPS/Test/run_lp_suite.sh: every update method on one pivot log
#-------------------------------------------------------------------------------

# Usage (from SPEX/SPEX_SDIPPS, after "make Test/sdipps_lp_test
# Test/sdipps_lp_refactor"):
#
#   Test/run_lp_suite.sh pivot_log [N] [Y] [out_dir]
#
# Runs, one after another (so the timings do not interfere), on the last N
# updates of pivot_log (default 50), forward solves only (-f):
#
#   spex.csv             SPEX_update_lu_colrep
#   refactor.csv         SPEX LU from scratch every Y updates (default 10)
#   sdipps_<p>[_lazy][_exit].csv
#                        SDIPPS_LU_col_rep with every combination of the
#                        options: <p> = diag or col (push priority), lazy
#                        scaling on/off, early exit on/off
#
# Each run's printed output goes to a .log file next to its .csv.  Results go
# to out_dir (default Results/<name>, <name> = the log's file name without
# "_bases.txt").

set -u
LOG=$1
N=${2:-50}
Y=${3:-10}
NAME=$(basename "$LOG" .txt)
NAME=${NAME%_bases}
OUT=${4:-Results/$NAME}
mkdir -p "$OUT"

run ()      # run <name> <program> <arguments...>
{
    local name=$1 ; shift
    echo "$(date +%H:%M:%S)  $NAME: $name"
    "$@" -o "$OUT/$name.csv" > "$OUT/$name.log" 2>&1
    local status=$?
    tail -1 "$OUT/$name.log"
    return $status
}

failed=0
run spex     Test/sdipps_lp_test     "$LOG" -n "$N" -f -s       || failed=1
run refactor Test/sdipps_lp_refactor "$LOG" -n "$N" -f -y "$Y"  || failed=1
for p in diag col ; do
    for lazy in "" _lazy ; do
        for ex in "" _exit ; do
            opts=""
            [ $p = col ] && opts="$opts -c"
            [ -n "$lazy" ] && opts="$opts -l"
            [ -n "$ex" ] && opts="$opts -e"
            run sdipps_$p$lazy$ex Test/sdipps_lp_test "$LOG" -n "$N" -f $opts \
                || failed=1
        done
    done
done
echo "$(date +%H:%M:%S)  $NAME: done"
exit $failed
