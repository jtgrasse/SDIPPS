#!/bin/bash
#-------------------------------------------------------------------------------
# SPEX/SPEX_SDIPPS/Test/run_lp_suite.sh: every update method on the pivot logs
#-------------------------------------------------------------------------------

# Usage (from SPEX/SPEX_SDIPPS, after "make Test/sdipps_lp_test
# Test/sdipps_lp_refactor"):
#
#   Test/run_lp_suite.sh [-n N] [-y Y] [-o out_root] [-s "scalings"]
#                        pivot_log ...
#
#       -n N     the last N updates of each log (default 50)
#       -y Y     refactorize every Y updates in the baseline (default 10)
#       -o DIR   results go to DIR/<name>/ (default Test/Results); <name> is
#                the log's file name without "_bases.txt"
#       -s LIST  the SDIPPS scalings to run (default "regular lazy eager")
#
# With no pivot_log given, every Test/Data/*_bases.txt is used.
#
# Runs, one after another (so the timings do not interfere), for every log:
# SPEX's update, the refactorization baseline, and SDIPPS with every
# combination of push priority (diag, col), scaling (from -s) and early exit
# (on, off).  All use forward solves only (-f).  Each writes <case>.csv and
# <case>.log:
#
#   spex                        SPEX_update_lu_colrep
#   refactor                    SPEX LU from scratch every Y updates
#   sdipps_<p>[_lazy|_eager][_exit]
#
# run_lp_parallel.sh runs the same test cases in parallel.

set -u

N=50
Y=10
OUT=Test/Results
SCALINGS="regular lazy eager"
while getopts "n:y:o:s:" opt ; do
    case $opt in
        n) N=$OPTARG ;;
        y) Y=$OPTARG ;;
        o) OUT=$OPTARG ;;
        s) SCALINGS=$OPTARG ;;
        *) echo "usage: $0 [-n N] [-y Y] [-o out_root] [-s scalings] log ..." ;
           exit 1 ;;
    esac
done
shift $((OPTIND - 1))
LOGS=("$@")
if [ ${#LOGS[@]} -eq 0 ] ; then
    LOGS=(Test/Data/*_bases.txt)
fi

for prog in Test/sdipps_lp_test Test/sdipps_lp_refactor ; do
    if [ ! -x $prog ] ; then
        echo "$prog not found: run make Test/sdipps_lp_test" \
             "Test/sdipps_lp_refactor first"
        exit 1
    fi
done

run ()      # run <name> <program> <arguments...>
{
    local name=$1 ; shift
    echo "$(date +%H:%M:%S)  $NAME: $name"
    "$@" -o "$DIR/$name.csv" > "$DIR/$name.log" 2>&1
    local status=$?
    tail -1 "$DIR/$name.log"
    return $status
}

failed=0
for LOG in "${LOGS[@]}" ; do
    NAME=$(basename "$LOG" .txt)
    NAME=${NAME%_bases}
    DIR="$OUT/$NAME"
    mkdir -p "$DIR"
    run spex     Test/sdipps_lp_test     "$LOG" -n "$N" -f -s       || failed=1
    run refactor Test/sdipps_lp_refactor "$LOG" -n "$N" -f -y "$Y"  || failed=1
    for p in diag col ; do
        for scaling in $SCALINGS ; do
            for ex in "" _exit ; do
                opts=""
                suffix=""
                [ $p = col ] && opts="$opts -c"
                case $scaling in
                    lazy)  opts="$opts -l" ; suffix=_lazy ;;
                    eager) opts="$opts -g" ; suffix=_eager ;;
                esac
                [ -n "$ex" ] && opts="$opts -e"
                run sdipps_$p$suffix$ex Test/sdipps_lp_test "$LOG" -n "$N" \
                    -f $opts || failed=1
            done
        done
    done
    echo "$(date +%H:%M:%S)  $NAME: done"
done

echo "$(date +%H:%M:%S)  all test cases done" \
     "($( [ $failed -eq 0 ] && echo "all passed" || echo "SOME FAILED" ))"
exit $failed
