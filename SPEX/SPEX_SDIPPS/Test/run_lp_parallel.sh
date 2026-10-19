#!/bin/bash
#-------------------------------------------------------------------------------
# SPEX/SPEX_SDIPPS/Test/run_lp_parallel.sh: the LP test cases, one per L3 cache
#-------------------------------------------------------------------------------

# Usage (from SPEX/SPEX_SDIPPS, after "make Test/sdipps_lp_test
# Test/sdipps_lp_refactor"):
#
#   Test/run_lp_parallel.sh [-n N] [-y Y] [-o out_root] [-s "scalings"]
#                           pivot_log ...
#
#       -n N     the last N updates of each log (default 50)
#       -y Y     refactorize every Y updates in the baseline (default 10)
#       -o DIR   results go to DIR/<name>/ (default Results); <name> is the
#                log's file name without "_bases.txt"
#       -s LIST  the SDIPPS scalings to run (default "regular lazy eager")
#
# With no pivot_log given, every Test/Data/*_bases.txt is used.
#
# The test cases are those of run_lp_suite.sh, for every log: SPEX's update,
# the refactorization baseline, and SDIPPS with every combination of push
# priority (diag, col), scaling (from -s) and early exit (on, off).  All use
# forward solves only (-f).  Each writes <case>.csv and <case>.log:
#
#   spex                        SPEX_update_lu_colrep
#   refactor                    SPEX LU from scratch every Y updates
#   sdipps_<p>[_lazy|_eager][_exit]
#
# Parallel runs: one worker per L3 cache, found from
# /sys/devices/system/cpu/*/cache/index3.  Each worker is pinned (taskset) to
# the CPUs that share its L3, and runs one test case at a time, so no two
# test cases share an L3 cache.  The workers take test cases from one queue
# (largest log first), so a long case on one L3 does not leave the other
# idle.  The test cases still share the memory bandwidth, so their timings
# can differ somewhat from runs made one at a time.
#
# Progress (start and end of each case, its L3 and its result) goes to
# <out_root>/parallel_progress.log.

set -u

N=50
Y=10
OUT=Results
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

#-------------------------------------------------------------------------------
# the L3 caches: one CPU list per cache, e.g. "0-7" and "8-15"
#-------------------------------------------------------------------------------

L3_GROUPS=($(cat /sys/devices/system/cpu/cpu*/cache/index3/shared_cpu_list \
    2>/dev/null | sort -u -t- -k1,1n))
if [ ${#L3_GROUPS[@]} -eq 0 ] ; then
    echo "could not find the L3 caches; running on one worker"
    L3_GROUPS=("0-$(( $(nproc) - 1 ))")
fi

#-------------------------------------------------------------------------------
# the queue of test cases: one line per case, "name|out_dir|command ..."
#-------------------------------------------------------------------------------

mkdir -p "$OUT"
QUEUE="$OUT/.parallel_queue"
NEXT="$OUT/.parallel_next"         # index of the next case to take
LOCK="$OUT/.parallel_lock"
PROGRESS="$OUT/parallel_progress.log"
: > "$QUEUE"
echo 0 > "$NEXT"

# largest log first: its cases take longest
SORTED_LOGS=($(ls -S "${LOGS[@]}"))

for LOG in "${SORTED_LOGS[@]}" ; do
    NAME=$(basename "$LOG" .txt)
    NAME=${NAME%_bases}
    DIR="$OUT/$NAME"
    mkdir -p "$DIR"
    echo "spex|$DIR|Test/sdipps_lp_test $LOG -n $N -f -s" >> "$QUEUE"
    echo "refactor|$DIR|Test/sdipps_lp_refactor $LOG -n $N -f -y $Y" >> "$QUEUE"
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
                echo "sdipps_$p$suffix$ex|$DIR|Test/sdipps_lp_test $LOG -n $N -f$opts" \
                    >> "$QUEUE"
            done
        done
    done
done

NCASES=$(wc -l < "$QUEUE")
echo "$(date +%H:%M:%S)  $NCASES test cases, ${#L3_GROUPS[@]} workers" \
     "(L3 caches: ${L3_GROUPS[*]})" | tee "$PROGRESS"

#-------------------------------------------------------------------------------
# a worker: take the next case from the queue until it is empty
#-------------------------------------------------------------------------------

worker ()       # worker <L3 number> <CPU list>
{
    local l3=$1 cpus=$2
    local failed=0
    while true ; do
        # take the next case (the lock keeps two workers from taking the
        # same one)
        local index
        index=$(
            flock "$LOCK" bash -c '
                i=$(cat "'"$NEXT"'")
                echo $((i + 1)) > "'"$NEXT"'"
                echo $i'
        )
        if [ "$index" -ge "$NCASES" ] ; then
            break
        fi
        local line name dir cmd
        line=$(sed -n "$((index + 1))p" "$QUEUE")
        name=${line%%|*}
        line=${line#*|}
        dir=${line%%|*}
        cmd=${line#*|}

        echo "$(date +%H:%M:%S)  L3 $l3 (CPUs $cpus): start $dir/$name" \
            >> "$PROGRESS"
        taskset -c "$cpus" $cmd -o "$dir/$name.csv" > "$dir/$name.log" 2>&1
        local status=$?
        [ $status -ne 0 ] && failed=1
        echo "$(date +%H:%M:%S)  L3 $l3 (CPUs $cpus): done  $dir/$name:" \
            "$(tail -1 "$dir/$name.log")" >> "$PROGRESS"
    done
    return $failed
}

#-------------------------------------------------------------------------------
# one worker per L3 cache
#-------------------------------------------------------------------------------

pids=()
for l3 in "${!L3_GROUPS[@]}" ; do
    worker "$l3" "${L3_GROUPS[$l3]}" &
    pids+=($!)
done

failed=0
for pid in "${pids[@]}" ; do
    wait "$pid" || failed=1
done
rm -f "$QUEUE" "$NEXT" "$LOCK"

echo "$(date +%H:%M:%S)  all test cases done" \
     "($( [ $failed -eq 0 ] && echo "all passed" || echo "SOME FAILED" ))" \
     | tee -a "$PROGRESS"
exit $failed
