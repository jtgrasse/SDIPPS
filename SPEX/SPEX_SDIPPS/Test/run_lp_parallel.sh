#!/bin/bash
#-------------------------------------------------------------------------------
# SPEX/SPEX_SDIPPS/Test/run_lp_parallel.sh: the LP test cases, one per L3 cache
#-------------------------------------------------------------------------------

# Usage (from SPEX/SPEX_SDIPPS, after "make Test/sdipps_lp_test
# Test/sdipps_lp_refactor"):
#
#   Test/run_lp_parallel.sh [-n N] [-y Y] [-u on|off] [-d data_dir]
#                           [-o out_root] [-p "priorities"] [-s "scalings"]
#                           [-e "early exits"] [-r "rwsop modes"]
#                           pivot_log ...
#
#       -n N     the last N updates of each log (default 50)
#       -y Y     refactorize every Y updates in the baseline (default 10);
#                0 skips the baseline
#       -u SW    run SPEX's update: on or off (default on)
#       -d DIR   with no pivot_log given, every DIR/*_bases.txt is used
#                (default Test/Data)
#       -o DIR   results go to DIR/<name>/ (default Test/Results); <name> is
#                the log's file name without "_bases.txt"
#       -p LIST  the SDIPPS push priorities to run (default "diag col")
#       -s LIST  the SDIPPS scalings to run (default "regular lazy eager")
#       -e LIST  the SDIPPS early exit settings to run (default "off on")
#       -r LIST  the SDIPPS RwSOP modes to run (default "history"):
#                history uses the seed history; eager does every RwSOP in
#                full, without it (sdipps_lp_test -r)
#
# The test cases are those of run_lp_suite.sh, for every log: SPEX's update,
# the refactorization baseline, and SDIPPS with every combination of push
# priority (from -p), scaling (from -s), early exit (from -e) and RwSOP mode
# (from -r).  All use forward solves only (-f).  Each writes <case>.csv and
# <case>.log:
#
#   spex                        SPEX_update_lu_colrep
#   refactor                    SPEX LU from scratch every Y updates
#   sdipps_<p>[_lazy|_eager][_exit][_rwsop]
#
# Parallel runs: one worker per L3 cache, found from
# /sys/devices/system/cpu/*/cache/index3, except that the last L3 is left
# open for other work when there is more than one.  Each worker is pinned
# (taskset) to the CPUs that share its L3, and runs one test case at a time,
# so no two test cases share an L3 cache.  The workers take test cases from
# one queue (largest log first), so a long case on one L3 does not leave the
# others idle.  The test cases still share the memory bandwidth, so their
# timings can differ somewhat from runs made one at a time.
#
# Progress (start and end of each case, its L3 and its result) goes to
# <out_root>/parallel_progress.log.

set -u

N=50
Y=10
SPEX=on
DATA=Test/Data
OUT=Test/Results
PRIORITIES="diag col"
SCALINGS="regular lazy eager"
EXITS="off on"
RWSOPS="history"
usage ()
{
    echo "usage: $0 [-n N] [-y Y] [-u on|off] [-d data_dir] [-o out_root]" \
         "[-p priorities] [-s scalings] [-e early_exits] [-r rwsop_modes]" \
         "log ..."
    exit 1
}
while getopts "n:y:u:d:o:p:s:e:r:" opt ; do
    case $opt in
        n) N=$OPTARG ;;
        y) Y=$OPTARG ;;
        u) SPEX=$OPTARG ;;
        d) DATA=$OPTARG ;;
        o) OUT=$OPTARG ;;
        p) PRIORITIES=$OPTARG ;;
        s) SCALINGS=$OPTARG ;;
        e) EXITS=$OPTARG ;;
        r) RWSOPS=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

case $SPEX in on|off) ;; *) echo "-u: must be on or off" ; usage ;; esac
case $Y in ''|*[!0-9]*) echo "-y: not a number: $Y" ; usage ;; esac
for p in $PRIORITIES ; do
    case $p in diag|col) ;; *) echo "-p: unknown priority $p" ; usage ;; esac
done
for scaling in $SCALINGS ; do
    case $scaling in regular|lazy|eager) ;;
        *) echo "-s: unknown scaling $scaling" ; usage ;; esac
done
for ex in $EXITS ; do
    case $ex in off|on) ;; *) echo "-e: unknown setting $ex" ; usage ;; esac
done
for rw in $RWSOPS ; do
    case $rw in history|eager) ;;
        *) echo "-r: unknown RwSOP mode $rw" ; usage ;; esac
done

LOGS=("$@")
if [ ${#LOGS[@]} -eq 0 ] ; then
    LOGS=("$DATA"/*_bases.txt)
    if [ ! -e "${LOGS[0]}" ] ; then
        echo "no *_bases.txt in $DATA"
        exit 1
    fi
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

# leave the last L3 cache open for other work, unless it is the only one
if [ ${#L3_GROUPS[@]} -gt 1 ] ; then
    L3_FREE=${L3_GROUPS[-1]}
    unset 'L3_GROUPS[-1]'
    echo "leaving the L3 cache of CPUs $L3_FREE open"
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
    if [ $SPEX = on ] ; then
        echo "spex|$DIR|Test/sdipps_lp_test $LOG -n $N -f -s" >> "$QUEUE"
    fi
    if [ $Y -gt 0 ] ; then
        echo "refactor|$DIR|Test/sdipps_lp_refactor $LOG -n $N -f -y $Y" \
            >> "$QUEUE"
    fi
    for p in $PRIORITIES ; do
        for scaling in $SCALINGS ; do
            for ex in $EXITS ; do
                for rw in $RWSOPS ; do
                    opts=""
                    suffix=""
                    [ $p = col ] && opts="$opts -c"
                    case $scaling in
                        lazy)  opts="$opts -l" ; suffix=_lazy ;;
                        eager) opts="$opts -g" ; suffix=_eager ;;
                    esac
                    if [ $ex = on ] ; then
                        opts="$opts -e"
                        suffix=${suffix}_exit
                    fi
                    if [ $rw = eager ] ; then
                        opts="$opts -r"
                        suffix=${suffix}_rwsop
                    fi
                    echo "sdipps_$p$suffix|$DIR|Test/sdipps_lp_test $LOG" \
                        "-n $N -f$opts" >> "$QUEUE"
                done
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
