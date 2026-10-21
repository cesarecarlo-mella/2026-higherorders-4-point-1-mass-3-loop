#!/bin/bash
# Iterated-integral solution of all families, no scheduler: runs the (family, crossing) jobs of
# run_iterint.wls as background processes, NJOBS at a time, then assembles every family whose six
# crossings all succeeded (assemble_II.wls -> II/<fam>_II.m).
#
#   nohup ./run_all.sh 8 > logs/run_all.log 2>&1 &          # all 22 families, 8 jobs at a time
#   nohup ./run_all.sh 4 H I3 > logs/run_all.log 2>&1 &     # only some families
#
# NJOBS must not exceed the Mathematica kernels the license allows, and the jobs must fit in the
# machine's memory (H at weight 8: ~2 GB per job; I1, I2, I3 not measured).
# Environment: WMAX (default 8); VALIDATE_DIR=<copy of all.ints/iterint_solutions> also compares
# Sum_{k<=6} of every crossing with the old file. Resumable: rerunning reloads the weights already
# in out/<fam>/. Per-job logs: logs/<fam>_<cross>.log; summary: logs/run_all.log.
# WL is the command that runs a Wolfram script (default "wolframscript -file"); on a machine
# without wolframscript use e.g. WL="math -script" or WL="/path/to/WolframKernel -script".
# A kernel that does not come up within STARTUP seconds (default 180; e.g. no free license or a
# failed link) is killed and the job retried, up to 3 times.
set -u
cd "$(dirname "$0")"
NJOBS=${1:-4}; shift || true
export WMAX=${WMAX:-8} VALIDATE_DIR=${VALIDATE_DIR:-} STARTUP=${STARTUP:-180} WL=${WL:-wolframscript -file}
if ! command -v ${WL%% *} > /dev/null; then echo "'${WL%% *}' not found: load Mathematica (module load ...) or set WL" >&2; exit 1; fi
mkdir -p logs out II
# largest families first, so the long jobs start early
if [ $# -gt 0 ]; then FAMS="$*"; else
  FAMS=$(for f in $(cat families.txt); do echo "$(wc -c < input/atilde/${f}_atilde.m) $f"; done | sort -rn | awk '{print $2}'); fi
echo "$(date '+%F %T')  families: $(echo $FAMS | tr '\n' ' ')  NJOBS=$NJOBS WMAX=$WMAX WL=$WL"

# wl <log> <script> <args...>: run wolframscript with a startup watchdog and retries
wl() {
  local log=$1; shift
  local attempt pid s
  for attempt in 1 2 3; do
    sleep $((RANDOM % 10))                       # stagger kernel launches
    $WL "$@" > "$log" 2>&1 &
    pid=$!
    for ((s = 0; s < STARTUP; s++)); do
      grep -q "kernel started" "$log" 2>/dev/null && break
      kill -0 $pid 2>/dev/null || break
      sleep 1
    done
    if grep -q "kernel started" "$log" 2>/dev/null; then wait $pid; return $?; fi
    if ! kill -0 $pid 2>/dev/null; then wait $pid; return $?; fi
    echo "$(date '+%F %T')  kernel did not start ($log), attempt $attempt: killing" >&2
    pkill -9 -P $pid 2>/dev/null; kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null
    sleep 30
  done
  return 1
}

run_one() {
  local fam=$1 cross=$2 log="logs/$1_$2.log" extra=""
  if [ -n "$VALIDATE_DIR" ]; then
    if [ "$fam" = H ] || [ "$fam" = I3 ]; then extra="validate=$VALIDATE_DIR/$fam.wl vmax=6"
    else extra="validate=$VALIDATE_DIR/${fam}_II.m vmax=6"; fi
  fi
  if wl "$log" run_iterint.wls family="$fam" atilde="input/atilde/${fam}_atilde.m" \
       bc="../families/results/${fam}_bc.m" cross="$cross" wmax="$WMAX" out="out/$fam" $extra \
     && grep -q "^w$WMAX: " "$log" && ! grep -q "differing keys" "$log"; then
    echo "$(date '+%F %T')  done   $fam $cross  $(grep "^w$WMAX: " "$log")"; touch "out/$fam/.ok_$cross"
  else
    echo "$(date '+%F %T')  FAILED $fam $cross ($log)"
  fi
}

assemble_one() {
  local fam=$1 log="logs/$1_assemble.log"
  if wl "$log" assemble_II.wls family="$fam" in="out/$fam" outdir=II wmax="$WMAX" && grep -q "^wrote" "$log"; then
    echo "$(date '+%F %T')  $(grep "^wrote" "$log")"
  else
    echo "$(date '+%F %T')  FAILED assembly $fam ($log)"
  fi
}
export -f wl run_one assemble_one

for f in $FAMS; do mkdir -p "out/$f"; rm -f out/$f/.ok_*; for c in x x12 x13 x23 x123 x132; do echo "$f $c"; done; done \
  | xargs -P "$NJOBS" -n 2 bash -c 'run_one "$0" "$1"'

for f in $FAMS; do
  if [ "$(ls out/$f/.ok_* 2>/dev/null | wc -l)" -eq 6 ]; then echo "$f"
  else echo "$(date '+%F %T')  SKIP assembly of $f: not all crossings succeeded" >&2; fi
done | xargs -P "$NJOBS" -n 1 bash -c 'assemble_one "$0"'
echo "$(date '+%F %T')  finished"
