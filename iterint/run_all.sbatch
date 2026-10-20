#!/bin/bash
# Iterated-integral solution of all families, all crossings: one array task per
# (family, crossing), 22 x 6 = 132 tasks. Task id t -> family line t/6 + 1 of families.txt,
# crossing t%6 of (x x12 x13 x23 x123 x132). Resumable: a resubmitted task reloads the weights
# already in out/<fam>/.
#
#   git clone <repo>; cd <repo>/iterint && mkdir -p logs && sbatch run_all.sbatch   # everything
#   sbatch --array=90-95 run_all.sbatch                    # only H (line 16 -> tasks 90..95)
#   VALIDATE_DIR=/path/iterint_solutions sbatch --export=ALL run_all.sbatch   # + check w<=6
#
# Then, per family:  wolframscript -file assemble_II.wls family=H
#
#SBATCH --job-name=iterint
#SBATCH --array=0-131
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -e
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs
WMAX=${WMAX:-8}
C=(x x12 x13 x23 x123 x132)
t=${SLURM_ARRAY_TASK_ID}
fam=$(sed -n "$((t / 6 + 1))p" families.txt)
cross=${C[$((t % 6))]}
echo "task $t: family $fam crossing $cross wmax $WMAX on $(hostname)"
# optional: VALIDATE_DIR=<copy of all.ints/iterint_solutions> compares Sum_{k<=6} with the old file
extra=""
if [ -n "$VALIDATE_DIR" ]; then
  if [ "$fam" = H ] || [ "$fam" = I3 ]; then old="$VALIDATE_DIR/${fam}.wl"; else old="$VALIDATE_DIR/${fam}_II.m"; fi
  extra="validate=$old vmax=6"
fi
wolframscript -file run_iterint.wls family="$fam" atilde="input/atilde/${fam}_atilde.m" \
  bc="../families/results/${fam}_bc.m" cross="$cross" wmax="$WMAX" out="out/$fam" $extra
