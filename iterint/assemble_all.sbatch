#!/bin/bash
# After run_all.sbatch: one task per family (line t+1 of families.txt) writing II/<fam>_II.m
# (II/H.wl, II/I3.wl), the all.ints/iterint_solutions format, weights 0..WMAX.
#
#   cd iterint && sbatch --dependency=afterok:<run_all job id> assemble_all.sbatch
#
#SBATCH --job-name=assemble_II
#SBATCH --array=0-21
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -e
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
WMAX=${WMAX:-8}
fam=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" families.txt)
wolframscript -file assemble_II.wls family="$fam" in="out/$fam" outdir=II wmax="$WMAX"
