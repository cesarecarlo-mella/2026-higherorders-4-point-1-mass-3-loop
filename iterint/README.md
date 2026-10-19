# iterint — iterated-integral solution (cluster step)

Builds the iterated-integral solution of one family in one crossing, as at the end of
`all.ints/reorganize.nb`:

    w_0 = BC_0 II[]
    w_k = BC_k II[] + Atilde_c . w_{k-1}        II[a..] l[b] -> II[a.., b]
    Atilde_c = atilde /. <c>lettersrule

The boundary conditions `BC` come from the local boundary pipeline (for H: `H_opt/run_fast.wls`
then `H_opt/export_fast.wls`, which writes `H_opt/results/H_bc.m`). This folder depends on
nothing else: no PolyLogTools and no repo paths. Copy `iterint/` plus the two input files to
the cluster.

## Files

- `run_iterint.wls` — one family, one crossing, weights 0..wmax.
- `load_iterint.wl` — read the output back (`IIWeight`, `IISum`).
- `input/lettersrules.m` — the letter permutations for x12, x13, x23, x123, x132, copied from
  reorganize.nb. They are the same for every family because all families share the 20-letter
  alphabet in `input/letters.dict.m`.

## Run

    wolframscript -file run_iterint.wls family=H atilde=H_atilde.m bc=H_bc.m \
        cross=x12 wmax=8 out=out/H [text=yes]

- `cross` is one of `x` (uncrossed), `x12`, `x13`, `x23`, `x123`, `x132`; run one job per crossing.
- Output is `out/H/H_<cross>_w<k>.wdx`, one file per weight, holding `{keys, S}`. WDX is
  platform-independent, so files written on Linux can be read on the Mac (MX files cannot).
  A rerun resumes from the weights already on disk.
- `text=yes` also writes `H_<cross>_w<k>.m` = `{BASISLC[H, c, i] -> w_k[[i]]}`. The old
  `<fam>_II.m` files are the sum over k of these.
- `validate=<old II file> vmax=6` compares `Sum_{k<=6} w_k` with an existing `<fam>_II.m`.

SLURM job array (6 crossings):

    #SBATCH --array=0-5
    #SBATCH --mem=16G
    C=(x x12 x13 x23 x123 x132)
    wolframscript -file run_iterint.wls family=H atilde=H_atilde.m bc=H_bc.m \
        cross=${C[$SLURM_ARRAY_TASK_ID]} wmax=8 out=out/H

## Cost (family H, 371 masters, measured on the laptop)

| weight | terms per crossing | time   | memory                                  |
|--------|--------------------|--------|-----------------------------------------|
| 6      | 1.5 M              | 20 s   | 1.3 GB peak                             |
| 7      | 9 M                | 90 s   | 0.4 GB in use after the step            |
| 8      | 52 M               | 8 min  | 2 GB in use after the step; peak not measured |

The weight-8 peak was not recorded; requesting 16 GB per job leaves plenty of headroom.
Only `w_{k-1}` is kept in memory. As text, weight 8 is roughly 2 GB per crossing.

## Read back

    Get["load_iterint.wl"];
    IIWeight["out/H", "H", "x12", 8, -20 ;;]   (* weight-8 piece, last 20 masters *)
    IISum["out/H", "H", "x", 8, {1, 2}]        (* Sum_{k<=8}, masters 1 and 2 *)

## Checks done

The sum through weight 6 is identical to `all.ints/iterint_solutions/H.wl` for all six
crossings, comparing coefficient by coefficient after expansion.
