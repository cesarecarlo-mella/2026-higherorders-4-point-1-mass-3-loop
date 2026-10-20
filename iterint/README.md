# iterint — iterated-integral solution (cluster step)

Builds the iterated-integral solution of every family in every crossing to weight 8, in the
format of `all.ints/iterint_solutions`. The construction is the one at the end of
`all.ints/reorganize.nb`:

    w_0 = BC_0 II[]
    w_k = BC_k II[] + Atilde_c . w_{k-1}        II[a..] l[b] -> II[a.., b]
    Atilde_c = atilde /. <c>lettersrule

The boundary conditions `BC` (ep^0..8) come from the local pipeline in `families/` and are read
from `families/results/<fam>_bc.m`; the Atildes are committed in `input/atilde/`. A clone of the
repo therefore has every input. No PolyLogTools is needed.

> ⚠️ w7 and w8 of the NPL families have not been checked independently
> (see `families/README.cesare.md`).

## Full run (all 22 families)

On the cluster, clone the repo (branch `pl-w8-crossings`), then:

    cd iterint      # logs/ is in the repo (SLURM needs it to exist before a task starts)
    sbatch run_all.sbatch                                          # 132 tasks = 22 families x 6 crossings
    sbatch --dependency=afterok:<jobid> assemble_all.sbatch        # 22 tasks, one per family

Result: `II/<fam>_II.m` (and `II/H.wl`, `II/I3.wl`, the names used in `all.ints`). Each file holds

    {BASISLC[fam, "", 1] -> Sum_{k<=8} w_k[[1]], ..., BASISLC[fam, "x132", n] -> ...}

with the crossings in the order `"", x12, x13, x23, x123, x132`, the same keys and order as the
`all.ints` files. There is no `ep`: the weight is the length of the II word plus the weight of
the constants.

- Task `t` of `run_all.sbatch` runs family line `t/6 + 1` of `families.txt` in crossing `t%6` of
  `(x x12 x13 x23 x123 x132)`. For example, H is line 16, so `sbatch --array=90-95 run_all.sbatch`
  runs only H.
- Jobs are resumable: a resubmitted task reloads the weights already in `out/<fam>/`.
- `WMAX=6 sbatch --export=ALL ...` stops at a lower weight (default 8).
- Optional check: `VALIDATE_DIR=<copy of all.ints/iterint_solutions> sbatch --export=ALL run_all.sbatch`
  compares `Sum_{k<=6} w_k` of each crossing with the old file (needs WMAX ≥ 6).
- Resources: 32 GB and 24 h per task are requested. H at weight 8 took 8 min and about 2 GB in
  use per crossing on the laptop (table below). The largest families (I1, I2, I3) have not been
  timed at weight 8.
- Size: text output at weight 8 is about 2 GB per crossing for H, so the assembled `H.wl` is
  roughly 12 GB. The assembled M3 file is 245 MB. `out/` and `II/` are gitignored (too big for
  GitHub); copy `II/` back with rsync.

## Files

- `run_iterint.wls` — one family, one crossing, weights 0..wmax → `out/<fam>/<fam>_<cross>_w<k>.wdx`.
- `assemble_II.wls` — one family, all crossings → `II/<fam>_II.m`. Writes entry by entry, so the
  whole expression is never held in memory.
- `run_all.sbatch`, `assemble_all.sbatch` — the SLURM arrays above.
- `families.txt` — the family list (task order).
- `input/atilde/<fam>_atilde.m` — copies of `all.ints/a_tildes` (`prepare_inputs.sh` refreshes them).
- `load_iterint.wl` — read the per-weight output back (`IIWeight`, `IISum`).
- `input/lettersrules.m` — the letter permutations for x12, x13, x23, x123, x132, copied from
  reorganize.nb. They are the same for every family because all families share the 20-letter
  alphabet in `input/letters.dict.m`.

## Single job

    wolframscript -file run_iterint.wls family=H atilde=input/atilde/H_atilde.m bc=../families/results/H_bc.m \
        cross=x12 wmax=8 out=out/H [text=yes] [validate=<old II file> vmax=6]

- `cross` is one of `x` (uncrossed), `x12`, `x13`, `x23`, `x123`, `x132`.
- Output is `out/H/H_<cross>_w<k>.wdx`, one file per weight, holding `{keys, S}`. WDX is
  platform-independent, so files written on Linux can be read on the Mac (MX files cannot).
- `text=yes` also writes `H_<cross>_w<k>.m` = `{BASISLC[H, c, i] -> w_k[[i]]}`.

## Cost (family H, 371 masters, measured on the laptop)

| weight | terms per crossing | time   | memory                                  |
|--------|--------------------|--------|-----------------------------------------|
| 6      | 1.5 M              | 20 s   | 1.3 GB peak                             |
| 7      | 9 M                | 90 s   | 0.4 GB in use after the step            |
| 8      | 52 M               | 8 min  | 2 GB in use after the step; peak not measured |

Family M3 (80 masters): 10 s and 0.9 GB peak per crossing at weight 8; assembly 53 s.

## Read back

    Get["load_iterint.wl"];
    IIWeight["out/H", "H", "x12", 8, -20 ;;]   (* weight-8 piece, last 20 masters *)
    IISum["out/H", "H", "x", 8, {1, 2}]        (* Sum_{k<=8}, masters 1 and 2 *)

## Checks done

- **H:** the sum through weight 6 is identical to `all.ints/iterint_solutions/H.wl` in all six
  crossings, coefficient by coefficient after expansion.
- **M3 (full pipeline: array tasks to w8, then `assemble_II.wls`):** every crossing is identical
  to `M3_II.m` through weight 6. The assembled file has the same 480 keys in the same order as
  `M3_II.m`, and truncated at total weight 6 all 480 entries are identical.
- **D** (square-root letters; the name is `System`D`): all six crossings are identical to
  `D_II.m` through weight 6. The assembled file has the same keys and order.
