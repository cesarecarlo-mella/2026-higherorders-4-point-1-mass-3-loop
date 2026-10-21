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

## Status (2026-10-21)

**Done for all 22 families, weights 0..8, all six crossings.** The results are on the cluster
only (not in git):

    newton4:/scratch/ge45tuy/2026-higherorders-4-point-1-mass-3-loop/iterint/II/<fam>_II.m.gz   assembled, gzipped (~95 GB unzipped)
    newton4:/scratch/ge45tuy/2026-higherorders-4-point-1-mass-3-loop/iterint/out/<fam>/          per-weight .wdx (load_iterint.wl)

The assembled files from this run are named H.wl / I3.wl, renamed by hand to `H_II.m` / `I3_II.m`.

**What this run checked:** nothing against the old results. It ran without `VALIDATE_DIR`, so
the only checks were internal to each job:
- the Atilde is linear in the letters;
- every BC entry is present;
- every crossing reached w8;
- each assembled file has 6 × masters entries.

**Compared with `all.ints` through w6 earlier, locally:** H, M3, D and R2 (see "Checks done"
below). The other 18 families have not been compared at w≤6 on these files. w7 and w8 have no
check at all.

Weight-8 cost per crossing (newton4):

| family | memory peak | time |
|---|---|---|
| H | 38 GB | 16 min |
| I1 | 31 GB | 1.8 h |
| G1 | 22 GB | 15 min |

Assembly took up to 1 h each for H and I3.

On newton4, launch with `setsid` (`nohup` alone is not enough: the kernels die when the SSH
session closes):

    setsid nohup env WL="math -script" ./run_all.sh 8 > logs/run_all.log 2>&1 < /dev/null &

## Full run (all 22 families)

On the cluster, clone the repo (branch `pl-w8-crossings`), then, with no scheduler:

    cd iterint
    setsid nohup ./run_all.sh 8 > logs/run_all.log 2>&1 < /dev/null &   # 132 jobs (22 families x 6 crossings), 8 at a time
    tail -f logs/run_all.log                            # one line per finished job / assembled family

`run_all.sh NJOBS [families...]` runs the (family, crossing) jobs of `run_iterint.wls` as
background processes, NJOBS at a time, largest families first. When they are done it assembles
every family whose six crossings all succeeded. Result: `II/<fam>_II.m` for every family
(`all.ints` names H and I3 `H.wl`, `I3.wl`). Each file holds

    {BASISLC[fam, "", 1] -> Sum_{k<=8} w_k[[1]], ..., BASISLC[fam, "x132", n] -> ...}

with the crossings in the order `"", x12, x13, x23, x123, x132`, the same keys and order as the
`all.ints` files. There is no `ep`: the weight is the length of the II word plus the weight of
the constants.

- **Mathematica command:** `WL` (default `wolframscript -file`). Where there is no
  `wolframscript`, use the kernel directly, e.g. `WL="math -script" ./run_all.sh 8` or
  `WL="/path/to/WolframKernel -script" ./run_all.sh 8` (both tested), or `module load` Mathematica
  first.
- **NJOBS:** keep it within the number of Mathematica kernels the license allows and within the
  machine's memory.
- **Families:** `./run_all.sh 4 H I3` runs only those families.
- **Resumable:** a rerun reloads the weights already in `out/<fam>/`.
- **Kernel watchdog:** a kernel that does not start within `STARTUP` seconds (default 180), for
  example with no free license or a failed link, is killed and the job retried up to 3 times.
- **Logs:** `logs/<fam>_<cross>.log` per job, `logs/<fam>_assemble.log` per family, and the
  summary `logs/run_all.log` (`done` / `FAILED` per job).
- **Lower weight:** `WMAX=6 ./run_all.sh ...` stops at a lower weight (default 8).
- **Optional check:** `VALIDATE_DIR=<copy of all.ints/iterint_solutions> ./run_all.sh ...` also
  compares `Sum_{k<=6} w_k` of every crossing with the old file. A mismatch marks the job FAILED.
- **Resources:** H at weight 8 took 8 min and about 2 GB in use per crossing on the laptop
  (table below). The largest families (I1, I2, I3) have not been timed at weight 8.
- **Size:** text output at weight 8 is about 2 GB per crossing for H, so the assembled `H_II.m` is
  roughly 12 GB. The assembled M3 file is 245 MB. `out/` and `II/` are gitignored (too big for
  GitHub); copy `II/` back with rsync.

## Files

- `run_iterint.wls` — one family, one crossing, weights 0..wmax → `out/<fam>/<fam>_<cross>_w<k>.wdx`.
- `assemble_II.wls` — one family, all crossings → `II/<fam>_II.m`. Writes entry by entry, so the
  whole expression is never held in memory.
- `run_all.sh` — the launcher above (no scheduler).
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
- **M3 (to w8, then `assemble_II.wls`):** every crossing is identical
  to `M3_II.m` through weight 6. The assembled file has the same 480 keys in the same order as
  `M3_II.m`, and truncated at total weight 6 all 480 entries are identical.
- **D** (square-root letters; the name is `System`D`): all six crossings are identical to
  `D_II.m` through weight 6. The assembled file has the same keys and order.
- **`run_all.sh` with WMAX=6 on M3 and R2:** the assembled `II/M3_II.m` and `II/R2_II.m` are
  identical to the `all.ints` files (keys, order and every entry).
