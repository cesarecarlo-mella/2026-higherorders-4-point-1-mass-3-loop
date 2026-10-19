# families — boundary conditions to weight 8 for all families

Goal: produce `BC[fam, crossing, i]` to `ep^8` for the 22 families. The format is the same as
`all.ints/boundaries_all/<fam>_bc.m` (currently `ep^6`). The iterated-integral solution is
then built on the cluster with `iterint/` from those files.

Method: the H pipeline (`H_opt/run_fast.wls`), made generic. Everything is rebuilt from the
`all.ints` inputs (`a_tildes`, `bases`, shift maps), not from the old `SLC/` files.

## Folder

| file | purpose |
|------|---------|
| `config.wl` | family names (new labels ↔ legacy `SLC/` labels), input paths, family groups |
| `check_inputs.wls` → `check_inputs.log` | inventory of the inputs + coverage of the planar reductions |
| `extra_shifts/<fam>_missingPL.m` | PL integrals that need a new Kira reduction to BASISLC |
| `precompute.wls <fam>` (to do) | weight-independent data from `a_tildes/<fam>_atilde.m`: segment letter matrices, projectors at x→0, y→0, x+y→1 (every positive eigenvalue), regularity residues at x=1, y=1, x=−y, planar coefficient matrix, weight 0 |
| `solve_bc.wls <fam> <wmax>` (to do) | generic `H_opt/run_fast.wls`: segments, matching, exact boundary solve |
| `export_bc.wls <fam>` (to do) | `results/<fam>_bc.m` (ep^0..8) in the `boundaries_all` format: corners 6 ("", x23), 4 (x12, x132), 2 (x13, x123) |
| `work/<fam>/`, `results/` | generated, gitignored |

## Plan

0. **Inputs** — done, see below.
1. **Fully planar families** A, E1, E2, M3, R2. No differential equation is needed: the boundary
   conditions are the shift map × the PL crossings (`planar_stuff_w8`) at the three corners.
   Validate against `boundaries_all` through `ep^6`.
2. **Generic solver on H.** Rebuild every matrix from `a_tildes/H_atilde.m` instead of the
   `H/` files. It must reproduce `H_opt/run_fast.wls` exactly. This checks the generic
   projectors and the square-root conventions of the regularity residues.
3. **Rational alphabet** B2, C, F2, G2, G3 (letters 1–6 only, the same geometry as PL).
4. **Algebraic alphabets with complete planar data**: D, G4, I1, I3, L3 (H done in step 2).
   Before each run, check that the extra letters of the segments cancel before t→1, as they
   do for H.
5. **After the new Kira reductions / shifts**: B1, N3, P3, then F1, G1, I2.

Per-family checks: weights ≤ 6 equal to `boundaries_all`, zero residual of the
over-determined boundary system, projected equations independent of the segment variable,
only `{0,1}` words evaluated at t→1.

## Inputs (from `check_inputs.log`)

All 22 families have an atilde in `a_tildes` (letters `l[1..20]`, `letters.dict.m`), a
canonical basis in `bases`, and fully resolved `ep^6` boundary conditions in `boundaries_all`
(6 crossings each).

- **Names.** `all.ints` uses the new labels. The `SLC/` folders use the legacy labels
  (N3=B3, L3=B8, M3=A3, P3=B4, R2=A2), and their `shiftedBasis` files use their own integral
  labels (`Df` for D, `B8` for L3). The `all.ints/result` shift maps use the new labels and
  sometimes map to crossed copies of the family itself (`Cx12`, `N3x13`). Those entries are
  not planar.
- **Normalizations.** The bases agree between `all.ints` and `SLC/`; the only exception is L3
  entry 3, where the old basis has an extra `1/sqrt`. The `SLC/` Atilde agrees with `a_tildes`
  except for D and L3, which differ by ±I in the rows and columns of their square-root-normalized
  masters (`bases/rootpos.m`: D 74, L3 3). `all.ints` is therefore taken as the reference.
- **Square-root-normalized masters** (`rootpos.m`): D 74, F1 82, G1 24 71, H 156 157,
  I2 201, I3 15 203, L3 3, N3 10, P3 7 37.

### Planar reductions (task: which families need extra shifts)

Method (`check_inputs.wls`): the Kira shift map `all.ints/result/<fam>_shift.m` is applied
to `bases/<fam>_canon.m` with exact labels. The map's own label can be an alias (`DD` for D:
`D` is Mathematica's derivative, so `INT[D, ...]` was relabeled for Kira). Crossed self-copies
(`DDx12`, `Cx12`, `B1x12`, `P3x12`, `N3x13`) are kept as separate, non-planar families, and
chained rules are applied until nothing changes. Every PL integral that appears (crossing tag
dropped) is then looked up in the PL → BASISLC tables: 448 PL integrals are covered
(`result/kira_inputtoinvertPL`, `result/kira_inputextra`, `SLC/kira_inputPLtoBASISLC_5`,
`SLC/kira_inputtoPL_BASISLC`, `SLC/*/kira_inputextra`).

| family | planar entries | basis integrals not in the map | PL integrals | not in tables |
|--------|----------------|--------------------------------|--------------|---------------|
| A, E1, E2, M3, R2 | all | 0 | 92, 162, 110, 47, 97 | 0 |
| C | 76 / 173 | 80 | 66 | 0 |
| **B1** | 89 / 150 | 36 | 69 | **5** |
| **B2** | 72 / 114 | 62 | 49 | **3** |
| **D** | 54 / 121 | 43 | 34 | **1** |
| **L3** | 85 / 131 | 50 | 76 | **14** |
| **N3** | 73 / 90 | 8 | 52 | **16** |
| **P3** | 86 / 143 | 43 | 49 | **10** |
| F2, G2, G3, G4, H, I1, I3 | `SLC/Solutions/<fam>/shiftedBasis` (bases and Atilde identical to `all.ints`); planar entries 26, 25, 40, 57, 133, 77, 103 | – | – | 0 |
| **F1, G1, I2** | no planar data in `all.ints` or `SLC/Solutions` | | | |

The lists for Kira are in `extra_shifts/{B1,B2,D,L3,N3,P3}_missingPL.m`. They are a
fallback: the entries that are already completely shifted to BASISLC are probably enough
to fix the boundary system. The missing reductions only matter if a family turns out
rank-deficient.

Planar sources, decided: `all.ints/result` shift maps where they exist; otherwise
`SLC/Solutions/<fam>/shiftedBasis`. The `Leading_Color` folder is not used.
"Basis integrals not in the map" are left as the family's own integrals. For non-planar
families most of them sit in genuinely non-planar sectors, but this check cannot prove that
none of them could have been shifted to PL; only the shift finder can.
