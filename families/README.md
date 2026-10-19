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

## Status (2026-10-19)

Boundary conditions to ep^8 (`results/<fam>_bc.m`) and segment solutions
(`segments/<fam>/eps_sep/w0..w8/{sol,const}<w>_<seg>.m`, H's format, gitignored, ~1 GB):

- **done**: A, B1, B2, C, D, E1, E2, F1, F2, G1, G2, G3, G4, H, I1, I2, L3, M3, N3, P3, R2;
- **running**: I3 (weights 1-7 done).

Every family: ep^1..6 identical to `all.ints/boundaries_all` for all six crossings; the
over-determined projector + regularity system has full rank and zero residuals at every weight,
7 and 8 included. H additionally reproduces `H_opt` exactly at every weight 0-8. The planar
families go through the same solver with no unknowns, so all their equations are checks; their
BC equal the planar route (`planar_bc.wls`).

## Scripts

- `solve_bc.wls fam=<f> [wmax=8]`: the H algorithm for any family, built per letter from
  `a_tildes` and the 20-letter dictionary (H's sqrt-letter formulas and branches); projectors on
  all positive eigenvalues at x->0, y->0, x+y->1; regularity residues at x=1, y=1, x=-y; planar
  constants from the shifted basis x PL crossings (`PL/pl_gpls_7_8`). Checkpoints in `work/<f>/`.
- `export_segments.wls fam=<f,...>`: segment solutions in H's `eps_sep` format.
- `planar_bc.wls`: BC of the planar families directly from shifts x PL (no DE).
- `export_shifted.wls`: readable shifted bases `results/<fam>_shiftedBasis.m`.
- `check_inputs.wls`, `config.wl`: inputs, names, planar sources.

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
| **B1** | 89 / 150 | 36 | 69 | **3** |
| **B2** | 72 / 114 | 62 | 49 | **3** |
| **D** | 54 / 121 | 43 | 34 | **1** |
| **L3** | 85 / 131 | 50 | 76 | **12** |
| **N3** | 73 / 90 | 8 | 52 | **15** |
| **P3** | 86 / 143 | 43 | 49 | **9** |
| F2, G2, G3, G4, H, I1, I3 | `SLC/Solutions/<fam>/shiftedBasis` (bases and Atilde identical to `all.ints`); planar entries 26, 25, 40, 57, 133, 77, 103 | – | – | 0 |
| F1 | 53 / 214 | 167 | 63 | 0 |
| G1 | 114 / 254 | 161 | 119 | 0 |
| **I2** | 116 / 305 | 115 | 92 | **1** |

The lists for Kira are in `extra_shifts/{B1,B2,D,I2,L3,N3,P3}_missingPL.m`. "Planar entries" counts every entry whose shift reaches only PL integrals, reducible or not. They are a
fallback: the entries that are already completely shifted to BASISLC are probably enough
to fix the boundary system. The missing reductions only matter if a family turns out
rank-deficient.

Planar sources, decided: `all.ints/result` shift maps where they exist; otherwise
`SLC/Solutions/<fam>/shiftedBasis`. F1, G1 and I2 have no shift map in `all.ints`: their Kira
shifts (`toshift/list.m`) and PL -> BASISLC tables come from
`DEQ_solutions/Leading_Color/<fam>/*_rotated_boundaries/toshift/`. The propagators are the same
as in `all.ints`: identical for I2 (Kira `integralfamilies.yaml`); for F1 and G1 the
`Leading_Color` `old_basis` uses the same integral index vectors as `all.ints` with proportional
coefficients. Nothing else from `Leading_Color` is used. (The rotated `shiftedBasis` there gives
the same planar entries up to an overall eps^6 normalization.)
"Basis integrals not in the map" are left as the family's own integrals. For non-planar
families most of them sit in genuinely non-planar sectors, but this check cannot prove that
none of them could have been shifted to PL; only the shift finder can.
