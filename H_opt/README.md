# H_opt — optimized boundary-constant pipeline for family H

Goal: reproduce **exactly** what `H/solution.nb` (section *Boundary via
Expansion*) produces for family H, weights 0–6, but much faster and with a small
memory footprint. `H/` is the reference and must **not** be modified.

This file is a briefing for whoever (human or Claude Code) works in this folder.
It was written by reading the notebook; statements marked **(verify)** are
hypotheses that still need to be checked numerically/symbolically.

`reference/solution_nb_inputs.wl` = all *Input* cells of `H/solution.nb`
converted to plain text (auto-extracted from the box structure; read it instead of
the .nb, but if something looks odd compare with the notebook).

---

## 0. Ground rules (laptop safety + repo hygiene)

* **One kernel only.** No `LaunchKernels`, no `Parallel*`.
* Every script: `$HistoryLength = 0;` and wrap the main body in
  `MemoryConstrained[..., 8*1024^3]` (8 GB). Wrap expensive single steps in
  `TimeConstrained` and print timings (`AbsoluteTiming`) for every step.
* Run in background with low priority and a log, e.g.
  `nice -n 10 wolframscript -file H_opt/run_w3.wls > H_opt/logs/w3.log 2>&1 &`
* **Checkpoint** after every weight and every segment (`DumpSave` to
  `H_opt/work/…mx` or `Export` .m). A script must be restartable from the last
  checkpoint; never recompute something that is on disk.
* Scratch goes to `H_opt/work/` and `H_opt/logs/` (git-ignored). Before Cesare
  pushes: delete anything not needed; only scripts, this README and final small
  results should be committed.
* Never write inside `H/`. Read-only from there.

## 1. Physics / notation

* Family H: 3-loop, 4-point, 1 massive leg; 371 master integrals in canonical
  form, `dJ = ε dÃ J`, `Ã = Σ_i Axy[i] Log[letter_i]` (`H/input/Atilde.m`,
  variables `s→x, t→y, m2→1`), 14 letters `arglogs` incl. two square-root letters
  (`Sqrt[-x(1-x-y)y] → I Sqrt[x(1-x-y)y]`).
* Basis: `H/input/basis.m`. `onlyplanar` = positions (133 of them) whose
  integrals are fully expressible through planar ones
  (`Position[Map[FreeQ[#, INT[b__,a__]]&, shiftedBasis], True]` with
  `shiftedBasis = Import["H/shiftedBasis/shiftedBasis.m"]`).
* Planar input: `SolutionWithPlanar = Import["H/results/planarInsertedXY.m"]`
  (ε-expanded GPLs in x,y for the 133 planar entries). **Frozen input.**
* Unknowns at weight w: `cnst[w,i][seg]`, i = 1..371, seg ∈ {1..6, x1, y1, xmy}.
  After using planar info, 238 unknowns remain on segment 5.

### Paths ("segments") — weight-independent data
Defined in notebook subsection *segments*; `δ→0`, `t∈(0,1)`:

| seg | parametrisation | role |
|----|----|----|
| 1 | x→t, y→δ | boundary eq. from P_1(toexp2) (y→0) |
| 2 | y→δt², x→1−δ((1−t)²+t²) | transport |
| 3 | x→1−δ−t, y→t | boundary eq. from P_1(toexp7) (x+y→1); analytic cont. |
| 4 | x→δ(1−t)², y→1−δ((1−t)²+t²) | transport |
| 5 | x→δ, y→1−t | **base segment**: all constants transported here; P_1(toexp1) (x→0); planar input |
| 6 | x→δt², y→δ(1−t)² | transport (closes the loop) |
| x1 | x=1, integrate in y (`dyx1`) | regularity `regCondxy[1]` |
| y1 | y=1, integrate in x (`dxy1`) | regularity `regCondxy[2]` |
| xmy | x→−y, integrate in y (`dxmy`) | regularity `regCondxy[3]` |

`dt[j]` = pull-back connection on segment j (a 371×371 matrix of rational
functions of t), `todelta[j]` = residue matrix of the δ→0 direction (used for the
analytic continuation `δ^{εA} → e^{-iπεA}` in `impart`).

### Boundary matrices — weight-independent, frozen
* `boundaryMxy[1] = H/matrix_to_JD/toJD/P_1_toexp1.m` (x→0)
* `boundaryMxy[2] = …/P_1_toexp2.m` (y→0)
* `boundaryMxy[3] = …/P_2_toexp3.m` (x,y→0) — loaded but not used in the segment pipeline
* `boundaryMxy[7] = …/P_1_toexp7.m` (x+y→1)
* `regCondxy[k] = H/boundary_matrices/breg{k}.m`, k=1,2,3 (numerical, contain I)

The P's are spectral projectors computed **without Jordan decomposition**, see
`H/matrix_to_JD/projector_method.pdf` (P_λ = q_λ(A) g_λ(A), exact, checked by
Tr P = a_λ and P² = P). Reuse them, do not recompute.

## 2. What the notebook does per weight w (w = 1..6)

1. **Planar constants on seg 5** (`plrules[w][5]`): coefficient of ε^w of the
   planar entries, `x→0`, drop `G[0,…,0]` (`Log00`), `y→1−t`, apply `tto1Red`
   (change of fibration `G[…,1−t] → G[…,t]`, from `ReduceCG`), `ShuffleG`, then
   `_G→0`. → 133 constants fixed.
2. **Integrate on each segment j=1..6**:
   `tempSol[j] = GIntegrate[#, t]& /@ Apart[dt[j].sol[w-1][j], t] + CNST[w][j]`.
   Same for x1 (in y, `dyx1`), y1 (in x, `dxy1`), xmy (in y, `dxmy`), followed by
   `Collect[..., FullSimplify[# /. roots[xmy]]&]` (cube roots of unity).
3. **Matching**: `match_k` = value at t→1 of segment k−1 (via `t→1−t`, `tto1Red`,
   `ShuffleG`, `Logt0`) minus value at t→0 of segment k. `match00/01/10` connect
   x1, y1, xmy to segments 5/3 using `impart = Σ_k (−iπ)^k/k! todelta^k . sol[w−k]`.
4. **Boundary equations**: `eq1 = P1.tempSol[5]`, `eq2 = P2.tempSol[1]`,
   `eq7 = P7.tempSol[3]`, `eq4/5/6 = regCondxy[1/2/3].tempSol[x1/y1/xmy]`.
5. **Transport** all constants to segment 5 with repeated `Solve` on the 371-dim
   matching systems (`c1to5`, …, `xmyto5`), substitute into the equations.
6. **Solve** the stacked system for `cnst[w,i][5]`
   (`H/numerical_eq/solve_numerical_equations.wls`, RowReduce on
   `CoefficientArrays`), then back-substitute to all segments and export
   `H/solutions/eps_sep/w{w}/sol{w}_{seg}.m` and `const{w}_{seg}.m`.

## 3. Where the time goes, and what to do instead

Ordered by expected gain. Each item must be validated against the reference
before moving on.

### 3.1 Replace `Apart` + `GIntegrate` by a pure index operation (biggest win)
The segment solutions are **pure**: every entry of `sol[w][j]` is
`Σ (number) × G[a_1,…,a_w', t]` with number ∈ Q[π, ζ, I, …] (checked on the
reference files at w6: segments 1–6 only contain letters {0,1}; x1 and y1 use
{0, ±1, 1/2, 2, −cxmy, −cxmyc}; xmy uses {0, ±1, ±2, ±cxmy, ±cxmyc}, with
cxmy, cxmyc the cube roots of unity from `roots[xmy]`). On each segment the connection is
d-log in t with *constant* matrices:

    dt[j] = Σ_a  M_a^{(j)} / (t − a)        (a ∈ letters of segment j)

Compute the `M_a^{(j)}` **once** (weight-independent; e.g. `Apart` of the
matrix once, or residues), store them as `SparseArray`. Then, since
∫_0^t dt'/(t'−a) G(w⃗; t') = G(a, w⃗; t),

    tempSol[w][j] = Σ_a M_a^{(j)} . Prepend_a(sol[w−1][j]) + CNST[w][j]

where `Prepend_a` maps `G[w⃗,t] → G[a,w⃗,t]` and a constant `c → c G[a,t]`.
No `Apart`, no `GIntegrate`, no `Together`. **(verify)** that this reproduces
`GIntegrate` exactly at w=1,2 (including the trailing-zero / `G[0,…,0,t]`
conventions of PolyLogTools).

### 3.2 Use a linear-algebra representation instead of big expressions
Represent `sol[w][j]` as `{words, C}`: a list of distinct G-words and a sparse
matrix C (371 × #words) of coefficients (coefficients themselves expanded in a
basis of transcendental constants {1, π², ζ3, …} if useful). Then
* prepending a letter = relabel words, matrix product with M_a = sparse dot;
* evaluation at t→0 (regularised) = keep only the empty word;
* evaluation at t→1 = contract C with a precomputed vector of values of the words
  at t=1 (from `tto1Red`, see §4);
* Collect/Together/Simplify on huge expressions disappear.

### 3.3 Matching without `Solve`
Each matching equation has the form `c_k − c_l + Δ = 0` (identity coefficient on
the unknowns; also true for `match00/01/10`, where the unknowns only enter through
`tempSol`). So `c_l = c_k + Δ` is a vector subtraction. All the
`Flatten[Solve[...]]` calls (and the unused "inverse relations" `c5to*`) go away.
`eq3top5` in the notebook references an undefined `eq3` — drop it.

### 3.4 Boundary equations: evaluate first, multiply few rows, invert once
* The exported systems (`H/numerical_eq/numerical_equations_w.m`) are purely
  numerical: t-dependence cancels in `P.tempSol`. Since regularised evaluation is
  linear and acts as the identity on constants, `P.F(t) = P.Reg_{t→0}F`. So
  evaluate `tempSol` at the regularised point (a vector of numbers) **before**
  multiplying by P / regCond. **(verify)** at w=1,2 that the resulting system is
  identical to the reference one.
* The unknowns enter with the **same coefficient matrix at every weight** (rows of
  P1, P2, P7, regCond1–3 after transport-by-identity and removal of the 133
  planar columns). So: once, select a set of 238 independent rows (exact
  `RowReduce`/rank on rationals with I), compute its inverse; at each weight only
  compute the RHS **for those rows** and multiply. Check the remaining rows
  numerically (cheap) as a consistency test.
* P matrices are dense 371×371 but low rank — consider using only pivot rows
  from the start.

### 3.5 Powers of `todelta` and analytic continuation
`impart = Σ_k (−iπ)^k/k! todelta^k . sol[w−k]`: precompute the sparse powers
`todelta[j]^k`, k ≤ 6, once.

### 3.6 Algebraic numbers on x1/y1/xmy
Instead of `FullSimplify[# /. roots[xmy]]` on every coefficient, keep
coefficients in an explicit basis of Q(ω) (ω = cube root of unity) or use
`RootReduce` on the (few) distinct coefficient values.

## 4. Out of scope for now: change of fibration (`ReduceCG`)
`tto1Red` (rules `G[…,1−t] → G[…,t]`) is **frozen input**. In the notebook it is
loaded from `../../change_fibration/others/fibrtto1.m` (outside this repo) and
then *extended in-session* by `ReduceCG` at w=3,4,5,6 — those extensions were
**never exported**. Therefore:
* Cesare: copy `fibrtto1.m` to `H_opt/input/fibrtto1.m`.
* The pipeline must, at each weight, list the `G[…,1−t]` that survive
  `/. tto1Red` and report how many are missing. If any are missing, stop and ask
  Cesare (either he provides the extended file, or run `ReduceCG` only on the
  missing ones and cache the result to `H_opt/work/tto1Red_w{w}.m`).
* Better long-term: the t→1 values of the words only need to be computed once per
  word (vector of constants), see §3.2.

Also frozen: `SolutionWithPlanar`, `shiftedBasis`, the P/regCond matrices,
`H/input/*`.

## 5. Validation target ("reproduced" = exact match)
For each weight w = 0..6 and each seg ∈ {1..6, x1, y1, xmy}:
* `sol[w][seg]` equals `H/solutions/eps_sep/w{w}/sol{w}_{seg}.m`
  (difference simplifies to exactly 0 after `ShuffleG`/canonical G form; compare
  entry by entry, report the first mismatch).
* constants equal `H/solutions/eps_sep/w{w}/const{w}_{seg}.m` and
  `H/numerical_eq/solution_{w}.m`.
* Report a timing table per step and per weight vs. the original notebook.

Validate w=1,2 first (fast), then go up one weight at a time.

## 6. Suggested layout
    H_opt/
      README.md                 this file
      reference/solution_nb_inputs.wl
      input/fibrtto1.m          (copied by Cesare, frozen)
      lib/                      .wl package(s): segment data, word algebra, solver
      run_precompute.wls        weight-independent data (M_a, todelta^k, pivots, inverse)
      run_weight.wls            `wolframscript -file run_weight.wls <w>`
      check_weight.wls          comparison with H/solutions/eps_sep/w<w>
      work/  logs/              git-ignored
