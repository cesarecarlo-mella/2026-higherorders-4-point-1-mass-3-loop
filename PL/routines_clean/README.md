# PL pipeline — clean routines

This folder holds reorganized, lightly-commented copies of the PL
boundary-via-expansion pipeline, for inspection. The originals in `../` are
untouched and remain the ones actually run. Detailed "why" explanations that
used to live as long comment blocks in the originals are written out here
instead; the code itself only carries short comments.

Pipeline order: `precompute.wls` → `spectral_projection.wls` (once per
locus) → `word_table_builder.wls` (once, weight>6 only) →
`weight_solver.wls` (once per weight, sequentially) →
`constant_extraction.wls` (optional, post-hoc inspection).

---

## 1. `precompute.wls`

Builds every **weight-independent** matrix the rest of the pipeline needs,
from the canonical connection matrix `atilde` (`input/atilde.m`, in letters
`{x, y, 1-x, 1-y, x+y, 1-x-y}`). Output: `work/precompute_pl.mx`.

Routines:
- `loadAtilde[]` / `buildAxyTable[atilde]` — splits `atilde` into one
  matrix `Axy[letter]` per alphabet letter (the coefficient of `Log[letter]`).
- `buildDxDy[]` — full-(x,y)-plane partial-derivative connection matrices
  `dx`, `dy` (just `Sum[Axy[l] D[Log[l], var]]` over all 6 letters — no
  square-root letters in this alphabet, so no extra care needed here).
- `buildSegment1[]` / `buildSegmentGeneric[j, param, order, axyTable]` —
  for each of the 6 path segments used by the boundary-via-expansion
  algorithm, parametrizes `(x,y)` as a function of `t` (and a vanishing
  parameter `delta`), series-expands in `delta`, and reads off `dt[j] =
  d/dt` of the connection matrix along that path.
  - Segment 1 is a genuine special case, not a simplification: it keeps
    only the letters whose `delta->0` value is *literally* the symbol `t`,
    deliberately dropping letters whose value is `1-t`. This is the exact
    recipe used for the analogous family H; reconstructing it by summing
    over *all* t-dependent letters (including the `1-t` ones) is a
    different, wrong recipe and was the direct cause of an early weight-2
    mismatch.
  - Segments 2–6 share one recipe but segment 2 uses delta-series order 2
    (the others use order 1), and all of them need an inner `Apart` *before*
    substituting the quadratic root `1-2t+2t^2 -> 2(t-c2)(t-c2c)` — omitting
    either step silently gives a different (wrong) `dt[j]`.
- `buildTodelta[]` — `todelta[j]`, the weight-independent sum of `Axy`
  over whichever letters vanish as `delta->0` on segment `j`. Used later
  (in `weight_solver.wls`'s `impartSum`) for the `(-I*Pi)^k/k!` analytic
  continuation across the branch cut at each segment's degenerate endpoint.
- `buildRegularityMatrices[dx, dy]` — `dyx1`, `dxy1`, `dxmy`: the
  connection matrices for the three regularity segments (`x->1`, `y->1`,
  `x->-y`), built by direct substitution (again, no square-root letters,
  so no extra machinery needed).

---

## 2. `spectral_projection.wls`

Builds the spectral projector matrices used to pick out, from a connection
matrix's eigenspace, the specific combination of master integrals that a
boundary or regularity condition constrains. Consolidates what used to be
five near-duplicate scripts (`run_higher_eigenvalues.wls` for locus `x0`,
`run_higher_y0.wls` for `y0`, `run_higher_xy1.wls` for `xy1`, plus the
generic `run_one_locus.wls` and its debug variant `run_one_locus_inline.wls`)
into one parametrized script.

Usage: `wolframscript -file spectral_projection.wls <locus> <all|smallest>`

Routines:
- `getLocusMatrix[nm]` — extracts the single `Axy`-combination relevant to
  one locus (`x0`, `y0`, `x1`, `y1`, `xmy`, `xy1`, or `xy0`). Deliberately
  extracts **only** the one matrix needed: pulling extra, unrelated `Axy`
  matrices into the same kernel session was empirically observed to make
  the later Horner/projector computation hang, for reasons never diagnosed
  (a Mathematica session-state quirk, not a math bug).
- `getSpectrum[A]` — eigenvalues of `A`, which must come out integer
  (these are resonance/weight-shift eigenvalues); cross-checked against
  `Tr[A^n]` trace identities for `n=1..4` as a solver-independent sanity
  test.
- `projLambda[A, lam, spec]` — the spectral projector onto eigenvalue
  `lam`'s generalized eigenspace, via the standard Lagrange/Horner
  construction (`P = q(A).g(A)`), verified via `Tr[P] == multiplicity` and
  idempotency `P.P - P == 0`.
- **Mode `all`** (used for `x0`, `y0`, `xy1`): keeps *every* positive
  eigenvalue's projector as a separate list entry; `weight_solver.wls`
  later stacks these (`Join@@...`) into the full boundary-condition matrix
  (`allP_x0.m`, `allP_y0.m`, `allP_xy1.m`).
- **Mode `smallest`** (used for `x1`, `y1`, `xmy`): keeps only the
  smallest positive eigenvalue's projector as a single bare matrix — the
  regularity-condition matrix (`regX1.m`, `regY1.m`, `regXmY.m`).

> **Reproducibility note:** the `all` vs `smallest` split above is inferred
> from the *shapes* of the existing data files (`allP_x0.m` is a list of 2
> stacked 291×291 matrices; `regX1.m` is a single bare 291×291 matrix) and
> from which original scripts exist for which locus. The exact commands
> that produced `regX1.m`, `regY1.m`, `regXmY.m` were run interactively
> mid-session and were never saved as a committed script, so this
> consolidated version is a faithful reconstruction, not a byte-for-byte
> replay — if you need to regenerate these three files, verify the output
> against the existing ones (`Tr[P]`, idempotency, and ideally a direct
> comparison against the checked-in `work/reg*.m`) before trusting it.

---

## 3. `word_table_builder.wls`

Builds `work/word_table_0to8.m`: the lookup table that replaces
`ReduceCG`/`GetGs`-based fibration for evaluating `G[...,t]` at `t->1`,
needed because `GetGs`'s internal basis generator is hard-capped at
weight ≤ 6 and simply fails above it.

**The method:** plain substitution `G[...,t] /. t -> 1`, with
`$HPLAutoConvertToKnownFunctions = True` already in effect. PolyLogTools
then auto-separates the finite closed-form MZV value from an *explicit*
symbolic divergence marker `HPL[{1},1]` (representing the `log(1-t)`
singularity), which is dropped — the same regularization convention this
whole pipeline already uses elsewhere (`Log00`, the all-ones-word kill at
low weight). No PSLQ, no `ReduceCG`, no recursion needed; this is literally
the shuffle-regularization PolyLogTools does internally.

The table covers **every** sub-word length 0 through 8 (not just pure
weight-7 or weight-8 tuples), because `GIntegrate`'s output at weight 7/8
frequently factors into a *shorter* G-word (e.g. `G[1,t]`, `G[0,0,1,t]`)
times an explicit Pi/Zeta coefficient carrying the rest of the weight. A
table built only from pure-length-7 (or -8) tuples never matches these
shorter sub-words — this was the actual cause of every "Solve: 0 solutions"
failure first encountered at weight 7, not the fibration method itself.

This reconstruction was checked against the existing `work/word_table_0to8.m`
when written: **511 rules, 0 mismatches** — confirmed byte-for-byte
equivalent (as a rule set) to the table currently used in production.

---

## 4. `weight_solver.wls`

The heart of the pipeline: given `sol[w-1][seg]` for every segment, solves
for `sol[w][seg]` at weight `w`. Six sections, run in order:

**Section 1 — load inputs.** `loadPrecompute[]` pulls in the weight-
independent matrices from stage 1/2; `loadPreviousWeightSolutions[wgt]`
pulls in every `sol[k][seg]` for `k = 0..wgt-1` (weight 0 is a single
constant vector reused on every segment; weights ≥ 1 are real per-segment
files).

**Section 2 — `IntPart`.** For each of the 9 segments, integrates
`dt[seg].sol[w-1][seg]` over the path parameter via `GIntegrate`. This is
the single most expensive step in the whole pipeline. Two things protect
against losing hours of compute to a crash/kill:
- `integrateWithProgress` integrates the 291 vector components **one at a
  time** in an explicit loop (not a single opaque `Map`), printing
  progress every 20 components and flagging any component that takes
  > 5s — without this, a stuck-looking run gave zero visibility into
  *which* of the 291 components was actually slow.
- `getOrCompute` checks for an existing `work/intpart{w}_{seg}.m`
  checkpoint before recomputing, and exports immediately after computing
  (not batched at the end) — so a kill only loses the segment currently in
  progress, never ones already finished.

**Section 3 — evaluation operators.** `evalStart` (at `t->0`) is always a
direct substitution. `evalAt1` (at `t->1`) dispatches on weight, because
three genuinely different regularization methods were needed as the
alphabet got harder to handle symbolically:
- **weight ≤ 2:** kill only the literal all-ones word
  (`G[1,1,...,1,t]->0`), then substitute `t->1` directly — the preloaded
  MZV table handles every other word. Using the heavier fibration method
  (next bullet) already at this weight gives a *numerically different,
  wrong* answer, not just an unnecessary slowdown — it is not simply "more
  general."
- **3 ≤ weight ≤ 6:** full fibration-change: `t -> 1-t`, reduce via the
  `tto1Red` table (extended on demand by `extendTto1Red`, since PL's
  alphabet introduces words beyond H's frozen 126-rule table), `ShuffleG`
  to canonicalize, then kill both the literal-`t`-ending words (`G[__,t]
  -> 0`) **and** the degenerate all-zero constants that already lost their
  `t`-dependence mid-reduction (`Log00`) — missing this second kill was the
  weight-3 bug: those degenerate constants leaked into every downstream
  delta vector, not just one.
- **weight ≥ 7:** `ReduceCG`/`GetGs` simply fails above weight 6, so
  instead apply the precomputed `word_table_0to8.m` (section above)
  directly.

**Section 4 — matching.** Builds the delta vectors `d2..d6` (differences
between adjacent segments' endpoint values — all computable from
*known*, previous-weight-only data) and `d00, d01, d10` (connections to
the `x1`/`y1`/`xmy` regularity segments via `impartSum`, the
`(-I*Pi)^k/k!`-weighted analytic continuation across each segment's branch
point).

**Section 5 — solve.** Expresses every segment's unknown integration
constant symbolically in terms of the single unknown vector `CNST5`
(291 constants for the base segment), builds the boundary equations
(`stackX0`, `stackY0`, `stackXy1`) and regularity equations (`resX1`,
`resY1`, `resXmY`), adds the one external input this family has — the
exact value of the fully massless banana integral at position 1 — and
solves the full linear system for `CNST5`.

**Section 6 — export.** Back-substitutes `CNST5`'s value into every
segment's constant, then applies `ShuffleG` to each exported solution
**before** writing it to disk. This `ShuffleG` pass is not cosmetic: raw
exported solutions were never canonicalized before, letting redundant
shuffle-equivalent G-word representations compound weight-over-weight. By
weight 7, `sol7_y1.m` had ballooned to 2.33M `LeafCount` (975 distinct
G-words) versus 327 words at weight 6 — and this bloat, not `GIntegrate`
itself, was the real cause of weight 8's `IntParty1` step running 2+ hours
without finishing. `ShuffleG` is value-preserving and gave up to a 25×
size reduction when applied retroactively to the weight-1..7 files.

---

## 5. `constant_extraction.wls`

Post-hoc inspection utility: for every weight 0..8 and segments `{2,4,6}`,
imports the solution and sets every `G`-polylog to zero (`/. _G -> 0`),
keeping just the pure-constant (Pi/Zeta) part. Writes
`work/constant/const{w}_{i}.m`. Weight 0 gives the same file for all three
segments by construction (weight 0 has no segment dependence). Weight 1's
constants are all identically zero — the weight-1 solution is pure
log/G-words with no constant piece, which is expected at that order.

---

## Validation status (as of weight 8)

- Weights 1–6: validated against the independent reference `planar_stuff/pl.m`
  — 0/291 mismatches at every weight, both before and after the retroactive
  `ShuffleG` cleanup.
- Weights 7–8: no external ground truth exists beyond weight 6, so these
  are validated via internal consistency instead — segment 6's entry 1,
  evaluated at `t->1` (which kills the `G[0,...,0,t] = Log[t]^k/k!` terms
  since `Log[1]=0`), reduces to exactly the known banana boundary value at
  both weights.
