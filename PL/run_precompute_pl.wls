$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

Print["=== step 1: load atilde (x,y form, /ep stripped) ==="];
{t0, atilde} = AbsoluteTiming[(Import["input/atilde.m"] /. z -> 1 - x - y)/ep];
Print["atilde dims: ", Dimensions[atilde], "  t=", t0, "s"];

letters = {x, y, 1 - x, 1 - y, x + y, 1 - x - y};
Do[Axy[letters[[i]]] = Coefficient[atilde, Log[letters[[i]]]], {i, 1, Length[letters]}];
Print["Axy[letter] built for all 6 letters"];

(* dx, dy: the x,y-partial-derivative connection matrices. No sqrt letters, so this
   is just a direct sum -- no rrs/revertRRS dance needed (that was H/A-specific). *)
dx = Sum[Axy[l] D[Log[l], x], {l, letters}];
dy = Sum[Axy[l] D[Log[l], y], {l, letters}];
Print["dx, dy built"];

(* === segments 1-6, same parametrizations as H_opt === *)
Print["=== step 2: segment residue-decomposition (M_a matrices) ==="];

(* helper: given arglogs/param, find M_a^(j) = residue-decomposed dt[j] as {denom -> matrix} rules *)
buildSegment[paramRule_, tvar_] := Module[{argsub, dtfull, aparts, denoms, result},
   argsub = letters /. paramRule;
   dtfull = Sum[Axy[letters[[i]]] D[Log[argsub[[i]]], tvar], {i, Length[letters]}];
   dtfull = Together[dtfull];
   result = dtfull;
   result
   ];

(* segment 1: x->t, y->delta. EXACT reference recipe (Cesare's original H notebook
   code): ser1 = the delta->0 VALUE of each letter (via Series to order delta^1,
   Normal, THEN /.delta->0 -- NOT differentiate-then-limit). dt[1] then sums
   Axy[i].D[Log[ser1[[i]]],t] ONLY over letters whose delta->0 value is literally
   the symbol t (selected via Position[ser1,t]) -- this deliberately EXCLUDES
   letters whose delta->0 value is "1-t" (a different, real structural choice,
   not an approximation/limit-order detail). My earlier reconstruction summed
   over ALL t-dependent letters including the 1-t ones, which is NOT the same
   recipe and was the real source of the weight-2 discrepancy. *)
param[1] = {x -> t, y -> delta};
ser1 = (Series[letters /. param[1], {delta, 0, 1}, Assumptions -> 0 < t < 1] // Normal) /. delta -> 0;
(* project every match (direct or nested, e.g. "t" inside "1-t") down to its outer
   list index, matching the reference's {a_,b_,c_}->{a} trick but robust to any
   nesting depth *)
ser1Positions = Union[First /@ Position[ser1, t]];
Print["segment 1: letters whose delta->0 value is/contains t: positions ", ser1Positions, " (", letters[[#]] & /@ ser1Positions, ")"];
dt[1] = Sum[Axy[letters[[i]]] D[Log[ser1[[i]]], t], {i, ser1Positions}];
Print["dt[1] built"];

(* segments 2-6: exact reference recipe --
   dt[j] = Apart[Apart[(D[Log[#]&/@ser[j],t]//Simplify /.delta->0//Simplify)
                  /. 1-2t+2t^2 -> 2(t-c2)(t-c2c), t] . AxyTable, t]
   -- note the INNER Apart (before the quadratic substitution) and the
   quadratic-root substitution itself, BOTH omitted in an earlier version of
   this script. The quadratic 1-2t+2t^2 arises from the parametrization's own
   algebra (delta*t^2, delta*(1-t)^2 type substitutions), independent of
   whether the alphabet has sqrt letters -- PL needs this exactly like H does.
   Also: segment 2 uses Series order delta^2 (not ^1 like the others) in the
   reference -- replicated exactly. *)
c2 = (1 - I)/2; c2c = (1 + I)/2;
axyTable = Table[Axy[letters[[i]]], {i, Length[letters]}];

(* segment 2: y->delta t^2, x-> 1-delta((1-t)^2+t^2) *)
param[2] = {y -> delta t^2, x -> 1 - delta ((1 - t)^2 + t^2)};
ser2 = Series[letters /. param[2], {delta, 0, 2}, Assumptions -> {delta > 0, 0 < t < 1}] // Normal;
dlog2 = ((D[Log[#], t] & /@ ser2) // Simplify) /. delta -> 0 // Simplify;
dt[2] = Apart[Apart[dlog2 /. 1 - 2 t + 2 t^2 -> 2 (t - c2) (t - c2c), t] . axyTable, t];
Print["dt[2] built"];

(* segment 3: x-> 1-delta-t, y->t *)
param[3] = {x -> 1 - delta - t, y -> t};
ser3 = Series[letters /. param[3], {delta, 0, 1}, Assumptions -> {delta > 0, 0 < t < 1}] // Normal;
dlog3 = ((D[Log[#], t] & /@ ser3) // Simplify) /. delta -> 0 // Simplify;
dt[3] = Apart[Apart[dlog3 /. 1 - 2 t + 2 t^2 -> 2 (t - c2) (t - c2c), t] . axyTable, t];
Print["dt[3] built"];

(* segment 4: x-> delta(1-t)^2, y-> 1-delta((1-t)^2+t^2) *)
param[4] = {x -> delta (1 - t)^2, y -> 1 - delta ((1 - t)^2 + t^2)};
ser4 = Series[letters /. param[4], {delta, 0, 1}, Assumptions -> {delta > 0, 0 < t < 1}] // Normal;
dlog4 = ((D[Log[#], t] & /@ ser4) // Simplify) /. delta -> 0 // Simplify;
dt[4] = Apart[Apart[dlog4 /. 1 - 2 t + 2 t^2 -> 2 (t - c2) (t - c2c), t] . axyTable, t];
Print["dt[4] built"];

(* segment 5: x->delta, y->1-t (the base segment) -- same general recipe as 2,3,4,6. *)
param[5] = {x -> delta, y -> 1 - t};
ser5 = Series[letters /. param[5], {delta, 0, 1}, Assumptions -> {delta > 0, 0 < t < 1}] // Normal;
dlog5 = ((D[Log[#], t] & /@ ser5) // Simplify) /. delta -> 0 // Simplify;
dt[5] = Apart[Apart[dlog5 /. 1 - 2 t + 2 t^2 -> 2 (t - c2) (t - c2c), t] . axyTable, t];
Print["dt[5] built"];

(* segment 6: x->delta t^2, y->delta(1-t)^2 *)
param[6] = {x -> delta t^2, y -> delta (1 - t)^2};
ser6 = Series[letters /. param[6], {delta, 0, 1}, Assumptions -> {delta > 0, 0 < t < 1}] // Normal;
dlog6 = ((D[Log[#], t] & /@ ser6) // Simplify) /. delta -> 0 // Simplify;
dt[6] = Apart[Apart[dlog6 /. 1 - 2 t + 2 t^2 -> 2 (t - c2) (t - c2c), t] . axyTable, t];
Print["dt[6] built"];

Print["=== step 3: residue-decompose each dt[j] into per-letter M_a matrices ==="];
Do[
  den = Denominator[Together[dt[j]]];
  denoms = Union[Cases[den, _?(! FreeQ[#, t] &), Infinity]];
  Print["segment ", j, " nontrivial t-dependence check, denom factors present: ", Length[denoms] > 0];
  , {j, 1, 6}];

(* Build M_a^(j): residues of dt[j] at each pole in t, via Apart decomposition per-entry is
   expensive; instead use the SAME pattern as H_opt's prepend: since we only ever need
   dt[j] to multiply a G-word-free previous-weight solution (a 291-vector of G-word sums)
   and GIntegrate, we can keep dt[j] as the full rational-function-in-t matrix and apply
   GIntegrate directly at weight-build time (no M_a decomposition needed for correctness,
   just for speed; given PL's cheap per-weight cost so far, skip M_a decomposition for now
   and revisit only if weight-build becomes slow). *)

Print["=== step 4: todelta[j] for j=1..6 (weight-independent, for impart/analytic continuation) ==="];
Do[
  argsub = letters /. param[j];
  leading = Series[argsub, {delta, 0, 0}, Assumptions -> 0 < t < 1] // Normal;
  pos0 = Union[First /@ Position[leading, 0]];
  todelta[j] = Sum[Axy[letters[[i]]], {i, pos0}];
  Print["todelta[", j, "] built, ", Length[pos0], " vanishing letters"];
  , {j, 1, 6}];

Print["=== step 5: x1/y1/xmy connection matrices (no sqrt letters, so direct substitution suffices) ==="];
dyx1 = Apart[dy /. x -> 1, y];
dxy1 = Apart[dx /. y -> 1, x];
dxmy = Apart[(dy - dx) /. x -> -y, y];
Print["dyx1, dxy1, dxmy built"];

DumpSave["work/precompute_pl.mx", {Axy, letters, dx, dy, dt, param, todelta, dyx1, dxy1, dxmy}];
Print["checkpoint: work/precompute_pl.mx written"];
Print["=== PRECOMPUTE DONE ==="];
