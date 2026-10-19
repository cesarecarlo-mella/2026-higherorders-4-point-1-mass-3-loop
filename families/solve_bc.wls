#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Boundary conditions of one family to ep^wmax: the H_opt/run_fast.wls algorithm (segments
   1-6, x1, y1, xmy; matching; projector + regularity equations; planar constants on segment 5),
   with every ingredient derived from all.ints:
     a_tildes/<fam>_atilde.m (letters l[i], letters.dict.m), bases/<fam>_canon.m,
     result/<fam>_shift.m or SLC/Solutions/<fam>/shiftedBasis (planar entries),
     boundaries_all/<fam>_bc.m (weight 0, and validation of ep^1..6).
   Letters: the common 20-letter alphabet (a_tildes/letters.dict.m). Everything is built per
   letter, exactly as in H_opt/run_precompute.wls, then combined with the family's A_l:
   - segments 1-6: d/dt log of the delta-series of each letter along the path, delta -> 0
     (H's recipe); letters vanishing as delta -> 0 give todelta (analytic continuation);
   - x1, y1, xmy: d log along the line; the sqrt letters l[7], l[8] use H's derivative formulas
     with H's branches sqrt(x y (1-x-y)) -> I y (x=1), I x (y=1), -I y (x=-y);
   - boundary equations: projectors on every positive eigenvalue of the residue at x->0
     (segment 5), y->0 (segment 1), x+y->1 (segment 3) (sum of A_l over the letters vanishing
     there, as in H); regularity: residues of the connection at x=1, y=1, x=-y (H's breg) on the
     x1, y1, xmy segments;
   - planar constants: planar entries (fully reduced to BASISLC) x PL crossings
     (PL/pl_gpls_7_8), regularized at x->0, then y->1.
   BC[fam, X, i] = sol0 + Sum_w ep^w const_w(corner): "",x23 -> segment 6; x12,x132 -> 4; x13,x123 -> 2.

   Usage: wolframscript -file solve_bc.wls fam=C [wmax=8]
   Writes families/work/<fam>/ (checkpoints) and, if ep^1..6 agree with boundaries_all,
   families/results/<fam>_bc.m. *)

$HistoryLength = 0;
SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "config.wl"}]];
Get[FileNameJoin[{here, "..", "PL", "load_gpl_xy.wl"}]];
args = Association[Cases[$ScriptCommandLine[[2 ;;]], s_String /; StringContainsQ[s, "="] :> Rule @@ StringSplit[s, "=", 2]]];
f = args["fam"]; fsym = ToExpression[f];
wmax = ToExpression[Lookup[args, "wmax", "8"]];
work = FileNameJoin[{here, "work", f}]; Quiet[CreateDirectory[work, CreateIntermediateDirectories -> True]];
outDir = FileNameJoin[{here, "results"}]; Quiet[CreateDirectory[outDir]];
ep = Global`ep; eps = Global`eps;
debug = MemberQ[$ScriptCommandLine, "debug"];
mark[s_] := If[debug, Print["  [", s, "]"]];
tic[label_, expr_] := Module[{t, r}, {t, r} = AbsoluteTiming[expr]; Print["  ", label, ": ", Round[t, 0.01], " s"]; r];
SetAttributes[tic, HoldRest];

(* ---------- inputs ---------- *)
atilde = Import[atildeFile[f]]; n = Length[atilde];
lets = Union@Cases[atilde, l[i_] :> i, Infinity];
Alet = Association[Table[i -> SparseArray[Coefficient[atilde, l[i]]], {i, lets}]];
If[Expand[atilde - Sum[Normal[Alet[i]] l[i], {i, lets}]] =!= ConstantArray[0, {n, n}], Print["atilde not linear in l[i]"]; Exit[1]];
arg = Association[(#[[1, 1]] -> (#[[2]] /. Log -> Identity)) & /@ Import[lettersDictFile]];
(* derivatives; sqrt = Sqrt[x y (1-x-y)] kept as a symbol for the two sqrt letters (H's formulas) *)
dX[i_] := Switch[i, 7, I y (-1 + 2 x + y)/(sqrt (-1 + x) (x + y)), 8, I y/(sqrt (1 - x)), _, D[Log[arg[i]], x]];
dY[i_] := Switch[i, 7, -I x/(sqrt (x + y)), 8, I x/(sqrt (1 - y)), _, D[Log[arg[i]], y]];
bcOld = Association[Import[bcFile[f]]];
sol0 = Table[Coefficient[bcOld[BC[fsym, "", i]], ep, 0], {i, n}];
Print[f, ": ", n, " masters, letters ", lets];

canon[v_] := Expand[ComplexExpand[v]];
constPat = Pi | _Zeta | _MZV;
normMZV = s_Symbol /; SymbolName[s] === "MZV" :> Global`MZV;

(* ---------- sparse {keys, S} algebra (as in H_opt/run_fast.wls) ---------- *)
emptyF = {{}, SparseArray[{}, {n, 0}]};
prune[{keys_, S_}] := Module[{p, S2, cols},
  If[keys === {}, Return[{{}, SparseArray[{}, {Length[S], 0}]}]];
  p = S["NonzeroPositions"]; If[p === {}, Return[{{}, SparseArray[{}, {Length[S], 0}]}]];
  S2 = SparseArray[Thread[p -> Expand[S["NonzeroValues"]]], Dimensions[S]];
  cols = Union[S2["NonzeroPositions"][[All, 2]]];
  If[cols === {}, {{}, SparseArray[{}, {Length[S], 0}]}, {keys[[cols]], S2[[All, cols]]}]];
mergeBlocks[blocks0_List] := Module[{blocks = Select[blocks0, #[[1]] =!= {} &], keys, mats, uk, idx, P},
  If[blocks === {}, Return[If[blocks0 === {}, emptyF, {{}, SparseArray[{}, {Length[blocks0[[1, 2]]], 0}]}]]];
  keys = Join @@ blocks[[All, 1]]; mats = Join[Sequence @@ blocks[[All, 2]], 2];
  uk = Union[keys]; idx = AssociationThread[uk -> Range[Length[uk]]];
  P = SparseArray[Table[{k, idx[keys[[k]]]} -> 1, {k, Length[keys]}], {Length[keys], Length[uk]}];
  prune[{uk, mats . P}]];
add[fs__] := mergeBlocks[{fs}];
neg[{k_, S_}] := {k, -S};
sub[a_, b_] := add[a, neg[b]];
matMul[M_, {k_, S_}] := If[k === {}, {{}, SparseArray[{}, {Length[M], 0}]}, prune[{k, M . S}]];
scal[{k_, S_}, c_, m_] := If[k === {}, {k, S}, {{#[[1]], #[[2]] m} & /@ k, c S}];
ev0[{k_, S_}] := Module[{sel = Flatten@Position[k, {{}, _}, {1}]}, If[sel === {}, {{}, SparseArray[{}, {Length[S], 0}]}, {k[[sel]], S[[All, sel]]}]];
nonConst[{k_, S_}] := Count[k, {w_, _} /; w =!= {}];
align[{k_, S_}, allk_] := Module[{idx = AssociationThread[allk -> Range[Length[allk]]]},
  If[k === {}, SparseArray[{}, {Length[S], Length[allk]}],
   S . SparseArray[Table[{c, idx[k[[c]]]} -> 1, {c, Length[k]}], {Length[k], Length[allk]}]]];
constF[v_List] := Module[{rows, keys, idx, terms},   (* list of constants -> {keys {{}, mono}, S} *)
  rows = Flatten[Table[With[{e = Expand[v[[i]]]}, If[e === 0, {}, terms = If[Head[e] === Plus, List @@ e, {e}];
       {i, {{}, # / (# /. constPat -> 1)}, # /. constPat -> 1} & /@ terms]], {i, Length[v]}], 1];
  If[rows === {}, Return[{{}, SparseArray[{}, {Length[v], 0}]}]];
  keys = Union[rows[[All, 2]]]; idx = AssociationThread[keys -> Range[Length[keys]]];
  prune[{keys, SparseArray[{#1, idx[#2]} -> #3 & @@@ rows, {Length[v], Length[keys]}]}]];
toVals[{k_, S_}] := If[k === {}, ConstantArray[0, Length[S]], Expand[S . (#[[2]] & /@ k)]];

(* t -> 1 for {0,1}-words *)
wordTable = Association[Get[FileNameJoin[{here, "..", "PL", "work", "word_table_0to8.m"}]] /. (G[a___, t] -> v_) :> ({a} -> v) /. (1 -> 1) -> ({} -> 1)];
splitMono[e_] := Module[{ex = Expand[e /. normMZV], terms}, If[ex === 0, Return[{}]];
  terms = If[Head[ex] === Plus, List @@ ex, {ex}]; ({# /. constPat -> 1, # / (# /. constPat -> 1)} & /@ terms)];
wordVal = Association[KeyValueMap[#1 -> splitMono[#2] &, wordTable]];
ev1[{k_, S_}] := Module[{bad, triples, nk, idx, T},
  If[k === {}, Return[{{}, SparseArray[{}, {Length[S], 0}]}]];
  bad = Select[k, ! SubsetQ[{0, 1}, #[[1]]] &];
  If[bad =!= {}, Print["ev1: words with letters outside {0,1}: ", Length[bad], " e.g. ", bad[[;; UpTo[3]]]]; Exit[1]];
  triples = Flatten[Table[{c, {{}, k[[c, 2]] #[[2]]}, #[[1]]} & /@ wordVal[k[[c, 1]]], {c, Length[k]}], 1];
  If[triples === {}, Return[{{}, SparseArray[{}, {Length[S], 0}]}]];
  nk = Union[triples[[All, 2]]]; idx = AssociationThread[nk -> Range[Length[nk]]];
  T = SparseArray[{#1, idx[#2]} -> #3 & @@@ triples, {Length[k], Length[nk]}];
  prune[{nk, S . T}]];

(* ---------- spectral projectors (PL/routines_clean/spectral_projection.wls) ---------- *)
getSpectrum[M_] := Module[{numEig, rounded, tally},
  numEig = Eigenvalues[N[Normal[M], 30]]; rounded = Round /@ numEig;
  If[Max[Abs[numEig - rounded]] > 10^-10, Print["non-integer spectrum"]; Exit[1]];
  tally = Tally[rounded];
  Do[If[Tr[MatrixPower[Normal[M], k]] =!= Sum[tt[[2]] tt[[1]]^k, {tt, tally}], Print["trace check failed"]; Exit[1]], {k, 1, 4}];
  tally];
projLambda[M0_, lam_, spec_] := Module[{M = Normal[M0], others, gScalar, gOfA, w, series, coeffs, B, qOfA, P, aLam, m},
  m = Length[M]; aLam = Cases[spec, {lam, mm_} :> mm][[1]]; others = DeleteCases[spec, {lam, _}];
  gOfA = If[others === {}, IdentityMatrix[m], Dot @@ (MatrixPower[M - #[[1]] IdentityMatrix[m], #[[2]]] & /@ others)];
  gScalar[zz_] := If[others === {}, 1, Times @@ ((zz - #[[1]])^#[[2]] & /@ others)];
  coeffs = CoefficientList[Normal[Series[1/gScalar[lam + w], {w, 0, aLam - 1}]], w];
  B = M - lam IdentityMatrix[m];
  qOfA = If[aLam === 1, coeffs[[1]] IdentityMatrix[m],
    Fold[(#1 . B + coeffs[[#2]] IdentityMatrix[m]) &, coeffs[[aLam]] IdentityMatrix[m], Reverse[Range[1, aLam - 1]]]];
  P = qOfA . gOfA;
  If[Tr[P] =!= aLam || Max[Abs[P . P - P]] =!= 0, Print["projector check failed for lambda=", lam]; Exit[1]];
  SparseArray[P]];
allPositive[M_, name_] := Module[{spec = getSpectrum[M], pos},
  pos = Select[spec[[All, 1]], # > 0 &];
  Print["  ", name, ": spectrum ", Sort[spec], ", projecting on ", pos];
  If[pos === {}, SparseArray[{}, {0, n}], Join @@ (projLambda[M, #, spec] & /@ pos)]];

(* ---------- weight-independent data ---------- *)
segs = {1, 2, 3, 4, 5, 6, "x1", "y1", "xmy"};
param = <|1 -> {x -> t, y -> delta}, 2 -> {y -> delta t^2, x -> 1 - delta ((1 - t)^2 + t^2)}, 3 -> {x -> 1 - delta - t, y -> t},
   4 -> {x -> delta (1 - t)^2, y -> 1 - delta ((1 - t)^2 + t^2)}, 5 -> {x -> delta, y -> 1 - t}, 6 -> {x -> delta t^2, y -> delta (1 - t)^2},
   "x1" -> {x -> 1}, "y1" -> {y -> 1}, "xmy" -> {x -> -y}|>;
segVar = <|1 -> t, 2 -> t, 3 -> t, 4 -> t, 5 -> t, 6 -> t, "x1" -> y, "y1" -> x, "xmy" -> y|>;
(* simple poles of a rational function of var: {root -> residue}; anything else is an error *)
poles[e_, var_] := Module[{tg = Together[e], num, den, sl, rts},
  If[tg === 0, Return[{}]];
  If[! FreeQ[tg, Complement[{x, y, t, delta, sqrt}, {var}]], Print["poles: unexpected variables in ", tg]; Exit[1]];
  num = Numerator[tg]; den = Denominator[tg];
  If[Exponent[num, var] >= Exponent[den, var], Print["poles: polynomial part in ", tg]; Exit[1]];
  sl = Solve[den == 0, var]; rts = Tally[If[sl === {}, {}, var /. sl]];
  If[AnyTrue[rts[[All, 2]], # > 1 &], Print["poles: double pole in ", tg]; Exit[1]];
  (canon[#] -> canon[Together[num/D[den, var] /. var -> #]]) & /@ rts[[All, 1]]];
addPoles[acc_, prs_, M_] := Fold[Function[{a, r}, Append[a, r[[1]] -> Lookup[a, r[[1]], 0] + r[[2]] M]], acc, prs];
(* segments 1-6: H's recipe per letter *)
segLetter[j_, i_] := Module[{ser, v0, dl},
  ser = Normal[Series[arg[i] /. param[j], {delta, 0, 2}, Assumptions -> {delta > 0, 0 < t < 1}]];
  v0 = Simplify[ser /. delta -> 0, 0 < t < 1];
  dl = Simplify[Simplify[D[Log[ser], t], {delta > 0, 0 < t < 1}] /. delta -> 0, 0 < t < 1];
  {v0 === 0, dl}];
(* lines: d log along the line, H's branches for sqrt *)
lineDlog[s_, i_] := Switch[s,
  "x1", Together[dY[i]] /. x -> 1 /. sqrt -> I y,
  "y1", Together[dX[i]] /. y -> 1 /. sqrt -> I x,
  "xmy", Together[dY[i] - dX[i]] /. x -> -y /. sqrt -> -I y];
(* regularity residues (H's breg): (x-1) dx at x=1, (y-1) dy at y=1, (x+y) dx at x=-y *)
regRes[r_, i_] := Module[{v = Switch[r,
     "R1", Together[(x - 1) dX[i]] /. x -> 1 /. sqrt -> I y,
     "R2", Together[(y - 1) dY[i]] /. y -> 1 /. sqrt -> I x,
     "R3", Together[(x + y) dX[i]] /. x -> -y /. sqrt -> -I y]},
  v = Simplify[v];
  If[! FreeQ[v, x | y], Print["regularity residue of l[", i, "] at ", r, " is not constant: ", v]; Exit[1]]; canon[v]];
(* letters vanishing at the boundary loci (H: Series of the letters, zero coefficient at order 0) *)
vanishAt["x0", i_] := Simplify[Coefficient[Normal[Series[arg[i], {x, 0, 1}, Assumptions -> {x > 0, y > 0, y < 1}]], x, 0]] === 0;
vanishAt["y0", i_] := Simplify[Coefficient[Normal[Series[arg[i], {y, 0, 1}, Assumptions -> {y > 0, x > 0, x < 1}]], y, 0]] === 0;
vanishAt["xy1", i_] := Simplify[Coefficient[Normal[Series[arg[i] /. {x -> (1 - delta) rho, y -> (1 - delta) (1 - rho)}, {delta, 0, 1},
       Assumptions -> {delta > 0, 0 < rho < 1}]], delta, 0]] === 0;
dataFile = FileNameJoin[{work, "precompute.mx"}];
If[FileExistsQ[dataFile], Get[dataFile],
  Print["=== precompute ==="];
  Do[L[j] = <||>; todelta[j] = SparseArray[{}, {n, n}];
    Do[{v0, dl} = segLetter[j, i];
      If[v0, todelta[j] += Alet[i]];
      L[j] = addPoles[L[j], poles[dl, t], Alet[i]], {i, lets}], {j, 1, 6}];
  Do[L[s] = <||>; todelta[s] = SparseArray[{}, {n, n}];
    Do[L[s] = addPoles[L[s], poles[lineDlog[s, i], segVar[s]], Alet[i]], {i, lets}], {s, {"x1", "y1", "xmy"}}];
  Do[L[s] = Select[L[s], Max[Abs[#]] =!= 0 &], {s, segs}];
  res = Association[Table[r -> Sum[regRes[r, i] Alet[i], {i, lets}], {r, {"R1", "R2", "R3"}}]];
  loc = Association[Table[lc -> Select[lets, vanishAt[lc, #] &], {lc, {"x0", "y0", "xy1"}}]];
  Print["  letters vanishing at x->0: ", loc["x0"], ", y->0: ", loc["y0"], ", x+y->1: ", loc["xy1"]];
  bm = <|"P1" -> allPositive[Sum[Alet[i], {i, loc["x0"]}], "x->0"], "P2" -> allPositive[Sum[Alet[i], {i, loc["y0"]}], "y->0"],
     "P7" -> allPositive[Sum[Alet[i], {i, loc["xy1"]}], "x+y->1"], "R1" -> res["R1"], "R2" -> res["R2"], "R3" -> res["R3"]|>;
  DumpSave[dataFile, {L, todelta, bm}]];
(* powers of todelta for the analytic continuation, computed on demand (any weight) *)
tdP[j_, 1] := todelta[j];
tdP[j_, k_Integer] := tdP[j, k] = tdP[j, k - 1] . todelta[j];
Do[Print["  letters on ", s, ": ", Keys[L[s]]], {s, segs}];
integrate[s_, F_] := If[F[[1]] === {}, F,
  mergeBlocks[Table[{{Prepend[#[[1]], a], #[[2]]} & /@ F[[1]], L[s][a] . F[[2]]}, {a, Keys[L[s]]}]]];

(* ---------- planar entries: shift map (or SLC shiftedBasis) -> BASISLC coefficient matrix ---------- *)
crossLabels = {"", "x12", "x13", "x23", "x123", "x132"};
cross = <|"" -> {y -> y, x -> x}, "x12" -> {y -> 1 - x - y, x -> x}, "x13" -> {y -> y, x -> 1 - x - y},
   "x23" -> {y -> x, x -> y}, "x123" -> {y -> x, x -> 1 - x - y}, "x132" -> {y -> 1 - x - y, x -> y}|>;
bl[c_] := ToExpression["BASISLC" <> c];
blVars = Flatten[Table[bl[c][k], {c, crossLabels}, {k, 291}]];
planarFile = FileNameJoin[{work, "planar.mx"}];
If[FileExistsQ[planarFile], Get[planarFile],
  Print["=== planar entries ==="];
  plTab = Flatten[Get /@ plTableFiles] /. {s23 -> 1 - y - x, s13 -> y};
  plRules = Dispatch[Flatten[Table[plTab /. PL[a__] :> INT[ToExpression["PL" <> c], a] /. BASISLC -> bl[c] /. cross[c], {c, crossLabels}]]];
  canonB = Import[canonFile[f]];
  shifted = If[! MissingQ[shiftMapFile[f]],
    With[{sh = Import[shiftMapFile[f]] /. INT[lb_String, __, {b__}] :> INT[lb, b]},
     With[{own = First[Commonest[Cases[sh[[All, 1]], INT[lb_String, ___] /; ! StringContainsQ[lb, "x"] && ! MemberQ[plHeads, lb] :> lb]]]},
      (canonB /. INT[fsym, b__] :> INT[own, b]) //. Dispatch[sh /. INT[lb_String, b__] /; MemberQ[plHeads, lb] :> INT[ToExpression[lb], b]]]],
    Import[slcShiftedBasisFile[f]]];
  shifted = shifted /. plRules /. {m2 -> 1, s12 -> x, s13 -> y, s23 -> 1 - x - y};
  op = Flatten@Position[shifted, e_ /; FreeQ[e, _INT], {1}, Heads -> False];
  ca = CoefficientArrays[shifted[[op]], blVars];
  plCoef = ca[[2]];
  Print["  ", Length[op], " planar entries (source: ", If[MissingQ[shiftMapFile[f]], "SLC shiftedBasis", "shift map"], ")"];
  DumpSave[planarFile, {op, plCoef}]];
unk = Complement[Range[n], op];
Print["  n = ", n, ", planar ", Length[op], " (distinct ", Length[Union[op]], ", range ", MinMax[op], "), unknowns ", Length[unk]];
w0PL = Get[FileNameJoin[{here, "..", "PL", "work", "weight0_solution.m"}]];
plF[m_, c_] := plF[m, c] = If[m == 0,
   {{{{}, {}, 1}}, SparseArray[List /@ w0PL]},
   With[{F = PLGPLSparse[m, If[c === "", "x", c]]}, {F[[1]] /. normMZV, F[[2]]}]];
(* planar constants on segment 5 at weight w: Sum_j R_j . PL_{w-j}, constant in x,y; then x->0, y->1 *)
planarConst[w_] := Module[{vals, dif, k, S, sel},
  vals = Table[Module[{R = Normal[plCoef] /. pt /. d -> 4 - 2 eps, ser, lo},
      ser = Normal[Series[R, {eps, 0, w}]];
      lo = Min[0, Min[Exponent[#, eps, Min] & /@ Flatten[ser]]];
      If[lo < 0, Print["planar coefficients have eps poles"]; Exit[1]];
      mergeBlocks[Flatten[Table[With[{Rj = SparseArray[Coefficient[ser, eps, j]]},
          Table[With[{F = plF[w - j, crossLabels[[ci]]]}, {F[[1]], Rj[[All, Range[291] + 291 (ci - 1)]] . F[[2]]}], {ci, 6}]], {j, 0, w}], 1]]],
    {pt, {{x -> 3/13, y -> 5/18}, {x -> 7/23, y -> 2/11}}}];
  dif = mergeBlocks[{vals[[1]], neg[vals[[2]]]}];
  If[dif[[1]] =!= {}, Print["planar combination depends on x,y at weight ", w]; Exit[1]];
  {k, S} = vals[[1]];
  sel = Flatten@Position[k, {{}, _, _}, {1}];
  If[sel === {}, {{}, SparseArray[{}, {Length[op], 0}]}, ev1[{{#[[2]], #[[3]]} & /@ k[[sel]], S[[All, sel]]}]]];

(* ---------- boundary system (weight independent) ---------- *)
eqs = {{"P1", 5}, {"P2", 1}, {"P7", 3}, {"R1", "x1"}, {"R2", "y1"}, {"R3", "xmy"}};
eqs = Select[eqs, Length[bm[#[[1]]]] > 0 &];   (* a locus without positive eigenvalues gives no equation *)
Amat = Join @@ (bm[#[[1]]] & /@ eqs);
Au = Amat[[All, unk]]; Ap = Amat[[All, op]];
pivFile = FileNameJoin[{work, "pivots.mx"}];
If[unk === {}, pivRows = {}; invA = {}; Print["  all entries planar: no unknowns, every equation is a check"]];
If[unk =!= {}, If[FileExistsQ[pivFile], Get[pivFile],
  rr = tic["RowReduce", RowReduce[Transpose[Normal[Au]]]];
  pivRows = Flatten[FirstPosition[#, Except[0], Missing[], {1}, Heads -> False] & /@ Select[rr, AnyTrue[#, # =!= 0 &] &]];
  Print["  boundary system: ", Length[Amat], " equations, rank ", Length[pivRows], " for ", Length[unk], " unknowns (", Length[op], " planar)"];
  If[Length[pivRows] < Length[unk], Print[f, ": boundary system is rank-deficient by ", Length[unk] - Length[pivRows], " -- more input needed"]; Exit[1]];
  invA = tic["inverse", Inverse[Normal[Au[[pivRows]]]]];
  DumpSave[pivFile, {pivRows, invA}]]];

(* ---------- weights ---------- *)
Do[sol[0][s] = constF[sol0], {s, segs}];
corner = <|"" -> 6, "x12" -> 4, "x13" -> 2, "x23" -> 6, "x123" -> 2, "x132" -> 4|>;
ok = True;
Do[
  solFile = FileNameJoin[{work, "sol" <> ToString[w] <> ".mx"}];
  If[FileExistsQ[solFile], With[{sw = Import[solFile, "MX"]}, Do[sol[w][s] = sw[s], {s, segs}]]; Continue[]];
  tw = AbsoluteTime[];
  mark["start"]; Do[IP[s] = integrate[s, sol[w - 1][s]], {s, segs}]; mark["integrate"];
  impart[j_] := add @@ Prepend[Table[matMul[tdP[j, k], scal[sol[w - k][j], (-I)^k/k!, Pi^k]], {k, 1, w}], IP[j]];
  dd = Table[sub[ev1[IP[k - 1]], ev0[IP[k]]], {k, 2, 5}]; dd = Prepend[dd, None];
  d6 = sub[ev1[IP[5]], ev0[IP[6]]];
  imp5 = impart[5]; imp3 = impart[3];
  mark["d2..d6, impart"]; d00 = sub[ev1[imp5], ev0[IP["xmy"]]]; d01 = sub[ev1[imp3], ev0[IP["y1"]]]; d10 = sub[ev0[imp3], ev0[IP["x1"]]];
  delta[5] = emptyF; delta[4] = neg[dd[[5]]]; delta[3] = sub[delta[4], dd[[4]]];
  delta[2] = sub[delta[3], dd[[3]]]; delta[1] = sub[delta[2], dd[[2]]];
  delta[6] = d6; delta["xmy"] = d00; delta["y1"] = add[delta[3], d01]; delta["x1"] = add[delta[3], d10];
  mark["matching"]; rhsB = Table[With[{full = matMul[bm[e[[1]]], IP[e[[2]]]]},
     If[nonConst[full] > 0, Print["  WARNING: ", e[[1]], ".tempSol[", e[[2]], "] depends on the variable (", nonConst[full], " keys)"]];
     neg[add[ev0[full], matMul[bm[e[[1]]], delta[e[[2]]]]]]], {e, eqs}];
  mark["rhs"]; pl = planarConst[w]; mark["planar"];
  monos = Union[Join @@ (#[[1]] & /@ Append[rhsB, pl])];
  If[monos === {}, monos = {{{}, 1}}];
  rhs = Join @@ (align[#, monos] & /@ rhsB);
  plV = If[pl[[1]] === {}, SparseArray[{}, {Length[op], Length[monos]}], align[pl, monos]];
  rhs2 = rhs - Ap . plV; mark["rhs2"];
  xs = If[unk === {}, {}, Expand[invA . Normal[rhs2[[pivRows]]]]];
  nres = Count[Flatten[Normal[Expand[If[unk === {}, -rhs2, Au . xs - rhs2]]]], Except[0]];
  c5M = Normal[SparseArray[{}, {n, Length[monos]}]]; c5M[[unk]] = xs; c5M[[op]] = Normal[plV];
  c5 = prune[{monos, SparseArray[c5M]}];
  Do[sol[w][s] = add[IP[s], c5, delta[s]], {s, segs}];
  Export[solFile, Association[Table[s -> sol[w][s], {s, segs}]], "MX"];
  (* validation against boundaries_all *)
  vmsg = If[w <= 6,
    bad = Sum[Count[Expand[toVals[ev0[sol[w][corner[X]]]] - Table[Coefficient[bcOld[BC[fsym, X, i]], ep, w], {i, n}]], Except[0]], {X, crossLabels}];
    If[bad > 0, ok = False]; "; vs boundaries_all: " <> If[bad == 0, "identical", ToString[bad] <> " entries differ"], ""];
  Print["w", w, ": ", Round[AbsoluteTime[] - tw, 0.1], " s, ", Length[monos], " transcendental columns, residuals ", nres, vmsg];
  If[nres > 0, Print["inconsistent boundary system"]; Exit[1]],
  {w, 1, wmax}];

(* validation of weights loaded from checkpoints is repeated here *)
Do[bad = Sum[Count[Expand[toVals[ev0[sol[w][corner[X]]]] - Table[Coefficient[bcOld[BC[fsym, X, i]], ep, w], {i, n}]], Except[0]], {X, crossLabels}];
  If[bad > 0, ok = False; Print["w", w, ": ", bad, " entries differ from boundaries_all"]], {w, 1, Min[6, wmax]}];
If[ok,
  bc = Flatten[Table[Table[BC[fsym, X, i] -> sol0[[i]] + Sum[ep^w toVals[ev0[sol[w][corner[X]]]][[i]], {w, 1, wmax}], {i, n}], {X, crossLabels}]];
  Export[FileNameJoin[{outDir, f <> "_bc.m"}], bc];
  Print[f, ": exported results/", f, "_bc.m (", Length[bc], " entries, ep^0..", wmax, ")"],
  Print[f, ": NOT exported -- disagreement with boundaries_all"]];
