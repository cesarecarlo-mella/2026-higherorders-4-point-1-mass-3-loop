#!/usr/bin/env wolframscript
(* ::Package:: *)

(* H_opt run_fast.wls -- family H, all weights 1..wmax with one script.

   Same algorithm as H/solution.nb "Boundary via Expansion" (segments 1-6, x1, y1, xmy;
   matching; projector + regularity equations; planar constants on segment 5), but:
   - every segment solution is a sparse matrix S (371 x #keys), key = {word, monomial},
     sol = Sum S[[i,k]] G[word; var] monomial, coefficients exact in Q(i, Sqrt[3]);
   - integration on a segment = prepending a letter: dt = Sum_a M_a/(var-a);
   - t->1 is only ever applied to {0,1}-words, via PL/work/word_table_0to8.m
     (direct G[w;1] with HPL[{1},1]->0 -- identical to the tto1Red reflection, no fibration);
   - x1/y1/xmy (cube-root letters) are only evaluated at their start point (all words -> 0);
   - P.tempSol and regCond.tempSol are checked to be variable-independent, then the stacked
     boundary system is solved exactly (it is not read from H/numerical_eq/solution_w.m);
   - planar constants: shiftedBasis (onlyplanar rows) x planar crossings from
     PL/pl_gpls_7_8 (run PL/run_gpl_xy_fast.wls 8 cross=<c> for all six crossings first),
     at regularized x->0, then y->1.

   Usage: wolframscript -file run_fast.wls <wmax> [validate]
   Writes H_opt/work/fast/sol<w>.mx (Association seg -> {keys, S}); restarts from the
   highest weight already on disk. "validate" compares w<=6 with H/solutions/eps_sep. *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
Global`x; Global`y; Global`t; Global`delta; Global`c2; Global`c2c; Global`cxmy; Global`cxmyc; Global`sqrt; Global`eps;
Global`Basis; Global`r1; Global`r2; Global`r3; Global`r4; Global`H;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];
SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];

root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
Hdir = FileNameJoin[{root, "H"}];
workdir = FileNameJoin[{root, "H_opt", "work"}];
fastdir = FileNameJoin[{workdir, "fast"}];
Quiet[CreateDirectory[fastdir]];
wmax = ToExpression[$ScriptCommandLine[[2]]];
validate = MemberQ[$ScriptCommandLine[[3 ;;]], "validate"];
tic[label_, expr_] := Module[{t, r}, {t, r} = AbsoluteTiming[expr]; Print["  ", label, ": ", Round[t, 0.01], " s"]; r];
SetAttributes[tic, HoldRest];

Get[FileNameJoin[{workdir, "precompute.mx"}]];
n = Length[Basis];
op = Sort[Flatten[onlyplanar]]; unk = Complement[Range[n], op];

(* ---------- exact numbers in Q(i, Sqrt[3]) ---------- *)
numRules = {c2 -> (1 - I)/2, c2c -> (1 + I)/2, cxmy -> -(-1)^(1/3), cxmyc -> (-1)^(2/3)};
canon[v_] := Expand[ComplexExpand[v /. numRules]];
constPat = Pi | _Zeta | _MZV;

(* ---------- sparse {keys, S} algebra ---------- *)
emptyF = {{}, SparseArray[{}, {n, 0}]};
prune[{keys_, S_}] := Module[{p, S2, cols},
  If[keys === {}, Return[emptyF]];
  p = S["NonzeroPositions"];
  If[p === {}, Return[emptyF]];
  S2 = SparseArray[Thread[p -> Expand[S["NonzeroValues"]]], Dimensions[S]];
  cols = Union[S2["NonzeroPositions"][[All, 2]]];
  If[cols === {}, emptyF, {keys[[cols]], S2[[All, cols]]}]];
mergeBlocks[blocks0_List] := Module[{blocks = Select[blocks0, #[[1]] =!= {} &], keys, mats, uk, idx, P},
  If[blocks === {}, Return[emptyF]];
  keys = Join @@ blocks[[All, 1]]; mats = Join[Sequence @@ blocks[[All, 2]], 2];
  uk = Union[keys]; idx = AssociationThread[uk -> Range[Length[uk]]];
  P = SparseArray[Table[{k, idx[keys[[k]]]} -> 1, {k, Length[keys]}], {Length[keys], Length[uk]}];
  prune[{uk, mats . P}]];
add[fs__] := mergeBlocks[{fs}];
neg[{k_, S_}] := {k, -S};
sub[a_, b_] := add[a, neg[b]];
matMul[M_, {k_, S_}] := If[k === {}, emptyF, prune[{k, M . S}]];
scal[{k_, S_}, c_, m_] := If[k === {}, emptyF, {{#[[1]], #[[2]] m} & /@ k, c S}];
(* regularized value at the start point: every non-empty word -> 0 *)
ev0[{k_, S_}] := Module[{sel = Flatten@Position[k, {{}, _}, {1}]}, If[sel === {}, emptyF, {k[[sel]], S[[All, sel]]}]];
nonConst[{k_, S_}] := Count[k, {w_, _} /; w =!= {}];
(* align a constant-only F to a given list of monomial keys *)
align[{k_, S_}, allk_] := Module[{idx = AssociationThread[allk -> Range[Length[allk]]]},
  If[k === {}, SparseArray[{}, {Length[S], Length[allk]}],
   S . SparseArray[Table[{c, idx[k[[c]]]} -> 1, {c, Length[k]}], {Length[k], Length[allk]}]]];
toExprs[{k_, S_}, var_] := Module[{vec = (If[#[[1]] === {}, 1, G[Sequence @@ #[[1]], var]] #[[2]]) & /@ k},
  If[k === {}, ConstantArray[0, Length[S]],
   Table[With[{r = S[[i]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]], {i, Length[S]}]]];

(* ---------- t -> 1 for {0,1}-words: word table ---------- *)
wordTable = Association[Import[FileNameJoin[{root, "PL", "work", "word_table_0to8.m"}]] /. (G[a___, t] -> v_) :> ({a} -> v) /. (1 -> 1) -> ({} -> 1)];
splitMono[e_] := Module[{ex = Expand[e], terms},
  If[ex === 0, Return[{}]];
  terms = If[Head[ex] === Plus, List @@ ex, {ex}];
  ({# /. constPat -> 1, # / (# /. constPat -> 1)} & /@ terms)];
wordVal = Association[KeyValueMap[#1 -> splitMono[#2] &, wordTable]];
ev1[{k_, S_}] := Module[{bad, triples, nk, idx, T},
  If[k === {}, Return[emptyF]];
  bad = Select[k, ! SubsetQ[{0, 1}, #[[1]]] &];
  If[bad =!= {}, Print["ev1: words with letters outside {0,1}: ", Length[bad], " e.g. ", bad[[;; UpTo[3]]]]; Abort[]];
  triples = Flatten[Table[{c, {{}, k[[c, 2]] #[[2]]}, #[[1]]} & /@ wordVal[k[[c, 1]]], {c, Length[k]}], 1];
  If[triples === {}, Return[emptyF]];
  nk = Union[triples[[All, 2]]]; idx = AssociationThread[nk -> Range[Length[nk]]];
  T = SparseArray[{#1, idx[#2]} -> #3 & @@@ triples, {Length[k], Length[nk]}];
  prune[{nk, S . T}]];

(* ---------- connection on each segment as Sum_a M_a/(var - a) ---------- *)
letterForm[Massoc_, var_] := Module[{acc = <||>, push, ar, dd, al, a, f, terms, den, num},
  push[aa_, rules_] := (acc[aa] = Join[Lookup[acc, aa, {}], rules]);
  Do[
    ar = Most[ArrayRules[SparseArray[Massoc[d]]]];
    dd = d /. numRules;
    If[PolynomialQ[dd, var] && Exponent[dd, var] == 1 && FreeQ[ar[[All, 2]], var],
      al = Coefficient[dd, var, 1]; a = canon[-(dd /. var -> 0)/al];
      push[a, (#[[1]] -> #[[2]]/al) & /@ ar],
      Do[
        f = Apart[Together[(r[[2]] /. numRules)/dd], var];
        terms = If[Head[f] === Plus, List @@ f, {f}];
        Do[den = Denominator[Together[tm]]; num = Numerator[Together[tm]];
          If[! (PolynomialQ[den, var] && Exponent[den, var] == 1 && FreeQ[num, var]),
            Print["letterForm: unexpected term ", tm]; Abort[]];
          al = Coefficient[den, var, 1]; a = canon[-(den /. var -> 0)/al];
          push[a, {r[[1]] -> num/al}], {tm, terms}], {r, ar}]],
    {d, Keys[Massoc]}];
  Association[KeyValueMap[#1 -> SparseArray[(#[[1]] -> canon[#[[2]]]) & /@ #2, {n, n}] &, acc]]];

lfFile = FileNameJoin[{fastdir, "letters.mx"}];
If[FileExistsQ[lfFile], Get[lfFile],
  Print["=== letter decomposition of segment connections ==="];
  L = <||>;
  Do[L[j] = tic["segment " <> ToString[j], letterForm[Mdt[j], t]], {j, 1, 6}];
  L["x1"] = tic["x1", letterForm[Mdyx1, y]];
  L["y1"] = tic["y1", letterForm[Mdxy1, x]];
  L["xmy"] = tic["xmy", letterForm[Mdxmy, y]];
  Do[tdP[j, 1] = SparseArray[canon /@ todelta[j]], {j, {3, 5}}];
  Do[tdP[j, k] = SparseArray[Expand[tdP[j, k - 1] . tdP[j, 1]]], {j, {3, 5}}, {k, 2, 8}];
  bm = SparseArray[canon /@ #] & /@ <|"P1" -> boundaryMxy1, "P2" -> boundaryMxy2, "P7" -> boundaryMxy7,
     "R1" -> regCondxy1, "R2" -> regCondxy2, "R3" -> regCondxy3|>;
  DumpSave[lfFile, {L, tdP, bm}]];
segs = {1, 2, 3, 4, 5, 6, "x1", "y1", "xmy"};
segVar = <|1 -> t, 2 -> t, 3 -> t, 4 -> t, 5 -> t, 6 -> t, "x1" -> y, "y1" -> x, "xmy" -> y|>;
Do[Print["letters ", s, ": ", Keys[L[s]]], {s, segs}];

integrate[s_, F_] := If[F[[1]] === {}, emptyF,
  mergeBlocks[Table[{{Prepend[#[[1]], a], #[[2]]} & /@ F[[1]], L[s][a] . F[[2]]}, {a, Keys[L[s]]}]]];

(* ---------- exact boundary system: fixed matrix, solved once ---------- *)
Amat = Join[bm["P1"], bm["P2"], bm["P7"], bm["R1"], bm["R2"], bm["R3"]];
Au = Amat[[All, unk]]; Ap = Amat[[All, op]];
pivFile = FileNameJoin[{fastdir, "pivots.mx"}];
If[FileExistsQ[pivFile], Get[pivFile],
  Print["=== selecting independent boundary rows ==="];
  rr = tic["RowReduce", RowReduce[Transpose[Normal[Au]]]];
  pivRows = Flatten[FirstPosition[#, Except[0], Missing[], {1}, Heads -> False] & /@ Select[rr, AnyTrue[#, # =!= 0 &] &]];
  Print["  rank ", Length[pivRows], " for ", Length[unk], " unknowns"];
  If[Length[pivRows] =!= Length[unk], Print["boundary system is not of full rank"]; Abort[]];
  invA = tic["inverse", Inverse[Normal[Au[[pivRows]]]]];
  DumpSave[pivFile, {pivRows, invA}]];

(* ---------- planar constants on segment 5 ---------- *)
Get[FileNameJoin[{root, "PL", "load_gpl_xy.wl"}]];
crossHead = <|"x" -> BASISLC, "x12" -> BASISLCx12, "x13" -> BASISLCx13, "x23" -> BASISLCx23,
   "x123" -> BASISLCx123, "x132" -> BASISLCx132|>;
plsExpr = Extract[shiftedBasis, onlyplanar[[Ordering[Flatten[onlyplanar]]]]];
plVars = Flatten[Table[crossHead[c][i], {c, Keys[crossHead]}, {i, 291}]];
plCoef = CoefficientArrays[plsExpr, plVars];
If[Length[plCoef] =!= 2 || Max[Abs[Normal[plCoef[[1]]]]] =!= 0, Print["shiftedBasis planar rows not linear"]; Abort[]];
plCoef = plCoef[[2]];
planarAt[w_, pt_] := Module[{R = SparseArray[Normal[plCoef] /. pt], blocks},
  blocks = Table[With[{F = PLGPLSparse[w, Keys[crossHead][[ci]]], cols = Range[291] + 291 (ci - 1)},
     (* PL's .mx may carry MZV in another context (written without PolyLogTools): unify *)
     {F[[1]] /. s_Symbol /; SymbolName[s] === "MZV" :> HPL`MZV, SparseArray[R[[All, cols]]] . F[[2]]}], {ci, Length[crossHead]}];
  mergeBlocks[blocks]];
planarConst[w_] := Module[{p1, p2, d, sel, k, S},
  p1 = planarAt[w, {x -> 3/13, y -> 5/18}]; p2 = planarAt[w, {x -> 7/23, y -> 2/11}];
  d = mergeBlocks[{p1, neg[p2]}];
  If[d[[1]] =!= {}, Print["planar combination is not constant in x,y: ", Length[d[[1]]], " keys differ"]; Abort[]];
  {k, S} = p1;
  sel = Flatten@Position[k, {{}, _, _}, {1}];      (* regularized x -> 0 *)
  ev1[{{#[[2]], #[[3]]} & /@ k[[sel]], S[[All, sel]]}]];  (* regularized y -> 1 *)

(* ---------- reference parsing (validation only) ---------- *)
refKeys[e_, var_] := Module[{ex, terms, f, gw, r, rat},
  ex = Expand[(e /. HoldPattern[G[__, 0]] -> 0) /. numRules];
  If[ex === 0, Return[{}]];
  If[Count[ex, G[__, var], {0, Infinity}] > 0 && ! FreeQ[ex, Power[G[__, var], _] | Times[___, G[__, var], ___, G[__, var], ___]],
    ex = Expand[ShuffleG[ex]]];
  terms = If[Head[ex] === Plus, List @@ ex, {ex}];
  Table[f = If[Head[tm] === Times, List @@ tm, {tm}];
    gw = Cases[f, G[a__, var] :> canon /@ {a}];
    r = Times @@ DeleteCases[f, G[__, var]];
    rat = r /. constPat -> 1;
    If[Length[gw] > 1 || ! NumericQ[rat], Print["refKeys: cannot parse ", tm]; Abort[]];
    {{First[gw, {}], r/rat}, canon[rat]}, {tm, terms}]];
toF[exprs_, var_] := Module[{rows, keys, idx},
  rows = Flatten[Table[{i, ##} & @@@ refKeys[exprs[[i]], var], {i, Length[exprs]}], 1];
  If[rows === {}, Return[emptyF]];
  keys = Union[rows[[All, 2]]]; idx = AssociationThread[keys -> Range[Length[keys]]];
  prune[{keys, SparseArray[{#1, idx[#2]} -> #3 & @@@ rows, {Length[exprs], Length[keys]}]}]];
segName[s_] := If[IntegerQ[s], ToString[s], s];

(* ---------- weight loop ---------- *)
sol0 = Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w0", "sol0.m"}]];
Do[sol[0][s] = toF[sol0, t], {s, segs}];
w0 = 1;
Do[If[FileExistsQ[FileNameJoin[{fastdir, "sol" <> ToString[w] <> ".mx"}]],
   solw = Import[FileNameJoin[{fastdir, "sol" <> ToString[w] <> ".mx"}], "MX"];
   Do[sol[w][s] = solw[s], {s, segs}]; w0 = w + 1, Break[]], {w, 1, wmax}];
If[w0 > 1, Print["loaded weights 1..", w0 - 1, " from ", fastdir]];

Do[
  Print["=== weight ", w, " ==="];
  tw = AbsoluteTime[];
  Do[IP[s] = tic["integrate " <> segName[s], integrate[s, sol[w - 1][s]]], {s, segs}];
  Print["  sizes (keys, nonzeros): ", Table[{Length[IP[s][[1]]], Length[IP[s][[2]]["NonzeroValues"]]}, {s, segs}]];

  impart[j_] := add @@ Prepend[Table[matMul[tdP[j, k], scal[sol[w - k][j], (-I)^k/k!, Pi^k]], {k, 1, w}], IP[j]];
  tic["matching",
    d = Table[ev1[IP[k - 1]] ~sub~ ev0[IP[k]], {k, 2, 5}]; d = Prepend[d, None]; (* d[[k]] = d_k for k=2..5 *)
    d6 = sub[ev1[IP[5]], ev0[IP[6]]];
    imp5 = impart[5]; imp3 = impart[3];
    d00 = sub[ev1[imp5], ev0[IP["xmy"]]];
    d01 = sub[ev1[imp3], ev0[IP["y1"]]];
    d10 = sub[ev0[imp3], ev0[IP["x1"]]]];
  delta[5] = emptyF; delta[4] = neg[d[[5]]]; delta[3] = sub[delta[4], d[[4]]];
  delta[2] = sub[delta[3], d[[3]]]; delta[1] = sub[delta[2], d[[2]]];
  delta[6] = d6; delta["xmy"] = d00; delta["y1"] = add[delta[3], d01]; delta["x1"] = add[delta[3], d10];

  (* boundary equations: Mat.(IP[seg] + c5 + delta[seg]) = 0 *)
  eqs = {{"P1", 5}, {"P2", 1}, {"P7", 3}, {"R1", "x1"}, {"R2", "y1"}, {"R3", "xmy"}};
  rhsB = tic["boundary rhs", Table[With[{full = matMul[bm[e[[1]]], IP[e[[2]]]]},
      If[nonConst[full] > 0, Print["  WARNING: ", e[[1]], ".tempSol[", e[[2]], "] depends on the variable (", nonConst[full], " keys)"]];
      neg[add[ev0[full], matMul[bm[e[[1]]], delta[e[[2]]]]]]], {e, eqs}]];
  pl = tic["planar constants", planarConst[w]];
  monos = Union[Join @@ (#[[1]] & /@ Append[rhsB, pl])];
  If[monos === {}, monos = {{{}, 1}}];
  Print["  transcendental columns: ", monos[[All, 2]], "; per block: ", (#[[1]][[All, 2]] & /@ Append[rhsB, pl])];  (* homogeneous system: keep one (zero) column *)
  rhs = Join @@ (align[#, monos] & /@ rhsB);
  plV = If[pl[[1]] === {}, SparseArray[{}, {Length[op], Length[monos]}], align[pl, monos]];
  rhs2 = rhs - Ap . plV;
  xs = tic["solve", Expand[invA . Normal[rhs2[[pivRows]]]]];
  resid = Expand[Au . xs - rhs2];
  nres = Count[Flatten[Normal[resid]], Except[0]];
  Print["  boundary system: ", Length[Amat], " equations, ", Length[unk], " unknowns, ", Length[monos],
   " transcendental columns; nonzero residuals: ", nres];
  If[nres > 0 && w <= 6,
    Print["  --- diagnosing against reference constants ---"];
    refC[seg_] := Module[{r = Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w" <> ToString[w], "const" <> ToString[w] <> "_" <> segName[seg] <> ".m"}]], vals},
      vals = Table[cnst[w, i][If[IntegerQ[seg], seg, ToExpression[seg]]], {i, n}] /. r;
      toF[vals /. _cnst -> 0, t]];
    rc5 = refC[5];
    plRef = With[{f = rc5}, {f[[1]], f[[2]][[op]]}];
    With[{dp = mergeBlocks[{pl, neg[plRef]}]},
      Print["  planar constants vs reference c5[onlyplanar]: ", Length[dp[[1]]], " differing keys"];
      If[dp[[1]] =!= {}, With[{rows = Union[dp[[2]]["NonzeroPositions"][[All, 1]]]},
        Print["  differing planar rows (basis index): ", op[[rows]]];
        Print["  mine: ", toExprs[pl, t][[rows[[;; UpTo[6]]]]], "  ref: ", toExprs[plRef, t][[rows[[;; UpTo[6]]]]]]]]];
    Do[dd = sub[sub[refC[s], rc5], delta[s]];
      Print["  delta[", segName[s], "] vs reference (c_seg - c5): ", If[dd[[1]] === {}, "identical", {dd[[1]], Count[Normal[dd[[2]]], Except[0], {2}]}]], {s, segs}];
    resRef = Table[sub[matMul[bm[e[[1]]], rc5], rhsB[[Position[eqs, e][[1, 1]]]]], {e, eqs}];
    Print["  residual with reference c5, per block (nonzero entries): ", Count[Normal[#[[2]]], Except[0], {2}] & /@ resRef]];
  If[nres > 0, Print["  inconsistent boundary system; xs dims ", Dimensions[xs], " resid sample ", Short[Normal[resid][[;; 5]], 3],
     " block residual counts ", Table[Count[Flatten[Normal[resid[[371 (b - 1) + 1 ;; 371 b]]]], Except[0]], {b, 6}]]; Abort[]];
  c5M = SparseArray[{}, {n, Length[monos]}];
  c5M = Normal[c5M]; c5M[[unk]] = xs; c5M[[op]] = Normal[plV];
  c5 = prune[{monos, SparseArray[c5M]}];
  Do[sol[w][s] = add[IP[s], c5, delta[s]], {s, segs}];
  Export[FileNameJoin[{fastdir, "sol" <> ToString[w] <> ".mx"}], Association[Table[s -> sol[w][s], {s, segs}]], "MX"];
  Print["  weight ", w, " done in ", Round[AbsoluteTime[] - tw, 0.1], " s; sizes: ",
   Table[Length[sol[w][s][[2]]["NonzeroValues"]], {s, segs}]];

  If[validate && w <= 6,
    Do[ref = toF[Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w" <> ToString[w], "sol" <> ToString[w] <> "_" <> segName[s] <> ".m"}]], segVar[s]];
      dd = sub[sol[w][s], ref];
      Print["  validate sol", w, "_", segName[s], ": ", If[dd[[1]] === {}, "identical", ToString[Length[dd[[1]]]] <> " differing keys, e.g. " <> ToString[dd[[1, ;; UpTo[3]]], InputForm]]],
      {s, segs}]],
  {w, w0, wmax}];
