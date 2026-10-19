#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Iterated-integral solution of family H, as in reorganize.nb (end of notebook):
     w_0 = BC_0 II[],   w_k = BC_k II[] + Atilde_c . w_{k-1}  with  II[a..] l[b] -> II[a.., b],
     solution(c) = Sum_k w_k,  for each crossing c, Atilde_c = H_atilde /. <c>lettersrule.
   Inputs: H_opt/input/H_atilde.m, H_opt/input/H_lettersrules.m (copied from reorganize.nb),
           H_opt/results/H_bc.m (export_fast.wls, eps^0..8).
   Done as sparse linear algebra: II-words are column keys, appending b = A_b . S.

   Usage: wolframscript -file run_iterint.wls <wmax> [validate] [m] [cross=x,x12,...]
   ("m": also write the text file; weight 8 is ~50M terms per crossing, so restrict crossings)
   Writes H_opt/work/fast/II_<c>.mx (Association k -> {keys, S}) and H_opt/results/H_II.m
   ({BASISLC[H, c, i] -> Sum_k w_k}, same format as all.ints/iterint_solutions/H.wl).
   "validate" compares Sum_{k<=6} with all.ints/iterint_solutions/H.wl. *)

$HistoryLength = 0;
Global`ep; Global`H;
SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];
root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
hopt = FileNameJoin[{root, "H_opt"}];
wmax = ToExpression[$ScriptCommandLine[[2]]];
validate = MemberQ[$ScriptCommandLine[[3 ;;]], "validate"];
writeM = MemberQ[$ScriptCommandLine[[3 ;;]], "m"];
oldII = "/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/all.ints/iterint_solutions/H.wl";

constPat = Pi | _Zeta | _MZV;
atilde = Import[FileNameJoin[{hopt, "input", "H_atilde.m"}]];
n = Length[atilde];
lrules = Append[Get[FileNameJoin[{hopt, "input", "H_lettersrules.m"}]], "" -> {}];
crossings = {"", "x12", "x13", "x23", "x123", "x132"};
(* cross=x,x12,... restricts the crossings ("x" = uncrossed) *)
With[{sel = StringCases[$ScriptCommandLine[[3 ;;]], "cross=" ~~ c__ :> c] // Flatten},
  If[sel =!= {}, crossings = StringSplit[First[sel], ","] /. "x" -> ""]];
bcA = Association[Import[FileNameJoin[{hopt, "results", "H_bc.m"}]]];

emptyF = {{}, SparseArray[{}, {n, 0}]};
prune[{keys_, S_}] := Module[{p, S2, cols},
  If[keys === {}, Return[emptyF]];
  p = S["NonzeroPositions"]; If[p === {}, Return[emptyF]];
  S2 = SparseArray[Thread[p -> Expand[S["NonzeroValues"]]], Dimensions[S]];
  cols = Union[S2["NonzeroPositions"][[All, 2]]];
  If[cols === {}, emptyF, {keys[[cols]], S2[[All, cols]]}]];
mergeBlocks[blocks0_List] := Module[{blocks = Select[blocks0, #[[1]] =!= {} &], keys, mats, uk, idx, P},
  If[blocks === {}, Return[emptyF]];
  keys = Join @@ blocks[[All, 1]]; mats = Join[Sequence @@ blocks[[All, 2]], 2];
  uk = Union[keys]; idx = AssociationThread[uk -> Range[Length[uk]]];
  P = SparseArray[Table[{k, idx[keys[[k]]]} -> 1, {k, Length[keys]}], {Length[keys], Length[uk]}];
  prune[{uk, mats . P}]];
(* expressions linear in II[...] with coefficients in Q(i)[Pi, Zeta, MZV] -> {keys, S} *)
toF[exprs_] := Module[{rows, keys, idx, terms, f, w, r, rat},
  rows = Flatten[Table[With[{e = Expand[exprs[[i]]]},
      If[e === 0, {}, terms = If[Head[e] === Plus, List @@ e, {e}];
       Table[f = If[Head[tm] === Times, List @@ tm, {tm}];
        w = Cases[f, II[a___] :> {a}]; r = Times @@ DeleteCases[f, _II]; rat = r /. constPat -> 1;
        If[Length[w] =!= 1 || ! NumericQ[rat], Print["cannot parse term ", tm]; Abort[]];
        {i, {w[[1]], r/rat}, rat}, {tm, terms}]]], {i, Length[exprs]}], 1];
  If[rows === {}, Return[emptyF]];
  keys = Union[rows[[All, 2]]]; idx = AssociationThread[keys -> Range[Length[keys]]];
  prune[{keys, SparseArray[{#1, idx[#2]} -> #3 & @@@ rows, {Length[exprs], Length[keys]}]}]];
toExprs[{k_, S_}] := Module[{vec = (II @@ #[[1]]) #[[2]] & /@ k},
  If[k === {}, ConstantArray[0, Length[S]],
   Table[With[{r = S[[i]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]], {i, Length[S]}]]];

Do[
  tc = AbsoluteTime[];
  effA = Expand[atilde /. lrules[c]];
  lets = Union@Cases[effA, l[b_] :> b, Infinity];
  Amat = Association[Table[b -> SparseArray[Coefficient[effA, l[b]]], {b, lets}]];
  If[Expand[effA - Sum[Normal[Amat[b]] l[b], {b, lets}]] =!= ConstantArray[0, {n, n}], Print["Atilde not linear in l[..]"]; Abort[]];
  bcv = Table[bcA[BC[H, c, i]], {i, n}];
  mxFile = FileNameJoin[{hopt, "work", "fast", "II_" <> If[c === "", "x", c] <> ".mx"}];
  wk = If[FileExistsQ[mxFile], Import[mxFile, "MX"], <|0 -> toF[Coefficient[bcv, ep, 0] II[]]|>];
  Do[If[KeyExistsQ[wk, k], Continue[]]; tk = AbsoluteTime[];
    wk[k] = mergeBlocks[Append[
       Table[{{Append[#[[1]], b], #[[2]]} & /@ wk[k - 1][[1]], Amat[b] . wk[k - 1][[2]]}, {b, lets}],
       toF[Coefficient[bcv, ep, k] II[]]]];
    Export[mxFile, wk, "MX"];
    Print["   weight ", k, ": ", Length[wk[k][[2]]["NonzeroValues"]], " nonzeros, ", Round[AbsoluteTime[] - tk, 0.1], " s, memory in use ",
     Round[MemoryInUse[]/10.^9, 0.1], " GB"],
    {k, 1, wmax}];
  Print["crossing \"", c, "\": ", Round[AbsoluteTime[] - tc, 0.1], " s; nonzeros per weight ",
   Table[Length[wk[k][[2]]["NonzeroValues"]], {k, 0, wmax}]];
  If[validate,
    If[! ValueQ[oldSol], oldSol = Association[Get[oldII]]];
    ref = toF[Table[oldSol[BASISLC[H, c, i]], {i, n}]];
    mine = mergeBlocks[Table[wk[k], {k, 0, Min[6, wmax]}]];
    d = mergeBlocks[{mine, {ref[[1]], -ref[[2]]}}];
    Print["   sum_{k<=6} vs H.wl: ", If[d[[1]] === {}, "identical", ToString[Length[d[[1]]]] <> " differing keys"]]];
  If[writeM, sols[c] = Total[Table[toExprs[wk[k]], {k, 0, wmax}]]];
  Clear[wk, Amat, effA],
  {c, crossings}];

If[! writeM, Exit[0]];
full = Flatten[Table[Table[BASISLC[H, c, i] -> sols[c][[i]], {i, n}], {c, crossings}]];
{t, ok} = AbsoluteTiming[Export[FileNameJoin[{hopt, "results", "H_II.m"}], full]];
Print["wrote H_opt/results/H_II.m (", Round[t, 0.1], " s, ", Round[FileByteCount[FileNameJoin[{hopt, "results", "H_II.m"}]]/10.^6], " MB)"];
