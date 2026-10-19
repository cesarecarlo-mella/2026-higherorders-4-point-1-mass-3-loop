#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Iterated-integral solution of one family, one crossing (reorganize.nb, end of notebook):
     w_0 = BC_0 II[],   w_k = BC_k II[] + Atilde_c . w_{k-1}   with  II[a..] l[b] -> II[a.., b],
   Atilde_c = atilde /. <c>lettersrule. Self-contained (no PolyLogTools, no repo paths), meant
   to run as one cluster job per crossing.

   Arguments (key=value):
     family=H            family symbol used in BC[family, crossing, i] and BASISLC[family, crossing, i]
     atilde=<file>       atilde in letters l[1..20] (all.ints/a_tildes/<fam>_atilde.m)
     bc=<file>           boundary conditions {BC[fam, crossing, i] -> series in ep} (H_opt/export_fast.wls)
     cross=x12           one of x (uncrossed), x12, x13, x23, x123, x132
     wmax=8
     out=<dir>           output directory (created if missing)
     letters=<file>      crossing letter rules (default: input/lettersrules.m next to this script)
     text=yes            also write <out>/<fam>_<cross>_w<k>.m, the weight-k piece as
                         {BASISLC[fam, crossing, i] -> w_k[[i]]}  (default: no)
     validate=<file>     optional: old II file ({BASISLC[..] -> Sum_k w_k}); compares Sum_{k<=vmax} w_k
     vmax=6

   Output: <out>/<fam>_<cross>_w<k>.wdx = {keys, S} for each k (platform-independent; the run
   resumes from the highest weight found). keys = {II-word, transcendental monomial},
   w_k[[i]] = Sum_j S[[i,j]] II @@ keys[[j,1]] keys[[j,2]]. Read with load_iterint.wl. *)

$HistoryLength = 0;
SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];
args = Association[Cases[$ScriptCommandLine[[2 ;;]], s_String /; StringContainsQ[s, "="] :>
     Rule @@ StringSplit[s, "=", 2]]];
need[k_] := If[KeyExistsQ[args, k], args[k], Print["missing argument ", k, "="]; Exit[1]];
(* resolve the family name exactly as the input files do: "D" is System`D (Mathematica's
   derivative), "H" is Global`H *)
famName = need["family"]; fam = ToExpression[famName];
cross = need["cross"]; crossLabel = If[cross === "x", "", cross];
wmax = ToExpression[need["wmax"]];
out = need["out"]; Quiet[CreateDirectory[out, CreateIntermediateDirectories -> True]];
lettersFile = Lookup[args, "letters", FileNameJoin[{DirectoryName[$InputFileName], "input", "lettersrules.m"}]];
writeText = Lookup[args, "text", "no"] === "yes";
ep = Global`ep;

constPat = Pi | _Zeta | _MZV;
atilde = Import[need["atilde"]];
n = Length[atilde];
lrules = Append[Get[lettersFile], "" -> {}];
If[! KeyExistsQ[lrules, crossLabel], Print["unknown crossing ", cross]; Exit[1]];
bcA = Association[Import[need["bc"]]];
bcv = Table[bcA[Global`BC[fam, crossLabel, i]], {i, n}];
If[! FreeQ[bcv, Global`BC | Missing], Print["bc file has no BC[", famName, ", \"", crossLabel, "\", i] for some i"]; Exit[1]];

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
        w = Cases[f, Global`II[a___] :> {a}]; r = Times @@ DeleteCases[f, _Global`II]; rat = r /. constPat -> 1;
        If[Length[w] =!= 1 || ! NumericQ[rat], Print["cannot parse term ", tm]; Exit[1]];
        {i, {w[[1]], r/rat}, rat}, {tm, terms}]]], {i, Length[exprs]}], 1];
  If[rows === {}, Return[emptyF]];
  keys = Union[rows[[All, 2]]]; idx = AssociationThread[keys -> Range[Length[keys]]];
  prune[{keys, SparseArray[{#1, idx[#2]} -> #3 & @@@ rows, {Length[exprs], Length[keys]}]}]];
toExprs[{k_, S_}] := Module[{vec = (Global`II @@ #[[1]]) #[[2]] & /@ k},
  If[k === {}, ConstantArray[0, Length[S]],
   Table[With[{r = S[[i]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]], {i, Length[S]}]]];

effA = Expand[atilde /. lrules[crossLabel]];
lets = Union@Cases[effA, Global`l[b_] :> b, Infinity];
Amat = Association[Table[b -> SparseArray[Coefficient[effA, Global`l[b]]], {b, lets}]];
If[Expand[effA - Sum[Normal[Amat[b]] Global`l[b], {b, lets}]] =!= ConstantArray[0, {n, n}],
  Print["atilde is not linear in l[..]"]; Exit[1]];
Print[famName, " crossing ", cross, ": ", n, " masters, letters ", lets];

file[k_, ext_] := FileNameJoin[{out, famName <> "_" <> cross <> "_w" <> ToString[k] <> "." <> ext}];
wk = <||>;
Do[
  If[FileExistsQ[file[k, "wdx"]], wk[k] = Import[file[k, "wdx"]]; Print["w", k, ": loaded"]; Continue[]];
  tk = AbsoluteTime[];
  wk[k] = If[k == 0, toF[Coefficient[bcv, ep, 0] Global`II[]],
    mergeBlocks[Append[
      Table[{{Append[#[[1]], b], #[[2]]} & /@ wk[k - 1][[1]], Amat[b] . wk[k - 1][[2]]}, {b, lets}],
      toF[Coefficient[bcv, ep, k] Global`II[]]]]];
  Export[file[k, "wdx"], wk[k]];
  If[writeText, Export[file[k, "m"], Thread[Table[Global`BASISLC[fam, crossLabel, i], {i, n}] -> toExprs[wk[k]]]]];
  Print["w", k, ": ", Length[wk[k][[2]]["NonzeroValues"]], " terms, ", Round[AbsoluteTime[] - tk, 0.1], " s, ",
   Round[MaxMemoryUsed[]/10.^9, 0.1], " GB peak"];
  If[k > 1, KeyDropFrom[wk, k - 2]],   (* only w_{k-1} is needed for the next step *)
  {k, 0, wmax}];

If[KeyExistsQ[args, "validate"],
  vmax = ToExpression[Lookup[args, "vmax", "6"]];
  old = Association[Get[args["validate"]]];
  ref = toF[Table[old[Global`BASISLC[fam, crossLabel, i]], {i, n}]];
  mine = mergeBlocks[Table[If[KeyExistsQ[wk, k], wk[k], Import[file[k, "wdx"]]], {k, 0, vmax}]];
  d = mergeBlocks[{mine, {ref[[1]], -ref[[2]]}}];
  Print["Sum_{k<=", vmax, "} vs ", args["validate"], ": ", If[d[[1]] === {}, "identical", ToString[Length[d[[1]]]] <> " differing keys"]]];
