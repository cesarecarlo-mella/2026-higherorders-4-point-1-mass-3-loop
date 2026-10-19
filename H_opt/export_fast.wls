#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Exports the output of run_fast.wls (H_opt/work/fast/sol<w>.mx) in the formats used so far:
   - H_opt/results/eps_sep/w<w>/sol<w>_<seg>.m   : 371 expressions in G[...,t] (x1, xmy: y; y1: x)
   - H_opt/results/eps_sep/w<w>/const<w>_<seg>.m : {cnst[w,i][seg] -> value, i = 1..371}
   - H_opt/results/H_bc.m : {BC[H, crossing, i] -> sol0 + Sum_w ep^w const_w(corner)},
     corners as in reorganize.nb: "",x23 -> 6; x12,x132 -> 4; x13,x123 -> 2.
   Cube-root / c2 letters are written back as the symbols cxmy, cxmyc, c2, c2c (as in H/).

   Usage: wolframscript -file export_fast.wls <wmax> [validate]
   "validate" compares const and BC files with H/solutions/eps_sep and the existing H_bc.m (w <= 6). *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
Global`x; Global`y; Global`t; Global`c2; Global`c2c; Global`cxmy; Global`cxmyc; Global`ep; Global`H;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];

root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
fastdir = FileNameJoin[{root, "H_opt", "work", "fast"}];
outdir = FileNameJoin[{root, "H_opt", "results"}];
wmax = ToExpression[$ScriptCommandLine[[2]]];
validate = MemberQ[$ScriptCommandLine[[3 ;;]], "validate"];

numRules = {c2 -> (1 - I)/2, c2c -> (1 + I)/2, cxmy -> -(-1)^(1/3), cxmyc -> (-1)^(2/3)};
canon[v_] := Expand[ComplexExpand[v /. numRules]];
backLetter = Association[Table[canon[s] -> s, {s, {c2, c2c, cxmy, -cxmy, cxmyc, -cxmyc}}]];
segs = {1, 2, 3, 4, 5, 6, "x1", "y1", "xmy"};
segVar = <|1 -> t, 2 -> t, 3 -> t, 4 -> t, 5 -> t, 6 -> t, "x1" -> y, "y1" -> x, "xmy" -> y|>;
segName[s_] := If[IntegerQ[s], ToString[s], s];
segSym[s_] := If[IntegerQ[s], s, Symbol[s]];

toExprs[{k_, S_}, var_] := Module[{vec},
  If[k === {}, Return[ConstantArray[0, Length[S]]]];
  vec = (If[#[[1]] === {}, 1, G[Sequence @@ (Lookup[backLetter, #, #] & /@ #[[1]]), var]] #[[2]]) & /@ k;
  Table[With[{r = S[[i]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]], {i, Length[S]}]];
constPart[{k_, S_}] := Module[{sel = Flatten@Position[k, {{}, _}, {1}]}, toExprs[{k[[sel]], S[[All, sel]]}, t]];

sol0 = Import[FileNameJoin[{root, "H", "solutions", "eps_sep", "w0", "sol0.m"}]];
n = Length[sol0];
Do[
  sw = Import[FileNameJoin[{fastdir, "sol" <> ToString[w] <> ".mx"}], "MX"];
  dir = FileNameJoin[{outdir, "eps_sep", "w" <> ToString[w]}]; Quiet[CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  Do[
    Export[FileNameJoin[{dir, "sol" <> ToString[w] <> "_" <> segName[s] <> ".m"}], toExprs[sw[s], segVar[s]]];
    cst[w, s] = constPart[sw[s]];
    Export[FileNameJoin[{dir, "const" <> ToString[w] <> "_" <> segName[s] <> ".m"}],
      Table[cnst[w, i][segSym[s]] -> cst[w, s][[i]], {i, n}]],
    {s, segs}];
  Print["w", w, ": exported sol/const for all segments to ", dir],
  {w, 1, wmax}];

corner = <|"" -> 6, "x12" -> 4, "x13" -> 2, "x23" -> 6, "x123" -> 2, "x132" -> 4|>;
bc = Flatten[Table[Table[BC[H, c, i] -> sol0[[i]] + Sum[ep^w cst[w, corner[c]][[i]], {w, 1, wmax}], {i, n}], {c, Keys[corner]}]];
Export[FileNameJoin[{outdir, "H_bc.m"}], bc];
Print["exported ", FileNameJoin[{outdir, "H_bc.m"}], " (", Length[bc], " entries, eps^0..", wmax, ")"];

If[validate,
  Do[With[{ref = Import[FileNameJoin[{root, "H", "solutions", "eps_sep", "w" <> ToString[w], "const" <> ToString[w] <> "_" <> segName[s] <> ".m"}]]},
      bad = Count[Table[Expand[canon[(cnst[w, i][segSym[s]] /. ref /. _cnst -> 0) /. HoldPattern[G[__, 0]] -> 0] - canon[cst[w, s][[i]]]], {i, n}], Except[0]];
      If[bad > 0, Print["  const", w, "_", segName[s], ": ", bad, " entries differ from H/"]]],
    {w, 1, Min[6, wmax]}, {s, segs}];
  Print["const files w<=6 compared with H/solutions/eps_sep (only differences are printed)"];
  old = Import["/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/all.ints/boundary_conditions/H/H_bc.m"];
  newA = Association[bc];
  bad = Select[old, Expand[(#[[2]] - Normal[Series[newA[#[[1]]], {ep, 0, 6}]]) /. HoldPattern[G[__, 0]] -> 0] =!= 0 &];
  Print["H_bc.m vs existing boundary_conditions/H/H_bc.m through ep^6: ", Length[old], " entries, ", Length[bad], " differ",
   If[bad =!= {}, ", e.g. " <> ToString[bad[[1, 1]], InputForm], ""]]];
