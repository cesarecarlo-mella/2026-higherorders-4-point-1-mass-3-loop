#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Writes ../planar_stuff_w8/pl<cross>.m in exactly the format of ../planar_stuff/:
     {BASISLC<cross>[i] -> sum_{w=0}^{8} eps^w sol[w][[i]], i = 1..291}
   from the sparse output of run_gpl_xy_fast.wls (run that first, to weight 8, for
   every crossing). Then re-imports the written file and checks it:
   eps^0..6 against ../planar_stuff/pl<cross>.m, eps^7..8 against the sparse data.

   Usage: wolframscript -file export_planar_w8.wls <cross>   (x, x12, x13, x23, x123, x132) *)

SetDirectory[DirectoryName[$InputFileName]];
Get["./load_gpl_xy.wl"];
cross = $ScriptCommandLine[[2]];
suffix = If[cross === "x", "", cross];
head = Symbol["BASISLC" <> suffix];
outDir = "../planar_stuff_w8"; Quiet[CreateDirectory[outDir]];
outFile = FileNameJoin[{outDir, "pl" <> suffix <> ".m"}];
wmax = 8;

w0 = Import["./work/weight0_solution.m"];
{t, sol} = AbsoluteTiming[Table[PLGPLExprs[w, cross], {w, 1, wmax}]];
Print[cross, ": rebuilt weights 1..", wmax, " (", Round[t, 0.1], " s)"];
full = Table[head[i] -> w0[[i]] + Sum[eps^w sol[[w, i]], {w, 1, wmax}], {i, 291}];
{t, ok} = AbsoluteTiming[Export[outFile, full]];
Print["wrote ", outFile, " (", Round[t, 0.1], " s, ", Round[FileByteCount[outFile]/10.^6, 0.1], " MB)"];

(* ---- checks on the file as written ---- *)
back = Import[outFile];
Print["heads/order as in planar_stuff: ", back[[All, 1]] === Table[head[i], {i, 291}]];
ref = Import["../planar_stuff/pl" <> suffix <> ".m"];
Print["same heads as reference: ", back[[All, 1]] === ref[[All, 1]]];
back = back[[All, 2]]; ref = ref[[All, 2]];
Print["max eps power: ", Max[Exponent[#, eps] & /@ back]];
Do[Print["eps^", w, " vs planar_stuff: ",
   If[Union[Expand[Coefficient[back, eps, w] - Coefficient[ref, eps, w]]] === {0}, "identical", "DIFFERENT"]], {w, 0, 6}];
Do[Print["eps^", w, " vs sparse data: ",
   If[Union[Expand[Coefficient[back, eps, w] - sol[[w]]]] === {0}, "identical", "DIFFERENT"]], {w, 7, wmax}];
