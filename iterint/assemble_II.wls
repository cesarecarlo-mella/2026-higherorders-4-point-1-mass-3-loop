#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Assembles the run_iterint.wls output of one family into the all.ints/iterint_solutions format:
     {BASISLC[fam, "", 1] -> Sum_{k<=wmax} w_k[[1]], ..., BASISLC[fam, "x132", n] -> ...}
   crossings in the order "", x12, x13, x23, x123, x132 (as in all.ints), no ep (the weight is
   the length of the II word). File name as in all.ints: H.wl, I3.wl, otherwise <fam>_II.m.
   The file is written entry by entry, so the full expression is never held in memory.

   Arguments (key=value):
     family=H
     in=<dir>      run_iterint.wls output (default out/<fam>)
     outdir=<dir>  (default II)
     wmax=8 *)

Print["kernel started"];   (* run_all.sh watches for this line *)
$HistoryLength = 0;
args = Association[Cases[$ScriptCommandLine[[2 ;;]], s_String /; StringContainsQ[s, "="] :>
     Rule @@ StringSplit[s, "=", 2]]];
need[k_] := If[KeyExistsQ[args, k], args[k], Print["missing argument ", k, "="]; Exit[1]];
famName = need["family"]; fam = ToExpression[famName];   (* "D" is System`D, as in the input files *)
in = Lookup[args, "in", FileNameJoin[{"out", famName}]];
outDir = Lookup[args, "outdir", "II"]; Quiet[CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
wmax = ToExpression[Lookup[args, "wmax", "8"]];
outName = If[MemberQ[{"H", "I3"}, famName], famName <> ".wl", famName <> "_II.m"];
outFile = FileNameJoin[{outDir, outName}];

crosses = {"x", "x12", "x13", "x23", "x123", "x132"};
label[c_] := If[c === "x", "", c];
file[c_, k_] := FileNameJoin[{in, famName <> "_" <> c <> "_w" <> ToString[k] <> ".wdx"}];
missing = Select[Flatten[Table[file[c, k], {c, crosses}, {k, 0, wmax}]], ! FileExistsQ[#] &];
If[missing =!= {}, Print[Length[missing], " weight files missing, e.g. ", missing[[1]]]; Exit[1]];

(* row i of the weight-k piece, built only from the keys that row uses *)
rowExpr[{keys_, S_}, i_] := If[keys === {}, 0,
   With[{r = S[[i]]}, r["NonzeroValues"] . ((Global`II @@ #[[1]]) #[[2]] & /@ keys[[Flatten[r["NonzeroPositions"]]]])]];

t0 = AbsoluteTime[];
str = OpenWrite[outFile, PageWidth -> Infinity];
WriteString[str, "(* Created with the Wolfram Language : www.wolfram.com *)\n{"];
first = True; nTot = 0;
Do[
  tc = AbsoluteTime[];
  data = Table[Import[file[c, k]], {k, 0, wmax}];
  n = Length[data[[1, 2]]];
  If[Union[Length[#[[2]]] & /@ data] =!= {n}, Print[c, ": weights have different numbers of masters"]; Exit[1]];
  Do[
    If[! first, WriteString[str, ", \n "]]; first = False;
    WriteString[str, ToString[Global`BASISLC[fam, label[c], i] -> Sum[rowExpr[data[[k + 1]], i], {k, 0, wmax}], InputForm]],
    {i, n}];
  nTot += n;
  Print[famName, " ", c, ": ", n, " entries written (", Round[AbsoluteTime[] - tc], " s)"];
  Clear[data],
  {c, crosses}];
WriteString[str, "}\n"];
Close[str];
Print["wrote ", outFile, ": ", nTot, " entries, weights 0..", wmax, ", ",
  Round[FileByteCount[outFile]/10.^6], " MB, ", Round[AbsoluteTime[] - t0], " s"];
