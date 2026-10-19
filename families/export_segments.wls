#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Segment solutions of solve_bc.wls in the format of H/solutions/eps_sep:
     segments/<fam>/eps_sep/w0/sol0.m
     segments/<fam>/eps_sep/w<w>/sol<w>_<seg>.m   : 371-style list of expressions in G[..., t]
                                                     (x1, xmy: G[..., y]; y1: G[..., x])
     segments/<fam>/eps_sep/w<w>/const<w>_<seg>.m : {cnst[w, i][seg] -> constant part}
   seg in {1..6, x1, y1, xmy}. Letters (1-+I)/2 and the cube roots of unity are written back as
   c2, c2c, cxmy, cxmyc (as in H/); other letters stay explicit numbers.
   Input: work/<fam>/sol<w>.mx (Association seg -> {keys, S}), weight 0 from boundaries_all.
   Usage: wolframscript -file export_segments.wls fam=C,D,... [wmax=8] *)

$HistoryLength = 0;
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "config.wl"}]];
args = Association[Cases[$ScriptCommandLine[[2 ;;]], s_String /; StringContainsQ[s, "="] :> Rule @@ StringSplit[s, "=", 2]]];
fams = StringSplit[args["fam"], ","];
wmax = ToExpression[Lookup[args, "wmax", "8"]];
ep = Global`ep;

canon[v_] := Expand[ComplexExpand[v]];
numRules = {c2 -> (1 - I)/2, c2c -> (1 + I)/2, cxmy -> -(-1)^(1/3), cxmyc -> (-1)^(2/3)};
backLetter = Association[Table[canon[s /. numRules] -> s, {s, {c2, c2c, cxmy, -cxmy, cxmyc, -cxmyc}}]];
segs = {1, 2, 3, 4, 5, 6, "x1", "y1", "xmy"};
segVar = <|1 -> t, 2 -> t, 3 -> t, 4 -> t, 5 -> t, 6 -> t, "x1" -> y, "y1" -> x, "xmy" -> y|>;
segName[s_] := If[IntegerQ[s], ToString[s], s];
segSym[s_] := If[IntegerQ[s], s, Symbol[s]];
toExprs[{k_, S_}, var_] := Module[{vec},
  If[k === {}, Return[ConstantArray[0, Length[S]]]];
  vec = (If[#[[1]] === {}, 1, G[Sequence @@ (Lookup[backLetter, #, #] & /@ #[[1]]), var]] #[[2]]) & /@ k;
  Table[With[{r = S[[i]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]], {i, Length[S]}]];
constPart[{k_, S_}] := Module[{sel = Flatten@Position[k, {{}, _}, {1}]},
  If[sel === {}, ConstantArray[0, Length[S]], S[[All, sel]] . (k[[sel, 2]])]];

Do[
  tf = AbsoluteTime[];
  fsym = ToExpression[f];
  n = Length[Import[canonFile[f]]];
  base = FileNameJoin[{here, "segments", f, "eps_sep"}];
  bcOld = Association[Import[bcFile[f]]];
  Quiet[CreateDirectory[FileNameJoin[{base, "w0"}], CreateIntermediateDirectories -> True]];
  Export[FileNameJoin[{base, "w0", "sol0.m"}], Table[Coefficient[bcOld[BC[fsym, "", i]], ep, 0], {i, n}]];
  done = 0;
  Do[
    file = FileNameJoin[{here, "work", f, "sol" <> ToString[w] <> ".mx"}];
    If[! FileExistsQ[file], Print[f, ": no weight ", w, " checkpoint, stopping"]; Break[]];
    sw = Import[file, "MX"];
    dir = FileNameJoin[{base, "w" <> ToString[w]}]; Quiet[CreateDirectory[dir]];
    Do[
      Export[FileNameJoin[{dir, "sol" <> ToString[w] <> "_" <> segName[s] <> ".m"}], toExprs[sw[s], segVar[s]]];
      Export[FileNameJoin[{dir, "const" <> ToString[w] <> "_" <> segName[s] <> ".m"}],
        Thread[Table[cnst[w, i][segSym[s]], {i, n}] -> Expand[constPart[sw[s]]]]],
      {s, segs}];
    done = w,
    {w, 1, wmax}];
  Print[f, ": exported weights 0..", done, " to segments/", f, "/eps_sep (", Round[AbsoluteTime[] - tf, 0.1], " s)"],
  {f, fams}];
