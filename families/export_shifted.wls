#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Readable shifted bases: results/<fam>_shiftedBasis.m = list over the canonical basis
   (bases/<fam>_canon.m order) of each entry written as
     Sum coeff(x, y, d) BASISLC<c>[k]        (c = "", x12, x13, x23, x123, x132)
   plus, for non-planar entries, the family's own integrals INT[fam, ...], crossed copies of
   the family (INT[famx12, ...]) or PL integrals with no BASISLC reduction (INT[PLx13, ...]).
   Same construction as planar_bc.wls / solve_bc.wls: shift map all.ints/result/<fam>_shift.m
   if present, else SLC/Solutions/<fam>/shiftedBasis; then the PL -> BASISLC tables; then
   m2 -> 1, s12 -> x, s13 -> y, s23 -> 1-x-y.

   Usage: wolframscript -file export_shifted.wls [fam=A,C,...]   (default: every family with planar data) *)

$HistoryLength = 0;
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "config.wl"}]];
args = Association[Cases[$ScriptCommandLine[[2 ;;]], s_String /; StringContainsQ[s, "="] :> Rule @@ StringSplit[s, "=", 2]]];
fams = If[KeyExistsQ[args, "fam"], StringSplit[args["fam"], ","],
   Select[families, ! MissingQ[shiftMapFile[#]] || ! MissingQ[slcShiftedBasisFile[#]] &]];
outDir = FileNameJoin[{here, "results"}]; Quiet[CreateDirectory[outDir]];

crossLabels = {"", "x12", "x13", "x23", "x123", "x132"};
cross = <|"" -> {y -> y, x -> x}, "x12" -> {y -> 1 - x - y, x -> x}, "x13" -> {y -> y, x -> 1 - x - y},
   "x23" -> {y -> x, x -> y}, "x123" -> {y -> x, x -> 1 - x - y}, "x132" -> {y -> 1 - x - y, x -> y}|>;
bl[c_] := ToExpression["BASISLC" <> c];
plTab = Flatten[Get /@ plTableFiles] /. {s23 -> 1 - y - x, s13 -> y};
plRules = Dispatch[Flatten[Table[plTab /. PL[a__] :> INT[ToExpression["PL" <> c], a] /. BASISLC -> bl[c] /. cross[c], {c, crossLabels}]]];
blPat = h_[_] /; StringStartsQ[SymbolName[h], "BASISLC"];

Do[
  fsym = ToExpression[f];
  canonB = Import[canonFile[f]];
  {src, shifted} = If[! MissingQ[shiftMapFile[f]],
    Module[{sh = Import[shiftMapFile[f]] /. INT[lb_String, __, {b__}] :> INT[lb, b], own},
     own = First[Commonest[Cases[sh[[All, 1]], INT[lb_String, ___] /; ! StringContainsQ[lb, "x"] && ! MemberQ[plHeads, lb] :> lb]]];
     {"shift map",
      ((canonB /. INT[fsym, b__] :> INT[own, b]) //. Dispatch[sh /. INT[lb_String, b__] /; MemberQ[plHeads, lb] :> INT[ToExpression[lb], b]])
        (* remaining non-PL labels: the family itself and its crossed copies, with the Kira alias (DD for D) undone *)
        /. INT[lb_String, b__] :> INT[ToExpression[StringReplace[lb, StartOfString ~~ own -> f]], b]}],
    {"SLC/Solutions shiftedBasis", Import[slcShiftedBasisFile[f]]}];
  shifted = shifted /. plRules /. {m2 -> 1, s12 -> x, s13 -> y, s23 -> 1 - x - y};
  shifted = Collect[#, {blPat, _INT}, Together] & /@ shifted;
  planar = Count[shifted, e_ /; FreeQ[e, _INT], {1}];
  Export[FileNameJoin[{outDir, f <> "_shiftedBasis.m"}], shifted];
  Print[StringPadRight[f, 3], ": ", Length[shifted], " entries, ", planar, " planar (only BASISLC), source ", src,
   "; other integrals left: ", Union@Cases[shifted, INT[h_, ___] :> h, Infinity]],
  {f, fams}];
