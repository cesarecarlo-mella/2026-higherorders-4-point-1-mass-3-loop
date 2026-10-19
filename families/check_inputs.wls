#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Inventory of the inputs of every family + coverage of the planar (PL) reductions.
   Prints one line per family and writes, for every family that needs new Kira reductions,
   families/extra_shifts/<fam>_missingPL.m = list of PL[...] integrals not covered by any
   PL -> BASISLC table (the crossing label is dropped: the reduction of PL[...] serves all
   crossings).
   Usage: wolframscript -file check_inputs.wls *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "config.wl"}]];
outDir = FileNameJoin[{DirectoryName[$InputFileName], "extra_shifts"}];
Quiet[CreateDirectory[outDir]];

known = <||>;
Do[Scan[(known[#] = True) &, Cases[Flatten[{Quiet@Get[f]}], (PL[a__] -> _) :> {a}]], {f, plTableFiles}];
Print["PL -> BASISLC tables: ", Length[plTableFiles], " files, ", Length[known], " PL integrals covered"];
dict = Import[lettersDictFile];
pt = {x -> 13/97, y -> 29/101};

Do[
  info = {};
  a = Import[atildeFile[f]]; n = Length[a];
  AppendTo[info, "n=" <> ToString[n] <> " letters " <> ToString[Union@Cases[a, l[i_] :> i, Infinity]]];
  canon = Import[canonFile[f]];
  If[Length[canon] =!= n, AppendTo[info, "CANON LENGTH " <> ToString[Length[canon]]]];
  bc = If[FileExistsQ[bcFile[f]], Import[bcFile[f]], {}];
  AppendTo[info, If[bc === {}, "NO BC", "bc " <> ToString[Length[bc]] <> " (6 crossings x " <> ToString[Length[bc]/6] <> ") to ep^" <>
      ToString[Max[Exponent[#[[2]], Global`ep] & /@ bc]] <> If[FreeQ[bc[[All, 2]], Global`cnst], "", " UNRESOLVED"]]];
  If[KeyExistsQ[rootPos, f], AppendTo[info, "sqrt-normalized " <> ToString[rootPos[f]]]];
  (* old SLC Atilde vs all.ints: numerically, d/dx at a point *)
  oldA = Select[{$SLC <> legacyName[f] <> "/input/Atilde.m", $SLC <> legacyName[f] <> "/input/" <> legacyName[f] <> "_Atilde.m"}, FileExistsQ];
  If[oldA =!= {},
    dn = N[D[a /. dict, x] /. pt, 30]; do = N[D[Import[First[oldA]] /. {s -> x, t -> y, m2 -> 1, eps -> 1, ep -> 1}, x] /. pt, 30];
    nd = Length[Position[Chop[dn - do, 10^-15], v_ /; v != 0, {2}]];
    AppendTo[info, "SLC Atilde " <> If[nd == 0, "identical", ToString[nd] <> " entries differ (+-I)"]]];
  (* planar reduction *)
  If[! MissingQ[slcShiftedBasisFile[f]],
    sb = Import[slcShiftedBasisFile[f]];
    left = Union@Cases[sb, INT[h_, b__] /; MemberQ[plHeads, ToString[h]] :> {b}, Infinity];
    miss1 = Select[left, ! KeyExistsQ[known, #] &];
    AppendTo[info, "SLC shiftedBasis: " <> ToString[Count[sb, e_ /; FreeQ[e, INT], {1}]] <> " planar entries, own label " <>
       ToString[Complement[Union@Cases[sb, INT[h_, ___] :> ToString[h], Infinity], plHeads]] <> ", unreduced PL " <> ToString[Length[left]] <>
       " (" <> ToString[Length[miss1]] <> " not in tables)"],
    miss1 = {}];
  If[! MissingQ[shiftMapFile[f]],
    (* exact labels: the family's own integrals are INT[f, ...] in the basis; the map may use an
       alias for the own label (e.g. "DD" for D), crossed self-copies ("Cx12", "DDx12", ...) and
       chained rules (crossed copies or crossed PL integrals shifted again) *)
    sh = Import[shiftMapFile[f]] /. INT[lab_String, __, {b__}] :> INT[lab, b];
    own = First[Commonest[Cases[sh[[All, 1]], INT[lab_String, ___] /; ! StringContainsQ[lab, "x"] && ! MemberQ[plHeads, lab] :> lab]]];
    shifted = (canon /. INT[Symbol[f], b__] :> INT[own, b]) //. Dispatch[sh];
    shifted = Expand /@ shifted;
    plInts = Union@Cases[shifted, INT[lab_String, b__] /; MemberQ[plHeads, lab] :> {b}, Infinity];
    miss2 = Select[plInts, ! KeyExistsQ[known, #] &];
    nonPL = Complement[Union@Cases[shifted, INT[lab_String, ___] :> lab, Infinity], plHeads];
    AppendTo[info, "shift map (own label \"" <> own <> "\"): " <>
       ToString[Count[shifted, e_ /; FreeQ[e, INT[lab_String, ___] /; ! MemberQ[plHeads, lab]], {1}]] <>
       " planar entries, non-PL integrals left " <> ToString[nonPL] <> ", " <>
       ToString[Length[plInts]] <> " PL integrals (" <> ToString[Length[miss2]] <> " not in tables), basis integrals not in the map " <>
       ToString[Length[Complement[Union@Cases[canon, INT[_, b__] :> {b}, Infinity], Union@Cases[sh[[All, 1]], INT[own, b__] :> {b}]]]]],
    miss2 = {}];
  If[MissingQ[slcShiftedBasisFile[f]] && MissingQ[shiftMapFile[f]], AppendTo[info, "NO PLANAR REDUCTION DATA"]];
  (* what still needs Kira: the route that will be used (SLC shiftedBasis if present, else shift map) *)
  need = If[! MissingQ[shiftMapFile[f]], miss2, miss1];
  If[need =!= {}, Export[FileNameJoin[{outDir, f <> "_missingPL.m"}], PL @@@ need]; AppendTo[info, "-> " <> ToString[Length[need]] <> " PL integrals to reduce"]];
  Print[StringPadRight[f, 3], " (", legacyName[f], "): ", StringRiffle[info, "; "]],
  {f, families}];
