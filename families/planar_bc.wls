#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Boundary conditions of the fully planar families (A, E1, E2, M3, R2) to ep^wmax, as in
   reorganize.nb: basis /. shift map /. PL->BASISLC (all crossings), then for crossing X
   (x,y) -> cross[X] and BASISLC labels composed with X, PL crossings inserted, d -> 4-2eps,
   BC = solution /. _G -> 0.
   Only the G-free part of the PL crossings is needed: BC_w = Sum_j R_j . C_{w-j}, with
   R = sum_j eps^j R_j the shift/reduction coefficients and C_m the PL constants at weight m
   (PL/pl_gpls_7_8, run PL/run_gpl_xy_fast.wls 8 cross=<c> first). R depends on x,y; the result
   must not: it is evaluated at two points and the two must agree exactly.

   Usage: wolframscript -file planar_bc.wls [wmax=8] [fam=A,E1,...]
   Writes families/results/<fam>_bc.m ({BC[fam, X, i] -> Sum_w ep^w BC_w}) and compares
   ep^0..6 with all.ints/boundaries_all/<fam>_bc.m. *)

$HistoryLength = 0;
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "config.wl"}]];
Get[FileNameJoin[{here, "..", "PL", "load_gpl_xy.wl"}]];
args = Association[Cases[$ScriptCommandLine[[2 ;;]], s_String /; StringContainsQ[s, "="] :> Rule @@ StringSplit[s, "=", 2]]];
wmax = ToExpression[Lookup[args, "wmax", "8"]];
fams = If[KeyExistsQ[args, "fam"], StringSplit[args["fam"], ","], planarFamilies];
outDir = FileNameJoin[{here, "results"}]; Quiet[CreateDirectory[outDir]];
ep = Global`ep; eps = Global`eps;

crossLabels = {"", "x12", "x13", "x23", "x123", "x132"};
cross = <|"" -> {y -> y, x -> x}, "x12" -> {y -> 1 - x - y, x -> x}, "x13" -> {y -> y, x -> 1 - x - y},
   "x23" -> {y -> x, x -> y}, "x123" -> {y -> x, x -> 1 - x - y}, "x132" -> {y -> 1 - x - y, x -> y}|>;
crossmatch = {{y, x} -> "", {1 - x - y, x} -> "x12", {y, 1 - x - y} -> "x13", {x, y} -> "x23", {x, 1 - x - y} -> "x123", {1 - x - y, y} -> "x132"};
bl[c_] := ToExpression["BASISLC" <> c];
plh[c_] := ToExpression["PL" <> c];
(* crossing X applied to a BASISLC^c label (reorganize.nb crossBASISLC) *)
compose[c_, X_] := (cross[c][[All, 2]] /. cross[X]) /. crossmatch;
lab[c_] := If[c === "", "x", c];   (* label used by the PL sparse files *)

(* PL -> BASISLC in every crossing *)
plTab = Flatten[Get /@ plTableFiles] /. {s23 -> 1 - y - x, s13 -> y};
plRules = Dispatch[Flatten[Table[plTab /. PL[a__] :> INT[plh[c], a] /. BASISLC -> bl[c] /. cross[c], {c, crossLabels}]]];
blVars = Flatten[Table[bl[c][k], {c, crossLabels}, {k, 291}]];
blIndex = AssociationThread[blVars -> Range[Length[blVars]]];

(* PL constants: C[m] = list over (crossing, k) of the G-free part at weight m *)
w0PL = Get[FileNameJoin[{here, "..", "PL", "work", "weight0_solution.m"}]];
normMZV = s_Symbol /; SymbolName[s] === "MZV" :> Global`MZV;
plConst[m_] := plConst[m] = If[m == 0, Flatten[Table[w0PL, {6}]],
   Flatten[Table[With[{F = PLGPLSparse[m, lab[c]]},
      With[{sel = Flatten@Position[F[[1]], {{}, {}, _}, {1}]},
       If[sel === {}, ConstantArray[0, 291], F[[2]][[All, sel]] . (F[[1]][[sel, 3]] /. normMZV)]]], {c, crossLabels}]]];
pts = {{x -> 3/13, y -> 5/18}, {x -> 7/23, y -> 2/11}};

Do[
  tf = AbsoluteTime[];
  fsym = ToExpression[f];
  canon = Import[canonFile[f]];
  n = Length[canon];
  sh = Import[shiftMapFile[f]] /. INT[l_String, __, {b__}] :> INT[l, b];
  own = First[Commonest[Cases[sh[[All, 1]], INT[l_String, ___] /; ! StringContainsQ[l, "x"] && ! MemberQ[plHeads, l] :> l]]];
  shRules = Dispatch[sh /. INT[l_String, b__] /; MemberQ[plHeads, l] :> INT[ToExpression[l], b]];
  shifted = (canon /. INT[fsym, b__] :> INT[own, b]) //. shRules;
  shifted = shifted /. plRules /. {m2 -> 1, s12 -> x, s13 -> y, s23 -> 1 - x - y};
  left = Union@Cases[shifted, _INT, Infinity];
  If[left =!= {}, Print[f, ": ", Length[left], " integrals not reduced to BASISLC, e.g. ", left[[;; UpTo[3]]]]; Continue[]];
  ca = CoefficientArrays[shifted, blVars];
  If[Length[ca] =!= 2 || Max[Abs[Normal[ca[[1]]] /. {x -> 1/3, y -> 1/5, d -> 7/3}]] =!= 0, Print[f, ": not linear in BASISLC"]; Continue[]];
  R = ca[[2]];
  bc = {};
  Do[
    (* columns of R for crossing X: BASISLC^c[k] -> BASISLC^(c o X)[k] *)
    perm = Flatten[Table[Range[291] + 291 (Position[crossLabels, compose[c, X]][[1, 1]] - 1), {c, crossLabels}]];
    vals = Table[
      RX = SparseArray[Normal[R] /. cross[X] /. pt /. d -> 4 - 2 eps];
      ser = Normal[Series[Normal[RX], {eps, 0, wmax}]];
      lo = Min[0, Min[Exponent[#, eps, Min] & /@ Flatten[ser]]];
      If[lo < 0, Print[f, " ", X, ": coefficients have poles eps^", lo, " -- PL needed beyond weight ", wmax]; Abort[]];
      Table[Expand[Sum[Coefficient[ser, eps, j] . plConst[w - j][[perm]], {j, 0, w}]], {w, 0, wmax}],
      {pt, pts}];
    nd = Count[Flatten[Expand[vals[[1]] - vals[[2]]]], Except[0]];
    If[nd > 0, Print[f, " ", X, ": result depends on x,y in ", nd, " entries"]; Abort[]];
    bc = Join[bc, Table[Global`BC[fsym, X, i] -> Sum[ep^w vals[[1, w + 1, i]], {w, 0, wmax}], {i, n}]],
    {X, crossLabels}];
  old = Association[Import[bcFile[f]]];
  bc = bc[[Ordering[Lookup[PositionIndex[Keys[old]], bc[[All, 1]]][[All, 1]]]]];   (* entry order of boundaries_all *)
  Export[FileNameJoin[{outDir, f <> "_bc.m"}], bc];
  diff = Select[bc, Expand[Normal[Series[#[[2]], {ep, 0, 6}]] - Lookup[old, #[[1]], Missing[]]] =!= 0 &];
  Print[f, ": ", n, " masters, ", Length[bc], " BC to ep^", wmax, " in ", Round[AbsoluteTime[] - tf, 0.1], " s; vs boundaries_all through ep^6: ",
   If[diff === {}, "identical", ToString[Length[diff]] <> " differ, e.g. " <> ToString[diff[[1, 1]], InputForm]]],
  {f, fams}];
