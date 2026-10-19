#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Fast weight-by-weight construction of the PL solution in GPLs of (x,y),
   G[wx; x] G[wy; y] with x-letters {0,1,-y,1-y} and y-letters {0,1}.

   Same mathematics as the solution.nb loop
     tofix = GIntegrate[D[atilde,x].sol[w-1], x]
     toadd = GIntegrate[D[atilde,y].sol[w-1] - DG[tofix,y], y]
     sol[w] = tofix + toadd + cnst[w]
   but with both integrations done as pure letter-prepending:
   - D[atilde,x] = Sum_l Ax[l]/(x-l), l in {0,1,-y,1-y}, so the x-integral of
     G[w;x] is G[l,w;x] (and of an x-free c, c G[l;x]).
   - D[atilde,y].sol - D[tofix,y] is x-independent, so it equals its regularized
     x->0 value. There tofix -> 0 (every G[..;x] has vanishing regularized limit,
     and D[atilde,y] is analytic at x=0), and D[atilde,y] -> Ay0[0]/y + Ay0[1]/(y-1).
     Hence toadd = Sum_m Ay0[m].prepend_m( sol[w-1] /. G[..;x] -> 0 ).
   Everything is then sparse rational linear algebra on the coefficient matrix
   S (291 x #keys), key = {x-word, y-word, transcendental monomial}.

   Usage: wolframscript -file run_gpl_xy_fast.wls <wmax> [validate] [m] [cross=x12|x13|x23|x123|x132]
   Writes pl_gpls_7_8/w<w>_sparse.mx = {keys, S} for w = 1..wmax (seconds);
   load with load_gpl_xy.wl. With "m", also writes the G-expression text file
   pl_gpls_7_8/w<w>_fast.m for w >= 7 (slow: ~15 min at w=8). *)

SetDirectory[DirectoryName[$InputFileName]];
wmax = ToExpression[$ScriptCommandLine[[2]]];
(* crossing: I^c(s12,s13,s23) = I(T[c]) with (s12,s13,s23) = (x,y,1-x-y); the crossed
   solution's constants are the PL constants at the corner where T[c] sits for x=y=0 *)
cross = First[StringCases[$ScriptCommandLine[[3 ;;]], "cross=" ~~ c__ :> c] // Flatten, "x"];
crossT = <|"x" -> {x, y, zz}, "x12" -> {x, zz, y}, "x13" -> {zz, y, x}, "x23" -> {y, x, zz},
   "x123" -> {zz, x, y}, "x132" -> {y, zz, x}|>[cross] /. zz -> 1 - x - y;
corner = <|"x" -> 6, "x23" -> 6, "x12" -> 4, "x132" -> 4, "x13" -> 2, "x123" -> 2|>[cross];
outDir = "./pl_gpls_7_8/" <> If[cross === "x", "", cross <> "/"];
Quiet[CreateDirectory[outDir]];
Print["crossing ", cross, ": (s12,s13,s23) -> ", crossT, ", constants from corner ", corner];
validate = MemberQ[$ScriptCommandLine[[3 ;;]], "validate"];
writeM = MemberQ[$ScriptCommandLine[[3 ;;]], "m"];
SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}];

(* ---------- connection matrices ---------- *)
atilde = Import["./input/atilde.m"] /. {y -> crossT[[2]], z -> crossT[[3]]} /. ep -> 1;
lg = {Log[x], Log[1 - x], Log[x + y], Log[1 - x - y], Log[y], Log[1 - y]};
A = AssociationMap[SparseArray[Coefficient[atilde, #]] &, lg];
rest = Expand[atilde - Sum[A[l] l, {l, lg}]];
If[! FreeQ[rest, x | y], Print["atilde is not of the expected dlog form"]; Abort[]];
If[! VectorQ[Flatten[Normal /@ Values[A]], NumericQ], Print["non-numeric connection entries"]; Abort[]];
Ax = <|0 -> A[Log[x]], 1 -> A[Log[1 - x]], -y -> A[Log[x + y]], 1 - y -> A[Log[1 - x - y]]|>;
Ay0 = <|0 -> A[Log[y]] + A[Log[x + y]], 1 -> A[Log[1 - y]] + A[Log[1 - x - y]]|>;

(* ---------- expression <-> coefficient matrix ---------- *)
constPat = Pi | _Zeta | _MZV;
termKey[t_] := Module[{f, xw, yw, r, rat},
  f = If[Head[t] === Times, List @@ t, {t}];
  xw = Cases[f, G[a__, x] :> {a}]; yw = Cases[f, G[a__, y] :> {a}];
  If[Length[xw] > 1 || Length[yw] > 1 || ! FreeQ[DeleteCases[f, G[__, x | y]], x | y | _G],
    Print["unexpected term: ", t]; Abort[]];
  r = Times @@ DeleteCases[f, G[__, x | y]];
  rat = r /. constPat -> 1;
  If[! NumericQ[rat], Print["bad coefficient: ", r]; Abort[]];
  {{First[xw, {}], First[yw, {}], r/rat}, rat}];

toMatrix[exprs_List] := Module[{rows, keys, idx},
  rows = Table[With[{e = Expand[exprs[[i]]]},
     If[e === 0, {}, {i, ##} & @@@ (termKey /@ If[Head[e] === Plus, List @@ e, {e}])]], {i, Length[exprs]}];
  rows = Flatten[rows, 1];
  keys = Union[rows[[All, 2]]];
  idx = AssociationThread[keys -> Range[Length[keys]]];
  {keys, SparseArray[{#1, idx[#2]} -> #3 & @@@ rows, {Length[exprs], Max[1, Length[keys]]}]}];

(* merge column blocks with possibly repeated keys *)
mergeBlocks[blocks0_List] := Module[{blocks = Select[blocks0, #[[1]] =!= {} &], keys, mats, uk, idx, P},
  keys = Join @@ blocks[[All, 1]]; mats = Join[Sequence @@ blocks[[All, 2]], 2];
  uk = Union[keys]; idx = AssociationThread[uk -> Range[Length[uk]]];
  P = SparseArray[Table[{k, idx[keys[[k]]]} -> 1, {k, Length[keys]}], {Length[keys], Length[uk]}];
  {uk, mats.P}];

(* row by row: a SparseArray . (symbolic vector) is ~300x slower *)
toExprs[{keys_, S_}] := Module[{gx, gy, vec},
  gx[{}] = 1; gx[w_] := G[Sequence @@ w, x]; gy[{}] = 1; gy[w_] := G[Sequence @@ w, y];
  vec = gx[#[[1]]] gy[#[[2]]] #[[3]] & /@ keys;
  Table[With[{r = S[[i]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]], {i, Length[S]}]];

(* ---------- one weight step ---------- *)
step[{keys_, S_}, cnst_] := Module[{xb, yb, sel},
  xb = Table[{{Prepend[#[[1]], l], #[[2]], #[[3]]} & /@ keys, Ax[l].S}, {l, Keys[Ax]}];
  sel = Flatten@Position[keys, {{}, _, _}, {1}];
  yb = If[sel === {}, {},
    Table[{{{}, Prepend[#[[2]], m], #[[3]]} & /@ keys[[sel]], Ay0[m].S[[All, sel]]}, {m, Keys[Ay0]}]];
  mergeBlocks[Join[xb, yb, {toMatrix[cnst]}]]];

(* drop numerically-zero columns *)
prune[{keys_, S_}] := Module[{nz = Flatten@Position[Total[Abs[S]], _?(# != 0 &), {1}]}, {keys[[nz]], S[[All, nz]]}];

(* ---------- run ---------- *)
cur = toMatrix[Import["./work/weight0_solution.m"]];
Print["w0: ", Length[cur[[1]]], " keys"];
If[validate,
  plm = Import["../planar_stuff/pl" <> If[cross === "x", "", cross] <> ".m"][[All, 2]];
  (* number of (entry, key) coefficients that differ *)
  nDiff[a_, b_] := Module[{d = prune[mergeBlocks[{a, {b[[1]], -b[[2]]}}]]},
    If[d[[1]] =!= {}, Print["   sample differing keys: ", d[[1, ;; UpTo[4]]]]];
    Length[d[[1]]]];
];
Do[
  cn = Import["./work/constant/const" <> ToString[w] <> "_" <> ToString[corner] <> ".m"];
  {t, cur} = AbsoluteTiming[prune[step[cur, cn]]];
  Print["w", w, ": ", Length[cur[[1]]], " keys, ", Length[cur[[2]]["NonzeroValues"]], " nonzeros, ", Round[t, 0.1], " s"];
  If[validate && w <= 6,
    Print["   differing coefficients vs planar_stuff/pl" <> If[cross === "x", "", cross] <> ".m at eps^", w, ": ", nDiff[cur, toMatrix[Coefficient[plm, eps, w]]]]];
  If[validate && w == 7 && cross === "x",
    Print["   differing coefficients vs pl_gpls_7_8/w7.m: ", nDiff[cur, toMatrix[Import["./pl_gpls_7_8/w7.m"]]]]];
  {t, ex} = AbsoluteTiming[Export[outDir <> "w" <> ToString[w] <> "_sparse.mx", cur, "MX"]];
  Print["   exported ", outDir, "w", w, "_sparse.mx (", Round[t, 0.1], " s)"];
  If[writeM && w >= 7,
    {t, ex} = AbsoluteTiming[toExprs[cur]];
    Export[outDir <> "w" <> ToString[w] <> "_fast.m", ex];
    Print["   exported ", outDir, "w", w, "_fast.m (", Round[t, 0.1], " s)"]],
  {w, 1, wmax}];
