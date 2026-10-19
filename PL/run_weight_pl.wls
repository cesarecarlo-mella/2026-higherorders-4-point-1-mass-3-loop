(* PL run_weight_pl.wls <w>
   Builds sol[w][seg] for seg in {1..6,x1,y1,xmy} by actually solving PL's
   boundary system at each weight (unlike H_opt, which reused H's
   already-solved frozen numbers -- PL has no such frozen solution for
   w>=1, only the single known banana boundary value).
   Usage: wolframscript -file run_weight_pl.wls <w> *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
Global`x; Global`y; Global`t; Global`delta; Global`ep; Global`eps;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];

SetDirectory[DirectoryName[$InputFileName]];
wgt = ToExpression[$ScriptCommandLine[[2]]];
Print["=== PL weight ", wgt, " ==="];

Print["=== loading precompute + boundary/regularity matrices ==="];
Get["work/precompute_pl.mx"];
allPx0 = Import["work/allP_x0.m"];
allPy0 = Import["work/allP_y0.m"];
allPxy1 = Import["work/allP_xy1.m"];
resX1 = Import["work/regX1.m"];
resY1 = Import["work/regY1.m"];
resXmY = Import["work/regXmY.m"];
stackX0 = Join @@ allPx0;
stackY0 = Join @@ allPy0;
stackXy1 = Join @@ allPxy1;

Log00 = {G[0, 0] -> 0, G[0, 0, 0] -> 0, G[0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0, 0] -> 0};
endPt1 = Table[G[Sequence @@ ConstantArray[1, k], t] -> 0, {k, 1, 8}];

n = 291;
Print["=== loading sol[", wgt - 1, "] ==="];
(* weight 0 is a single constant vector (same on every segment), not per-segment files *)
sol0 = Import["work/weight0_solution.m"];
Do[sol[0][j] = sol0, {j, 1, 6}];
sol[0][x1] = sol0; sol[0][y1] = sol0; sol[0][xmy] = sol0;
(* weights 1..wgt-1 are real per-segment files *)
Do[sol[k][j] = Import["work/sol" <> ToString[k] <> "_" <> ToString[j] <> ".m"], {j, 1, 6}, {k, 1, wgt - 1}];
Do[
  sol[k][x1] = Import["work/sol" <> ToString[k] <> "_x1.m"];
  sol[k][y1] = Import["work/sol" <> ToString[k] <> "_y1.m"];
  sol[k][xmy] = Import["work/sol" <> ToString[k] <> "_xmy.m"];
  , {k, 1, wgt - 1}];

Print["=== building IntPart via GIntegrate (correctness-first, no M_a optimization yet) ==="];
(* per-component progress: integrate each of the 291 vector entries individually (instead of
   a single opaque Map) so we can see WHICH component we're on and how long each one takes --
   previously a step like IntParty1 could run for 2+ hours with zero visibility into progress.
   Prints every `step` components and whenever a single component takes >5s (flagging the
   expensive ones specifically). *)
integrateWithProgress[label_, rawVec_, var_] := Module[{n, result, t0, tComp, step},
   n = Length[rawVec];
   step = 20;
   result = Table[Null, {n}];
   t0 = AbsoluteTime[];
   Do[
     {tComp, result[[i]]} = AbsoluteTiming[GIntegrate[rawVec[[i]], var]];
     If[Mod[i, step] == 0 || i == n || tComp > 5,
       Print["    ", label, " [", i, "/", n, "]  this-component=", tComp, "s  elapsed=", AbsoluteTime[] - t0, "s"]
       ];
     , {i, n}];
   result
   ];
(* checkpoint + resume: if work/intpart{wgt}_{seg}.m already exists on disk, load it directly
   instead of recomputing (GIntegrate is the single most expensive step in this pipeline, and
   runs have needed killing/restarting more than once -- losing a 20+ minute IntPartx1/y1 to a
   kill was the direct motivation). Otherwise compute it AND export it immediately (not batched
   at the end), so a future kill only loses the segment currently in progress, never ones
   already finished. *)
getOrCompute[label_, filename_, rawVec_, var_] := Module[{tC, result},
   If[FileExistsQ[filename],
     Print["  ", label, ": found existing checkpoint ", filename, ", loading (skipping recompute)"];
     {0., Import[filename]},
     {tC, result} = AbsoluteTiming[integrateWithProgress[label, rawVec, var]];
     Export[filename, result];
     Print["  ", label, ": computed and exported to ", filename];
     {tC, result}
     ]
   ];
ipFile[seg_] := "work/intpart" <> ToString[wgt] <> "_" <> ToString[seg] <> ".m";
{tI1, IntPart1} = getOrCompute["IntPart1", ipFile[1], Apart[dt[1].sol[wgt - 1][1], t], t]; Print["  [1/9] IntPart1 done, t=", tI1, "s"];
{tI2, IntPart2} = getOrCompute["IntPart2", ipFile[2], Apart[dt[2].sol[wgt - 1][2], t], t]; Print["  [2/9] IntPart2 done, t=", tI2, "s"];
{tI3, IntPart3} = getOrCompute["IntPart3", ipFile[3], Apart[dt[3].sol[wgt - 1][3], t], t]; Print["  [3/9] IntPart3 done, t=", tI3, "s"];
{tI4, IntPart4} = getOrCompute["IntPart4", ipFile[4], Apart[dt[4].sol[wgt - 1][4], t], t]; Print["  [4/9] IntPart4 done, t=", tI4, "s"];
{tI5, IntPart5} = getOrCompute["IntPart5", ipFile[5], Apart[dt[5].sol[wgt - 1][5], t], t]; Print["  [5/9] IntPart5 done, t=", tI5, "s"];
{tI6, IntPart6} = getOrCompute["IntPart6", ipFile[6], Apart[dt[6].sol[wgt - 1][6], t], t]; Print["  [6/9] IntPart6 done, t=", tI6, "s"];
{tIx1, IntPartx1} = getOrCompute["IntPartx1", ipFile["x1"], Apart[dyx1.sol[wgt - 1][x1], y], y]; Print["  [7/9] IntPartx1 done, t=", tIx1, "s"];
{tIy1, IntParty1} = getOrCompute["IntParty1", ipFile["y1"], Apart[dxy1.sol[wgt - 1][y1], x], x]; Print["  [8/9] IntParty1 done, t=", tIy1, "s"];
{tIxmy, IntPartxmy} = getOrCompute["IntPartxmy", ipFile["xmy"], Apart[dxmy.sol[wgt - 1][xmy], y], y]; Print["  [9/9] IntPartxmy done, t=", tIxmy, "s"];
Print["IntPart timings: ", {tI1, tI2, tI3, tI4, tI5, tI6, tIx1, tIy1, tIxmy}];

CNST5 = Table[cS5[i], {i, 1, n}];

Print["=== loading fibration-change table (needed for t->1 evaluation beyond weight 1) ==="];
tto1Red = Import["../H_opt/input/fibrtto1.m"];
Print["tto1Red: ", Length[tto1Red], " rules"];

(* "start" evaluation (var->0): direct, no fibration change needed (matches H's weight-1 pattern). *)
evalStart[expr_, var_] := expr /. G[__, var] -> 0 /. var -> 0;
(* "end" evaluation (var->1), weight<=2 (H's OWN recipe, confirmed from reference lines 973-978
   and 807-812): only kill the literal ALL-ONES word of the current weight (mirrors Logt0's
   all-zero pattern at t->0), then substitute t->1 directly -- PolyLogTools/HPL auto-evaluates
   every other {0,1}-letter word to its MZV value via the preloaded MZV table (startup banner:
   "Table of MZVs loaded up to weight 6"), no ReduceCG/fibration needed. H itself only switches
   to the heavier fibration (t->1-t/.tto1Red//ShuffleG) starting at weight 3 (reference line
   1139) -- using it already at weight <=2 was the likely bug: it does NOT just regularize, it
   can give a numerically different (wrong) answer even though "fibration" sounds more general. *)
evalEndSimple[expr_] := expr /. endPt1 /. t -> 1;
(* "end" evaluation (var->1), weight>=3: full fibration-change, needed once literal "1-t"-type
   mixed arguments can appear (per H reference, starting weight 3). *)
evalEndFibration[expr_] := ((expr /. t -> 1 - t /. tto1Red // ShuffleG) /. t -> 0 /. G[__, t] -> 0 /. HPL[a__] -> 0);
evalAt0 = evalStart[#, t] &;
(* H's own notebook EXTENDS tto1Red at every weight>=3 (reference lines 1133, 1159, 1338,
   1525, 1757): tto1Red=Join[tto1Red, ReduceCG[GetGs[expr/.t->1-t],{t},Point->{t->3/14}]]//Union
   -- before reducing each expression that needs it. The frozen 126-rule table alone is only
   guaranteed complete for H's own weight<=2-worth of words; PL's own alphabet needs its own
   extension once genuinely new weight>=3 words appear. Mutates tto1Red as a side effect. *)
extendTto1Red[expr_] := Module[{new, oldLhs},
   new = ReduceCG[GetGs[expr /. t -> 1 - t], {t}, Point -> {t -> 3/14}];
   oldLhs = tto1Red[[All, 1]];
   new = Select[new, ! MemberQ[oldLhs, #[[1]]] &];
   tto1Red = Join[tto1Red, new];
   ];
(* evalEndFibration's final step (G[__,t]->0) only kills words still ending in literal t;
   degenerate all-zero constants like G[0,0,0] that lost their t-dependence earlier in the
   t->1-t/.tto1Red//ShuffleG reduction sail through un-killed and need an explicit /.Log00
   pass too (this was the actual weight>=3 bug: these leaked into EVERY d2..d6, not just d01). *)
(* weight>=7: ReduceCG's internal basis-generator (GetGs) is hard-capped at weight<=6, so the
   symbolic fibration-change approach above stops working entirely. First attempt (PSLQ-fitting
   numeric values against a small known MZV basis) was WRONG: direct numeric evaluation at t=1
   genuinely diverges for any word with a LEADING RUN OF 1's (not just the fully-all-ones word --
   e.g. G[1,1,0,t], G[1,0,1,t] both fail), confirmed empirically and by Cesare's counterexample.
   The actual fix, found with Cesare: plain substitution G[...,t]/.t->1, WITH
   $HPLAutoConvertToKnownFunctions=True already in effect, automatically separates the finite
   closed-form part (in terms of the correct weight-7/8 MZV basis, e.g. Zeta[7], or even
   MZV[{5,3}] for the one genuinely-irreducible weight-8 piece) from an EXPLICIT symbolic
   divergence marker HPL[{1},1] (representing the log(1-t) singularity), to ANY power. This is
   literally the shuffle-regularization PolyLogTools already does internally -- no PSLQ, no
   ReduceCG, no recursion needed. Dropping HPL[{1},1]->0 afterward is the SAME regularization
   convention this whole pipeline already uses (identical in spirit to /.HPL[a__]->0 at the end
   of evalEndFibration above, and to Log00/endPt1) -- it is NOT assuming/requiring the divergent
   pieces to algebraically cancel across the full combination; dropping them IS the prescription.
   Precomputed once as a lookup table over all 2^7/2^8 {0,1}-letter words (PL's segment alphabet
   is pure {0,1} after parametrization, same as weight<=6) since re-substituting inline for every
   call is unnecessary work. *)
(* CRITICAL FIX: GIntegrate, with $HPLAutoConvertToKnownFunctions=True in effect, does NOT
   always produce pure weight-7/8 (7-/8-letter) G-words -- it frequently factors a weight-7/8
   expression into a LOWER-weight G-word (e.g. G[1,t], G[0,0,1,t], G[0,0,0,0,1,t]) multiplied by
   an explicit Pi/Zeta coefficient carrying the rest of the weight. A table built only from pure
   7-letter (or 8-letter) tuples NEVER matches these shorter sub-words, so they sailed through
   evalAt1 completely unsubstituted -- confirmed by direct inspection of d2's own entries, which
   were full of literal G[1,t], G[0,0,1,t], etc. after the "fix". This was the actual bug behind
   every "Solve: 0 solutions" result at weight 7, not the fibration method, not Expand, not
   Log00. Fix: use the UNIFIED table covering every {0,1}-word from weight 0 through 8 (511
   rules total, already built and independently validated via two methods -- direct
   substitution and from-scratch shuffle algebra, both agreeing exactly) so ANY sub-word of
   ANY length that appears gets correctly matched and replaced. Expand the RHS once at load
   time for the same reason as before (un-expanded nested arithmetic can block Solve's internal
   cancellation). *)
wordTableAll = Import["work/word_table_0to8.m"] /. (lhs_ -> rhs_) :> (lhs -> Expand[rhs]);
wordTable7 = wordTableAll;
wordTable8 = wordTableAll;
evalEndTable[expr_, table_] := expr /. table;
evalAt1 = Which[
   wgt <= 2, evalEndSimple[#],
   wgt <= 6, (extendTto1Red[#]; evalEndFibration[#] /. Log00),
   wgt == 7, evalEndTable[#, wordTable7],
   True, evalEndTable[#, wordTable8]
   ] &;

Print["=== matching: build delta vectors d2..d6, d00,d01,d10 (all KNOWN, previous-weight data only) ==="];
{td2, d2} = AbsoluteTiming[
   tempSol1v = evalAt1[IntPart1]; tempSol2v0 = evalAt0[IntPart2]; tempSol1v - tempSol2v0];
Print["  d2 done, t=", td2, "s"];
{td3, d3} = AbsoluteTiming[
   tempSol2v = evalAt1[IntPart2]; tempSol3v0 = evalAt0[IntPart3]; tempSol2v - tempSol3v0];
Print["  d3 done, t=", td3, "s"];
{td4, d4} = AbsoluteTiming[
   tempSol3v = evalAt1[IntPart3]; tempSol4v0 = evalAt0[IntPart4]; tempSol3v - tempSol4v0];
Print["  d4 done, t=", td4, "s"];
{td5, d5} = AbsoluteTiming[
   tempSol4v = evalAt1[IntPart4]; tempSol5v0 = evalAt0[IntPart5]; tempSol4v - tempSol5v0];
Print["  d5 done, t=", td5, "s"];
{td6, d6} = AbsoluteTiming[
   tempSol5v = evalAt1[IntPart5]; tempSol6v0 = evalAt0[IntPart6]; tempSol5v - tempSol6v0];
Print["  d6 done, t=", td6, "s"];

(* impart[j] = Sum_{k=1}^{wgt} (-I Pi)^k/k! * todelta[j]^k . sol[wgt-k][j]
   -- weight 1 only needs k=1; weight >=2 needs the higher-k terms too. *)
impartSum[j_] := Sum[
   ((-I Pi)^k/k!) MatrixPower[todelta[j], k] . sol[wgt - k][j],
   {k, 1, wgt}];
impart5sol = impartSum[5];
(* H's own reference: match00 uses t->0 ONLY at weight 1 (reference lines 826-829); starting
   at weight 2 it switches to t->1 (lines 990-993, reconfirmed at weight 3, lines 1158-1162).
   Using evalAt0 unconditionally here (as this script did before) is wrong for weight>=2. *)
{td00, d00} = AbsoluteTiming[If[wgt == 1,
   (evalAt0[impart5sol + IntPart5] /. HPL[a__] -> 0) - evalStart[IntPartxmy, y],
   (evalAt1[impart5sol + IntPart5] /. Log00 /. HPL[a__] -> 0) - (evalStart[IntPartxmy, y] /. Log00)
   ]];
Print["  d00 done, t=", td00, "s"];
impart3sol = impartSum[3];
(* H's literal reference code for match01/match10 (lines 830-837 at weight 1, 994-1000 at
   weight 2) never includes the all-ones-word kill that the segment chain (match2..match6)
   needs -- only t->1/.Log00/.HPL->0 (or at weight2: t->1/.G[0,0]->0/.HPL->0). Replicating
   literally rather than going through the generic evalAt1 wrapper (which bakes in the
   all-ones kill via endPt1, which is NOT what H does here). *)
{td01, d01} = AbsoluteTiming[
   tmp01 = If[wgt <= 2,
      (impart3sol + IntPart3) /. t -> 1 /. Log00 /. HPL[a__] -> 0,
      evalAt1[impart3sol + IntPart3]
      ];
   tmp01 - evalStart[IntParty1, x]];
Print["  d01 done, t=", td01, "s"];
{td10, d10} = AbsoluteTiming[
   tmp10 = If[wgt <= 2,
      (impart3sol + IntPart3) /. t -> 0 /. Log00 /. HPL[a__] -> 0,
      evalStart[impart3sol + IntPart3, t] /. HPL[a__] -> 0
      ];
   tmp10 - evalStart[IntPartx1, y]];
Print["  d10 done, t=", td10, "s"];
Print["delta vectors built"];

Print["=== express CNST1..6,x1,y1,xmy symbolically in terms of unknown CNST5 ==="];
CNST4sym = CNST5 - d5;
CNST3sym = CNST4sym - d4;
CNST2sym = CNST3sym - d3;
CNST1sym = CNST2sym - d2;
CNST6sym = CNST5 + d6;
CNSTxmysym = CNST5 + d00;
CNSTy1sym = CNST3sym + d01;
CNSTx1sym = CNST3sym + d10;

tempSol5 = IntPart5 + CNST5;
tempSol1 = IntPart1 + CNST1sym;
tempSol3 = IntPart3 + CNST3sym;
tempSolx1 = IntPartx1 + CNSTx1sym;
tempSoly1 = IntParty1 + CNSTy1sym;
tempSolxmy = IntPartxmy + CNSTxmysym;

Print["=== boundary + regularity equations ==="];
eq1raw = stackX0.tempSol5;
eq2raw = stackY0.tempSol1;
eq7raw = stackXy1.tempSol3;
eq4raw = resX1.tempSolx1;
eq5raw = resY1.tempSoly1;
eq6raw = resXmY.tempSolxmy;

knownBanana = -Coefficient[
   1 - Pi^2/4 eps^2 - 29 Zeta[3] eps^3 - 71 Pi^4/160 eps^4 + (29 Pi^2 Zeta[3]/4 - 1263 Zeta[5]/5) eps^5
     + (841 Zeta[3]^2/2 - 11539 Pi^6/24192) eps^6
     + (2059 Pi^4 Zeta[3]/160 + 1263 Pi^2 Zeta[5]/20 - 18567 Zeta[7]/7) eps^7
     + (-841 Pi^2 Zeta[3]^2/8 + 36627 Zeta[3] Zeta[5]/5 - 737687 Pi^8/1382400) eps^8,
   eps, wgt];
Print["known banana boundary value at weight ", wgt, " (position 1): ", knownBanana];

allEqs = Join[
   Thread[eq1raw == 0], Thread[eq2raw == 0], Thread[eq7raw == 0],
   Thread[eq4raw == 0], Thread[eq5raw == 0], Thread[eq6raw == 0],
   {CNST5[[1]] == knownBanana}
   ] // Union;
Print["total equations: ", Length[allEqs]];

{tSolve, solveResult} = AbsoluteTiming[Solve[allEqs, CNST5]];
Print["Solve time: ", tSolve, "s, num solutions: ", Length[solveResult]];
s1 = solveResult[[1]];
CNST5val = CNST5 /. s1;
If[! FreeQ[CNST5val, cS5], Print["WARNING: unresolved cS5 symbols remain!"]; Print[Select[CNST5val, ! FreeQ[#, cS5] &]]];

Print["=== back-substituting to all segments ==="];
CNST4val = CNST4sym /. s1;
CNST3val = CNST3sym /. s1;
CNST2val = CNST2sym /. s1;
CNST1val = CNST1sym /. s1;
CNST6val = CNST6sym /. s1;
CNSTxmyval = CNSTxmysym /. s1;
CNSTy1val = CNSTy1sym /. s1;
CNSTx1val = CNSTx1sym /. s1;

(* CRITICAL FIX: exported solutions were NEVER simplified before saving, letting redundant
   shuffle-equivalent G-word representations compound weight over weight. Found via direct
   measurement: sol7_y1.m had ballooned to 2.3M leaves (975 distinct words, one single entry
   at 54806 leaves) by weight 7, versus 327 words at weight 6 -- this is almost certainly why
   IntParty1 at weight 8 ran for 2+ hours without even finishing Apart/GIntegrate on the first
   component: the matrix-vector product + Apart had to process an already-bloated input.
   ShuffleG cuts this vector's LeafCount in half (2.3M->1.07M) in 7.4s and is verified to
   preserve the value exactly -- applying it here keeps every future weight's input compact. *)
Print["=== simplifying solutions via ShuffleG before export (prevents weight-over-weight bloat) ==="];
{tShuf1, solOut1} = AbsoluteTiming[ShuffleG /@ (IntPart1 + CNST1val)];
{tShuf2, solOut2} = AbsoluteTiming[ShuffleG /@ (IntPart2 + CNST2val)];
{tShuf3, solOut3} = AbsoluteTiming[ShuffleG /@ (IntPart3 + CNST3val)];
{tShuf4, solOut4} = AbsoluteTiming[ShuffleG /@ (IntPart4 + CNST4val)];
{tShuf5, solOut5} = AbsoluteTiming[ShuffleG /@ (IntPart5 + CNST5val)];
{tShuf6, solOut6} = AbsoluteTiming[ShuffleG /@ (IntPart6 + CNST6val)];
{tShufXmy, solOutXmy} = AbsoluteTiming[ShuffleG /@ (IntPartxmy + CNSTxmyval)];
{tShufX1, solOutX1} = AbsoluteTiming[ShuffleG /@ (IntPartx1 + CNSTx1val)];
{tShufY1, solOutY1} = AbsoluteTiming[ShuffleG /@ (IntParty1 + CNSTy1val)];
Print["ShuffleG timings: ", {tShuf1, tShuf2, tShuf3, tShuf4, tShuf5, tShuf6, tShufXmy, tShufX1, tShufY1}];
solOut[1] = solOut1; solOut[2] = solOut2; solOut[3] = solOut3; solOut[4] = solOut4;
solOut[5] = solOut5; solOut[6] = solOut6; solOut[xmy] = solOutXmy; solOut[x1] = solOutX1; solOut[y1] = solOutY1;

Do[Export["work/sol" <> ToString[wgt] <> "_" <> ToString[j] <> ".m", solOut[j]], {j, 1, 6}];
Export["work/sol" <> ToString[wgt] <> "_xmy.m", solOut[xmy]];
Export["work/sol" <> ToString[wgt] <> "_x1.m", solOut[x1]];
Export["work/sol" <> ToString[wgt] <> "_y1.m", solOut[y1]];
Print["=== DONE: exported sol", wgt, "_{1..6,x1,y1,xmy}.m ==="];
