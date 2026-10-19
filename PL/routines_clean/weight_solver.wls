(* ============================================================================
   weight_solver.wls  --  stage 3 of the PL pipeline: solve for sol[w][seg] at
   a single weight w, given sol[w-1][*] already on disk.
   Clean reorganized copy of ../run_weight_pl.wls. See README.md for the full
   narrative (why each piece exists, which bugs it fixes); comments here are
   kept short.
   Usage: wolframscript -file weight_solver.wls <w>
   ============================================================================ *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
Global`x; Global`y; Global`t; Global`delta; Global`ep; Global`eps;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];
SetDirectory[DirectoryName[$InputFileName]];

wgt = ToExpression[$ScriptCommandLine[[2]]];
n = 291;
Print["=== PL weight ", wgt, " ==="];

(* ---------------------------------------------------------------------
   section 1: load weight-independent matrices + all sol[w-1][*] inputs
   --------------------------------------------------------------------- *)
loadPrecompute[] := (
   Get["../work/precompute_pl.mx"];
   {allPx0, allPy0, allPxy1} = Import /@ {"../work/allP_x0.m", "../work/allP_y0.m", "../work/allP_xy1.m"};
   {resX1, resY1, resXmY} = Import /@ {"../work/regX1.m", "../work/regY1.m", "../work/regXmY.m"};
   {stackX0, stackY0, stackXy1} = Join @@ # & /@ {allPx0, allPy0, allPxy1};
   );

loadPreviousWeightSolutions[wgt_] := (
   sol0 = Import["../work/weight0_solution.m"];
   Do[sol[0][j] = sol0, {j, 1, 6}];
   sol[0][x1] = sol0; sol[0][y1] = sol0; sol[0][xmy] = sol0;
   Do[sol[k][j] = Import["../work/sol" <> ToString[k] <> "_" <> ToString[j] <> ".m"], {j, 1, 6}, {k, 1, wgt - 1}];
   Do[
     sol[k][x1] = Import["../work/sol" <> ToString[k] <> "_x1.m"];
     sol[k][y1] = Import["../work/sol" <> ToString[k] <> "_y1.m"];
     sol[k][xmy] = Import["../work/sol" <> ToString[k] <> "_xmy.m"];
     , {k, 1, wgt - 1}];
   );

(* ---------------------------------------------------------------------
   section 2: IntPart[seg] = GIntegrate[ dt[seg] . sol[w-1][seg] ], with
   per-component progress logging and disk checkpointing (GIntegrate is the
   single most expensive step; checkpointing means a kill only loses the
   segment currently in progress).
   --------------------------------------------------------------------- *)
integrateWithProgress[label_, rawVec_, var_] := Module[{len, result, t0, tComp, step},
   len = Length[rawVec]; step = 20;
   result = Table[Null, {len}];
   t0 = AbsoluteTime[];
   Do[
     {tComp, result[[i]]} = AbsoluteTiming[GIntegrate[rawVec[[i]], var]];
     If[Mod[i, step] == 0 || i == len || tComp > 5,
       Print["    ", label, " [", i, "/", len, "]  this-component=", tComp, "s  elapsed=", AbsoluteTime[] - t0, "s"]
       ];
     , {i, len}];
   result
   ];

getOrCompute[label_, filename_, rawVec_, var_] := Module[{tC, result},
   If[FileExistsQ[filename],
     Print["  ", label, ": found existing checkpoint ", filename, ", loading"];
     {0., Import[filename]},
     {tC, result} = AbsoluteTiming[integrateWithProgress[label, rawVec, var]];
     Export[filename, result];
     Print["  ", label, ": computed and exported to ", filename];
     {tC, result}
     ]
   ];

ipFile[wgt_, seg_] := "../work/intpart" <> ToString[wgt] <> "_" <> ToString[seg] <> ".m";

buildAllIntParts[wgt_] := Module[{t},
   {t[1], IntPart1} = getOrCompute["IntPart1", ipFile[wgt, 1], Apart[dt[1].sol[wgt - 1][1], t], t];
   {t[2], IntPart2} = getOrCompute["IntPart2", ipFile[wgt, 2], Apart[dt[2].sol[wgt - 1][2], t], t];
   {t[3], IntPart3} = getOrCompute["IntPart3", ipFile[wgt, 3], Apart[dt[3].sol[wgt - 1][3], t], t];
   {t[4], IntPart4} = getOrCompute["IntPart4", ipFile[wgt, 4], Apart[dt[4].sol[wgt - 1][4], t], t];
   {t[5], IntPart5} = getOrCompute["IntPart5", ipFile[wgt, 5], Apart[dt[5].sol[wgt - 1][5], t], t];
   {t[6], IntPart6} = getOrCompute["IntPart6", ipFile[wgt, 6], Apart[dt[6].sol[wgt - 1][6], t], t];
   {t[x1], IntPartx1} = getOrCompute["IntPartx1", ipFile[wgt, "x1"], Apart[dyx1.sol[wgt - 1][x1], y], y];
   {t[y1], IntParty1} = getOrCompute["IntParty1", ipFile[wgt, "y1"], Apart[dxy1.sol[wgt - 1][y1], x], x];
   {t[xmy], IntPartxmy} = getOrCompute["IntPartxmy", ipFile[wgt, "xmy"], Apart[dxmy.sol[wgt - 1][xmy], y], y];
   Print["IntPart timings: ", {t[1], t[2], t[3], t[4], t[5], t[6], t[x1], t[y1], t[xmy]}];
   ];

(* ---------------------------------------------------------------------
   section 3: t->0 / t->1 evaluation operators.
   evalStart: t->0, always just a direct substitution (no regularization
   needed at the lower endpoint).
   evalAt1 dispatches on weight:
     w<=2  : kill only the literal all-ones word, then substitute t->1
             directly (MZV table auto-evaluates the rest).
     3<=w<=6: full fibration-change (t->1-t, reduce via tto1Red, ShuffleG),
             then kill degenerate all-zero words via Log00.
     w>=7  : ReduceCG/GetGs hard-caps at weight<=6, so instead use the
             precomputed word_table_0to8.m (see word_table_builder.wls):
             plain t->1 substitution auto-separates finite MZV content from
             an explicit HPL[{1},1] divergence marker, which is dropped.
   --------------------------------------------------------------------- *)
(* kills degenerate all-zero words (no "t") that lose t-dependence mid-reduction *)
Log00 = {G[0, 0] -> 0, G[0, 0, 0] -> 0, G[0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0] -> 0,
   G[0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0, 0] -> 0};
endPt1 = Table[G[Sequence @@ ConstantArray[1, k], t] -> 0, {k, 1, 8}];

evalStart[expr_, var_] := expr /. G[__, var] -> 0 /. var -> 0;
evalAt0 = evalStart[#, t] &;

evalEndSimple[expr_] := expr /. endPt1 /. t -> 1;
evalEndFibration[expr_] := ((expr /. t -> 1 - t /. tto1Red // ShuffleG) /. t -> 0 /. G[__, t] -> 0 /. HPL[a__] -> 0);

(* extends tto1Red at weight>=3 with any genuinely new words PL's own
   alphabet introduces beyond H's frozen 126-rule table; mutates tto1Red. *)
extendTto1Red[expr_] := Module[{new, oldLhs},
   new = ReduceCG[GetGs[expr /. t -> 1 - t], {t}, Point -> {t -> 3/14}];
   oldLhs = tto1Red[[All, 1]];
   new = Select[new, ! MemberQ[oldLhs, #[[1]]] &];
   tto1Red = Join[tto1Red, new];
   ];

evalEndTable[expr_, table_] := expr /. table;

setupEvalAt1[wgt_] := (
   tto1Red = Import["../../H_opt/input/fibrtto1.m"];
   wordTableAll = Import["../work/word_table_0to8.m"] /. (lhs_ -> rhs_) :> (lhs -> Expand[rhs]);
   evalAt1 = Which[
      wgt <= 2, evalEndSimple[#],
      wgt <= 6, (extendTto1Red[#]; evalEndFibration[#] /. Log00),
      True, evalEndTable[#, wordTableAll]
      ] &;
   );

(* ---------------------------------------------------------------------
   section 4: matching -- build delta vectors d2..d6 (adjacent segments)
   and d00,d01,d10 (x1/y1/xmy regularity connections via impartSum, the
   analytic continuation across the branch cut at each segment's delta->0
   endpoint).
   --------------------------------------------------------------------- *)
impartSum[j_, wgt_] := Sum[((-I Pi)^k/k!) MatrixPower[todelta[j], k] . sol[wgt - k][j], {k, 1, wgt}];

buildDeltaVectors[wgt_] := Module[{},
   d2 = evalAt1[IntPart1] - evalAt0[IntPart2];
   d3 = evalAt1[IntPart2] - evalAt0[IntPart3];
   d4 = evalAt1[IntPart3] - evalAt0[IntPart4];
   d5 = evalAt1[IntPart4] - evalAt0[IntPart5];
   d6 = evalAt1[IntPart5] - evalAt0[IntPart6];

   impart5sol = impartSum[5, wgt];
   d00 = If[wgt == 1,
      (evalAt0[impart5sol + IntPart5] /. HPL[a__] -> 0) - evalStart[IntPartxmy, y],
      (evalAt1[impart5sol + IntPart5] /. Log00 /. HPL[a__] -> 0) - (evalStart[IntPartxmy, y] /. Log00)
      ];

   impart3sol = impartSum[3, wgt];
   tmp01 = If[wgt <= 2, (impart3sol + IntPart3) /. t -> 1 /. Log00 /. HPL[a__] -> 0, evalAt1[impart3sol + IntPart3]];
   d01 = tmp01 - evalStart[IntParty1, x];
   tmp10 = If[wgt <= 2, (impart3sol + IntPart3) /. t -> 0 /. Log00 /. HPL[a__] -> 0,
      evalStart[impart3sol + IntPart3, t] /. HPL[a__] -> 0];
   d10 = tmp10 - evalStart[IntPartx1, y];
   ];

(* ---------------------------------------------------------------------
   section 5: boundary + regularity equations, solve for the single unknown
   vector CNST5 (291 constants of integration for the base segment), using
   the one known external input: the fully massless banana's exact value.
   --------------------------------------------------------------------- *)
knownBananaValue[wgt_] := -Coefficient[
   1 - Pi^2/4 eps^2 - 29 Zeta[3] eps^3 - 71 Pi^4/160 eps^4 + (29 Pi^2 Zeta[3]/4 - 1263 Zeta[5]/5) eps^5
     + (841 Zeta[3]^2/2 - 11539 Pi^6/24192) eps^6
     + (2059 Pi^4 Zeta[3]/160 + 1263 Pi^2 Zeta[5]/20 - 18567 Zeta[7]/7) eps^7
     + (-841 Pi^2 Zeta[3]^2/8 + 36627 Zeta[3] Zeta[5]/5 - 737687 Pi^8/1382400) eps^8,
   eps, wgt];

solveForCNST5[wgt_] := Module[{eq1, eq2, eq7, eq4, eq5, eq6, allEqs, sr},
   CNST5 = Table[cS5[i], {i, 1, n}];
   CNST4sym = CNST5 - d5; CNST3sym = CNST4sym - d4; CNST2sym = CNST3sym - d3; CNST1sym = CNST2sym - d2;
   CNST6sym = CNST5 + d6; CNSTxmysym = CNST5 + d00; CNSTy1sym = CNST3sym + d01; CNSTx1sym = CNST3sym + d10;

   tempSol5 = IntPart5 + CNST5;
   tempSol1 = IntPart1 + CNST1sym;
   tempSol3 = IntPart3 + CNST3sym;
   tempSolx1 = IntPartx1 + CNSTx1sym;
   tempSoly1 = IntParty1 + CNSTy1sym;
   tempSolxmy = IntPartxmy + CNSTxmysym;

   eq1 = stackX0.tempSol5; eq2 = stackY0.tempSol1; eq7 = stackXy1.tempSol3;
   eq4 = resX1.tempSolx1; eq5 = resY1.tempSoly1; eq6 = resXmY.tempSolxmy;

   knownBanana = knownBananaValue[wgt];
   Print["known banana boundary value at weight ", wgt, " (position 1): ", knownBanana];

   allEqs = Join[Thread[eq1 == 0], Thread[eq2 == 0], Thread[eq7 == 0],
      Thread[eq4 == 0], Thread[eq5 == 0], Thread[eq6 == 0], {CNST5[[1]] == knownBanana}] // Union;
   Print["total equations: ", Length[allEqs]];

   {tSolve, sr} = AbsoluteTiming[Solve[allEqs, CNST5]];
   Print["Solve time: ", tSolve, "s, num solutions: ", Length[sr]];
   sr[[1]]
   ];

(* ---------------------------------------------------------------------
   section 6: back-substitute, ShuffleG-simplify (prevents weight-over-
   weight G-word bloat -- this was the real cause of the weight-8
   IntParty1 2h+ slowdown, not GIntegrate itself), and export.
   --------------------------------------------------------------------- *)
backSubstituteAndExport[wgt_, s1_] := Module[{solOut},
   CNST4val = CNST4sym /. s1; CNST3val = CNST3sym /. s1; CNST2val = CNST2sym /. s1; CNST1val = CNST1sym /. s1;
   CNST6val = CNST6sym /. s1; CNSTxmyval = CNSTxmysym /. s1; CNSTy1val = CNSTy1sym /. s1; CNSTx1val = CNSTx1sym /. s1;

   solOut[1] = ShuffleG /@ (IntPart1 + CNST1val);
   solOut[2] = ShuffleG /@ (IntPart2 + CNST2val);
   solOut[3] = ShuffleG /@ (IntPart3 + CNST3val);
   solOut[4] = ShuffleG /@ (IntPart4 + CNST4val);
   solOut[5] = ShuffleG /@ (IntPart5 + CNST5val);
   solOut[6] = ShuffleG /@ (IntPart6 + CNST6val);
   solOut[xmy] = ShuffleG /@ (IntPartxmy + CNSTxmyval);
   solOut[x1] = ShuffleG /@ (IntPartx1 + CNSTx1val);
   solOut[y1] = ShuffleG /@ (IntParty1 + CNSTy1val);

   Do[Export["../work/sol" <> ToString[wgt] <> "_" <> ToString[j] <> ".m", solOut[j]], {j, 1, 6}];
   Export["../work/sol" <> ToString[wgt] <> "_xmy.m", solOut[xmy]];
   Export["../work/sol" <> ToString[wgt] <> "_x1.m", solOut[x1]];
   Export["../work/sol" <> ToString[wgt] <> "_y1.m", solOut[y1]];
   Print["=== DONE: exported sol", wgt, "_{1..6,x1,y1,xmy}.m ==="];
   ];

(* ============================ run ============================ *)
loadPrecompute[];
loadPreviousWeightSolutions[wgt];
buildAllIntParts[wgt];
setupEvalAt1[wgt];
buildDeltaVectors[wgt];
s1 = solveForCNST5[wgt];
backSubstituteAndExport[wgt, s1];
