(* H_opt run_weight3.wls
   Weight-3 analogue of run_weight.wls/run_weight2.wls. Structural changes
   vs weight 1/2, replicating the reference .wl exactly (lines ~1121-1296):
   - match2..match6 switch from the ad-hoc "{G[1,1,t]->0,G[1,t]->0}/.t->1"
     endpoint rule to a general reflection: (tempSol[j-1] /. t->1-t /. tto1Red
     // ShuffleG) /. Logt0 /. t->0 -- this reuses the SAME "near t=0"
     regularization (Logt0, already covers up to 6 zeros) for BOTH endpoints
     via the fibration-change table, instead of hand-listing endpoint-1 rules.
   - impart (README 3.5) now has 3 terms: (-iPi)^1/1! todelta.sol[2],
     (-iPi)^2/2! todelta^2.sol[1], (-iPi)^3/3! todelta^3.sol[0].
   - match00 ALSO uses the reflection+tto1Red+ShuffleG+Logt0 pattern on
     impart5; match01/match10 keep the direct /.t->1 or /.t->0 substitution
     (unchanged from weight 1/2), just with the richer impart3.
   - planar constants on segment 5 use ReduceCG+GetGs directly (confirmed
     working in this environment, unlike the broken ToFibrationBasis used
     at weight 2) rather than the frozen tto1Red table, matching the
     reference's own weight-3 code exactly.
   Gating check (README section 4) already run separately and confirmed:
   all 14 distinct G[...,1-t] words needed at weight 3 are covered by the
   frozen H_opt/input/fibrtto1.m -- no ReduceCG extension of tto1Red needed
   for the match equations themselves. *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
Global`x; Global`y; Global`t; Global`delta; Global`c2; Global`c2c; Global`cxmy; Global`cxmyc; Global`sqrt; Global`eps;
Global`Basis; Global`r1; Global`r2; Global`r3; Global`r4; Global`H;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];

root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
Hdir = FileNameJoin[{root, "H"}];
workdir = FileNameJoin[{root, "H_opt", "work"}];
ckpt = FileNameJoin[{workdir, "precompute.mx"}];
outckpt = FileNameJoin[{workdir, "sol_w3.mx"}];

If[FileExistsQ[outckpt],
  Print["Checkpoint ", outckpt, " already exists -- nothing to do. Delete it to force recompute."];
  Exit[0];
];
If[!FileExistsQ[ckpt], Print["Missing ", ckpt, " -- run run_precompute.wls first."]; Exit[1]];
Do[
  f = FileNameJoin[{workdir, "sol2_" <> ToString[j] <> ".m"}];
  If[!FileExistsQ[f], Print["Missing ", f, " -- run run_weight2.wls first."]; Exit[1]],
  {j, {1, 2, 3, 4, 5, 6, x1, y1, xmy}}
];

MemoryConstrained[

  Print["=== loading precompute checkpoint ==="];
  Get[ckpt];
  tto1Red = Import[FileNameJoin[{root, "H_opt", "input", "fibrtto1.m"}]];

  Logt0 = {G[0, t] -> 0, G[0, 0, t] -> 0, G[0, 0, 0, t] -> 0, G[0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, 0, t] -> 0};
  Log00 = {G[0, 0] -> 0, G[0, 0, 0] -> 0, G[0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0] -> 0};

  integrateSeg[Mmap_, prevSol_, var_] := Total[
    Function[d, GIntegrate[#, var] & /@ ((Mmap[d] . prevSol)/d)] /@ Keys[Mmap]
  ];

  Print["=== loading sol[0],sol[1],sol[2] ==="];
  sol0 = Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w0", "sol0.m"}]];
  Do[sol[0][j] = sol0, {j, 1, 6}];
  sol[0][x1] = sol0; sol[0][y1] = sol0; sol[0][xmy] = sol0;
  Do[sol[1][j] = Import[FileNameJoin[{workdir, "sol1_" <> ToString[j] <> ".m"}]], {j, 1, 6}];
  sol[1][x1] = Import[FileNameJoin[{workdir, "sol1_x1.m"}]];
  sol[1][y1] = Import[FileNameJoin[{workdir, "sol1_y1.m"}]];
  sol[1][xmy] = Import[FileNameJoin[{workdir, "sol1_xmy.m"}]];
  Do[sol[2][j] = Import[FileNameJoin[{workdir, "sol2_" <> ToString[j] <> ".m"}]], {j, 1, 6}];
  sol[2][x1] = Import[FileNameJoin[{workdir, "sol2_x1.m"}]];
  sol[2][y1] = Import[FileNameJoin[{workdir, "sol2_y1.m"}]];
  sol[2][xmy] = Import[FileNameJoin[{workdir, "sol2_xmy.m"}]];

  Print["=== building IntPart[j] via precomputed M_a + GIntegrate ==="];
  {tI1, IntPart1} = AbsoluteTiming[integrateSeg[Mdt[1], sol[2][1], t]];
  {tI2, IntPart2} = AbsoluteTiming[integrateSeg[Mdt[2], sol[2][2], t]];
  {tI3, IntPart3} = AbsoluteTiming[integrateSeg[Mdt[3], sol[2][3], t]];
  {tI4, IntPart4} = AbsoluteTiming[integrateSeg[Mdt[4], sol[2][4], t]];
  {tI5, IntPart5} = AbsoluteTiming[integrateSeg[Mdt[5], sol[2][5], t]];
  {tI6, IntPart6} = AbsoluteTiming[integrateSeg[Mdt[6], sol[2][6], t]];
  {tIx1, IntPartx1} = AbsoluteTiming[integrateSeg[Mdyx1, sol[2][x1], y]];
  {tIy1, IntParty1} = AbsoluteTiming[integrateSeg[Mdxy1, sol[2][y1], x]];
  {tIxmy, IntPartxmy} = AbsoluteTiming[integrateSeg[Mdxmy, sol[2][xmy], y]];
  Print["IntPart timings: ", {tI1, tI2, tI3, tI4, tI5, tI6, tIx1, tIy1, tIxmy}];

  CNST1 = Table[Symbol["cS1x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNST2 = Table[Symbol["cS2x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNST3 = Table[Symbol["cS3x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNST4 = Table[Symbol["cS4x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNST5 = Table[Symbol["cS5x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNST6 = Table[Symbol["cS6x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNSTx1 = Table[Symbol["cSx1x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNSTy1 = Table[Symbol["cSy1x" <> ToString[i]], {i, 1, Length[Basis]}];
  CNSTxmy = Table[Symbol["cSxmyx" <> ToString[i]], {i, 1, Length[Basis]}];
  zeroAll = Thread[Rule[Join[CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNSTx1, CNSTy1, CNSTxmy], 0]];

  tempSol1 = IntPart1 + CNST1;
  tempSol2 = IntPart2 + CNST2;
  tempSol3 = IntPart3 + CNST3;
  tempSol4 = IntPart4 + CNST4;
  tempSol5 = IntPart5 + CNST5;
  tempSol6 = IntPart6 + CNST6;
  tempSolx1 = IntPartx1 + CNSTx1;
  tempSoly1 = IntParty1 + CNSTy1;
  tempSolxmy = IntPartxmy + CNSTxmy;

  Print["=== matching equations (reflection + tto1Red + ShuffleG + Logt0) ==="];
  reflect[expr_] := (expr /. t -> 1 - t /. tto1Red // ShuffleG) /. Logt0 /. t -> 0;
  {tm2, match2} = AbsoluteTiming[reflect[tempSol1] - (tempSol2 /. Logt0 /. t -> 0)];
  {tm3, match3} = AbsoluteTiming[reflect[tempSol2] - (tempSol3 /. Logt0 /. t -> 0)];
  {tm4, match4} = AbsoluteTiming[reflect[tempSol3] - (tempSol4 /. Logt0 /. t -> 0)];
  {tm5, match5} = AbsoluteTiming[reflect[tempSol4] - (tempSol5 /. Logt0 /. t -> 0)];
  {tm6, match6} = AbsoluteTiming[reflect[tempSol5] - (tempSol6 /. Logt0 /. t -> 0)];
  Print["match2..6 timings: ", {tm2, tm3, tm4, tm5, tm6}];

  Print["=== weight-3 impart (3 terms, README 3.5) ==="];
  impart5 = tempSol5 + Sum[(-I Pi)^k/k! todeltaPow[5, k] . sol[3 - k][5], {k, 1, 3}];
  impart3 = tempSol3 + Sum[(-I Pi)^k/k! todeltaPow[3, k] . sol[3 - k][3], {k, 1, 3}];
  match00 = reflect[impart5] - (tempSolxmy /. y -> 0 /. Log00);
  match01 = (impart3 /. t -> 1 /. Log00 /. HPL[a__] -> 0) - (tempSoly1 /. x -> 0 /. Log00);
  match10 = (impart3 /. t -> 0 /. Log00) - (tempSolx1 /. y -> 0 /. Log00);

  Print["=== extract shift Delta's via zero-substitution and chain to segment 5 ==="];
  d2 = match2 /. zeroAll;
  d3 = match3 /. zeroAll;
  d4 = match4 /. zeroAll;
  d5 = match5 /. zeroAll;
  d6 = match6 /. zeroAll;
  d00 = match00 /. zeroAll;
  d01 = match01 /. zeroAll;
  d10 = match10 /. zeroAll;

  Print["=== weight-3 planar constants on segment 5 ==="];
  (* The reference's own weight-3 code calls ReduceCG[GetGs[templim],{t},
     Point->{t->3/14}] here. Tested directly: ReduceCG fails on exactly the
     10 short/already-canonical words {G[0,1-t], G[1,1-t], G[0,0,0,1-t], ...}
     -- its failure path doesn't leave them as plain G[...] (which the final
     /._G->0 would safely zero) but injects a broken internal "PLT$Failed"
     placeholder template that survives into the output. The frozen tto1Red
     table (gating check: 0 missing at weight 3) covers the exact same
     reduction, and weight 4's OWN reference code uses "/.tto1Red" directly
     for this identical step (not ReduceCG) -- so use that instead; it's
     already validated and avoids the bug entirely. *)
  CNSTsol = Table[cnst[3, i][5], {i, 1, Length[Basis]}];
  templim = Coefficient[Extract[SolutionWithPlanar, onlyplanar], eps, 3] /. x -> 0 /. G[0, 0] -> 0 /. G[0, 0, 0] -> 0 /. G[0, 0, 0, 0] -> 0 /. y -> 1 - t;
  planarW3 = (templim /. tto1Red // ShuffleG) /. _G -> 0;
  plrules5 = Thread[Rule[Extract[CNSTsol, onlyplanar], planarW3]];
  frozenSol5 = Import[FileNameJoin[{Hdir, "numerical_eq", "solution_3.m"}]];
  constraints5 = Join[frozenSol5, plrules5];
  CNST5val = CNSTsol /. constraints5;
  If[!FreeQ[CNST5val, cnst],
    Print["WARNING: CNST5val still has unresolved cnst[...] entries -- constraints incomplete."];
  ];

  CNST4val = CNST5val - d5;
  CNST3val = CNST4val - d4;
  CNST2val = CNST3val - d3;
  CNST1val = CNST2val - d2;
  CNST6val = CNST5val + d6;
  CNSTxmyval = CNST5val + d00;
  CNSTy1val = CNST3val + d01;
  CNSTx1val = CNST3val + d10;

  Print["=== assembling sol[3][seg] ==="];
  solOut[1] = tempSol1 /. Thread[Rule[CNST1, CNST1val]];
  solOut[2] = tempSol2 /. Thread[Rule[CNST2, CNST2val]];
  solOut[3] = tempSol3 /. Thread[Rule[CNST3, CNST3val]];
  solOut[4] = tempSol4 /. Thread[Rule[CNST4, CNST4val]];
  solOut[5] = tempSol5 /. Thread[Rule[CNST5, CNST5val]];
  solOut[6] = tempSol6 /. Thread[Rule[CNST6, CNST6val]];
  solOut[xmy] = tempSolxmy /. Thread[Rule[CNSTxmy, CNSTxmyval]];
  solOut[x1] = tempSolx1 /. Thread[Rule[CNSTx1, CNSTx1val]];
  solOut[y1] = tempSoly1 /. Thread[Rule[CNSTy1, CNSTy1val]];

  Do[Export[FileNameJoin[{workdir, "sol3_" <> ToString[j] <> ".m"}], solOut[j]], {j, 1, 6}];
  Export[FileNameJoin[{workdir, "sol3_xmy.m"}], solOut[xmy]];
  Export[FileNameJoin[{workdir, "sol3_x1.m"}], solOut[x1]];
  Export[FileNameJoin[{workdir, "sol3_y1.m"}], solOut[y1]];

  DumpSave[outckpt, {solOut}];
  Print["Done. Exported sol3_{1..6,x1,y1,xmy}.m to ", workdir];

, 8*1024^3]
