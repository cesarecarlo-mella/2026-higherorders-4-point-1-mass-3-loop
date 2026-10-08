(* H_opt run_weight2.wls
   Weight-2 analogue of run_weight.wls. Differs from weight 1 in two ways
   that the reference notebook itself requires (reference lines 957-1119):
   - the analytic-continuation "impart" now has two terms, order (-i*Pi)^1
     (todelta.sol[1]) and order (-i*Pi)^2/2 (todelta^2.sol[0]) -- README 3.5;
   - plrules[2][5] needs ToFibrationBasis (G[...,1-t] -> G[...,t]) + ShuffleG
     before the final _G->0 cleanup, because at weight 2 the eps^2 coefficient
     of the planar solution is no longer already letter-free once y->1-t is
     substituted (unlike weight 1, where /._G->0 alone sufficed).
   Everything else (IntPart via precomputed M_a + GIntegrate, matching via
   plain vector subtraction, reuse of the frozen H/numerical_eq/solution_2.m)
   is identical in structure to run_weight.wls. *)

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
outckpt = FileNameJoin[{workdir, "sol_w2.mx"}];

If[FileExistsQ[outckpt],
  Print["Checkpoint ", outckpt, " already exists -- nothing to do. Delete it to force recompute."];
  Exit[0];
];
If[!FileExistsQ[ckpt], Print["Missing ", ckpt, " -- run run_precompute.wls first."]; Exit[1]];
Do[
  f = FileNameJoin[{workdir, "sol1_" <> ToString[j] <> ".m"}];
  If[!FileExistsQ[f], Print["Missing ", f, " -- run run_weight.wls 1 first."]; Exit[1]],
  {j, {1, 2, 3, 4, 5, 6, x1, y1, xmy}}
];

MemoryConstrained[

  Print["=== loading precompute checkpoint ==="];
  Get[ckpt];

  Logt0 = {G[0, t] -> 0, G[0, 0, t] -> 0, G[0, 0, 0, t] -> 0, G[0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, 0, t] -> 0};
  Log00 = {G[0, 0] -> 0, G[0, 0, 0] -> 0, G[0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0] -> 0};

  integrateSeg[Mmap_, prevSol_, var_] := Total[
    Function[d, GIntegrate[#, var] & /@ ((Mmap[d] . prevSol)/d)] /@ Keys[Mmap]
  ];

  Print["=== loading sol[0] (frozen) and sol[1] (my own validated weight-1 output) ==="];
  sol0 = Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w0", "sol0.m"}]];
  Do[sol[1][j] = Import[FileNameJoin[{workdir, "sol1_" <> ToString[j] <> ".m"}]], {j, 1, 6}];
  sol[1][x1] = Import[FileNameJoin[{workdir, "sol1_x1.m"}]];
  sol[1][y1] = Import[FileNameJoin[{workdir, "sol1_y1.m"}]];
  sol[1][xmy] = Import[FileNameJoin[{workdir, "sol1_xmy.m"}]];

  Print["=== building IntPart[j] via precomputed M_a + GIntegrate ==="];
  {tI1, IntPart1} = AbsoluteTiming[integrateSeg[Mdt[1], sol[1][1], t]];
  {tI2, IntPart2} = AbsoluteTiming[integrateSeg[Mdt[2], sol[1][2], t]];
  {tI3, IntPart3} = AbsoluteTiming[integrateSeg[Mdt[3], sol[1][3], t]];
  {tI4, IntPart4} = AbsoluteTiming[integrateSeg[Mdt[4], sol[1][4], t]];
  {tI5, IntPart5} = AbsoluteTiming[integrateSeg[Mdt[5], sol[1][5], t]];
  {tI6, IntPart6} = AbsoluteTiming[integrateSeg[Mdt[6], sol[1][6], t]];
  {tIx1, IntPartx1} = AbsoluteTiming[integrateSeg[Mdyx1, sol[1][x1], y]];
  {tIy1, IntParty1} = AbsoluteTiming[integrateSeg[Mdxy1, sol[1][y1], x]];
  {tIxmy, IntPartxmy} = AbsoluteTiming[integrateSeg[Mdxmy, sol[1][xmy], y]];
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

  Print["=== matching equations (plain vector subtraction, no Solve) ==="];
  match2 = (tempSol1 /. {G[1, 1, t] -> 0, G[1, t] -> 0} /. t -> 1) - (tempSol2 /. {G[0, 0, t] -> 0, G[0, t] -> 0} /. t -> 0);
  match3 = (tempSol2 /. {G[1, 1, t] -> 0, G[1, t] -> 0} /. t -> 1) - (tempSol3 /. {G[0, 0, t] -> 0, G[0, t] -> 0} /. t -> 0);
  match4 = (tempSol3 /. {G[1, 1, t] -> 0, G[1, t] -> 0} /. t -> 1) - (tempSol4 /. {G[0, 0, t] -> 0, G[0, t] -> 0} /. t -> 0);
  match5 = (tempSol4 /. {G[1, 1, t] -> 0, G[1, t] -> 0} /. t -> 1) - (tempSol5 /. {G[0, 0, t] -> 0, G[0, t] -> 0} /. t -> 0);
  match6 = (tempSol5 /. {G[1, 1, t] -> 0, G[1, t] -> 0} /. t -> 1) - (tempSol6 /. {G[0, 0, t] -> 0, G[0, t] -> 0} /. t -> 0);

  Print["=== weight-2 impart (todelta^1 and todelta^2 terms, README 3.5) ==="];
  impart5 = tempSol5 + todelta[5] . sol[1][5] (tag I Pi) + (tag I Pi)^2 (todelta[5] . todelta[5]) . sol0/2 /. tag -> -1;
  match00 = (impart5 /. G[1, t] -> 0 /. G[1, 1, t] -> 0 /. t -> 1 /. Log00) - (tempSolxmy /. y -> 0 /. Log00);
  impart3 = tempSol3 + todelta[3] . sol[1][3] (tag I Pi) + (tag I Pi)^2 (todelta[3] . todelta[3]) . sol0/2 /. tag -> -1;
  match01 = (impart3 /. t -> 1 /. G[0, 0] -> 0 /. HPL[a__] -> 0) - (tempSoly1 /. x -> 0 /. Log00);
  match10 = (impart3 /. t -> 0 /. HPL[a__] -> 0 /. Log00) - (tempSolx1 /. y -> 0 /. Log00);

  Print["=== extract shift Delta's via zero-substitution and chain to segment 5 ==="];
  d2 = match2 /. zeroAll;
  d3 = match3 /. zeroAll;
  d4 = match4 /. zeroAll;
  d5 = match5 /. zeroAll;
  d6 = match6 /. zeroAll;
  d00 = match00 /. zeroAll;
  d01 = match01 /. zeroAll;
  d10 = match10 /. zeroAll;

  Print["=== weight-2 planar constants on segment 5 (G[...,1-t]->G[...,t] + ShuffleG) ==="];
  (* ToFibrationBasis itself is broken in this environment -- it shells out to
     ginsh for a custom GiNaC extension ('PolylogTools_H') that this machine's
     ginsh build doesn't have, confirmed hanging/erroring on even a single
     trivial G[1,1,1-t] term. H_opt/input/fibrtto1.m (tto1Red) already IS the
     frozen G[...,1-t]->G[...,t] reduction table ToFibrationBasis would have
     produced (README section 4), so use it directly instead. *)
  tto1Red = Import[FileNameJoin[{root, "H_opt", "input", "fibrtto1.m"}]];
  CNSTsol = Table[cnst[2, i][5], {i, 1, Length[Basis]}];
  preRed = Coefficient[Extract[SolutionWithPlanar, onlyplanar], eps, 2] /. x -> 0 /. G[0, 0] -> 0 /. G[0, 0, 0] -> 0 /. y -> 1 - t;
  survivingOneMinusT = Union[Cases[preRed /. tto1Red, G[a___, 1 - t] :> G[a, 1 - t], Infinity]];
  If[survivingOneMinusT =!= {},
    Print["WARNING: ", Length[survivingOneMinusT], " G[...,1-t] terms are NOT covered by tto1Red: ", survivingOneMinusT];
  ];
  {tFib, planarW2} = AbsoluteTiming[(preRed /. tto1Red // ShuffleG) /. _G -> 0];
  Print["tto1Red+ShuffleG time: ", tFib, "s"];
  plrules5 = Thread[Rule[Extract[CNSTsol, onlyplanar], planarW2]];
  frozenSol5 = Import[FileNameJoin[{Hdir, "numerical_eq", "solution_2.m"}]];
  constraints5 = Join[frozenSol5, plrules5];
  CNST5val = CNSTsol /. constraints5;
  If[!FreeQ[CNST5val, cnst],
    Print["WARNING: CNST5val still has unresolved cnst[...] entries -- constraints incomplete."];
  ];

  CNST4val = CNST5val - d5;
  CNST3val = CNST4val - d4;
  CNST2val = CNST3val - d3;
  CNST1val = CNST2val - d2;
  (* CNST6 sits in the MINUS position of match6 = tempSol5|t1 - tempSol6|t0
     (unlike CNST1..4, each of which is the PLUS/first-term variable in ITS
     OWN match equation, correctly giving "-d"); solving match6=0 for CNST6
     gives CNST5 + d6, not CNST5 - d6. Confirmed by re-deriving the algebra
     and checking row 7: CNST5val[[7]]=-Pi^2/9, d6[[7]]=Pi^2/9, needs to sum
     to 0 (the reference's sol2_6[[7]]), which only "+d6" gives. *)
  CNST6val = CNST5val + d6;
  CNSTxmyval = CNST5val + d00;
  CNSTy1val = CNST3val + d01;
  CNSTx1val = CNST3val + d10;

  Print["=== assembling sol[2][seg] ==="];
  solOut[1] = tempSol1 /. Thread[Rule[CNST1, CNST1val]];
  solOut[2] = tempSol2 /. Thread[Rule[CNST2, CNST2val]];
  solOut[3] = tempSol3 /. Thread[Rule[CNST3, CNST3val]];
  solOut[4] = tempSol4 /. Thread[Rule[CNST4, CNST4val]];
  solOut[5] = tempSol5 /. Thread[Rule[CNST5, CNST5val]];
  solOut[6] = tempSol6 /. Thread[Rule[CNST6, CNST6val]];
  solOut[xmy] = tempSolxmy /. Thread[Rule[CNSTxmy, CNSTxmyval]];
  solOut[x1] = tempSolx1 /. Thread[Rule[CNSTx1, CNSTx1val]];
  solOut[y1] = tempSoly1 /. Thread[Rule[CNSTy1, CNSTy1val]];

  Do[Export[FileNameJoin[{workdir, "sol2_" <> ToString[j] <> ".m"}], solOut[j]], {j, 1, 6}];
  Export[FileNameJoin[{workdir, "sol2_xmy.m"}], solOut[xmy]];
  Export[FileNameJoin[{workdir, "sol2_x1.m"}], solOut[x1]];
  Export[FileNameJoin[{workdir, "sol2_y1.m"}], solOut[y1]];

  DumpSave[outckpt, {solOut}];
  Print["Done. Exported sol2_{1..6,x1,y1,xmy}.m to ", workdir];

, 8*1024^3]
