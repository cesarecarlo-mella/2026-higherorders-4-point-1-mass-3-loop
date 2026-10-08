(* H_opt run_weight4.wls
   Weight-4 analogue of run_weight3.wls (same reflection+tto1Red+ShuffleG+
   Logt0 pattern for match2-6/match00, same direct-substitution pattern for
   match01/match10), with:
   - impart now has 4 terms (README 3.5): k=1..4, (-iPi)^k/k! todelta^k.sol[4-k].
   - planar constants on segment 5 via "/.tto1Red" directly (matches the
     reference's OWN weight-4 code exactly -- no ReduceCG at all at this
     weight, unlike weight 3's reference code; gating check confirmed 0
     missing for both the match-equations and this planar step).
   - x1/y1/xmy get an extra Collect[...,{_G},FullSimplify[#/.rootsXmy]&]
     algebraic-reduction pass (README 3.6), first needed starting here per
     the reference's weight-4 code.
   Gating check (separate run, after sol[3] existed) confirmed: 30 distinct
   G[...,1-t] words needed for match-equations and 22 for the planar step,
   both 0 missing against the frozen H_opt/input/fibrtto1.m. *)

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
outckpt = FileNameJoin[{workdir, "sol_w4.mx"}];

If[FileExistsQ[outckpt],
  Print["Checkpoint ", outckpt, " already exists -- nothing to do. Delete it to force recompute."];
  Exit[0];
];
If[!FileExistsQ[ckpt], Print["Missing ", ckpt, " -- run run_precompute.wls first."]; Exit[1]];
Do[
  f = FileNameJoin[{workdir, "sol3_" <> ToString[j] <> ".m"}];
  If[!FileExistsQ[f], Print["Missing ", f, " -- run run_weight3.wls first."]; Exit[1]],
  {j, {1, 2, 3, 4, 5, 6, x1, y1, xmy}}
];

MemoryConstrained[

  Print["=== loading precompute checkpoint ==="];
  Get[ckpt];
  tto1Red = Import[FileNameJoin[{root, "H_opt", "input", "fibrtto1.m"}]];
  tto1RedD = Dispatch[tto1Red];

  Logt0 = {G[0, t] -> 0, G[0, 0, t] -> 0, G[0, 0, 0, t] -> 0, G[0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, 0, t] -> 0};
  Log00 = {G[0, 0] -> 0, G[0, 0, 0] -> 0, G[0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0] -> 0};

  integrateSeg[Mmap_, prevSol_, var_] := Total[
    Function[d, GIntegrate[#, var] & /@ ((Mmap[d] . prevSol)/d)] /@ Keys[Mmap]
  ];

  (* Structural fix (measured: 53s -> ~0.4s on match01 at weight 4; this
     scales WORSE at higher weight since impart's k-sum literally grows by
     one term each weight, so apply it everywhere a vector containing
     repeated G-words gets a substitution that can trigger non-trivial
     per-occurrence work (an endpoint evaluation like /.t->1 that invokes
     $HPLAutoConvertToKnownFunctions's MZV reduction, or a /.tto1Red lookup).
     Plain ReplaceAll re-derives that reduction once per OCCURRENCE across
     the 371-dim vector, not once per DISTINCT word -- GetGs finds the
     distinct words (31 of them for impart3 at weight 4, each trivial to
     reduce alone), Dispatch turns the per-word lookup into O(1) hashing
     instead of a linear rule scan, and substituting through that table
     costs the reduction exactly once per distinct word. *)
  (* No trailing "/.rule" after "/.disp": measured that it's pure, costly
     (11.5s) dead weight. GetGs[expr] finds every distinct G-word including
     ones whose argument is just the bare symbol "t" itself (e.g. the
     tempSol[j] additive CNST pieces have no G/t at all, and every actual
     t-dependent piece is INSIDE some G[...]), so substituting via Dispatch
     over those words already produces the exact correct result on its own
     (verified against the direct "/.rule" result). *)
  fastSubst[expr_, rule_] := Module[{gs, conv, disp},
    gs = GetGs[expr];
    conv = (# -> (# /. rule)) & /@ gs;
    disp = Dispatch[conv];
    expr /. disp
  ];

  Print["=== loading sol[0],sol[1],sol[2],sol[3] ==="];
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
  Do[sol[3][j] = Import[FileNameJoin[{workdir, "sol3_" <> ToString[j] <> ".m"}]], {j, 1, 6}];
  sol[3][x1] = Import[FileNameJoin[{workdir, "sol3_x1.m"}]];
  sol[3][y1] = Import[FileNameJoin[{workdir, "sol3_y1.m"}]];
  sol[3][xmy] = Import[FileNameJoin[{workdir, "sol3_xmy.m"}]];

  Print["=== building IntPart[j] via precomputed M_a + GIntegrate ==="];
  {tI1, IntPart1} = AbsoluteTiming[integrateSeg[Mdt[1], sol[3][1], t]];
  {tI2, IntPart2} = AbsoluteTiming[integrateSeg[Mdt[2], sol[3][2], t]];
  {tI3, IntPart3} = AbsoluteTiming[integrateSeg[Mdt[3], sol[3][3], t]];
  {tI4, IntPart4} = AbsoluteTiming[integrateSeg[Mdt[4], sol[3][4], t]];
  {tI5, IntPart5} = AbsoluteTiming[integrateSeg[Mdt[5], sol[3][5], t]];
  {tI6, IntPart6} = AbsoluteTiming[integrateSeg[Mdt[6], sol[3][6], t]];
  {tIx1, IntPartx1} = AbsoluteTiming[integrateSeg[Mdyx1, sol[3][x1], y]];
  {tIy1, IntParty1} = AbsoluteTiming[integrateSeg[Mdxy1, sol[3][y1], x]];
  {tIxmy, IntPartxmy} = AbsoluteTiming[integrateSeg[Mdxmy, sol[3][xmy], y]];
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
  tempSolx1 = Collect[IntPartx1 + CNSTx1, {_G}, FullSimplify[# /. rootsXmy] &];
  tempSoly1 = Collect[IntParty1 + CNSTy1, {_G}, FullSimplify[# /. rootsXmy] &];
  tempSolxmy = Collect[IntPartxmy + CNSTxmy, {_G}, FullSimplify[# /. rootsXmy] &];

  Print["=== matching equations (reflection + tto1Red + ShuffleG + Logt0) ==="];
  reflect[expr_] := fastSubst[(fastSubst[expr, t -> 1 - t] /. tto1RedD // ShuffleG) /. Logt0, t -> 0];
  {tm2, match2} = AbsoluteTiming[reflect[tempSol1] - (tempSol2 /. Logt0 /. t -> 0)];
  {tm3, match3} = AbsoluteTiming[reflect[tempSol2] - (tempSol3 /. Logt0 /. t -> 0)];
  {tm4, match4} = AbsoluteTiming[reflect[tempSol3] - (tempSol4 /. Logt0 /. t -> 0)];
  {tm5, match5} = AbsoluteTiming[reflect[tempSol4] - (tempSol5 /. Logt0 /. t -> 0)];
  {tm6, match6} = AbsoluteTiming[reflect[tempSol5] - (tempSol6 /. Logt0 /. t -> 0)];
  Print["match2..6 timings: ", {tm2, tm3, tm4, tm5, tm6}];

  Print["=== weight-4 impart (4 terms, README 3.5) ==="];
  {tImp5, impart5} = AbsoluteTiming[tempSol5 + Sum[(-I Pi)^k/k! todeltaPow[5, k] . sol[4 - k][5], {k, 1, 4}]];
  {tImp3, impart3} = AbsoluteTiming[tempSol3 + Sum[(-I Pi)^k/k! todeltaPow[3, k] . sol[4 - k][3], {k, 1, 4}]];
  Print["impart5, impart3 build: ", tImp5, " ", tImp3, "s"];
  {tM00, match00} = AbsoluteTiming[reflect[impart5] - (tempSolxmy /. y -> 0 /. Log00)];
  {tM01, match01} = AbsoluteTiming[(fastSubst[impart3, t -> 1] /. Log00 /. HPL[a__] -> 0) - (tempSoly1 /. x -> 0 /. Log00)];
  {tM10, match10} = AbsoluteTiming[(fastSubst[impart3, t -> 0] /. Log00) - (tempSolx1 /. y -> 0 /. Log00)];
  Print["match00, match01, match10 build: ", tM00, " ", tM01, " ", tM10, "s"];

  Print["=== extract shift Delta's via zero-substitution and chain to segment 5 ==="];
  {tD, {d2, d3, d4, d5, d6, d00, d01, d10}} = AbsoluteTiming[{
    match2 /. zeroAll, match3 /. zeroAll, match4 /. zeroAll, match5 /. zeroAll,
    match6 /. zeroAll, match00 /. zeroAll, match01 /. zeroAll, match10 /. zeroAll
  }];
  Print["zero-substitution extraction (d2..d10): ", tD, "s"];

  Print["=== weight-4 planar constants on segment 5 (tto1Red directly, matches reference) ==="];
  CNSTsol = Table[cnst[4, i][5], {i, 1, Length[Basis]}];
  {tTL, templim} = AbsoluteTiming[Coefficient[Extract[SolutionWithPlanar, onlyplanar], eps, 4] /. x -> 0 /. Log00 /. y -> 1 - t];
  {tPW, planarW4} = AbsoluteTiming[(templim /. tto1RedD // ShuffleG) /. _G -> 0];
  Print["templim build: ", tTL, "s; planarW4 (tto1Red+ShuffleG): ", tPW, "s"];
  plrules5 = Thread[Rule[Extract[CNSTsol, onlyplanar], planarW4]];
  frozenSol5 = Import[FileNameJoin[{Hdir, "numerical_eq", "solution_4.m"}]];
  constraints5 = Join[frozenSol5, plrules5];
  {tC5, CNST5val} = AbsoluteTiming[CNSTsol /. constraints5];
  Print["CNST5val substitution: ", tC5, "s"];
  If[!FreeQ[CNST5val, cnst],
    Print["WARNING: CNST5val still has unresolved cnst[...] entries -- constraints incomplete."];
  ];

  {tChain, {CNST4val, CNST3val, CNST2val, CNST1val, CNST6val, CNSTxmyval, CNSTy1val, CNSTx1val}} = AbsoluteTiming[
    Module[{c4, c3, c2, c1, c6, cxm, cy1, cx1},
      c4 = CNST5val - d5; c3 = c4 - d4; c2 = c3 - d3; c1 = c2 - d2;
      c6 = CNST5val + d6; cxm = CNST5val + d00; cy1 = c3 + d01; cx1 = c3 + d10;
      {c4, c3, c2, c1, c6, cxm, cy1, cx1}
    ]
  ];
  Print["chain to segment 5 (CNST1..6,x1,y1,xmy val): ", tChain, "s"];

  Print["=== assembling sol[4][seg] ==="];
  {tAsm, dummyAsm} = AbsoluteTiming[
    solOut[1] = tempSol1 /. Thread[Rule[CNST1, CNST1val]];
    solOut[2] = tempSol2 /. Thread[Rule[CNST2, CNST2val]];
    solOut[3] = tempSol3 /. Thread[Rule[CNST3, CNST3val]];
    solOut[4] = tempSol4 /. Thread[Rule[CNST4, CNST4val]];
    solOut[5] = tempSol5 /. Thread[Rule[CNST5, CNST5val]];
    solOut[6] = tempSol6 /. Thread[Rule[CNST6, CNST6val]];
    solOut[xmy] = tempSolxmy /. Thread[Rule[CNSTxmy, CNSTxmyval]];
    solOut[x1] = tempSolx1 /. Thread[Rule[CNSTx1, CNSTx1val]];
    solOut[y1] = tempSoly1 /. Thread[Rule[CNSTy1, CNSTy1val]];
  ];
  Print["assembly (substituting CNSTval into tempSol): ", tAsm, "s"];

  {tExp, dummyExp} = AbsoluteTiming[
    Do[Export[FileNameJoin[{workdir, "sol4_" <> ToString[j] <> ".m"}], solOut[j]], {j, 1, 6}];
    Export[FileNameJoin[{workdir, "sol4_xmy.m"}], solOut[xmy]];
    Export[FileNameJoin[{workdir, "sol4_x1.m"}], solOut[x1]];
    Export[FileNameJoin[{workdir, "sol4_y1.m"}], solOut[y1]];
  ];
  Print["export: ", tExp, "s"];

  DumpSave[outckpt, {solOut}];
  Print["Done. Exported sol4_{1..6,x1,y1,xmy}.m to ", workdir];

, 8*1024^3]
