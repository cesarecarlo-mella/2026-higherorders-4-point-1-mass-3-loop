(* H_opt run_weight.wls <w>
   Builds sol[w][seg] for seg in {1..6,x1,y1,xmy} using the precomputed
   M_a residue matrices (H_opt/work/precompute.mx) instead of
   Apart[dt[j].sol[w-1][j],t] + GIntegrate (README section 3.1), and
   plain vector subtraction instead of Solve for the matching/transport
   (README section 3.3). Reuses the already-solved numeric constants
   H/numerical_eq/solution_{w}.m (frozen input -- the reference notebook
   itself does exactly this, it does not re-derive them inside the
   notebook either) rather than re-solving the 238x238 linear system,
   since the goal here is to validate the integration/matching
   machinery reproduces the reference exactly.
   Usage: wolframscript -file run_weight.wls <w>   (w = 1 or 2 only) *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
(* PolyLogTools defines its own x,y,t,delta,c2,c2c,cxmy,cxmyc,sqrt,eps symbols and
   prepends its context ahead of Global` on $ContextPath; without this fix every
   bare occurrence of these letters typed below would silently bind to the
   PolyLogTools` symbol instead of the Global` one baked into precompute.mx,
   producing nonsense like "y/Global`y" that doesn't cancel. *)
(* Confirmed colliding names (Names["PolyLogTools`*"] intersected with everything
   this script references): x,y,t,delta,c2,Basis,r1,r2,H. Force all of them --
   plus every other symbol this script or the checkpoint uses -- into Global`
   before any bare reference, and keep Global` ahead of PolyLogTools` on the
   search path so that priority holds for the rest of the script. *)
Global`x; Global`y; Global`t; Global`delta; Global`c2; Global`c2c; Global`cxmy; Global`cxmyc; Global`sqrt; Global`eps;
Global`Basis; Global`r1; Global`r2; Global`r3; Global`r4; Global`H;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];

wgt = ToExpression[$ScriptCommandLine[[2]]];
If[wgt != 1 && wgt != 2,
  Print["This script is scoped to w=1,2 only (see H_opt/README.md section 5 validation order). Got w=", wgt];
  Exit[1];
];

root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
Hdir = FileNameJoin[{root, "H"}];
workdir = FileNameJoin[{root, "H_opt", "work"}];
ckpt = FileNameJoin[{workdir, "precompute.mx"}];
outckpt = FileNameJoin[{workdir, "sol_w" <> ToString[wgt] <> ".mx"}];

If[FileExistsQ[outckpt],
  Print["Checkpoint ", outckpt, " already exists -- nothing to do. Delete it to force recompute."];
  Exit[0];
];
If[!FileExistsQ[ckpt], Print["Missing ", ckpt, " -- run run_precompute.wls first."]; Exit[1]];

MemoryConstrained[

  Print["=== loading precompute checkpoint ==="];
  Get[ckpt];

  Logt0 = {G[0, t] -> 0, G[0, 0, t] -> 0, G[0, 0, 0, t] -> 0, G[0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, t] -> 0, G[0, 0, 0, 0, 0, 0, t] -> 0};
  Log00 = {G[0, 0] -> 0, G[0, 0, 0] -> 0, G[0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0] -> 0, G[0, 0, 0, 0, 0, 0, 0] -> 0};

  (* prepend-and-sum: replaces GIntegrate[#,t]&/@Apart[dt[j].prevSol,t] *)
  integrateSeg[Mmap_, prevSol_, var_] := Total[
    Function[d, GIntegrate[#, var] & /@ ((Mmap[d] . prevSol)/d)] /@ Keys[Mmap]
  ];

  Print["=== loading sol[", wgt - 1, "] ==="];
  If[wgt == 1,
    sol0 = Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w0", "sol0.m"}]];
    Do[sol[0][j] = sol0, {j, 1, 6}];
    sol[0][x1] = sol0; sol[0][y1] = sol0; sol[0][xmy] = sol0;
    ,
    Do[sol[1][j] = Import[FileNameJoin[{workdir, "sol1_" <> ToString[j] <> ".m"}]], {j, 1, 6}];
    sol[1][x1] = Import[FileNameJoin[{workdir, "sol1_x1.m"}]];
    sol[1][y1] = Import[FileNameJoin[{workdir, "sol1_y1.m"}]];
    sol[1][xmy] = Import[FileNameJoin[{workdir, "sol1_xmy.m"}]];
  ];

  Print["=== building IntPart[j] via precomputed M_a + GIntegrate (no Apart on dt.sol) ==="];
  {tI1, IntPart1} = AbsoluteTiming[integrateSeg[Mdt[1], sol[wgt - 1][1], t]];
  {tI2, IntPart2} = AbsoluteTiming[integrateSeg[Mdt[2], sol[wgt - 1][2], t]];
  {tI3, IntPart3} = AbsoluteTiming[integrateSeg[Mdt[3], sol[wgt - 1][3], t]];
  {tI4, IntPart4} = AbsoluteTiming[integrateSeg[Mdt[4], sol[wgt - 1][4], t]];
  {tI5, IntPart5} = AbsoluteTiming[integrateSeg[Mdt[5], sol[wgt - 1][5], t]];
  {tI6, IntPart6} = AbsoluteTiming[integrateSeg[Mdt[6], sol[wgt - 1][6], t]];
  {tIx1, IntPartx1} = AbsoluteTiming[integrateSeg[Mdyx1, sol[wgt - 1][x1], y]];
  {tIy1, IntParty1} = AbsoluteTiming[integrateSeg[Mdxy1, sol[wgt - 1][y1], x]];
  {tIxmy, IntPartxmy} = AbsoluteTiming[integrateSeg[Mdxmy, sol[wgt - 1][xmy], y]];
  Print["IntPart timings: ", {tI1, tI2, tI3, tI4, tI5, tI6, tIx1, tIy1, tIxmy}];

  (* each segment's unknown gets its own symbol head so zero-substitution is unambiguous *)
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

  If[wgt != 1,
    Print["wgt 2's impart/match00/01/10 use a different (todelta^2) formula -- not yet implemented here, see run_weight2.wls once w=1 validates."];
    Exit[1];
  ];
  impart5 = -I Pi todelta[5];
  match00 = ((impart5 . sol[0][5] + tempSol5) /. t -> 0 /. Log00) - (tempSolxmy /. y -> 0 /. Log00);
  impart3 = -I Pi todelta[3];
  match01 = ((impart3 . sol[0][3] + tempSol3) /. t -> 1 /. Log00 /. HPL[a__] -> 0) - (tempSoly1 /. x -> 0 /. Log00);
  match10 = ((impart3 . sol[0][3] + tempSol3) /. t -> 0 /. Log00 /. HPL[a__] -> 0) - (tempSolx1 /. y -> 0 /. Log00);

  Print["=== extract shift Delta's via zero-substitution (README 3.3) and chain to segment 5 ==="];
  d2 = match2 /. zeroAll; (* CNST1 - CNST2 + d2 = 0  =>  CNST1 = CNST2 - d2 *)
  d3 = match3 /. zeroAll;
  d4 = match4 /. zeroAll;
  d5 = match5 /. zeroAll;
  d6 = match6 /. zeroAll;
  d00 = match00 /. zeroAll; (* CNST5-part + d00 - CNSTxmy = 0 => CNSTxmy = CNST5 + d00 *)
  d01 = match01 /. zeroAll;
  d10 = match10 /. zeroAll;

  Print["=== known planar + frozen numeric constants on segment 5 ==="];
  CNSTsol = Table[cnst[wgt, i][5], {i, 1, Length[Basis]}];
  plrules5 = Thread[Rule[Extract[CNSTsol, onlyplanar], Coefficient[Extract[SolutionWithPlanar, onlyplanar], eps, wgt]]] /. _G -> 0;
  frozenSol5 = Import[FileNameJoin[{Hdir, "numerical_eq", "solution_" <> ToString[wgt] <> ".m"}]];
  constraints5 = Join[frozenSol5, plrules5];
  CNST5val = CNSTsol /. constraints5;
  If[!FreeQ[CNST5val, cnst],
    Print["WARNING: CNST5val still has unresolved cnst[...] entries -- constraints incomplete."];
  ];

  (* chain: CNST4 = CNST5 - d5 ; CNST3 = CNST4 - d4 = CNST5 - d5 - d4 ; CNST2 = CNST3-d3 ; CNST1 = CNST2-d2 ; CNST6 = CNST5 - d6 *)
  CNST4val = CNST5val - d5;
  CNST3val = CNST4val - d4;
  CNST2val = CNST3val - d3;
  CNST1val = CNST2val - d2;
  (* sign fix: CNST6 is the MINUS-position variable in match6 = tempSol5|t1
     - tempSol6|t0, so solving match6=0 for it gives CNST5+d6, not CNST5-d6
     (unlike CNST1..4, each the PLUS/first-term variable in its own match
     equation). Masked at weight 1 since length-1 G-words fully vanish at
     both endpoints there; confirmed wrong by a weight-2 row-7 bisection. *)
  CNST6val = CNST5val + d6;
  CNSTxmyval = CNST5val + d00;
  (* match01: (impart3.sol0_3 + tempSol3)|t1 - tempSoly1|x0 = 0 => CNST3part + d01 - CNSTy1 = 0 => CNSTy1 = CNST3 + d01 *)
  CNSTy1val = CNST3val + d01;
  CNSTx1val = CNST3val + d10;

  Print["=== assembling sol[", wgt, "][seg] ==="];
  solOut[1] = tempSol1 /. Thread[Rule[CNST1, CNST1val]];
  solOut[2] = tempSol2 /. Thread[Rule[CNST2, CNST2val]];
  solOut[3] = tempSol3 /. Thread[Rule[CNST3, CNST3val]];
  solOut[4] = tempSol4 /. Thread[Rule[CNST4, CNST4val]];
  solOut[5] = tempSol5 /. Thread[Rule[CNST5, CNST5val]];
  solOut[6] = tempSol6 /. Thread[Rule[CNST6, CNST6val]];
  solOut[xmy] = tempSolxmy /. Thread[Rule[CNSTxmy, CNSTxmyval]];
  solOut[x1] = tempSolx1 /. Thread[Rule[CNSTx1, CNSTx1val]];
  solOut[y1] = tempSoly1 /. Thread[Rule[CNSTy1, CNSTy1val]];

  Do[Export[FileNameJoin[{workdir, "sol" <> ToString[wgt] <> "_" <> ToString[j] <> ".m"}], solOut[j]], {j, 1, 6}];
  Export[FileNameJoin[{workdir, "sol" <> ToString[wgt] <> "_xmy.m"}], solOut[xmy]];
  Export[FileNameJoin[{workdir, "sol" <> ToString[wgt] <> "_x1.m"}], solOut[x1]];
  Export[FileNameJoin[{workdir, "sol" <> ToString[wgt] <> "_y1.m"}], solOut[y1]];

  DumpSave[outckpt, {solOut}];
  Print["Done. Exported sol", wgt, "_{1..6,x1,y1,xmy}.m to ", workdir];

, 8*1024^3]
