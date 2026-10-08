(* H_opt run_precompute.wls
   Weight-independent data: Axy[i] constant matrices, dx/dy connection,
   segment dt[j] (j=1..6) and dxy1/dyx1/dxmy decomposed into residue
   matrices M_a (keyed by the exact denominator factor Apart produced,
   e.g. "t", "-1+t", "t-c2" -- usable directly as GIntegrate[.../denom,t]).
   Residues are extracted by substitution+Coefficient (fast, syntactic),
   NOT by Residue/Limit on the full 371x371 matrix (empirically ~30x+
   slower and the cause of an earlier runaway kernel).
   Checkpoints to H_opt/work/precompute.mx. Restartable: delete the
   checkpoint to force a recompute; otherwise this script just exits. *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
(* Load PolyLogTools here too, even though this script never calls GIntegrate/G
   itself: H/results/planarInsertedXY.m (SolutionWithPlanar) contains literal
   G[...] text. If this script imported it WITHOUT PolyLogTools loaded, "G"
   would get created fresh as Global`G, incompatible with the PolyLogTools`G
   that GIntegrate/ShuffleG produce downstream in run_weight.wls -- the two
   G-functions would silently fail to combine. Load it first so "G" (and
   every other colliding name below) resolves the same way everywhere. *)
Global`x; Global`y; Global`t; Global`delta; Global`c2; Global`c2c; Global`cxmy; Global`cxmyc; Global`sqrt; Global`eps;
Global`Basis; Global`r1; Global`r2; Global`r3; Global`r4; Global`H;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];

root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
Hdir = FileNameJoin[{root, "H"}];
workdir = FileNameJoin[{root, "H_opt", "work"}];
ckpt = FileNameJoin[{workdir, "precompute.mx"}];

If[FileExistsQ[ckpt],
  Print["Checkpoint already exists at ", ckpt, " -- nothing to do. Delete it to force recompute."];
  Exit[0];
];

fastResidues[M_, var_] := Module[{dens},
  dens = Union[Cases[M, Power[d_, -1] /; !FreeQ[d, var] :> d, Infinity]];
  Association[# -> Coefficient[M /. Power[#, -1] -> dummy, dummy] & /@ dens]
];

MemoryConstrained[

  Print["=== step 1: load frozen inputs ==="];
  {t1, Basis} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "input", "basis.m"}]]];
  Print["Basis: ", Length[Basis], " entries, ", t1, "s"];

  {t2, Atilde} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "input", "Atilde.m"}]] /. s -> x /. t -> y /. m2 -> 1];
  Print["Atilde: ", Dimensions[Atilde], ", ", t2, "s"];

  {t3, shiftedBasis} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "shiftedBasis", "shiftedBasis.m"}]]];
  {t3b, SolutionWithPlanar} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "results", "planarInsertedXY.m"}]]];
  {t3c, onlyplanar} = AbsoluteTiming[Position[Map[FreeQ[#, INT[b__, a__]] &, shiftedBasis], True]];
  Print["shiftedBasis/SolutionWithPlanar/onlyplanar: ", t3, " ", t3b, " ", t3c, "s; Length[onlyplanar]=", Length[onlyplanar]];

  {t4, boundaryMxy1} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "matrix_to_JD", "toJD", "P_1_toexp1.m"}]]];
  {t5, boundaryMxy2} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "matrix_to_JD", "toJD", "P_1_toexp2.m"}]]];
  {t6, boundaryMxy3} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "matrix_to_JD", "toJD", "P_2_toexp3.m"}]]];
  {t7, boundaryMxy7} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "matrix_to_JD", "toJD", "P_1_toexp7.m"}]]];
  {t8, regCondxy1} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "boundary_matrices", "breg1.m"}]]];
  {t9, regCondxy2} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "boundary_matrices", "breg2.m"}]]];
  {t10, regCondxy3} = AbsoluteTiming[Import[FileNameJoin[{Hdir, "boundary_matrices", "breg3.m"}]]];
  Print["boundary/reg matrices loaded: ", {t4, t5, t6, t7, t8, t9, t10}];

  Print["=== step 2: alphabet / Axy[i] (constant matrices) ==="];
  {t11, alphabet} = AbsoluteTiming[Cases[Atilde, _Log, Infinity] // Union];
  alphabetI = alphabet /. Sqrt[-x (1 - x - y) y] -> I Sqrt[x (1 - x - y) y];
  arglogs = alphabetI /. Log -> Identity;
  LLs = Table[LL[i], {i, 1, Length[alphabetI]}];
  alphabetNaming = Thread[Rule[alphabetI, LLs]];
  Print["Length[alphabet] = ", Length[alphabet], ", ", t11, "s"];

  {t12, AtildeNamed} = AbsoluteTiming[Atilde /. alphabetNaming];
  Do[Axy[i] = Coefficient[AtildeNamed, LLs[[i]]], {i, 1, Length[LLs]}];
  Do[
    If[Variables[Flatten[Axy[i]]] =!= {},
      Print["WARNING: Axy[", i, "] is NOT constant! vars = ", Variables[Flatten[Axy[i]]]]
    ],
    {i, 1, Length[LLs]}
  ];
  Print["Axy[1..", Length[LLs], "] built, ", t12, "s"];

  Print["=== step 3: dx, dy (x,y connection, with sqrt-letter care) ==="];
  sqrtL1 = Log[((x (1 - x - y) - I Sqrt[x (1 - x - y) y])/(x (1 - x - y) + I Sqrt[x (1 - x - y) y]))];
  sqrtL2 = Log[((x y - I Sqrt[x (1 - x - y) y])/(x y + I Sqrt[x (1 - x - y) y]))];
  DSrules = {sqrtL1 -> ((I y (-1 + 2 x + y))/(sqrt (-1 + x) (x + y))), sqrtL2 -> ((I y)/(sqrt - sqrt x))};
  DTrules = {sqrtL1 -> -((I x)/(sqrt (x + y))), sqrtL2 -> ((I x)/(sqrt - sqrt y))};
  sqrtpart = Coefficient[Atilde, sqrtL1] sqrtL1 + Coefficient[Atilde, sqrtL2] sqrtL2;

  {t13, dx} = AbsoluteTiming[D[Atilde /. {sqrtL1 -> 0, sqrtL2 -> 0}, x] + (sqrtpart /. DSrules)];
  {t14, dy} = AbsoluteTiming[D[Atilde /. {sqrtL1 -> 0, sqrtL2 -> 0}, y] + (sqrtpart /. DTrules)];

  rrs = {((y)/((-1 + x) sqrt)) -> 1/r1, ((y)/(sqrt (x + y))) -> 1/r2, ((x)/((-1 + y) sqrt)) -> 1/r3, ((x)/(sqrt (x + y))) -> 1/r4};
  revertRRS = Thread[Rule[rrs[[;; , 2]], rrs[[;; , 1]]]];

  {t15, dx} = AbsoluteTiming[Apart[dx, x] /. rrs];
  {t16, dy} = AbsoluteTiming[Apart[dy, y] /. rrs];
  Print["dx,dy: deriv ", t13, " ", t14, "s; Apart ", t15, " ", t16, "s"];

  Print["=== step 4: segment parametrizations, dt[j] for j=1..6, FAST residue extraction ==="];
  Clear[param, ser, dt, todelta, Mdt];
  param[1] = {x -> t, y -> \[Delta]};
  posz = Position[Series[arglogs /. param[1], {\[Delta], 0, 0}, Assumptions -> 0 < t < 1] // Normal, 0] // Flatten;
  todelta[1] = Sum[Axy[i], {i, posz}];
  ser1 = (Series[arglogs /. {x -> t, y -> \[Delta]}, {\[Delta], 0, 1}, Assumptions -> 0 < t < 1] // Normal) /. \[Delta] -> 0;
  posT = (Position[ser1, t] /. {{a_, b_, c_} -> {a}} // Flatten) // Union;
  {td1, dt[1]} = AbsoluteTiming[Sum[Axy[i] D[Log[ser1[[i]]], t], {i, posT}]];
  {tr1, Mdt[1]} = AbsoluteTiming[fastResidues[dt[1], t]];
  Print["seg1: build ", td1, "s, residues ", tr1, "s, denoms=", Keys[Mdt[1]]];

  roots2 = {c2 -> (1 - I)/2, c2c -> (1 + I)/2};
  Do[
    param[j] = Switch[j,
      2, {y -> \[Delta] (t)^2, x -> 1 - \[Delta] ((1 - t)^2 + t^2)},
      3, {x -> 1 - \[Delta] - t, y -> t},
      4, {x -> \[Delta] (1 - t)^2, y -> 1 - \[Delta] ((1 - t)^2 + t^2)},
      5, {x -> \[Delta], y -> 1 - t},
      6, {x -> \[Delta] (t)^2, y -> \[Delta] (1 - t)^2}
    ];
    {tsa, serj} = AbsoluteTiming[Series[arglogs /. param[j], {\[Delta], 0, 2}, Assumptions -> {\[Delta] > 0, 0 < t < 1}] // Normal];
    pos0 = Position[serj /. \[Delta] -> 0, 0] // Flatten;
    todelta[j] = Sum[Axy[i], {i, pos0}];
    {tsb, dlogvec} = AbsoluteTiming[((D[Log[#] & /@ serj, t] // Simplify) /. \[Delta] -> 0 // Simplify) /. 1 - 2 t + 2 (t)^(2) -> 2 (t - c2) (t - c2c)];
    {tsc, dt[j]} = AbsoluteTiming[Apart[dlogvec . Table[Axy[i], {i, 1, Length[arglogs]}], t]];
    {tsd, Mdt[j]} = AbsoluteTiming[fastResidues[dt[j], t]];
    Print["seg", j, ": series ", tsa, "s, dlogvec ", tsb, "s, dt ", tsc, "s, residues ", tsd, "s, denoms=", Keys[Mdt[j]]],
    {j, 2, 6}
  ];

  Print["=== step 4b: checkpoint A (segments done) ==="];
  DumpSave[FileNameJoin[{workdir, "precompute_A.mx"}], {Basis, Atilde, shiftedBasis, SolutionWithPlanar, onlyplanar,
    boundaryMxy1, boundaryMxy2, boundaryMxy3, boundaryMxy7, regCondxy1, regCondxy2, regCondxy3,
    alphabet, alphabetI, arglogs, LLs, alphabetNaming,
    dt, todelta, Mdt, dx, dy, rrs, revertRRS}];
  Print["checkpoint A written."];

  Print["=== step 5: dxy1, dyx1, dxmy and their residue decompositions ==="];
  yto1root = Thread[Rule[{1/r1, 1/r2}, {1/r1 /. revertRRS /. 1/sqrt -> 1/I/x /. y -> 1, 1/r2 /. revertRRS /. 1/sqrt -> 1/I/x /. y -> 1}]];
  xto1root = Thread[Rule[{1/r3, 1/r4}, {1/r3 /. revertRRS /. 1/sqrt -> 1/I/y /. x -> 1, 1/r4 /. revertRRS /. 1/sqrt -> 1/I/y /. x -> 1}]];
  xtomyroot = Thread[Rule[{1/r1, 1/r2, 1/r3, 1/r4}, {1/r1 /. revertRRS /. 1/sqrt -> -1/I/y /. x -> -y, (Series[(Series[1/r2 /. revertRRS /. sqrt -> Sqrt[x y (1 - x - y)], {x, -y, 1}] // Normal) /. 1/Sqrt[-y^2] -> -1/(I y) /. x -> -y (1 - delta), {delta, 0, 0}] // Normal), 1/r3 /. revertRRS /. 1/sqrt -> -1/I/y /. x -> -y, (Series[(Series[1/r4 /. revertRRS /. sqrt -> Sqrt[x y (1 - x - y)], {x, -y, 1}] // Normal) /. 1/Sqrt[-y^2] -> -1/(I y) /. x -> -y (1 - delta), {delta, 0, 0}] // Normal)}]];

  Off[Power::infy]; Off[Infinity::indet];

  {t17, dxy1} = AbsoluteTiming[Apart[(dx //. yto1root /. y -> 1) /. 1/(1 - x + x^2) -> 1/(x + cxmy)/(x + cxmyc), x]];
  {t18, dyx1} = AbsoluteTiming[Apart[(dy //. xto1root /. x -> 1) /. 1/(1 - y + y^2) -> 1/(y + cxmy)/(y + cxmyc), y]];
  tmpxmy = Apart[((dy - dx) /. xtomyroot /. x -> -y), y] /. delta -> 0;
  {t19, dxmy} = AbsoluteTiming[Apart[tmpxmy /. 1/(y^2 + y + 1) -> 1/(y - cxmy)/(y - cxmyc) /. 1/(1 - y + y^2) -> 1/(y + cxmy)/(y + cxmyc), y]];
  Print["dxy1,dyx1,dxmy build: ", t17, " ", t18, " ", t19, "s"];

  {t20, Mdxy1} = AbsoluteTiming[fastResidues[dxy1, x]];
  {t21, Mdyx1} = AbsoluteTiming[fastResidues[dyx1, y]];
  {t22, Mdxmy} = AbsoluteTiming[fastResidues[dxmy, y]];
  Print["residue extraction: ", t20, " ", t21, " ", t22, "s"];
  Print["dxy1 denoms: ", Keys[Mdxy1]];
  Print["dyx1 denoms: ", Keys[Mdyx1]];
  Print["dxmy denoms: ", Keys[Mdxmy]];

  Print["=== step 6b: todelta powers for segments 3,5 (weight-independent, README 3.5) ==="];
  (* impart at weight w needs todelta[j]^k for k=1..w-1 (only j=3,5 are ever
     used in the impart formulas, per the reference's match00/01/10). Cache
     up to k=6 once here so weights 3-6 never recompute them. *)
  todeltaPow[5, 1] = todelta[5];
  todeltaPow[3, 1] = todelta[3];
  Do[todeltaPow[5, k] = todeltaPow[5, k - 1] . todelta[5], {k, 2, 6}];
  Do[todeltaPow[3, k] = todeltaPow[3, k - 1] . todelta[3], {k, 2, 6}];
  Print["todelta powers built for j=3,5, k=1..6."];

  Print["=== step 6c: roots[xmy] (algebraic reduction target for x1/y1/xmy, README 3.6) ==="];
  rootsXmy = {cxmy -> -(-1)^(1/3), cxmyc -> (-1)^(2/3)};

  Print["=== step 7: checkpoint ==="];
  DumpSave[ckpt, {Basis, Atilde, shiftedBasis, SolutionWithPlanar, onlyplanar,
    boundaryMxy1, boundaryMxy2, boundaryMxy3, boundaryMxy7, regCondxy1, regCondxy2, regCondxy3,
    alphabet, alphabetI, arglogs, LLs, alphabetNaming,
    dt, todelta, Mdt, dx, dy, dxy1, dyx1, dxmy, Mdxy1, Mdyx1, Mdxmy,
    todeltaPow, rootsXmy}];
  Print["Checkpoint written to ", ckpt];

, 8*1024^3]
