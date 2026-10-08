$HistoryLength = 0;

SetDirectory[DirectoryName[$InputFileName]];

Print["=== step 1: load atilde, substitute z -> 1-x-y, strip the ep prefactor ==="];
{t0, atilde} = AbsoluteTiming[(Import["input/atilde.m"] /. z -> 1 - x - y)/ep];
Print["atilde dims: ", Dimensions[atilde], "  t=", t0, "s"];
Print["atilde[[1,1]] after /ep: ", atilde[[1,1]], "  (should have no ep left)"];

letters = {x, y, 1 - x, 1 - y, x + y, 1 - x - y};
Print["letters: ", letters];
strayVars = Complement[Variables[atilde /. _Log -> 0], {}];
Print["stray non-Log variables remaining (should be {}): ", strayVars];
If[strayVars =!= {}, Print["!!! STRAY VARIABLES PRESENT, stopping."]; Abort[]];

Print["=== step 2: extract Axy[letter] via Coefficient, verify reconstruction ==="];
{t1, axyList} = AbsoluteTiming[
   Table[l -> Coefficient[atilde, Log[l]], {l, letters}]
   ];
Print["Axy extraction t=", t1, "s"];
Do[Axy[axyList[[i, 1]]] = axyList[[i, 2]], {i, 1, Length[axyList]}];

recon = atilde - Sum[Axy[l] Log[l], {l, letters}];
reconPos = Position[recon, _?(! PossibleZeroQ[#] &), {2}, Heads -> False];
Print["Reconstruction residual nonzero entries: ", Length[reconPos],
  " (should be 0)"];
If[Length[reconPos] > 0,
  Print["!!! RECONSTRUCTION FAILED, stopping. Sample: ",
    Extract[recon, Take[reconPos, Min[3, Length[reconPos]]]]];
  Abort[]];

DumpSave["work/axy.mx", {Axy, letters}];
Print["checkpoint: work/axy.mx written"];

Print["=== step 3: single-letter residue matrices (6 loci) + combined x,y->0 corner ==="];
locus["x0"] = Axy[x];
locus["y0"] = Axy[y];
locus["x1"] = Axy[1 - x];
locus["y1"] = Axy[1 - y];
locus["xmy"] = Axy[x + y];
locus["xy1"] = Axy[1 - x - y];
locus["xy0"] = Axy[x] + Axy[y] + Axy[x + y];

locusNames = {"x0"};

Print["=== step 4: spectra (numeric eigenvalues, rounded, trace-verified) ==="];
getSpectrum[A_] := Module[{numEig, rounded, tally, ok = True, n, tr, pred},
   numEig = Eigenvalues[N[A, 30]];
   rounded = Round /@ numEig;
   If[Max[Abs[numEig - rounded]] > 10^-10,
    Return[{$Failed, "non-integer eigenvalues", numEig}]];
   tally = Tally[rounded];
   Do[
     tr = Tr[MatrixPower[A, n]];
     pred = Sum[t[[2]] t[[1]]^n, {t, tally}];
     If[tr =!= pred,
      ok = False;
      Print["  trace mismatch at n=", n, ": Tr(A^n)=", tr, " predicted=", pred]],
     {n, 1, 4}];
   If[!ok, Return[{$Failed, "trace mismatch", tally}]];
   tally
   ];

Do[
  {tS, spec} = AbsoluteTiming[getSpectrum[locus[nm]]];
  Print["locus ", nm, ": spectrum=", spec, "  t=", tS, "s"];
  specOf[nm] = spec;
  , {nm, locusNames}];

DumpSave["work/spectra.mx", {locus, specOf, locusNames}];
Print["checkpoint: work/spectra.mx written"];

Print["=== step 5: projector method (lambda=1 at each locus, Jordan-free, Horner's rule) ==="];
projLambda[A_, lam_, specArg_] := Module[
   {others, gScalar, gOfA, w, series, coeffs, B, qOfA, P, aLam, tr, idem, idemPos, n},
   n = Length[A];
   aLam = Cases[specArg, {lam, m_} :> m][[1]];
   others = DeleteCases[specArg, {lam, _}];
   gOfA = If[others === {}, IdentityMatrix[n],
     Dot @@ (MatrixPower[A - #[[1]] IdentityMatrix[n], #[[2]]] & /@ others)];
   (* qLam(z) = Taylor series of 1/g(z) about z=lam, truncated at order aLam-1 *)
   gScalar[zz_] := If[others === {}, 1, Times @@ ((zz - #[[1]])^#[[2]] & /@ others)];
   series = Series[1/gScalar[lam + w], {w, 0, aLam - 1}] // Normal;
   coeffs = CoefficientList[series, w];
   (* qOfA via Horner's rule: incremental, never materializes a standalone high power of B *)
   B = A - lam IdentityMatrix[n];
   qOfA = If[aLam === 1, coeffs[[1]] IdentityMatrix[n],
     Fold[(#1 . B + coeffs[[#2]] IdentityMatrix[n]) &,
      coeffs[[aLam]] IdentityMatrix[n], Reverse[Range[1, aLam - 1]]]];
   P = qOfA . gOfA;
   tr = Tr[P];
   idem = P.P - P;
   idemPos = Position[idem, _?(! PossibleZeroQ[#] &), {2}, Heads -> False];
   {P, tr === aLam, Length[idemPos] === 0, aLam}
   ];

Do[
  spec = specOf[nm];
  If[Head[spec] === List && FreeQ[spec, $Failed],
   posEig = Select[spec[[All, 1]], # > 0 &];
   If[posEig =!= {},
    lamUse = Min[posEig];
    {tP, res} = AbsoluteTiming[projLambda[locus[nm], lamUse, spec]];
    {P1[nm], traceOK[nm], idemOK[nm], mult1[nm]} = res;
    lamUsed[nm] = lamUse;
    Print["locus ", nm, " lambda=", lamUse, " (smallest positive eigenvalue): mult=", mult1[nm],
     " traceOK=", traceOK[nm], " idemOK=", idemOK[nm], "  t=", tP, "s"];
    If[! traceOK[nm] || ! idemOK[nm],
     Print["!!! PROJECTOR VERIFICATION FAILED at locus ", nm, " -- flagging, not stopping."];
     P1[nm] = $Failed;
     ];
    ,
    Print["locus ", nm, ": NO positive eigenvalue present at all! spectrum=", spec];
    P1[nm] = $Failed;
    ],
   Print["locus ", nm, ": spectrum computation failed: ", spec];
   P1[nm] = $Failed;
   ];
  , {nm, locusNames}];

DumpSave["work/projectors.mx", {P1, traceOK, idemOK, mult1}];
Print["checkpoint: work/projectors.mx written"];

Print["=== step 6: ranks ==="];
sixLoci = {"x0", "y0", "x1", "y1", "xmy", "xy1"};
If[AllTrue[sixLoci, P1[#] =!= $Failed &],
  stack6 = Join @@ (P1 /@ sixLoci);
  {tR6, rank6} = AbsoluteTiming[MatrixRank[stack6]];
  Print["Rank of 6 single-letter lambda=1 projectors stacked: ", rank6, "  (291-rank=", 291 - rank6, ")  t=", tR6, "s"];
  ,
  Print["Cannot compute 6-locus rank, some projector(s) failed/missing."];
  ];

allSeven = Append[sixLoci, "xy0"];
If[AllTrue[allSeven, P1[#] =!= $Failed &],
  stack7 = Join @@ (P1 /@ allSeven);
  {tR7, rank7} = AbsoluteTiming[MatrixRank[stack7]];
  Print["Rank including combined x,y->0 corner: ", rank7, "  (291-rank=", 291 - rank7, ")  t=", tR7, "s"];
  ,
  Print["Cannot compute 7-locus rank, some projector(s) failed/missing."];
  ];

Print["=== DONE ==="];
