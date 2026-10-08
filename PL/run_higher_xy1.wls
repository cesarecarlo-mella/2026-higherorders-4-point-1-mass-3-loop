$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

atilde = (Import["input/atilde.m"] /. z -> 1 - x - y)/ep;
AxyXY1 = Coefficient[atilde, Log[1 - x - y]];

getSpectrum[A_] := Module[{numEig, rounded, tally, ok = True, n, tr, pred},
   numEig = Eigenvalues[N[A, 30]];
   rounded = Round /@ numEig;
   If[Max[Abs[numEig - rounded]] > 10^-10, Return[{$Failed, "non-integer", numEig}]];
   tally = Tally[rounded];
   Do[tr = Tr[MatrixPower[A, n]]; pred = Sum[t[[2]] t[[1]]^n, {t, tally}];
     If[tr =!= pred, ok = False], {n, 1, 4}];
   If[! ok, Return[{$Failed, "trace mismatch", tally}]];
   tally];

projLambda[A_, lam_, specArg_] := Module[
   {others, gScalar, gOfA, w, series, coeffs, B, qOfA, P, aLam, idemPos, n},
   n = Length[A];
   aLam = Cases[specArg, {lam, m_} :> m][[1]];
   others = DeleteCases[specArg, {lam, _}];
   gOfA = If[others === {}, IdentityMatrix[n],
     Dot @@ (MatrixPower[A - #[[1]] IdentityMatrix[n], #[[2]]] & /@ others)];
   gScalar[zz_] := If[others === {}, 1, Times @@ ((zz - #[[1]])^#[[2]] & /@ others)];
   series = Series[1/gScalar[lam + w], {w, 0, aLam - 1}] // Normal;
   coeffs = CoefficientList[series, w];
   B = A - lam IdentityMatrix[n];
   qOfA = If[aLam === 1, coeffs[[1]] IdentityMatrix[n],
     Fold[(#1 . B + coeffs[[#2]] IdentityMatrix[n]) &, coeffs[[aLam]] IdentityMatrix[n], Reverse[Range[1, aLam - 1]]]];
   P = qOfA . gOfA;
   idemPos = Position[P.P - P, _?(! PossibleZeroQ[#] &), {2}, Heads -> False];
   {P, Tr[P] === aLam, Length[idemPos] === 0, aLam}
   ];

Print["=== xy1 ==="];
{tS, specXY1} = AbsoluteTiming[getSpectrum[AxyXY1]];
Print["spectrum: ", specXY1, "  t=", tS];
posXY1 = Select[specXY1[[All, 1]], # > 0 &];
Print["positive eigenvalues: ", posXY1];
allPxy1 = {};
Do[
  {tP, res} = AbsoluteTiming[projLambda[AxyXY1, lam, specXY1]];
  {P, tOK, iOK, mult} = res;
  Print["  lambda=", lam, " mult=", mult, " traceOK=", tOK, " idemOK=", iOK, " t=", tP];
  AppendTo[allPxy1, P];
  , {lam, posXY1}];
Export["work/allP_xy1.m", allPxy1];
Print["xy1 done, exported work/allP_xy1.m"];
