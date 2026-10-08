$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

atilde = (Import["input/atilde.m"] /. z -> 1 - x - y)/ep;
AxyY = Coefficient[atilde, Log[y]];

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

Print["=== y0 ==="];
{tS, specY0} = AbsoluteTiming[getSpectrum[AxyY]];
Print["spectrum: ", specY0, "  t=", tS];
posY0 = Select[specY0[[All, 1]], # > 0 &];
Print["positive eigenvalues: ", posY0];
allPy0 = {};
Do[
  {tP, res} = AbsoluteTiming[projLambda[AxyY, lam, specY0]];
  {P, tOK, iOK, mult} = res;
  Print["  lambda=", lam, " mult=", mult, " traceOK=", tOK, " idemOK=", iOK, " t=", tP];
  AppendTo[allPy0, P];
  , {lam, posY0}];
Export["work/allP_y0.m", allPy0];
Print["y0 done, exported work/allP_y0.m"];
