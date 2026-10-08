(* Compute spectrum + lambda=1(or smallest positive) projector for ONE locus, in its own fresh kernel.
   Usage: wolframscript -file run_one_locus.wls <locusName>
   locusName in {x0,y0,x1,y1,xmy,xy1,xy0} *)
$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

nm = $ScriptCommandLine[[2]];
Print["=== locus ", nm, " ==="];

atilde = (Import["input/atilde.m"] /. z -> 1 - x - y)/ep;

(* Only extract the single letter(s) this specific locus actually needs -- extracting
   extra unrelated Axy matrices in the same session was empirically observed to make the
   later Horner computation hang (unexplained Mathematica quirk), so keep this minimal. *)
locusMx = Switch[nm,
   "x0", Coefficient[atilde, Log[x]],
   "y0", Coefficient[atilde, Log[y]],
   "x1", Coefficient[atilde, Log[1 - x]],
   "y1", Coefficient[atilde, Log[1 - y]],
   "xmy", Coefficient[atilde, Log[x + y]],
   "xy1", Coefficient[atilde, Log[1 - x - y]],
   "xy0", Coefficient[atilde, Log[x]] + Coefficient[atilde, Log[y]] + Coefficient[atilde, Log[x + y]]];

getSpectrum[A_] := Module[{numEig, rounded, tally, ok = True, n, tr, pred},
   numEig = Eigenvalues[N[A, 30]];
   rounded = Round /@ numEig;
   If[Max[Abs[numEig - rounded]] > 10^-10, Return[{$Failed, "non-integer eigenvalues", numEig}]];
   tally = Tally[rounded];
   Do[
     tr = Tr[MatrixPower[A, n]];
     pred = Sum[t[[2]] t[[1]]^n, {t, tally}];
     If[tr =!= pred, ok = False;
      Print["  trace mismatch n=", n, ": ", tr, " vs ", pred]],
     {n, 1, 4}];
   If[! ok, Return[{$Failed, "trace mismatch", tally}]];
   tally
   ];

{tS, spec} = AbsoluteTiming[getSpectrum[locusMx]];
Print["spectrum=", spec, "  t=", tS, "s"];
If[Head[spec] =!= List || ! FreeQ[spec, $Failed],
  Print["!!! spectrum computation failed."]; Export["work/result_" <> nm <> ".m", {"FAILED_SPECTRUM", spec}]; Exit[1]];

posEig = Select[spec[[All, 1]], # > 0 &];
If[posEig === {},
  Print["!!! no positive eigenvalue present."]; Export["work/result_" <> nm <> ".m", {"NO_POSITIVE_EIGENVALUE", spec}]; Exit[1]];
lamUse = Min[posEig];

projLambda[A_, lam_, specArg_] := Module[
   {others, gScalar, gOfA, w, series, coeffs, B, qOfA, P, aLam, tr, idem, idemPos, n},
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
     Fold[(#1 . B + coeffs[[#2]] IdentityMatrix[n]) &,
      coeffs[[aLam]] IdentityMatrix[n], Reverse[Range[1, aLam - 1]]]];
   P = qOfA . gOfA;
   tr = Tr[P];
   idem = P.P - P;
   idemPos = Position[idem, _?(! PossibleZeroQ[#] &), {2}, Heads -> False];
   {P, tr === aLam, Length[idemPos] === 0, aLam}
   ];

Print["lambda=", lamUse, " (smallest positive eigenvalue); computing projector..."];
{tP, res} = AbsoluteTiming[projLambda[locusMx, lamUse, spec]];
{P, traceOK, idemOK, mult} = res;
Print["traceOK=", traceOK, " idemOK=", idemOK, " mult=", mult, "  t=", tP, "s"];

If[! traceOK || ! idemOK,
  Print["!!! PROJECTOR VERIFICATION FAILED."];
  Export["work/result_" <> nm <> ".m", {"VERIFICATION_FAILED", lamUse, mult}];
  Exit[1]];

Export["work/result_" <> nm <> ".m", {"OK", lamUse, mult, P}];
Print["=== DONE, result exported to work/result_", nm, ".m ==="];
