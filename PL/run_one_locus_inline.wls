$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

nm = $ScriptCommandLine[[2]];
Print["=== locus ", nm, " (inline, no Module) ==="];

atilde = (Import["input/atilde.m"] /. z -> 1 - x - y)/ep;

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
     If[tr =!= pred, ok = False],
     {n, 1, 4}];
   If[! ok, Return[{$Failed, "trace mismatch", tally}]];
   tally
   ];

{tS, spec} = AbsoluteTiming[getSpectrum[locusMx]];
Print["spectrum=", spec, "  t=", tS, "s"];

posEig = Select[spec[[All, 1]], # > 0 &];
lamUse = Min[posEig];
Print["lambda=", lamUse];

(* === Horner computation INLINED at top level, no Module, no function === *)
A = locusMx;
lam = lamUse;
aLam = Cases[spec, {lam, m_} :> m][[1]];
others = DeleteCases[spec, {lam, _}];
n = Length[A];
{tG, gOfA} = AbsoluteTiming[
   If[others === {}, IdentityMatrix[n], Dot @@ (MatrixPower[A - #[[1]] IdentityMatrix[n], #[[2]]] & /@ others)]];
Print["gOfA done, t=", tG];

gScalar[zz_] := If[others === {}, 1, Times @@ ((zz - #[[1]])^#[[2]] & /@ others)];
{tSer, coeffs} = AbsoluteTiming[CoefficientList[Series[1/gScalar[lam + w], {w, 0, aLam - 1}] // Normal, w]];
Print["coeffs done, t=", tSer, " num=", Length[coeffs]];

B = A - lam IdentityMatrix[n];
{tH, qOfA} = AbsoluteTiming[
   If[aLam === 1, coeffs[[1]] IdentityMatrix[n],
    Fold[(#1 . B + coeffs[[#2]] IdentityMatrix[n]) &, coeffs[[aLam]] IdentityMatrix[n], Reverse[Range[1, aLam - 1]]]]];
Print["qOfA (Horner) done, t=", tH];

{tP, P} = AbsoluteTiming[qOfA . gOfA];
Print["P done, t=", tP, " Tr[P]=", Tr[P], " expected=", aLam];

idemPos = Position[P.P - P, _?(! PossibleZeroQ[#] &), {2}, Heads -> False];
Print["idempotency bad entries: ", Length[idemPos]];

Export["work/result_" <> nm <> ".m", {"OK", lamUse, aLam, P}];
Print["=== DONE, exported work/result_", nm, ".m ==="];
