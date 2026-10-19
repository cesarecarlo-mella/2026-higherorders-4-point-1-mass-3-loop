(* ============================================================================
   spectral_projection.wls  --  stage 2 of the PL pipeline (boundary/regularity
   projector matrices, one locus at a time)
   Clean consolidation of run_higher_eigenvalues.wls (x0), run_higher_y0.wls
   (y0), run_higher_xy1.wls (xy1), and the generic run_one_locus.wls -- all
   four were near-identical copies of the same two routines (getSpectrum,
   projLambda) differing only in which letter's connection matrix they used
   and whether they kept ALL positive-eigenvalue projectors (x0/y0/xy1, used
   as stacked boundary-condition matrices) or just the smallest one
   (x1/y1/xmy/xy1/xy0, used as single regularity matrices).
   Usage: wolframscript -file spectral_projection.wls <locus> <all|smallest>
   locus in {x0,y0,x1,y1,xmy,xy1,xy0}
   See README.md for what each output file is used for downstream.
   ============================================================================ *)

$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

nm = $ScriptCommandLine[[2]];
mode = If[Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]], "all"];

(* Only extract the ONE Axy matrix this locus needs. Extracting extra,
   unrelated Axy matrices in the same kernel session was empirically observed
   to make the later Horner/projector computation hang -- keep this minimal. *)
getLocusMatrix[nm_] := Module[{atilde},
   atilde = (Import["../input/atilde.m"] /. z -> 1 - x - y)/ep;
   Switch[nm,
     "x0", Coefficient[atilde, Log[x]],
     "y0", Coefficient[atilde, Log[y]],
     "x1", Coefficient[atilde, Log[1 - x]],
     "y1", Coefficient[atilde, Log[1 - y]],
     "xmy", Coefficient[atilde, Log[x + y]],
     "xy1", Coefficient[atilde, Log[1 - x - y]],
     "xy0", Coefficient[atilde, Log[x]] + Coefficient[atilde, Log[y]] + Coefficient[atilde, Log[x + y]]]
   ];

(* Spectrum of a connection matrix must be integer-valued (eigenvalues are
   the matrix's resonance/weight-shift structure); cross-checked via trace
   identities Tr[A^n] = Sum[mult*lambda^n] for n=1..4 as a cheap sanity test
   independent of the eigenvalue solver itself. *)
getSpectrum[A_] := Module[{numEig, rounded, tally, ok = True, tr, pred},
   numEig = Eigenvalues[N[A, 30]];
   rounded = Round /@ numEig;
   If[Max[Abs[numEig - rounded]] > 10^-10, Return[{$Failed, "non-integer", numEig}]];
   tally = Tally[rounded];
   Do[
     tr = Tr[MatrixPower[A, n]];
     pred = Sum[t[[2]] t[[1]]^n, {t, tally}];
     If[tr =!= pred, ok = False],
     {n, 1, 4}];
   If[! ok, Return[{$Failed, "trace mismatch", tally}]];
   tally
   ];

(* Spectral projector onto eigenvalue lam, built via the standard
   Lagrange/Horner construction P = q(A).g(A) where g(A) is the product over
   OTHER eigenvalues' generalized-eigenspace factors and q is the Taylor
   inverse of g needed to normalize P to be idempotent on lam's eigenspace.
   Verified via Tr[P]==multiplicity and P.P-P==0. *)
projLambda[A_, lam_, spec_] := Module[
   {others, gScalar, gOfA, w, series, coeffs, B, qOfA, P, aLam, idemPos, n},
   n = Length[A];
   aLam = Cases[spec, {lam, m_} :> m][[1]];
   others = DeleteCases[spec, {lam, _}];
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

(* ============================ run ============================ *)
locusMx = getLocusMatrix[nm];
{tS, spec} = AbsoluteTiming[getSpectrum[locusMx]];
Print[nm, ": spectrum=", spec, "  t=", tS, "s"];
If[Head[spec] =!= List || ! FreeQ[spec, $Failed], Print["!!! spectrum computation failed."]; Exit[1]];

posEig = Select[spec[[All, 1]], # > 0 &];
If[posEig === {}, Print["!!! no positive eigenvalue present."]; Exit[1]];

Which[
  mode == "all",
  (* x0 / y0 / xy1 pattern: keep ALL positive-eigenvalue projectors, stacked
     later (Join@@...) into the full boundary-condition matrix. *)
  allP = {};
  Do[
    {tP, res} = AbsoluteTiming[projLambda[locusMx, lam, spec]];
    {P, tOK, iOK, mult} = res;
    Print["  lambda=", lam, " mult=", mult, " traceOK=", tOK, " idemOK=", iOK, " t=", tP];
    AppendTo[allP, P];
    , {lam, posEig}];
  Export["../work/allP_" <> nm <> ".m", allP];
  Print[nm, " done, exported work/allP_", nm, ".m (", Length[allP], " projector(s))"],

  True,
  (* x1 / y1 / xmy pattern: keep only the smallest positive eigenvalue's
     projector, as a single bare matrix (the regularity-condition matrix). *)
  lamUse = Min[posEig];
  {tP, res} = AbsoluteTiming[projLambda[locusMx, lamUse, spec]];
  {P, tOK, iOK, mult} = res;
  Print["  lambda=", lamUse, " mult=", mult, " traceOK=", tOK, " idemOK=", iOK, " t=", tP];
  If[! tOK || ! iOK, Print["!!! PROJECTOR VERIFICATION FAILED."]; Exit[1]];
  regFileName = Switch[nm, "x1", "regX1.m", "y1", "regY1.m", "xmy", "regXmY.m", _, "reg" <> nm <> ".m"];
  Export["../work/" <> regFileName, P];
  Print[nm, " done, exported work/", regFileName];
  ];
