(* ::Package::

   Numerical evaluation of the 371x371 canonical system along the line y = y0,
   from the regular singular point x=0 to a target point x1, weight-by-weight
   in eps (order-by-order transcendentality), following the strategy of
   Czakon & Tancredi (arXiv:2606.30354) adapted to this system as derived in
   regularized_vector.tex and numerical_procedure.tex.

   HOW TO USE: paste this entire file into a Mathematica notebook cell (or
   Get/Import it), edit the three lines in the PARAMETERS block below, then
   evaluate. The final answer is in the variable JAtX1 (see bottom): a list
   {J^(0), J^(1), ..., J^(maxWeight)} of numeric 371-vectors, the Taylor
   coefficients in eps of J(x1,y0).
*)

(* ======================= PARAMETERS (edit these) ======================= *)

baseDir = "/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/NPL/DEQ_solutions/SLC/Solutions/H";

y0 = 3/15;              (* fixed y along the integration path *)
x1 = 1/15;               (* TARGET POINT x1 -- edit freely *)
maxWeight = 6;            (* highest weight (eps order) wanted, 0 <= maxWeight <= 6 -- edit freely *)

xStart = 10^-9;           (* regulator distance from the x=0 singularity: J' is
                              analytic at x=0, so J'(xStart) = J'(0) to excellent
                              approximation; xStart avoids Log[0] in the numerics *)
workingPrecisionND = MachinePrecision;  (* switch to e.g. 30 for higher precision (slower) *)
precisionGoalND = 8;

useSqrtVariable = True;   (* OPTION: integrate in u = Sqrt[x] instead of x directly.
                              B(x,y0) contains genuine 1/Sqrt[x]-type terms (from the
                              algebraic letter r(x,y0)=Sqrt[x y0(1-x-y0)] ~ Sqrt[x]
                              having its branch point exactly at x=0) that the residue
                              subtraction never removes -- only the plain 1/x pole was
                              handled. Substituting x=u^2 (dx=2u du) cancels this
                              1/Sqrt[x]=1/u divergence exactly via the Jacobian, giving
                              a genuinely smooth RHS for the ODE solver. Set to False to
                              integrate directly in x (the original, less robust method,
                              kept for comparison). *)

useAuxiliaryRootODE = True; (* OPTION (Czakon-Tancredi strategy): enlarge the system by
                              +1 dimension, treating r(x,y0)=Sqrt[x y0(1-x-y0)] as its
                              own dynamical variable satisfying dr/dx = r(1-2x-y0)/(2x(1-x-y0)),
                              instead of substituting the closed-form Sqrt[...] expression
                              directly. "sqrt" is kept as a free symbol in B(x,y0) and
                              replaced by the current NUMERIC value of the dynamical r at
                              each RHS evaluation. Verified to reproduce the closed-form
                              result (agrees to ~1e-10 relative on r(x1) itself); kept as
                              an option for cross-checking. Set to False to substitute the
                              closed form directly instead (skips the extra dimension). *)

(* ========================================================================= *)

If[maxWeight < 0 || maxWeight > 6, Print["ERROR: maxWeight must be between 0 and 6."]; Abort[]];

SetDirectory[baseDir];
Print["Loading exact data (projectors, nilpotents, residue matrix, connection matrix) ..."];

spectrum = {-4, -3, -2, -1, 0, 1, 2, 3};
Npow = <||>;
Do[
  Npow[lam] = Get["matrix_to_JD/toJD/Npowers_" <> ToString[lam] <> "_toexp1.m"];
  , {lam, spectrum}];
Ax = Get["matrix_to_JD/toexp1.m"];
n = Length[Ax];
id = IdentityMatrix[n];
Print["  n = ", n, ", spectrum = ", spectrum];

(* ----------------------------------------------------------------------- *)
(* Q_j: exact, x-INDEPENDENT matrices such that the eps^j Taylor coefficient
   of M(x,eps) = exp(eps Ax Log x) is simply  M_j(x) = (Log x)^j * Q_j.
   Derivation: numerical_procedure.tex, Section 3.1.  Q_0 = Identity (sum of
   projectors), Q_1 = Ax (Jordan-Chevalley S+N), verified below. *)
(* ----------------------------------------------------------------------- *)

(* Need Q_0,...,Q_maxWeight: the ODE hierarchy uses Q_0..Q_{maxWeight-1},
   and the final reconstruction M(x1,eps).J'(x1,eps) additionally needs
   Q_maxWeight (the top eps-order of the product). *)
Print["Building Q_j (j=0.." <> ToString[maxWeight] <> ") -- exact rational matrices ..."];
Qj = Table[
  Sum[
    If[lam == 0 && j - k == 0, 1, lam^(j - k)] / ((j - k)! k!) *
      If[k < Length[Npow[lam]], Npow[lam][[k + 1]], 0*id],
    {lam, spectrum}, {k, 0, Min[2, j]}
  ],
  {j, 0, maxWeight}
];
Print["  Q0 == Identity: ", Qj[[1]] === id];
If[maxWeight >= 1, Print["  Q1 == Ax: ", Qj[[2]] === Ax]];

MjExact[p_] := Qj[[p + 1]]; (* matrix such that M_p(x) = Log[x]^p * MjExact[p] *)

(* ----------------------------------------------------------------------- *)
(* B(x,y0) = dx.m(x,y0) - Ax/x : the regular (pole-free) part of the x-derivative
   of the connection one-form along y=y0. Evaluated by plain numeric substitution
   (no Compile) -- simpler, and works uniformly at any WorkingPrecision. *)
(* ----------------------------------------------------------------------- *)

Print["Loading numerical_evaluation/input/dx.m and specializing to y = ", y0, " ..."];

dxMat = Get["numerical_evaluation/input/dx.m"];
dxAtY0 = dxMat /. y -> y0 /. {r1 -> sqrt*(x - 1)/y0, r2 -> sqrt*(x + y0)/y0};
If[!useAuxiliaryRootODE,
  (* fallback: substitute the closed form directly, no extra ODE dimension *)
  dxAtY0 = dxAtY0 /. sqrt -> Sqrt[-x*y0*(-1 + x + y0)];
];

BbareAtY0 = dxAtY0 - Ax/x; (* symbolic in x, and (if useAuxiliaryRootODE) also in sqrt *)

precND = If[workingPrecisionND === MachinePrecision, MachinePrecision, workingPrecisionND];

(* Bfunc always takes (x, r): when useAuxiliaryRootODE is False, "sqrt" no longer
   appears in BbareAtY0, so r is simply ignored (the substitution is a no-op). *)
Bfunc[xIn_, rIn_] := N[BbareAtY0 /. {x -> xIn, sqrt -> rIn}, precND];

(* dr/dx for the auxiliary root ODE: r(x,y0)=Sqrt[x y0(1-x-y0)] satisfies
   dr = r/(2g) dg with g=x y0(1-x-y0), i.e. dr/dx = r(1-2x-y0)/(2x(1-x-y0)). *)
drdx[xIn_?NumericQ, rIn_?NumericQ] := rIn*(1 - 2*xIn - y0)/(2*xIn*(1 - xIn - y0));
rClosedForm[xIn_] := Sqrt[xIn*y0*(1 - xIn - y0)];

(* ----------------------------------------------------------------------- *)
(* Boundary vector J'(0,y0), weight by weight (already numerically evaluated
   from solutions/eps_sep/w<k>/sol<k>_5.m at t=1-y0). *)
(* ----------------------------------------------------------------------- *)

Print["Loading regularized boundary vector J'(0,y0), weights 0.." <> ToString[maxWeight] <> " ..."];
Jreg0 = Table[
  Get["numerical_evaluation/output/Jreg_boundary_w" <> ToString[w] <> "_num.m"],
  {w, 0, maxWeight}
]; (* Jreg0[[k+1]] = J'^{(k)}(0,y0), k=0..maxWeight *)

(* ----------------------------------------------------------------------- *)
(* RHS of the hierarchical ODE, Eq. (hierarchy) in numerical_procedure.tex.
   State vector Y = stack( J'^{(1)}(x), ..., J'^{(maxWeight)}(x) ), length n*maxWeight.
   (J'^{(0)} is the known constant Jreg0[[1]] and is NOT part of the dynamical state.)
*)
(* ----------------------------------------------------------------------- *)

(* C_j(x) = M(x,eps)^{-1} B(x,y0) M(x,eps) truncated to eps^j.
   Using M_p(x) = Log[x]^p Q_p (Q_p exact, x-independent -- see above) and
   M(x,eps)^{-1} = M(x,-eps), the L-powers factor out completely:
     C_j(x) = Log[x]^j * Sum_{p=0}^{j} (-1)^p  Q_p . B(x,y0) . Q_{j-p}
   so only ONE matrix (B(x,y0)) actually depends on x; the rest is a fixed,
   precomputable contraction of exact matrices against Bx.

   dJ'/dx = Sum_{k=0}^{m-1} C_{m-1-k}(x) . J'^{(k)}(x)  --  this is the plain
   x-derivative, coreDx[xVal,Yvec]. When integrating in u=Sqrt[x] instead
   (useSqrtVariable=True), dJ'/du = 2u * dJ'/dx|_{x=u^2} by the chain rule;
   the explicit 1/Sqrt[x]=1/u pieces hiding inside B(x,y0) then cancel exactly
   against this Jacobian factor 2u, instead of being evaluated as a genuine
   (numerically delicate) divergence. *)
coreDx[xVal_?NumericQ, rVal_?NumericQ, Yvec_?(VectorQ[#, NumericQ] &)] := Module[
  {L, Bx, QB, Cj, blocks, out, m, j, p, q, k},
  If[maxWeight == 0, Return[{}]];
  L = Log[xVal];
  Bx = Bfunc[xVal, rVal];
  QB = Table[Qj[[p + 1]] . Bx, {p, 0, maxWeight - 1}]; (* Q_p . B(x) *)
  Cj = Table[
    L^j * Sum[(-1)^p * QB[[p + 1]] . Qj[[j - p + 1]], {p, 0, j}],
    {j, 0, maxWeight - 1}
    ];
  blocks = Partition[Yvec, n];
  out = Table[
    Sum[Cj[[m - 1 - k + 1]] . If[k == 0, Jreg0[[1]], blocks[[k]]], {k, 0, m - 1}],
    {m, 1, maxWeight}
    ];
  Flatten[out]
  ];

(* The dynamical state is stack(J'^(1),...,J'^(maxWeight)) and, when
   useAuxiliaryRootODE, one extra trailing component r(x). *)
rhsX[xVal_?NumericQ, YvecExt_?(VectorQ[#, NumericQ] &)] := Module[{Yvec, rVal, dY, dr},
  If[useAuxiliaryRootODE,
   Yvec = Most[YvecExt]; rVal = Last[YvecExt];
   dY = coreDx[xVal, rVal, Yvec];
   dr = drdx[xVal, rVal];
   Append[dY, dr]
   ,
   coreDx[xVal, 0., YvecExt]
   ]
  ];

rhsU[uVal_?NumericQ, YvecExt_?(VectorQ[#, NumericQ] &)] := Module[{Yvec, rVal, xVal, dY, dr},
  xVal = uVal^2;
  If[useAuxiliaryRootODE,
   Yvec = Most[YvecExt]; rVal = Last[YvecExt];
   dY = 2*uVal*coreDx[xVal, rVal, Yvec];
   dr = 2*uVal*drdx[xVal, rVal];
   Append[dY, dr]
   ,
   2*uVal*coreDx[xVal, 0., YvecExt]
   ]
  ];

Print["\n*** Ready. Running the integration now. ***"];

(* ----------------------------------------------------------------------- *)
(* Integrate from xStart to x1. State vector Y[x] stacks J'^{(1)}(x),...,
   J'^{(maxWeight)}(x). Initial condition: J'^{(k)}(xStart) ~= J'^{(k)}(0,y0),
   the known boundary vector (valid to excellent approximation since J' is
   analytic/continuous at x=0 and xStart is tiny). *)
(* ----------------------------------------------------------------------- *)

rFinalNumeric = Null; (* populated below if useAuxiliaryRootODE *)

JregAtX1 = If[maxWeight == 0,
   {},
   Module[{Y0, sol, Yend, uStart, u1, rStart},
    Y0 = Flatten[Table[Jreg0[[k + 1]], {k, 1, maxWeight}]];
    If[useAuxiliaryRootODE,
     rStart = N[rClosedForm[xStart], precND]; (* closed form used only to SEED the ODE at xStart *)
     Y0 = Append[Y0, rStart];
    ];
    t0 = AbsoluteTime[];
    If[useSqrtVariable,
     uStart = Sqrt[xStart]; u1 = Sqrt[x1];
     Print["Integrating ODE (state dimension ", Length[Y0], ") in u=Sqrt[x], from u=",
      N[uStart], " to u=", N[u1],
      If[useAuxiliaryRootODE, "  [with auxiliary root ODE, +1 dim]", ""], " ..."];
     sol = NDSolve[
       {Y'[u] == rhsU[u, Y[u]], Y[uStart] == Y0},
       Y, {u, uStart, u1},
       WorkingPrecision -> workingPrecisionND,
       PrecisionGoal -> precisionGoalND,
       MaxSteps -> 100000
       ];
     Yend = Y[u1] /. First[sol];
     ,
     Print["Integrating ODE (state dimension ", Length[Y0], ") from x=", N[xStart],
      " to x=", N[x1],
      If[useAuxiliaryRootODE, "  [with auxiliary root ODE, +1 dim]", ""], " ..."];
     sol = NDSolve[
       {Y'[x] == rhsX[x, Y[x]], Y[xStart] == Y0},
       Y, {x, xStart, x1},
       WorkingPrecision -> workingPrecisionND,
       PrecisionGoal -> precisionGoalND,
       MaxSteps -> 100000
       ];
     Yend = Y[x1] /. First[sol];
     ];
    Print["  integration time: ", AbsoluteTime[] - t0];
    If[useAuxiliaryRootODE,
     rFinalNumeric = Last[Yend];
     Print["  r(x1) from auxiliary ODE: ", rFinalNumeric,
      "   (closed form: ", N[rClosedForm[x1], precND], ")"];
     Partition[Most[Yend], n]
     ,
     Partition[Yend, n]
    ]
    ]
   ]; (* JregAtX1[[k]] = J'^{(k)}(x1,y0), k=1..maxWeight *)

(* ----------------------------------------------------------------------- *)
(* Recover the physical vector: J(x1,y0) = M(x1,eps) . J'(x1,eps), truncated
   to eps^maxWeight.   M_{p}(x1) = Log[x1]^p * Q_p. *)
(* ----------------------------------------------------------------------- *)

JregFull = Prepend[JregAtX1, Jreg0[[1]]]; (* JregFull[[k+1]] = J'^{(k)}(x1), k=0..maxWeight *)
Lx1 = N[Log[x1], If[workingPrecisionND === MachinePrecision, MachinePrecision, workingPrecisionND]];

JAtX1 = Table[
   Sum[(Lx1^(m - k) * Qj[[m - k + 1]]) . JregFull[[k + 1]], {k, 0, m}],
   {m, 0, maxWeight}
   ]; (* JAtX1[[m+1]] = J^{(m)}(x1,y0), i.e. the coefficient of eps^m *)

Print["\n=== DONE ==="];
Print["J(x1,y0) Taylor coefficients in eps computed up to weight ", maxWeight, "."];
Print["Result stored in JAtX1 = {J^(0), J^(1), ..., J^(", maxWeight, ")}, each a numeric length-", n, " vector."];
Print["Example -- component 1 at each weight:"];
Do[Print["  eps^", m, " : ", JAtX1[[m + 1]][[1]]], {m, 0, maxWeight}];

