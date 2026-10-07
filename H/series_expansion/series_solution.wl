(* ::Package::

   Analytic series-expansion solution of the regularized canonical system,
   along the line y=y0, weight-by-weight in eps AND order-by-order in x
   around the regular singular point x=0. Companion to
   numerical_evaluation/integrate_canonical.wl (which solves the same
   regularized system numerically via NDSolve); this package instead
   Taylor/Puiseux-expands around x=0 and integrates term-by-term
   analytically, exploiting Integrate[x^n,x]=x^(n+1)/(n+1).

   METHOD (see numerical_evaluation/notes/numerical_evaluation_writeup.pdf
   for the derivation of the regularization itself):
     - B(x,y0) = dx.m(x,y0) - Ax/x  is expanded as a Puiseux series in x
       (half-integer powers appear from the algebraic letter
       r=Sqrt[x y0(1-x-y0)], which has its branch point at x=0).
     - R_j(x) = Sum_{p=0}^j (-1)^p Q_p . B(x,y0) . Q_{j-p}  is built as a
       PER-POWER NUMERIC matrix (fast: plain matrix multiplication, no
       symbolic overhead) -- this is the dominant cost, ~O(maxWeight^2)
       matrix triple-products.
     - The regularized solution J'^(m)(x) = Sum_l Log[x]^l S_{m,l}(x)
       satisfies the exact recursive hierarchy (same as the numeric script):
         dJ'^(0)/dx = 0
         dJ'^(m)/dx = Sum_{k=0}^{m-1} C_{m-1-k}(x) . J'^(k)(x)
       Each S_{m,l}(x) is obtained by a NUMERIC-matrix . (x-polynomial)-vector
       contraction (cheap: O(n^2) scalar-times-expression operations, not the
       O(n^3) dense-symbolic-matrix blowup a naive "multiply everything
       symbolically" approach would hit), followed by Expand[] to keep the
       polynomial compact, then trivial term-by-term integration.
     - IMPORTANT: the final reassembly J^(m)(x) = Sum_k Log[x]^(m-k) Q_{m-k}.J'^(k)(x)
       MUST be Expand[]-ed after construction -- the raw Sum/Dot introduces a
       large number of not-yet-combined terms (e.g. the weight-6 result was
       ~592MB unexpanded vs ~3.4MB after Expand, a 173x difference) even
       though the underlying content is compact.

   HOW TO USE: edit the PARAMETERS block, then evaluate. Results are saved,
   per weight, to series_expansion/output/{Jreg,J_full}_series_symbolic_w<m>.m
   (already Expand[]-ed) and series_expansion/output/J_series_w<m>_at_x1.m
   (numeric value at the target point x1).
*)

(* ======================= PARAMETERS (edit these) ======================= *)

baseDir = "/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/NPL/DEQ_solutions/SLC/Solutions/H";

y0 = 3/15;           (* fixed y along the path *)
x1 = 1/15;            (* target point x1 -- edit freely *)
maxWeight = 6;         (* highest weight (eps order) wanted, 0 <= maxWeight <= 6 *)
seriesOrder = 15;      (* x-power truncation; radius of convergence is 4/25=0.16,
                           order 15 gives ample precision for x1 well inside that *)

(* ========================================================================= *)

If[maxWeight < 0 || maxWeight > 6, Print["ERROR: maxWeight must be between 0 and 6."]; Abort[]];

SetDirectory[baseDir];
If[! DirectoryQ["series_expansion/output"], CreateDirectory["series_expansion/output"]];

Print["Loading exact data (projectors, nilpotents, residue matrix, boundary vectors) ..."];
spectrum = {-4, -3, -2, -1, 0, 1, 2, 3};
Npow = <||>;
Do[Npow[lam] = Get["matrix_to_JD/toJD/Npowers_" <> ToString[lam] <> "_toexp1.m"], {lam, spectrum}];
Ax = Get["matrix_to_JD/toexp1.m"];
n = Length[Ax]; id = IdentityMatrix[n];

(* Q_j: exact, x-independent matrices with M_j(x) = Log[x]^j Q_j -- see
   numerical_evaluation_writeup.pdf Sec. 5.1 for the derivation. *)
Qj = Table[
  Sum[If[lam == 0 && j - k == 0, 1, lam^(j - k)]/((j - k)! k!) *
     If[k < Length[Npow[lam]], Npow[lam][[k + 1]], 0*id],
    {lam, spectrum}, {k, 0, Min[2, j]}],
  {j, 0, maxWeight}];
Print["  Q0 == Identity: ", Qj[[1]] === id];
If[maxWeight >= 1, Print["  Q1 == Ax: ", Qj[[2]] === Ax]];

Jreg0 = Table[Get["numerical_evaluation/output/Jreg_boundary_w" <> ToString[w] <> "_num.m"], {w, 0, maxWeight}];

(* ----------------------------------------------------------------------- *)
(* B(x,y0) = dx.m(x,y0) - Ax/x, Puiseux-expanded around x=0. Closed-form
   Sqrt[] substitution is fine here (unlike the numeric ODE integrator):
   Series[] handles the branch point at x=0 exactly, no precision concerns. *)
(* ----------------------------------------------------------------------- *)

Print["Loading numerical_evaluation/input/dx.m, building B(x,y0) ..."];
dxMat = Get["numerical_evaluation/input/dx.m"];
dxAtY0 = dxMat /. y -> y0 /. {r1 -> sqrt*(x - 1)/y0, r2 -> sqrt*(x + y0)/y0} /.
   sqrt -> Sqrt[-x*y0*(-1 + x + y0)];
BbareAtY0 = dxAtY0 - Ax/x;

Print["Finding nonzero entries ..."];
nzPos = {};
Do[If[BbareAtY0[[i, j]] =!= 0, AppendTo[nzPos, {i, j}]], {i, 1, n}, {j, 1, n}];
Print["  nonzero entries: ", Length[nzPos]];

Print["Series-expanding each nonzero entry to order ", seriesOrder, " ..."];
t0 = AbsoluteTime[];
entrySeries = <||>;
Do[
  {i, j} = pos;
  s = Series[BbareAtY0[[i, j]], {x, 0, seriesOrder}];
  If[Head[s] === SeriesData,
   den = s[[6]]; nmin = s[[4]]; coeffs = s[[3]];
   pairs = Table[{(nmin + k - 1)/den, coeffs[[k]]}, {k, 1, Length[coeffs]}];
   entrySeries[pos] = Select[pairs, #[[2]] =!= 0 &];
   ];
  , {pos, nzPos}];
Print["  done, time=", AbsoluteTime[] - t0];

Print["Assembling Bseries[power] = numeric SparseArray (via GroupBy) ..."];
t0 = AbsoluteTime[];
flatTriples = Flatten[Table[{pos, #[[1]], #[[2]]} & /@ entrySeries[pos], {pos, nzPos}], 1];
groupedByPower = GroupBy[flatTriples, #[[2]] &, (#[[1]] -> #[[3]] &) /@ # &];
Bseries = <||>;
Do[Bseries[p] = SparseArray[groupedByPower[p], {n, n}], {p, Keys[groupedByPower]}];
Print["  done, time=", AbsoluteTime[] - t0, "   distinct x-powers: ", Length[Keys[Bseries]]];

Print["Computing R_j(x), j=0..", maxWeight - 1, " (per-power numeric matrices) ..."];
t0 = AbsoluteTime[];
Rseries = Table[
   assoc = <||>;
   Do[assoc[power] = Sum[(-1)^p*Qj[[p + 1]].Bseries[power].Qj[[j - p + 1]], {p, 0, j}],
     {power, Keys[Bseries]}];
   assoc
   , {j, 0, maxWeight - 1}];
Print["  done, time=", AbsoluteTime[] - t0];

(* ----------------------------------------------------------------------- *)
(* Weight-by-weight recursive series integration.
   S_{k,l}(x): 371-vector, an actual x-polynomial (numeric coefficients),
   for the coefficient of Log[x]^l in J'^(k)(x). Each recursion step is a
   NUMERIC-matrix . (numeric-or-polynomial)-vector contraction (cheap),
   Expand[]-ed to stay compact, then integrated via Integrate[.,x] (trivial
   power rule on an already-expanded polynomial). *)
(* ----------------------------------------------------------------------- *)

Print["\nRunning weight-by-weight recursion ..."];
Srep = <||>;
Srep[0] = <|0 -> Jreg0[[1]]|>;

Do[
  t0 = AbsoluteTime[];
  Tvec = <||>;
  Do[
    acc = ConstantArray[0, n];
    Do[
      lprime = l - (m - 1 - k);
      If[lprime >= 0 && KeyExistsQ[Srep[k], lprime],
       RjAssoc = Rseries[[(m - 1 - k) + 1]];
       Skl = Srep[k][lprime];
       Do[acc = acc + x^rpow*(RjAssoc[rpow].Skl), {rpow, Keys[RjAssoc]}];
       ];
      , {k, 0, m - 1}];
    Tvec[l] = Expand[acc];
    , {l, 0, m - 1}];

  Svec = <||>;
  Do[
    corrected = Tvec[l];
    If[KeyExistsQ[Svec, l + 1], corrected = Expand[corrected - (l + 1)*(Svec[l + 1]/x)]];
    integrated = Integrate[corrected, x];
    If[l == 0, integrated = integrated + Jreg0[[m + 1]]];
    Svec[l] = Expand[integrated];
    , {l, m - 1, 0, -1}];

  Srep[m] = Svec;
  Print["  weight ", m, " done, time=", AbsoluteTime[] - t0, "   log-levels: ", Keys[Svec]];
  , {m, 1, maxWeight}];

(* ----------------------------------------------------------------------- *)
(* Assemble the symbolic-in-x results (regularized J' and full physical J),
   Expand[]-ing explicitly -- this is the step that must not be skipped:
   the raw Sum/Dot is enormously bloated with un-combined terms even though
   the true content is compact (weight 6: ~592MB raw vs ~3.4MB expanded). *)
(* ----------------------------------------------------------------------- *)

Print["\nAssembling and expanding symbolic-in-x results ..."];
t0 = AbsoluteTime[];
JregSymbolic = Table[
   Expand[Sum[Log[x]^l*Srep[m][l], {l, Keys[Srep[m]]}]],
   {m, 0, maxWeight}];
JfullSymbolic = Table[
   Expand[Sum[Log[x]^(m - k)*(Qj[[m - k + 1]].JregSymbolic[[k + 1]]), {k, 0, m}]],
   {m, 0, maxWeight}];
Print["  done, time=", AbsoluteTime[] - t0];

Print["Saving results (already expanded/compact) ..."];
Do[
  Put[JregSymbolic[[m + 1]], "series_expansion/output/Jreg_series_symbolic_w" <> ToString[m] <> ".m"];
  Put[JfullSymbolic[[m + 1]], "series_expansion/output/J_full_series_symbolic_w" <> ToString[m] <> ".m"];
  , {m, 0, maxWeight}];

Print["Evaluating at x1=", x1, " and saving numeric result ..."];
JAtX1 = Table[JfullSymbolic[[m + 1]] /. x -> N[x1], {m, 0, maxWeight}];
Do[
  Put[JAtX1[[m + 1]], "series_expansion/output/J_series_w" <> ToString[m] <> "_at_x1.m"];
  , {m, 0, maxWeight}];

Print["\n=== DONE ==="];
Print["J(x1,y0) Taylor coefficients in eps computed up to weight ", maxWeight, "."];
Print["Result stored in JAtX1 = {J^(0), J^(1), ..., J^(", maxWeight, ")}, each a numeric length-", n, " vector."];
Print["Example -- component 3 at each weight:"];
Do[Print["  eps^", m, " : ", JAtX1[[m + 1]][[3]]], {m, 0, maxWeight}];
