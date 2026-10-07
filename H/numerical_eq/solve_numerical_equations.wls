(* Usage: wolframscript -file numerical_eq/solve_numerical_equations.wls N
   Solves numerical_equations_N.m (falling back to numerical_equations_N_exact.m
   if the plain file isn't found) and exports the result to solution_N.m.

   Method: instead of generic Solve[] (which chokes on symbol-laden elimination
   when coefficients carry Pi/Zeta[3] and can take 10+ minutes), extract the
   linear system's rational coefficient matrix via CoefficientArrays, confirm
   the unknowns (cnst[...]) never multiply a transcendental, split the
   right-hand side into its rational and transcendental parts (arbitrary
   powers of Pi, Zeta[3], etc. are auto-detected), and solve all parts in a
   single RowReduce pass. Verifies by direct substitution before exporting;
   refuses to export if verification fails. *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];

order = $ScriptCommandLine[[2]];

file = "numerical_eq/numerical_equations_" <> order <> ".m";
If[! FileExistsQ[file], file = "numerical_eq/numerical_equations_" <> order <> "_exact.m"];
Print["Reading ", file];
eqs = Import[file];
Print["Total entries: ", Length[eqs]];

vars = Union@Cases[eqs, cnst[_, _][_], Infinity];
If[vars === {}, vars = Union@Cases[eqs, cnst[_][_], Infinity]];
If[vars === {}, vars = Union@Cases[eqs, cnst[_], Infinity]];
Print["var count: ", Length[vars]];

real = Select[eqs, Head[#] === Equal &];
Print["Equal entries: ", Length[real]];

t0 = AbsoluteTime[];
{c0, c1} = CoefficientArrays[(#[[1]] - #[[2]] &) /@ real, vars];
c0 = Normal[c0]; c1 = Normal[c1];
Print["CoefficientArrays extraction time: ", AbsoluteTime[] - t0, "s"];
Print["matrix c1 purely rational (Pi/Zeta/I free)? ", FreeQ[c1, Pi] && FreeQ[c1, Zeta[_]] && FreeQ[c1, I]];
Print["dims: c1 = ", Dimensions[c1], "   c0 length = ", Length[c0]];

(* generic transcendental-basis detection: handles arbitrary powers of Pi, Zeta[3], etc, and combinations *)
gens = Union@Cases[c0, Pi | Zeta[_] | E | EulerGamma, Infinity];
Print["generators found: ", gens];

degs = Table[Max[Prepend[Cases[Exponent[#, g] & /@ c0, _Integer], 0]], {g, gens}];
Print["max degrees: ", Thread[gens -> degs]];

basisTuples = If[gens === {}, {{}}, Tuples[Range[0, #] & /@ degs]];
basis = If[gens === {}, {1}, (Times @@ MapThread[Power, {gens, #}]) & /@ basisTuples];
Print["basis (", Length[basis], " elements): ", basis];

extractCoeff[expr_, tuple_] := Module[{e = expr}, Do[e = Coefficient[e, gens[[i]], tuple[[i]]], {i, Length[gens]}]; e];
bs = Table[extractCoeff[#, tuple] & /@ c0, {tuple, basisTuples}];
recon = Total[MapThread[#1*#2 &, {bs, basis}]];
Print["decomposition check: max|c0 - recon| = ", Max[Abs[N[c0 - recon]]]];

t1 = AbsoluteTime[];
rr = RowReduce[Join[c1, Transpose[(-#) & /@ bs], 2]];
Print["RowReduce time: ", AbsoluteTime[] - t1, "s"];

nv = Length[vars];
Print["rank / pivots found: ", MatrixRank[c1], "  (need ", nv, ")"];

sols = Table[rr[[1 ;; nv, nv + k]], {k, Length[basis]}];
solvars = Total[MapThread[#1*#2 &, {sols, basis}]];
constraints = Thread[vars -> solvars];

checkVals = ((#[[1]] - #[[2]]) &) /@ real /. constraints // Simplify;
numBad = Count[checkVals, Except[0]];
Print["verification: nonzero residuals after substitution = ", numBad];

zeroCount = Count[solvars, 0];
Print["zero constants: ", zeroCount, " / ", nv];
Print["sample: ", Take[constraints, 5]];

outFile = "numerical_eq/solution_" <> order <> ".m";
If[numBad === 0,
  Export[outFile, constraints];
  Print["Exported to ", outFile];
  ,
  Print["NOT EXPORTED — verification failed with ", numBad, " bad residuals."];
]
