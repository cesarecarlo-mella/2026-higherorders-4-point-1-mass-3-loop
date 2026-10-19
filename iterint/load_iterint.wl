(* Loader for run_iterint.wls output (<dir>/<fam>_<cross>_w<k>.wdx = {keys, S}).
   IIWeight[dir, fam, cross, k, rows]  -> weight-k piece w_k for the given entries (rows default All)
   IISum[dir, fam, cross, kmax, rows]  -> Sum_{k<=kmax} w_k  (the old <fam>_II.m / H.wl form)
   fam is a string ("H"), cross is "x" (uncrossed), "x12", "x13", "x23", "x123" or "x132". *)

IIData[dir_String, fam_String, cross_String, k_Integer] := IIData[dir, fam, cross, k] =
  Import[FileNameJoin[{dir, fam <> "_" <> cross <> "_w" <> ToString[k] <> ".wdx"}]];

(* row by row: SparseArray . (symbolic vector) is ~300x slower *)
IIWeight[dir_String, fam_String, cross_String, k_Integer, rows_ : All] := Module[{keys, S, vec},
  {keys, S} = IIData[dir, fam, cross, k];
  If[keys === {}, Return[ConstantArray[0, Length[Range[Length[S]][[rows]]]]]];
  vec = (II @@ #[[1]]) #[[2]] & /@ keys;
  With[{r = S[[#]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]] & /@ Range[Length[S]][[rows]]];

IISum[dir_String, fam_String, cross_String, kmax_Integer, rows_ : All] :=
  Sum[IIWeight[dir, fam, cross, k, rows], {k, 0, kmax}];
