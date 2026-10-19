(* Loader for run_iterint.wls output: H_opt/work/fast/II_<c>.mx = Association k -> {keys, S},
   keys = {II-word, transcendental monomial}, w_k[[i]] = Sum S[[i,j]] II@@word monomial.
   HIIWeight[k, c, rows]  -> the weight-k piece w_k for the given entries
   HIISum[kmax, c, rows]  -> Sum_{k<=kmax} w_k  (the H.wl / H_II.m form)
   c = "x" (uncrossed), "x12", "x13", "x23", "x123", "x132"; rows defaults to All. *)
$HIIDir = FileNameJoin[{DirectoryName[$InputFileName], "work", "fast"}];
HIIData[c_String] := HIIData[c] = Import[FileNameJoin[{$HIIDir, "II_" <> c <> ".mx"}], "MX"];
HIIWeight[k_Integer, c_String : "x", rows_ : All] := Module[{keys, S, vec},
  {keys, S} = HIIData[c][k];
  If[keys === {}, Return[ConstantArray[0, Length[Range[Length[S]][[rows]]]]]];
  vec = (II @@ #[[1]]) #[[2]] & /@ keys;
  With[{r = S[[#]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]] & /@ Range[Length[S]][[rows]]];
HIISum[kmax_Integer, c_String : "x", rows_ : All] := Sum[HIIWeight[k, c, rows], {k, 0, kmax}];
