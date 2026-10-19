(* Loader for the output of run_gpl_xy_fast.wls.
   pl_gpls_7_8/[<cross>/]w<w>_sparse.mx holds {keys, S}: S is a 291 x Length[keys]
   SparseArray of rationals, keys[[k]] = {x-word, y-word, transcendental monomial}, and
     sol[w][[i]] = Sum_k S[[i,k]] G[x-word; x] G[y-word; y] monomial   (G[{};_] = 1).
   cross is one of "x" (uncrossed, default), "x12", "x13", "x23", "x123", "x132".

   PLGPLSparse[w, cross]       -> {keys, S}
   PLGPLExprs[w, cross]        -> list of 291 expressions in G[...,x], G[...,y]
   PLGPLExprs[w, cross, rows]  -> only the requested entries (much cheaper at w = 8) *)

$PLGPLDir = FileNameJoin[{DirectoryName[$InputFileName], "pl_gpls_7_8"}];

PLGPLSparse[w_Integer, cross_String : "x"] := PLGPLSparse[w, cross] =
  Import[FileNameJoin[Join[{$PLGPLDir}, If[cross === "x", {}, {cross}], {"w" <> ToString[w] <> "_sparse.mx"}]], "MX"];

(* row by row: a SparseArray . (symbolic vector) is ~300x slower *)
PLGPLExprs[w_Integer, cross_String : "x", rows_ : All] := Module[{keys, S, gx, gy, vec},
  {keys, S} = PLGPLSparse[w, cross];
  gx[{}] = 1; gx[l_] := G[Sequence @@ l, x]; gy[{}] = 1; gy[l_] := G[Sequence @@ l, y];
  vec = gx[#[[1]]] gy[#[[2]]] #[[3]] & /@ keys;
  With[{r = S[[#]]}, r["NonzeroValues"] . vec[[Flatten[r["NonzeroPositions"]]]]] & /@ Range[Length[S]][[rows]]];
