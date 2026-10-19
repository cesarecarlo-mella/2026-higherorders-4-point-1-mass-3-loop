(* ============================================================================
   constant_extraction.wls  --  extracts the "pure constant" part (all G-
   polylogs set to zero) of segments 2, 4, 6 at every weight 0..8, for
   external inspection.
   Clean copy of the one-off script used to build work/constant/const{w}_{i}.m.
   ============================================================================ *)

$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];
CreateDirectory["../work/constant", CreateIntermediateDirectories -> True];

sourceFile[w_, i_] := If[w == 0, "../work/weight0_solution.m", "../work/sol" <> ToString[w] <> "_" <> ToString[i] <> ".m"];

extractConstant[w_, i_] := Module[{sol, const, outFile, nonzero},
   sol = Import[sourceFile[w, i]];
   const = sol /. _G -> 0;
   outFile = "../work/constant/const" <> ToString[w] <> "_" <> ToString[i] <> ".m";
   Export[outFile, const];
   nonzero = Length[Select[const, # =!= 0 &]];
   Print["w=", w, " i=", i, "  -> ", outFile, "  (", nonzero, "/", Length[const], " nonzero)"];
   ];

Do[Do[extractConstant[w, i], {i, {2, 4, 6}}], {w, 0, 8}];
Print["=== DONE extracting constants ==="];
