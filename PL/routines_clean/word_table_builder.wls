(* ============================================================================
   word_table_builder.wls  --  builds work/word_table_0to8.m, the lookup table
   that replaces ReduceCG/GetGs-based fibration (which hard-caps at weight<=6)
   for evaluating G[...,t] at t->1 at weight 7 and 8.
   Method: plain substitution t->1 on every {0,1}-letter G-word, with
   $HPLAutoConvertToKnownFunctions=True already in effect. PolyLogTools
   auto-separates the finite closed-form MZV value from an explicit symbolic
   divergence marker HPL[{1},1] (the log(1-t) singularity), which is then
   dropped -- the same regularization convention used elsewhere in this
   pipeline (Log00, the all-ones-word kill at low weight).
   Covers every sub-word length 0..8 (not just pure weight-7/8 tuples),
   because GIntegrate's output frequently factors a weight-7/8 expression
   into a SHORTER G-word times an explicit Pi/Zeta coefficient.
   Verified: reproduces the existing work/word_table_0to8.m exactly (511
   rules, 0 mismatches) when this script was written.
   ============================================================================ *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
SetDirectory[DirectoryName[$InputFileName]];

buildWordTable[maxWeight_] := Module[{rules},
   rules = {1 -> 1};  (* weight-0 word (empty G) is the trivial constant 1 *)
   Do[
     Do[
       expr = G[Sequence @@ w, t];
       val = expr /. t -> 1 /. HPL[{1}, 1] -> 0;
       AppendTo[rules, expr -> val];
       , {w, Tuples[{0, 1}, k]}];
     , {k, 1, maxWeight}];
   rules
   ];

{tBuild, wordTableAll} = AbsoluteTiming[buildWordTable[8]];
Print["built ", Length[wordTableAll], " rules in ", tBuild, "s"];

Export["../work/word_table_0to8.m", wordTableAll];
Print["exported work/word_table_0to8.m"];
