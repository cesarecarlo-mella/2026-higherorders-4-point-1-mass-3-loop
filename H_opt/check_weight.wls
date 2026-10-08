(* H_opt check_weight.wls <w>
   Compares H_opt/work/sol{w}_{seg}.m against the frozen reference
   H/solutions/eps_sep/w{w}/sol{w}_{seg}.m for every seg, reducing the
   difference to canonical G-word form (ShuffleG) and reporting the
   first mismatch if the difference is not identically zero. *)

$HistoryLength = 0;
Needs["PolyLogTools`"];
$HPLAutoConvertToKnownFunctions = True;
Global`x; Global`y; Global`t; Global`delta; Global`c2; Global`c2c; Global`cxmy; Global`cxmyc; Global`sqrt; Global`eps;
Global`Basis; Global`r1; Global`r2; Global`r3; Global`r4; Global`H;
$ContextPath = Prepend[DeleteCases[$ContextPath, "Global`"], "Global`"];

wgt = ToExpression[$ScriptCommandLine[[2]]];

root = "/Users/cesarecarlomella/Documents/GitHub/2026-higherorders-4-point-1-mass-3-loop";
Hdir = FileNameJoin[{root, "H"}];
workdir = FileNameJoin[{root, "H_opt", "work"}];

segs = {1, 2, 3, 4, 5, 6, x1, y1, xmy};
segName[s_] := ToString[s];

MemoryConstrained[
  allOk = True;
  Do[
    mine = Import[FileNameJoin[{workdir, "sol" <> ToString[wgt] <> "_" <> segName[seg] <> ".m"}]];
    ref = Import[FileNameJoin[{Hdir, "solutions", "eps_sep", "w" <> ToString[wgt], "sol" <> ToString[wgt] <> "_" <> segName[seg] <> ".m"}]];
    (* G[0,0,...,0] (all letters zero) is 0 by the shuffle-regularization
       convention used throughout the reference (its Log00/Logt0 rule lists) --
       not an automatic PolyLogTools simplification, so apply it explicitly
       before comparing, or purely-cosmetic all-zero-word differences show up
       as false mismatches. *)
    (* FreeQ[{zeros},Except[0]] is WRONG here -- the list {zeros} itself (as a
       whole, at level 0) already matches Except[0] since it's not literally
       the integer 0, so FreeQ always returns False regardless of the actual
       elements. Union[{zeros}]==={0} correctly tests "every element is 0". *)
    zeroWord = HoldPattern[G[zeros__]] /; Union[{zeros}] === {0} -> 0;
    {tdiff, diff} = AbsoluteTiming[Collect[((mine - ref) /. zeroWord) // ShuffleG, _G, Together]];
    nz = Flatten[Position[diff, Except[0], {1}, Heads -> False]];
    If[nz === {},
      Print["seg ", seg, ": OK (", Length[mine], " entries), ", tdiff, "s"],
      allOk = False;
      Print["seg ", seg, ": MISMATCH at ", Length[nz], " entries, first at index ", First[nz], ":"];
      Print["  mine: ", mine[[First[nz]]]];
      Print["  ref:  ", ref[[First[nz]]]];
      Print["  diff: ", diff[[First[nz]]]];
    ],
    {seg, segs}
  ];
  Print["=== overall: ", If[allOk, "ALL SEGMENTS MATCH", "MISMATCHES FOUND"], " ==="];
, 4*1024^3]
