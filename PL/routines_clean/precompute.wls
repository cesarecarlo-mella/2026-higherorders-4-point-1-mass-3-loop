(* ============================================================================
   precompute.wls  --  stage 1 of the PL pipeline (weight-independent setup)
   Clean reorganized copy of ../run_precompute_pl.wls. See README.md for the
   full narrative explanation of each routine; comments here are kept short.
   Output: work/precompute_pl.mx (DumpSave of Axy, letters, dx, dy, dt, param,
   todelta, dyx1, dxy1, dxmy)
   ============================================================================ *)

$HistoryLength = 0;
SetDirectory[DirectoryName[$InputFileName]];

(* ---- letter alphabet + per-letter connection matrices Axy[letter] ---- *)
letters = {x, y, 1 - x, 1 - y, x + y, 1 - x - y};

loadAtilde[] := (Import["../input/atilde.m"] /. z -> 1 - x - y)/ep;

buildAxyTable[atilde_] := Do[
   Axy[letters[[i]]] = Coefficient[atilde, Log[letters[[i]]]],
   {i, 1, Length[letters]}];

(* full-plane x,y partial-derivative connection matrices *)
buildDxDy[] := {
   Sum[Axy[l] D[Log[l], x], {l, letters}],
   Sum[Axy[l] D[Log[l], y], {l, letters}]
   };

(* ---- segment parametrizations + path-derivative matrices dt[1..6] ----
   Each segment j maps (x,y) -> a path through the kinematic plane, t in
   [0,1], delta -> 0 the degenerate direction. dt[j] = d/dt of the connection
   matrix restricted to that path. Segments 2-6 share one recipe (series in
   delta, Simplify, delta->0, quadratic resolve, Apart); segment 1 is special
   (keeps only letters whose delta->0 value IS t, dropping "1-t" letters --
   a deliberate structural choice, not an approximation). *)
c2 = (1 - I)/2; c2c = (1 + I)/2;

buildSegment1[] := Module[{ser1, pos},
   param[1] = {x -> t, y -> delta};
   ser1 = (Series[letters /. param[1], {delta, 0, 1}, Assumptions -> 0 < t < 1] // Normal) /. delta -> 0;
   pos = Union[First /@ Position[ser1, t]];
   dt[1] = Sum[Axy[letters[[i]]] D[Log[ser1[[i]]], t], {i, pos}];
   ];

(* delta-series order is 2 for segment 2, 1 for segments 3-6 (reference recipe) *)
buildSegmentGeneric[j_, paramRule_, seriesOrder_, axyTable_] := Module[{ser, dlog},
   param[j] = paramRule;
   ser = Series[letters /. paramRule, {delta, 0, seriesOrder}, Assumptions -> {delta > 0, 0 < t < 1}] // Normal;
   dlog = ((D[Log[#], t] & /@ ser) // Simplify) /. delta -> 0 // Simplify;
   dt[j] = Apart[Apart[dlog /. 1 - 2 t + 2 t^2 -> 2 (t - c2) (t - c2c), t] . axyTable, t];
   ];

buildAllSegments[] := Module[{axyTable},
   axyTable = Table[Axy[letters[[i]]], {i, Length[letters]}];
   buildSegment1[];
   buildSegmentGeneric[2, {y -> delta t^2, x -> 1 - delta ((1 - t)^2 + t^2)}, 2, axyTable];
   buildSegmentGeneric[3, {x -> 1 - delta - t, y -> t}, 1, axyTable];
   buildSegmentGeneric[4, {x -> delta (1 - t)^2, y -> 1 - delta ((1 - t)^2 + t^2)}, 1, axyTable];
   buildSegmentGeneric[5, {x -> delta, y -> 1 - t}, 1, axyTable];
   buildSegmentGeneric[6, {x -> delta t^2, y -> delta (1 - t)^2}, 1, axyTable];
   ];

(* ---- todelta[j]: weight-independent "delta -> 0 leading letters" sum, used
   by impartSum in the weight solver for analytic continuation across the
   branch point at each segment's delta=0 endpoint. ---- *)
buildTodelta[] := Do[
   Module[{argsub, leading, pos0},
     argsub = letters /. param[j];
     leading = Series[argsub, {delta, 0, 0}, Assumptions -> 0 < t < 1] // Normal;
     pos0 = Union[First /@ Position[leading, 0]];
     todelta[j] = Sum[Axy[letters[[i]]], {i, pos0}];
     ],
   {j, 1, 6}];

(* ---- regularity-segment connection matrices (no sqrt letters => direct substitution) ---- *)
buildRegularityMatrices[dx_, dy_] := {
   Apart[dy /. x -> 1, y],   (* dyx1 *)
   Apart[dx /. y -> 1, x],   (* dxy1 *)
   Apart[(dy - dx) /. x -> -y, y]  (* dxmy *)
   };

(* ============================ run ============================ *)
atilde = loadAtilde[];
buildAxyTable[atilde];
{dx, dy} = buildDxDy[];
buildAllSegments[];
buildTodelta[];
{dyx1, dxy1, dxmy} = buildRegularityMatrices[dx, dy];

DumpSave["../work/precompute_pl.mx", {Axy, letters, dx, dy, dt, param, todelta, dyx1, dxy1, dxmy}];
Print["precompute done -> work/precompute_pl.mx"];
