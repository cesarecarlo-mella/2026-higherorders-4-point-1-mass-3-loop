(* ::Package:: *)
(* SingleScale3L.wl
   Closed forms for single-scale three-loop massless integrals.

   CONVENTIONS (all functions below)
     measure      : Prod_l  d^d k_l / (I Pi^(d/2))      (no e^(eps EulerGamma) included)
     propagators  : 1/(-k^2)^nu
     result       : coefficient C such that  integral = C * (-q^2)^(L d/2 - Sum nu)
   Family-H normalisation = eq. (3.1) of arXiv:2410.19088 (basis.m, planarInsertedXY.m): multiply by
   Exp[L eps EulerGamma] (-q^2)^(L eps); P_i = -(momentum)^2; Euclidean region
   -s12 -> x, -s13 -> y, -s23 -> 1-x-y, -m2 -> 1.  Same convention as FormFactorMasters3L.wl.

   Verified (numerically, 40-50 digits) against: iterated bubbles, a direct Feynman-parameter
   integration of the on-shell triangle, and the family-H boundary values w0..w6
   (H/solutions/eps_sep, H/results/planarInsertedXY.m).
   PROVENANCE: element 1-4 closed form = eq. (3.4) of arXiv:1306.2799 (literature); its eps^7, eps^8 terms,
   the closed form of 17/18 (derived here with iterated G-functions) and ReduceGamma were computed here.
   The form-factor values of 23/58, 67/74, 169/170 (w0-w6) come from the H pipeline.
   The Wolfram code itself has NOT been run yet (written without a kernel): run the
   self-test at the bottom of this file first. *)

BeginPackage["SingleScale3L`"];

GOne::usage = "GOne[a,b,d]: one-loop massless bubble coefficient, = C with (-q^2)^(d/2-a-b).";
Banana::usage = "Banana[{a1,..,an},d]: (n-1)-loop massless sunset/banana with n lines.";
OnShellTriangle::usage = "OnShellTriangle[a1,a2,a3,d]: one-loop massless vertex, p1^2=p2^2=0, q=p1+p2; a1: line q-vertex<->p1-vertex, a2: q-vertex<->p2-vertex, a3: p1-vertex<->p2-vertex.";
ReduceGamma::usage = "ReduceGamma[edges, ext, d] with edges = {{v1,v2,nu},...}. ext = {in,out} (2-point) or {q,p1,p2} (vertex with two on-shell legs). Reduces by bubble (parallel) and series merging; returns {C, exponent} with integral = C (-q^2)^exponent, or {$Failed, remaining edges} if the graph is not Gamma-reducible.";
HSingleScale::usage = "HSingleScale: Association with the single-scale masters of family H (basis.m positions), their integrals, kinematic scale, closed forms and expansions.";

Begin["`Private`"];

GOne[a_, b_, d_] := Gamma[a + b - d/2] Gamma[d/2 - a] Gamma[d/2 - b]/(Gamma[a] Gamma[b] Gamma[d - a - b]);

Banana[a_List, d_] := (Times @@ (Gamma[d/2 - #]/Gamma[#] & /@ a)) *
   Gamma[Total[a] - (Length[a] - 1) d/2]/Gamma[Length[a] d/2 - Total[a]];

OnShellTriangle[a1_, a2_, a3_, d_] :=
  Gamma[a1 + a2 + a3 - d/2] Gamma[d/2 - a1 - a3] Gamma[d/2 - a2 - a3]/(Gamma[a1] Gamma[a2] Gamma[d - a1 - a2 - a3]);

ReduceGamma[edges0_List, ext_List, d_] := Module[
  {es = edges0, fac = 1, changed = True, pair, v, pos, e1, e2, other, L, nu, verts, tri},
  L = Length[es] - Length[Union[Flatten[es[[All, {1, 2}]]]]] + 1;
  nu = Total[es[[All, 3]]];
  While[changed,
   changed = False;
   If[AnyTrue[es, #[[1]] === #[[2]] &], Return[{0, L d/2 - nu}, Module]]; (* tadpole: scaleless *)
   (* parallel edges = one-loop bubble *)
   pair = Catch[
     Do[If[Sort[es[[i, {1, 2}]]] === Sort[es[[j, {1, 2}]]], Throw[{i, j}]], {i, Length[es]}, {j, i + 1, Length[es]}];
     None];
   If[pair =!= None,
    fac = fac GOne[es[[pair[[1]], 3]], es[[pair[[2]], 3]], d];
    es = Append[Delete[es, List /@ pair],
      {es[[pair[[1]], 1]], es[[pair[[1]], 2]], es[[pair[[1]], 3]] + es[[pair[[2]], 3]] - d/2}];
    changed = True; Continue[]];
   (* internal 2-valent vertex = propagators in series *)
   verts = Complement[Union[Flatten[es[[All, {1, 2}]]]], ext];
   v = SelectFirst[verts, Count[Flatten[es[[All, {1, 2}]]], #] == 2 &, None];
   If[v =!= None,
    pos = Flatten[Position[es, {v, _, _} | {_, v, _}, {1}]];
    {e1, e2} = es[[pos]];
    other[e_] := If[e[[1]] === v, e[[2]], e[[1]]];
    es = Append[Delete[es, List /@ pos], {other[e1], other[e2], e1[[3]] + e2[[3]]}];
    changed = True]];
  Which[
   Length[ext] == 2 && Length[es] == 1 && Sort[es[[1, {1, 2}]]] === Sort[ext],
   {fac, L d/2 - nu},
   Length[ext] == 3 && Length[es] == 3,
   tri = Function[{a, b}, Total[Cases[es, {a, b, n_} | {b, a, n_} :> n]]];
   {fac OnShellTriangle[tri[ext[[1]], ext[[2]]], tri[ext[[1]], ext[[3]]], tri[ext[[2]], ext[[3]]], d], L d/2 - nu},
   True, {$Failed, es}]
  ];

(* ---------------------------------------------------------------------------
   Family H (H/input/basis.m, propagators H/input/props.m), d = 4 - 2 eps.
   Element n of the canonical basis (as in planarInsertedXY.m) = scale^(-3 eps) * C_n(eps).
   ------------------------------------------------------------------------- *)

With[{eps = Global`eps, x = Global`x, y = Global`y, d = 4 - 2 Global`eps},
 Module[{C1, C17, ser1, ser17},
  C1  = Exp[3 eps EulerGamma] Gamma[1 - eps]^4 Gamma[1 + 3 eps]/Gamma[1 - 4 eps];
  C17 = Exp[3 eps EulerGamma]/3 Gamma[1 + eps]^2 Gamma[1 - eps]^5 Gamma[1 + 3 eps] Gamma[1 - 3 eps]/
        (Gamma[1 - 2 eps]^2 Gamma[1 + 2 eps] Gamma[1 - 4 eps]);
  ser1 = 1 - Pi^2/4 eps^2 - 29 Zeta[3] eps^3 - 71 Pi^4/160 eps^4
    + (29 Pi^2 Zeta[3]/4 - 1263 Zeta[5]/5) eps^5 + (841 Zeta[3]^2/2 - 11539 Pi^6/24192) eps^6
    + (2059 Pi^4 Zeta[3]/160 + 1263 Pi^2 Zeta[5]/20 - 18567 Zeta[7]/7) eps^7
    + (-841 Pi^2 Zeta[3]^2/8 + 36627 Zeta[3] Zeta[5]/5 - 737687 Pi^8/1382400) eps^8;
  ser17 = 1/3 - Pi^2/12 eps^2 - 23 Zeta[3]/3 eps^3 - 11 Pi^4/96 eps^4
    + (23 Pi^2 Zeta[3]/12 - 351 Zeta[5]/5) eps^5 + (529 Zeta[3]^2/6 - 49199 Pi^6/362880) eps^6
    + (253 Pi^4 Zeta[3]/96 + 351 Pi^2 Zeta[5]/20 - 5503 Zeta[7]/7) eps^7
    + (-529 Pi^2 Zeta[3]^2/24 + 8073 Zeta[3] Zeta[5]/5 - 1006669 Pi^8/5806080) eps^8;
  HSingleScale = <|
    1 -> <|"type" -> "banana, 4 lines (2,2,2,1)", "scale" -> 1, "integral" -> Global`INT[Global`H, 2, 0, 2, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
           "prefactor" -> -eps^3, "closed" -> C1, "series" -> ser1|>,
    2 -> <|"type" -> "banana, 4 lines (2,2,1,2)", "scale" -> x, "integral" -> Global`INT[Global`H, 2, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           "prefactor" -> -eps^3 x, "closed" -> C1, "series" -> ser1|>,
    3 -> <|"type" -> "banana, 4 lines (2,1,2,2)", "scale" -> 1 - x - y, "integral" -> Global`INT[Global`H, 2, 0, 1, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
           "prefactor" -> -eps^3 (1 - x - y), "closed" -> C1, "series" -> ser1|>,
    4 -> <|"type" -> "banana, 4 lines (2,1,2,2)", "scale" -> y, "integral" -> Global`INT[Global`H, 2, 1, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
           "prefactor" -> -eps^3 y, "closed" -> C1, "series" -> ser1|>,
    17 -> <|"type" -> "two bubbles (2,1),(2,1) in series, parallel to a line (1)", "scale" -> x,
           "integral" -> Global`INT[Global`H, 2, 2, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
           "prefactor" -> -eps^3 (1 + 2 eps) x, "closed" -> C17, "series" -> ser17|>,
    18 -> <|"type" -> "as 17", "scale" -> y,
           "integral" -> Global`INT[Global`H, 0, 0, 1, 1, 0, 0, 2, 2, 0, 1, 0, 0, 0, 0, 0],
           "prefactor" -> -eps^3 (1 + 2 eps) y, "closed" -> C17, "series" -> ser17|>,
    (* massless 3-loop FORM-FACTOR integrals (vertex with two on-shell legs): no Gamma closed form.
       None of them is itself one of the form-factor masters in FormFactorMasters3L.wl:
       23/58 have the graph of A62 but with a dot; 67/74 the graph of A75 with other external legs;
       169/170 (non-planar) match none of the master graphs. Weight 7-8 values therefore need an IBP
       reduction onto FormFactorMasters3L.wl (or the H pipeline run to weight 8).
       Values = output of the H boundary pipeline (H/solutions/eps_sep, w0..w6).
       169: cross-checked against H/amflow_solution.m (point -s12 = 1) to 39 digits. *)
    23 -> <|"type" -> "form factor, K4 tetrahedron, 6 lines, one dot; legs p1,p3", "scale" -> y,
           "integral" -> Global`INT[Global`H, 0, 0, 1, 1, 0, 1, 0, 2, 1, 1, 0, 0, 0, 0, 0], "prefactor" -> -eps^5 y,
           "series" -> Pi^2/36 eps^2 - Zeta[3]/3 eps^3 - Pi^4/360 eps^4 + (-Pi^2 Zeta[3]/3 - Zeta[5]) eps^5 + (2393 Pi^6/362880 + 52 Zeta[3]^2/3) eps^6|>,
    58 -> <|"type" -> "as 23, legs p1,p2", "scale" -> x,
           "integral" -> Global`INT[Global`H, 1, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0], "prefactor" -> -eps^5 x,
           "series" -> Pi^2/36 eps^2 - Zeta[3]/3 eps^3 - Pi^4/360 eps^4 + (-Pi^2 Zeta[3]/3 - Zeta[5]) eps^5 + (2393 Pi^6/362880 + 52 Zeta[3]^2/3) eps^6|>,
    67 -> <|"type" -> "form factor, planar, 7 lines (K4 with one leg on a subdivided edge); legs p1,p3", "scale" -> y,
           "integral" -> Global`INT[Global`H, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0], "prefactor" -> eps^6 y,
           "series" -> -Zeta[3]/3 eps^3 - Pi^4/216 eps^4 + (-3 Zeta[5] - Pi^2 Zeta[3]/36) eps^5 + (-1181 Pi^6/90720 + 13 Zeta[3]^2/3) eps^6|>,
    74 -> <|"type" -> "as 67, legs p1,p2", "scale" -> x,
           "integral" -> Global`INT[Global`H, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0], "prefactor" -> eps^6 x,
           "series" -> -Zeta[3]/3 eps^3 - Pi^4/216 eps^4 + (-3 Zeta[5] - Pi^2 Zeta[3]/36) eps^5 + (-1181 Pi^6/90720 + 13 Zeta[3]^2/3) eps^6|>,
    169 -> <|"type" -> "form factor, NON-PLANAR, 8 lines (K4 core, both on-shell legs on edges adjacent at one vertex); legs p1,p2", "scale" -> x,
           "integral" -> Global`INT[Global`H, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0], "prefactor" -> eps^6 x^2,
           "series" -> 1/4 - 115 Pi^2/432 eps^2 - 305 Zeta[3]/36 eps^3 - 3149 Pi^4/17280 eps^4
              + (-6569 Zeta[5]/60 + 1273 Pi^2 Zeta[3]/144) eps^5 + (-516863 Pi^6/4354560 + 9817 Zeta[3]^2/72) eps^6|>,
    170 -> <|"type" -> "as 169, legs p1,p3", "scale" -> y,
           "integral" -> Global`INT[Global`H, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0], "prefactor" -> eps^6 y^2,
           "series" -> 1/4 - 115 Pi^2/432 eps^2 - 305 Zeta[3]/36 eps^3 - 3149 Pi^4/17280 eps^4
              + (-6569 Zeta[5]/60 + 1273 Pi^2 Zeta[3]/144) eps^5 + (-516863 Pi^6/4354560 + 9817 Zeta[3]^2/72) eps^6|>
  |>]];

End[];
EndPackage[];

(* ---------------------------------------------------------------------------
   SELF-TEST (evaluate after loading; every line should give True)

   d = 4 - 2 eps;
   (* banana via reducer == closed formula, element 1 == eq. (3.4) of arXiv:1306.2799 *)
   red1 = ReduceGamma[{{a, b, 2}, {a, b, 2}, {a, b, 2}, {a, b, 1}}, {a, b}, d];
   FullSimplify[red1[[1]]/Banana[{2, 2, 2, 1}, d]] === 1
   Normal@Series[-eps^3 Exp[3 eps EulerGamma] red1[[1]] - HSingleScale[1, "series"], {eps, 0, 8}] // FullSimplify
   (* element 17: A=q-in, B=q-out, C internal 4-valent vertex *)
   red17 = ReduceGamma[{{A, C, 2}, {A, C, 1}, {C, B, 2}, {C, B, 1}, {A, B, 1}}, {A, B}, d];
   Normal@Series[-eps^3 (1 + 2 eps) Exp[3 eps EulerGamma] red17[[1]] - HSingleScale[17, "series"], {eps, 0, 8}] // FullSimplify
   Normal@Series[HSingleScale[17, "closed"] - HSingleScale[17, "series"], {eps, 0, 8}] // FullSimplify
   (* a non-reducible graph must fail, e.g. the K4 form factor 58 *)
   ReduceGamma[{{g, f, 1}, {g, h, 1}, {g, c, 1}, {f, c, 2}, {f, h, 1}, {h, c, 1}}, {f, h, c}, d][[1]] === $Failed
   ------------------------------------------------------------------------- *)
