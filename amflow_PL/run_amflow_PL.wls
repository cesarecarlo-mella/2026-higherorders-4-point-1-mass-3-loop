#!/usr/bin/env wolframscript
(* ::Package:: *)

(* AMFlow evaluation pipeline for family PL (3-loop, 4-point, 1-mass).
   Adapted from Cesare's B13 template -- SAME propagator-count (15) and kinematic
   TYPE (3 massless legs + 1 massive leg), but PL's own propagator list and
   invariant names (m2, s13, s23), taken directly from Cesare's propagator
   definition (D1..D15, screenshot 2026-10-08):

     D1=k1^2            D6=(k3-p1)^2        D11=(k2-p123)^2
     D2=k2^2            D7=(k1-p12)^2       D12=(k3-p123)^2
     D3=k3^2            D8=(k2-p12)^2       D13=(k1-k2)^2
     D4=(k1-p1)^2       D9=(k3-p12)^2       D14=(k1-k3)^2
     D5=(k2-p1)^2       D10=(k1-p123)^2     D15=(k2-k3)^2

   with p12=p1+p2, p123=p1+p2+p3. PL[v1,...,v15] indices are assumed to follow
   this SAME D1..D15 ordering (standard convention, not independently verified
   against a separate index-ordering file -- flag/check this if results look
   off).

   Kinematics: p1^2=p2^2=p3^2=0 (three massless legs), p4^2=m2 (one massive
   leg), p4 = -(p1+p2+p3) by momentum conservation. s13=(p1+p3)^2=2 p1.p3,
   s23=(p2+p3)^2=2 p2.p3 (matching Canonical_basis.m's own s13/s23 naming and
   its recurring "m2-s13-s23" denominator pattern). By conservation,
   (p1+p2)^2 = m2-s13-s23, i.e. p1.p2 = (m2-s13-s23)/2.

   Usage: wolframscript -file run_amflow_PL.wls *)

Get["/Users/cesarecarlomella/Software/amflow/AMFlow.m"];
SetReductionOptions["IBPReducer" -> "Blade"];

AMFlowInfo["Family"] = PL;
AMFlowInfo["Loop"] = {k1, k2, k3};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
   p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> m2,
   p1 p2 -> (m2 - s13 - s23)/2, p1 p3 -> s13/2, p2 p3 -> s23/2
   };
AMFlowInfo["Propagator"] = {
   k1^2, k2^2, k3^2,
   (k1 - p1)^2, (k2 - p1)^2, (k3 - p1)^2,
   (k1 - p1 - p2)^2, (k2 - p1 - p2)^2, (k3 - p1 - p2)^2,
   (k1 - p1 - p2 - p3)^2, (k2 - p1 - p2 - p3)^2, (k3 - p1 - p2 - p3)^2,
   (k1 - k2)^2, (k1 - k3)^2, (k2 - k3)^2
   };
(* Euclidean numeric point -- arbitrary but non-degenerate (m2-s13-s23 != 0
   and well separated in scale), mirroring the B13 template's own choice
   style. Adjust if a different kinematic point is needed. *)
AMFlowInfo["Numeric"] = {s13 -> -1, s23 -> -1, m2 -> -3};
AMFlowInfo["NThread"] = 30;

(*norm=Exp[3 ep EulerGamma];*)

targetInts = Import["./tocompute_PL.m"];
Print["loaded ", Length[targetInts], " target integrals"];

amflowtarget = targetInts /. PL[a___] :> ToJ[{a}];

(* SolveIntegrals[jints,goal,order]'s raw epsilon-expansion target is
   order + leading, leading = -2*Length[Loop] = -6 for this 3-loop family.
   order=8 (raw eps^2) is enough for the DEEP top-sector term, whose own
   Canonical_basis.m coefficient already carries ~eps^6 (matching that
   term's own true leading pole). But each basis entry is a SUM of terms
   from many different sectors, and lower-propagator (milder-pole) masters
   carry a SMALLER explicit eps-power in their own coefficient (as low as
   eps^1, confirmed by direct inspection of Canonical_basis.m's last 10
   entries). To get THAT term's canonical eps^8 contribution, its bare
   master needs raw data up to eps^(8-1)=eps^7 -- order=8 never computed
   anywhere near that depth (uniformly stopped at raw eps^2 for everyone).
   order=13 reaches raw eps^7 exactly (order+leading = 13-6 = 7 = 8-1). *)
sol = SolveIntegrals[amflowtarget, 10, 13];
sol >> "amflow_solution_PL.m";
Print["done, exported amflow_solution_PL.m"];
