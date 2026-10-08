(* ::Package:: *)
(* FormFactorMasters3L.wl
   All master integrals of the massless three-loop form factor (q^2 off-shell, p1^2 = p2^2 = 0),
   i.e. all single-scale three-loop massless vertex and propagator integrals that are masters,
   expanded up to transcendental weight 8 (B81: weight 10).

   CONVENTION (= eq. (3.1) of arXiv:2410.19088, our family-H normalisation):
       I = (-q^2)^(3 eps) Exp[3 eps EulerGamma] Int prod_l d^dk_l/(I Pi^(d/2))  prod_i 1/(-D_i)^(nu_i),  d = 4 - 2 eps,
   i.e. Euclidean-sign propagators, evaluated at -q^2 = 1.  (Same convention as H/results/planarInsertedXY.m.)
   Conversion used (checked):  I = (-1)^N Exp[3 eps EulerGamma] I_LS = Exp[3 eps EulerGamma] X_GGHIS / Gamma[1-eps]^3,
   N = sum of propagator powers, I_LS = Lee-Smirnov arXiv:1010.1334 (propagators 1/(k^2+i0)),
   X_GGHIS = Gehrmann et al. arXiv:1004.3653 appendix A.

   NAMES: keys are the names of arXiv:1010.1334; "GGHIS" gives the name/figure of arXiv:1004.3653.
   Graph pictures: figures/<figure>.svg (q enters on the left, p1 up, p2 down; broken line = crossing).

   CONSTANTS: Zeta[n]; MZV[m1,m2] = Sum_{i>j>=1} sign(m1)^i sign(m2)^j / (i^|m1| j^|m2|)   (definition of 1010.1334)
       MZV[-6,-2] = -0.01476857327566665580645 (computed here, ~20 digits);
       MZV[5,3], MZV[8,2] only in B81 (definition of 1005.0362, same as above).

   PROVENANCE: every expansion in this file is taken from the literature (LS11 up to weight 8; B81 from LSS10).
   Only the conversion to the convention above (and the re-expansion it requires) was done here.
   VERIFICATION (in Python/mpmath, 30-40 digits):
     - transcription: byte-identical (whitespace aside) to the arXiv sources;
     - closed Gamma forms == their expansions, all orders;
     - Lee-Smirnov vs Gehrmann et al.: identical at every order both papers give (all 20 integrals);
     - B81 (Lee-Smirnov-Smirnov) vs Gehrmann et al. at eps^0..eps^2.
   The Wolfram code has not been run (no kernel available when written): load it once and evaluate
   FFMasterSeries["A62"] etc. as a smoke test.

   REFERENCES (short keys used below)
     GHHS06  T. Gehrmann, G. Heinrich, T. Huber, C. Studerus, Phys. Lett. B 640 (2006) 252 [hep-ph/0607185]
     HHM08   G. Heinrich, T. Huber, D. Maitre, Phys. Lett. B 662 (2008) 344 [arXiv:0711.3590]
     HHKS09  G. Heinrich, T. Huber, D.A. Kosower, V.A. Smirnov, Phys. Lett. B 678 (2009) 359 [arXiv:0902.3512]
     BCSSS09 P.A. Baikov, K.G. Chetyrkin, A.V. Smirnov, V.A. Smirnov, M. Steinhauser, PRL 102 (2009) 212002 [arXiv:0902.3519]
     LSS10a  R.N. Lee, A.V. Smirnov, V.A. Smirnov, JHEP 1004 (2010) 020 [arXiv:1001.2887]
     LSS10   R.N. Lee, A.V. Smirnov, V.A. Smirnov, Nucl. Phys. B Proc. Suppl. 205-206 (2010) 308 [arXiv:1005.0362]
     GGHIS10 T. Gehrmann, E.W.N. Glover, T. Huber, N. Ikizlerli, C. Studerus, JHEP 1006 (2010) 094 [arXiv:1004.3653]
     LS11    R.N. Lee, V.A. Smirnov, JHEP 02 (2011) 102 [arXiv:1010.1334]
     Chetyrkin-Tkachov81  K.G. Chetyrkin, F.V. Tkachov, Nucl. Phys. B 192 (1981) 159
     CKT80   K.G. Chetyrkin, A.L. Kataev, F.V. Tkachov, Nucl. Phys. B 174 (1980) 345
     Kazakov84  D.I. Kazakov, Teor. Mat. Fiz. 58 (1984) 343
     MINCER  S.G. Gorishnii, S.A. Larin, L.R. Surguladze, F.V. Tkachov, Comput. Phys. Commun. 55 (1989) 381;
             S.A. Larin, F.V. Tkachov, J.A.M. Vermaseren, NIKHEF-H-91-18
     Bekavac05  S. Bekavac, Comput. Phys. Commun. 175 (2006) 180 [hep-ph/0505174]
*)

BeginPackage["FormFactorMasters3L`"];
FFMasters::usage = "FFMasters: Association name -> <|\"GGHIS\", \"figure\", \"props\", \"planar\", \"description\", \"closed\", \"series\", \"orders\", \"refs\"|>.";
FFMasterSeries::usage = "FFMasterSeries[name] returns the eps-expansion (convention of arXiv:2410.19088 eq. 3.1, -q^2 = 1).";
MZV::usage = "MZV[m1,m2]: (alternating) double zeta value, Sum_{i>j} sign(m1)^i sign(m2)^j/(i^|m1| j^|m2|).";
Begin["`Private`"];
eps = Global`eps;
FFMasters = <|
  "A4" -> <|"GGHIS" -> "B41", "figure" -> "figures/GGHIS_B41.svg", "props" -> 4, "planar" -> "planar",
     "description" -> "three-loop banana (4 lines), propagator type",
     "refs" -> "Chetyrkin-Tkachov81, MINCER, Bekavac05, GGHIS10, LS11",
     "orders" -> {-1, 7},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[1 - eps]^4*Gamma[3*eps - 2]/Gamma[4 - 4*eps],
     "series" -> (1/36)*eps^-1 + (71/216)*eps^0 + (3115/1296 - 1/144*Pi^2)*eps^1 + (-29/36*Zeta[3] - 71/864*Pi^2 + 109403/7776)*eps^2 + (-2059/216*Zeta[3] - 3115/5184*Pi^2 - 71/5760*Pi^4 + 3386467/46656)*eps^3 + (-90335/1296*Zeta[3] - 109403/31104*Pi^2 - 5041/34560*Pi^4 - 421/60*Zeta[5] + (29/144)*Pi^2*Zeta[3] + 96885467/279936)*eps^4 + (-3172687/7776*Zeta[3] - 3386467/186624*Pi^2 - 44233/41472*Pi^4 - 29891/360*Zeta[5] - 11539/870912*Pi^6 + (841/72)*Zeta[3]^2 + (2059/864)*Pi^2*Zeta[3] + 2631913075/1679616)*eps^5 + (-98207543/46656*Zeta[3] - 96885467/1119744*Pi^2 - 262283/432*Zeta[5] - 7767613/1244160*Pi^4 - 819269/5225472*Pi^6 - 2063/28*Zeta[7] + (421/240)*Pi^2*Zeta[5] + (2059/5760)*Pi^4*Zeta[3] + (59711/432)*Zeta[3]^2 + (90335/5184)*Pi^2*Zeta[3] + 69020223371/10077696)*eps^6 + (-2809678543/279936*Zeta[3] - 2631913075/6718464*Pi^2 - 46058663/12960*Zeta[5] - 240439157/7464960*Pi^4 - 5134855/4478976*Pi^6 - 146473/168*Zeta[7] - 737687/49766400*Pi^8 - 841/288*Pi^2*Zeta[3]^2 + (29891/1440)*Pi^2*Zeta[5] + (12209/60)*Zeta[3]*Zeta[5] + (146189/34560)*Pi^4*Zeta[3] + (3172687/31104)*Pi^2*Zeta[3] + (2619715/2592)*Zeta[3]^2 + 1766128887523/60466176)*eps^7|>,
  "A53" -> <|"GGHIS" -> "B51", "figure" -> "figures/GGHIS_B51.svg", "props" -> 5, "planar" -> "planar",
     "description" -> "propagator type: bubble x sunset in series",
     "refs" -> "Chetyrkin-Tkachov81, MINCER, Bekavac05, GGHIS10, LS11",
     "orders" -> {-2, 6},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[eps]*Gamma[1 - eps]^5*Gamma[2*eps - 1]/(Gamma[2 - 2*eps]*Gamma[3 - 3*eps]),
     "series" -> (-1/4)*eps^-2 + (-17/8)*eps^-1 + (-183/16 + (1/16)*Pi^2)*eps^0 + (-1597/32 + (13/4)*Zeta[3] + (17/32)*Pi^2)*eps^1 + (-12359/64 + (17/384)*Pi^4 + (183/64)*Pi^2 + (221/8)*Zeta[3])*eps^2 + (-88629/128 - 13/16*Pi^2*Zeta[3] + (303/20)*Zeta[5] + (289/768)*Pi^4 + (1597/128)*Pi^2 + (2379/16)*Zeta[3])*eps^3 + (-603871/256 - 221/32*Pi^2*Zeta[3] - 169/8*Zeta[3]^2 + (9439/483840)*Pi^6 + (5151/40)*Zeta[5] + (1037/512)*Pi^4 + (12359/256)*Pi^2 + (20761/32)*Zeta[3])*eps^4 + (-3970589/512 - 2379/64*Pi^2*Zeta[3] - 2873/16*Zeta[3]^2 - 221/384*Pi^4*Zeta[3] - 303/80*Pi^2*Zeta[5] + (2439/28)*Zeta[7] + (160463/967680)*Pi^6 + (55449/80)*Zeta[5] + (27149/3072)*Pi^4 + (88629/512)*Pi^2 + (160667/64)*Zeta[3])*eps^5 + (-25455447/1024 - 20761/128*Pi^2*Zeta[3] - 30927/32*Zeta[3]^2 - 3757/768*Pi^4*Zeta[3] - 5151/160*Pi^2*Zeta[5] - 3939/20*Zeta[3]*Zeta[5] + (169/32)*Pi^2*Zeta[3]^2 + (76397/7741440)*Pi^8 + (41463/56)*Zeta[7] + (575779/645120)*Pi^6 + (483891/160)*Zeta[5] + (210103/6144)*Pi^4 + (603871/1024)*Pi^2 + (1152177/128)*Zeta[3])*eps^6|>,
  "A54" -> <|"GGHIS" -> "B52", "figure" -> "figures/GGHIS_B52.svg", "props" -> 5, "planar" -> "planar",
     "description" -> "propagator type: bubble with bubble-in-bubble (= family-H element 17/18 topology)",
     "refs" -> "Chetyrkin-Tkachov81, MINCER, Bekavac05, GGHIS10, LS11",
     "orders" -> {-2, 6},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[eps]^2*Gamma[1 - eps]^5*Gamma[2 - 3*eps]*Gamma[3*eps - 1]/(Gamma[2*eps]*Gamma[2 - 2*eps]^2*Gamma[3 - 4*eps]),
     "series" -> (-1/3)*eps^-2 + (-10/3)*eps^-1 + (-64/3 + (1/12)*Pi^2)*eps^0 + (-112 + (5/6)*Pi^2 + (23/3)*Zeta[3])*eps^1 + (-528 + (11/96)*Pi^4 + (16/3)*Pi^2 + (230/3)*Zeta[3])*eps^2 + (-2336 - 23/12*Pi^2*Zeta[3] + (351/5)*Zeta[5] + (55/48)*Pi^4 + 28*Pi^2 + (1472/3)*Zeta[3])*eps^3 + (-29824/3 - 115/6*Pi^2*Zeta[3] - 529/6*Zeta[3]^2 + (49199/362880)*Pi^6 + (22/3)*Pi^4 + 702*Zeta[5] + 132*Pi^2 + 2576*Zeta[3])*eps^4 + (-123904/3 - 368/3*Pi^2*Zeta[3] - 2645/3*Zeta[3]^2 - 253/96*Pi^4*Zeta[3] - 351/20*Pi^2*Zeta[5] + (5503/7)*Zeta[7] + (49199/36288)*Pi^6 + (77/2)*Pi^4 + (22464/5)*Zeta[5] + 584*Pi^2 + 12144*Zeta[3])*eps^5 + (-507136/3 - 16928/3*Zeta[3]^2 - 644*Pi^2*Zeta[3] - 1265/48*Pi^4*Zeta[3] - 8073/5*Zeta[3]*Zeta[5] - 351/2*Pi^2*Zeta[5] + (529/24)*Pi^2*Zeta[3]^2 + (1006669/5806080)*Pi^8 + (55030/7)*Zeta[7] + (49199/5670)*Pi^6 + (363/2)*Pi^4 + (117936/5)*Zeta[5] + (7456/3)*Pi^2 + 53728*Zeta[3])*eps^6|>,
  "A64" -> <|"GGHIS" -> "B62", "figure" -> "figures/GGHIS_B62.svg", "props" -> 6, "planar" -> "planar",
     "description" -> "propagator type B62",
     "refs" -> "Chetyrkin-Tkachov81, MINCER, Bekavac05, GGHIS10, LS11",
     "orders" -> {-3, 5},
     "closed" -> None,
     "series" -> (1/3)*eps^-3 + (7/3)*eps^-2 + (31/3 - 1/12*Pi^2)*eps^-1 + (-7/12*Pi^2 + (7/3)*Zeta[3] + 103/3)*eps^0 + (-31/12*Pi^2 + (5/96)*Pi^4 + (49/3)*Zeta[3] + 235/3)*eps^1 + (-103/12*Pi^2 - 7/12*Pi^2*Zeta[3] + 19/3 + (35/96)*Pi^4 + (289/3)*Zeta[3] + (599/5)*Zeta[5])*eps^2 + (-3953/3 - 235/12*Pi^2 - 599/6*Zeta[3]^2 - 49/12*Pi^2*Zeta[3] + (967/480)*Pi^4 + (108481/362880)*Pi^6 + (1729/3)*Zeta[3] + (4193/5)*Zeta[5])*eps^3 + (-31889/3 - 4193/6*Zeta[3]^2 - 1553/480*Pi^4*Zeta[3] - 599/20*Pi^2*Zeta[5] - 289/12*Pi^2*Zeta[3] - 19/12*Pi^2 + (5263/480)*Pi^4 + (13593/7)*Zeta[7] + (108481/51840)*Pi^6 + (10213/3)*Zeta[3] + (20609/5)*Zeta[5])*eps^4 + (576*MZV[-6, -2] - 188141/3 - 21449/6*Zeta[3]^2 - 10871/480*Pi^4*Zeta[3] - 9847/5*Zeta[3]*Zeta[5] - 4193/20*Pi^2*Zeta[5] - 1729/12*Pi^2*Zeta[3] + (599/24)*Pi^2*Zeta[3]^2 + (1913939/5806080)*Pi^8 + (3953/12)*Pi^2 + (28723/480)*Pi^4 + (3695263/362880)*Pi^6 + 13593*Zeta[7] + (90257/5)*Zeta[5] + (57445/3)*Zeta[3])*eps^5|>,
  "A64d" -> <|"GGHIS" -> "B62 (one squared propagator)", "figure" -> "figures/GGHIS_B62.svg", "props" -> 7, "planar" -> "planar",
     "description" -> "UT replacement of A64 (B62 with a squared denominator, Fig. 2 of 1010.1334)",
     "refs" -> "LS11",
     "orders" -> {-1, 4},
     "closed" -> None,
     "series" -> (-6*Zeta[3])*eps^-1 + (-12*Zeta[3] - 1/10*Pi^4)*eps^0 + (-102*Zeta[5] - 24*Zeta[3] - 1/5*Pi^4 + (3/2)*Pi^2*Zeta[3])*eps^1 + (-577/2520*Pi^6 - 204*Zeta[5] - 48*Zeta[3] - 2/5*Pi^4 + 3*Pi^2*Zeta[3] + 120*Zeta[3]^2)*eps^2 + (-1413*Zeta[7] - 577/1260*Pi^6 - 408*Zeta[5] - 96*Zeta[3] - 4/5*Pi^4 + 6*Pi^2*Zeta[3] + (51/2)*Pi^2*Zeta[5] + 240*Zeta[3]^2 + (301/80)*Pi^4*Zeta[3])*eps^3 + (-288*MZV[-6, -2] - 2826*Zeta[7] - 25579/100800*Pi^8 - 577/630*Pi^6 - 816*Zeta[5] - 30*Pi^2*Zeta[3]^2 - 192*Zeta[3] - 8/5*Pi^4 + 12*Pi^2*Zeta[3] + 51*Pi^2*Zeta[5] + 480*Zeta[3]^2 + (301/40)*Pi^4*Zeta[3] + (11028/5)*Zeta[3]*Zeta[5])*eps^4|>,
  "A66" -> <|"GGHIS" -> "B61", "figure" -> "figures/GGHIS_B61.svg", "props" -> 6, "planar" -> "planar",
     "description" -> "propagator type: three bubbles in series",
     "refs" -> "Chetyrkin-Tkachov81, MINCER, Bekavac05, GGHIS10, LS11",
     "orders" -> {-3, 5},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[eps]^3*Gamma[1 - eps]^6/Gamma[2 - 2*eps]^3,
     "series" -> (1)*eps^-3 + (6)*eps^-2 + (24 - 1/4*Pi^2)*eps^-1 + (-3/2*Pi^2 - 7*Zeta[3] + 80)*eps^0 + (-6*Pi^2 - 42*Zeta[3] - 37/480*Pi^4 + 240)*eps^1 + (-168*Zeta[3] - 20*Pi^2 - 37/80*Pi^4 - 93/5*Zeta[5] + (7/4)*Pi^2*Zeta[3] + 672)*eps^2 + (-560*Zeta[3] - 60*Pi^2 - 37/20*Pi^4 - 558/5*Zeta[5] - 943/120960*Pi^6 + (49/2)*Zeta[3]^2 + (21/2)*Pi^2*Zeta[3] + 1792)*eps^3 + (-1680*Zeta[3] - 168*Pi^2 - 37/6*Pi^4 - 2232/5*Zeta[5] - 381/7*Zeta[7] - 943/20160*Pi^6 + (93/20)*Pi^2*Zeta[5] + (259/480)*Pi^4*Zeta[3] + 147*Zeta[3]^2 + 42*Pi^2*Zeta[3] + 4608)*eps^4 + (-4704*Zeta[3] - 448*Pi^2 - 37/2*Pi^4 - 1488*Zeta[5] - 2286/7*Zeta[7] - 943/5040*Pi^6 - 49/8*Pi^2*Zeta[3]^2 + (6527/9676800)*Pi^8 + (651/5)*Zeta[3]*Zeta[5] + (279/10)*Pi^2*Zeta[5] + (259/80)*Pi^4*Zeta[3] + 588*Zeta[3]^2 + 140*Pi^2*Zeta[3] + 11520)*eps^5|>,
  "A65" -> <|"GGHIS" -> "C61", "figure" -> "figures/GGHIS_C61.svg", "props" -> 6, "planar" -> "planar",
     "description" -> "product: one-loop bubble x two-loop vertex with bubble",
     "refs" -> "GGHIS10, LS11",
     "orders" -> {-3, 5},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[eps]^2*Gamma[2*eps]*Gamma[1 - 2*eps]^2*Gamma[1 - eps]^4/(Gamma[2 - 3*eps]*Gamma[2 - 2*eps]^2),
     "series" -> (1/2)*eps^-3 + (7/2)*eps^-2 + ((1/24)*Pi^2 + 33/2)*eps^-1 + (-11/2*Zeta[3] + (7/24)*Pi^2 + 131/2)*eps^0 + (-77/2*Zeta[3] - 77/960*Pi^4 + (11/8)*Pi^2 + 473/2)*eps^1 + (-363/2*Zeta[3] - 539/960*Pi^4 - 273/10*Zeta[5] - 11/24*Pi^2*Zeta[3] + (131/24)*Pi^2 + 1611/2)*eps^2 + (-1441/2*Zeta[3] - 847/320*Pi^4 - 1911/10*Zeta[5] - 5233/80640*Pi^6 - 77/24*Pi^2*Zeta[3] + (121/4)*Zeta[3]^2 + (473/24)*Pi^2 + 5281/2)*eps^3 + (-5203/2*Zeta[3] - 10087/960*Pi^4 - 9009/10*Zeta[5] - 5233/11520*Pi^6 - 121/8*Pi^2*Zeta[3] - 2313/14*Zeta[7] - 91/40*Pi^2*Zeta[5] + (847/960)*Pi^4*Zeta[3] + (847/4)*Zeta[3]^2 + (537/8)*Pi^2 + 16867/2)*eps^4 + (-17721/2*Zeta[3] - 35763/10*Zeta[5] - 36421/960*Pi^4 - 57563/26880*Pi^6 - 2313/2*Zeta[7] - 1441/24*Pi^2*Zeta[3] - 2326579/58060800*Pi^8 - 637/40*Pi^2*Zeta[5] + (121/48)*Pi^2*Zeta[3]^2 + (3003/10)*Zeta[3]*Zeta[5] + (5929/960)*Pi^4*Zeta[3] + (3993/4)*Zeta[3]^2 + (5281/24)*Pi^2 + 52905/2)*eps^5|>,
  "A51" -> <|"GGHIS" -> "A51", "figure" -> "figures/GGHIS_A51.svg", "props" -> 5, "planar" -> "planar",
     "description" -> "one-loop vertex with two-loop sunset insertion",
     "refs" -> "GHHS06 (closed form, any indices), LS11",
     "orders" -> {-2, 6},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[1 - eps]^3*Gamma[2 - 3*eps]^2*Gamma[2*eps - 1]*Gamma[3*eps - 1]/(Gamma[3 - 4*eps]*Gamma[3 - 3*eps]),
     "series" -> (1/24)*eps^-2 + (19/48)*eps^-1 + ((1/32)*Pi^2 + 233/96)*eps^0 + (-23/24*Zeta[3] + (19/64)*Pi^2 + 2363/192)*eps^1 + (-437/48*Zeta[3] + (7/2304)*Pi^4 + (233/128)*Pi^2 + 7227/128)*eps^2 + (-5359/96*Zeta[3] - 351/40*Zeta[5] - 23/32*Pi^2*Zeta[3] + (133/4608)*Pi^4 + (2363/256)*Pi^2 + 62641/256)*eps^3 + (-54349/192*Zeta[3] - 6669/80*Zeta[5] - 437/64*Pi^2*Zeta[3] - 65243/2903040*Pi^6 + (529/48)*Zeta[3]^2 + (1631/9216)*Pi^4 + (21681/512)*Pi^2 + 1575481/1536)*eps^4 + (-166221/128*Zeta[3] - 81783/160*Zeta[5] - 5359/128*Pi^2*Zeta[3] - 1239617/5806080*Pi^6 - 5503/56*Zeta[7] - 1053/160*Pi^2*Zeta[5] - 161/2304*Pi^4*Zeta[3] + (16541/18432)*Pi^4 + (10051/96)*Zeta[3]^2 + (187923/1024)*Pi^2 + 12951211/3072)*eps^5 + (-1440743/256*Zeta[3] - 829413/320*Zeta[5] - 54349/256*Pi^2*Zeta[3] - 15201619/11612160*Pi^6 - 104557/112*Zeta[7] - 20007/320*Pi^2*Zeta[5] - 75527/1720320*Pi^8 - 3059/4608*Pi^4*Zeta[3] + (529/64)*Pi^2*Zeta[3]^2 + (8073/40)*Zeta[3]*Zeta[5] + (16863/4096)*Pi^4 + (123257/192)*Zeta[3]^2 + (1575481/2048)*Pi^2 + 105241601/6144)*eps^6|>,
  "A52" -> <|"GGHIS" -> "A52", "figure" -> "figures/GGHIS_A52.svg", "props" -> 5, "planar" -> "planar",
     "description" -> "one-loop vertex with two bubble insertions (lines q-p2 and p1-p2)",
     "refs" -> "GHHS06 (closed form, any indices), LS11",
     "orders" -> {-2, 6},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[eps]*Gamma[1 - 2*eps]*Gamma[1 - eps]^4*Gamma[2 - 3*eps]*Gamma[3*eps - 1]/(Gamma[2 - 2*eps]^2*Gamma[3 - 4*eps]),
     "series" -> (-1/6)*eps^-2 + (-5/3)*eps^-1 + (-32/3 - 1/24*Pi^2)*eps^0 + (-56 - 5/12*Pi^2 + (23/6)*Zeta[3])*eps^1 + (-264 - 8/3*Pi^2 + (25/576)*Pi^4 + (115/3)*Zeta[3])*eps^2 + (-1168 - 14*Pi^2 + (23/24)*Pi^2*Zeta[3] + (351/10)*Zeta[5] + (125/288)*Pi^4 + (736/3)*Zeta[3])*eps^3 + (-14912/3 - 66*Pi^2 - 529/12*Zeta[3]^2 + (66041/725760)*Pi^6 + (115/12)*Pi^2*Zeta[3] + (25/9)*Pi^4 + 351*Zeta[5] + 1288*Zeta[3])*eps^4 + (-61952/3 - 292*Pi^2 - 2645/6*Zeta[3]^2 - 575/576*Pi^4*Zeta[3] + (351/40)*Pi^2*Zeta[5] + (5503/14)*Zeta[7] + (184/3)*Pi^2*Zeta[3] + (66041/72576)*Pi^6 + (175/12)*Pi^4 + (11232/5)*Zeta[5] + 6072*Zeta[3])*eps^5 + (-253568/3 - 3728/3*Pi^2 - 8464/3*Zeta[3]^2 - 2875/288*Pi^4*Zeta[3] - 8073/10*Zeta[3]*Zeta[5] - 529/48*Pi^2*Zeta[3]^2 + (351/4)*Pi^2*Zeta[5] + (1513373/11612160)*Pi^8 + 322*Pi^2*Zeta[3] + (27515/7)*Zeta[7] + (66041/11340)*Pi^6 + (275/4)*Pi^4 + (58968/5)*Zeta[5] + 26864*Zeta[3])*eps^6|>,
  "A61" -> <|"GGHIS" -> "A61", "figure" -> "figures/GGHIS_A61.svg", "props" -> 6, "planar" -> "planar",
     "description" -> "one-loop vertex with two bubbles in series on line p1-p2",
     "refs" -> "GHHS06 (closed form, any indices), LS11",
     "orders" -> {-3, 5},
     "closed" -> Exp[3*EulerGamma*eps]*Gamma[eps]^2*Gamma[3*eps]*Gamma[1 - 3*eps]^2*Gamma[1 - eps]^4/(Gamma[2 - 4*eps]*Gamma[2 - 2*eps]^2),
     "series" -> (1/3)*eps^-3 + (8/3)*eps^-2 + ((1/4)*Pi^2 + 44/3)*eps^-1 + (-17/3*Zeta[3] + 2*Pi^2 + 208/3)*eps^0 + (-136/3*Zeta[3] + (83/1440)*Pi^4 + 11*Pi^2 + 304)*eps^1 + (-748/3*Zeta[3] - 281/5*Zeta[5] - 17/4*Pi^2*Zeta[3] + (83/180)*Pi^4 + 52*Pi^2 + 1280)*eps^2 + (-3536/3*Zeta[3] - 2248/5*Zeta[5] - 34*Pi^2*Zeta[3] - 44651/362880*Pi^6 + (289/6)*Zeta[3]^2 + (913/360)*Pi^4 + 228*Pi^2 + 15808/3)*eps^3 + (-5168*Zeta[3] - 12364/5*Zeta[5] - 187*Pi^2*Zeta[3] - 44651/45360*Pi^6 - 4817/7*Zeta[7] - 843/20*Pi^2*Zeta[5] - 1411/1440*Pi^4*Zeta[3] + (1156/3)*Zeta[3]^2 + (1079/90)*Pi^4 + 960*Pi^2 + 64256/3)*eps^4 + (-21760*Zeta[3] - 58448/5*Zeta[5] - 884*Pi^2*Zeta[3] - 38536/7*Zeta[7] - 491161/90720*Pi^6 - 1686/5*Pi^2*Zeta[5] - 960457/3225600*Pi^8 - 1411/180*Pi^4*Zeta[3] + (289/8)*Pi^2*Zeta[3]^2 + (4777/5)*Zeta[3]*Zeta[5] + (6358/3)*Zeta[3]^2 + (1577/30)*Pi^4 + 3952*Pi^2 + 259328/3)*eps^5|>,
  "A62" -> <|"GGHIS" -> "A62", "figure" -> "figures/GGHIS_A62.svg", "props" -> 6, "planar" -> "planar",
     "description" -> "K4 tetrahedron vertex (= family-H 23/58 topology, there with a dot)",
     "refs" -> "HHM08, GGHIS10, LS11",
     "orders" -> {-1, 4},
     "closed" -> None,
     "series" -> (2*Zeta[3])*eps^-1 + ((7/180)*Pi^4 + 18*Zeta[3])*eps^0 + (-7/6*Pi^2*Zeta[3] + 10*Zeta[5] + (7/20)*Pi^4 + 122*Zeta[3])*eps^1 + (-21/2*Pi^2*Zeta[3] - 78*Zeta[3]^2 - 473/15120*Pi^6 + 90*Zeta[5] + (427/180)*Pi^4 + 738*Zeta[3])*eps^2 + (-702*Zeta[3]^2 - 427/6*Pi^2*Zeta[3] - 415/144*Pi^4*Zeta[3] - 473/1680*Pi^6 - 445/2*Zeta[7] - 5/2*Pi^2*Zeta[5] + 610*Zeta[5] + (287/20)*Pi^4 + 4202*Zeta[3])*eps^3 + ((992/3)*MZV[-6, -2] - 4758*Zeta[3]^2 - 861/2*Pi^2*Zeta[3] - 415/16*Pi^4*Zeta[3] - 425689/1814400*Pi^8 - 4005/2*Zeta[7] - 28853/15120*Pi^6 - 2576/5*Zeta[3]*Zeta[5] - 45/2*Pi^2*Zeta[5] + (253/6)*Pi^2*Zeta[3]^2 + 3690*Zeta[5] + (14707/180)*Pi^4 + 23058*Zeta[3])*eps^4|>,
  "A63" -> <|"GGHIS" -> "A63", "figure" -> "figures/GGHIS_A63.svg", "props" -> 6, "planar" -> "planar",
     "description" -> "two-loop planar vertex with bubble insertion",
     "refs" -> "GHHS06 (3F2), LS11",
     "orders" -> {-3, 5},
     "closed" -> None,
     "series" -> (1/6)*eps^-3 + (3/2)*eps^-2 + ((1/8)*Pi^2 + 55/6)*eps^-1 + (-35/6*Zeta[3] + (9/8)*Pi^2 + 95/2)*eps^0 + (-105/2*Zeta[3] - 77/2880*Pi^4 + (55/8)*Pi^2 + 1351/6)*eps^1 + (-1925/6*Zeta[3] - 651/10*Zeta[5] - 49/24*Pi^2*Zeta[3] - 77/320*Pi^4 + (285/8)*Pi^2 + 2023/2)*eps^2 + (-3325/2*Zeta[3] - 5859/10*Zeta[5] - 147/8*Pi^2*Zeta[3] - 847/576*Pi^4 - 93451/725760*Pi^6 + (1141/12)*Zeta[3]^2 + (1351/8)*Pi^2 + 26335/6)*eps^3 + (-47285/6*Zeta[3] - 7161/2*Zeta[5] - 2695/24*Pi^2*Zeta[3] - 93451/80640*Pi^6 - 1463/192*Pi^4 - 9017/14*Zeta[7] - 713/40*Pi^2*Zeta[5] + (511/320)*Pi^4*Zeta[3] + (3423/4)*Zeta[3]^2 + (6069/8)*Pi^2 + 37215/2)*eps^4 + (-544/9*MZV[-6, -2] - 70805/2*Zeta[3] - 37107/2*Zeta[5] - 4655/8*Pi^2*Zeta[3] - 1027961/145152*Pi^6 - 81153/14*Zeta[7] - 104027/2880*Pi^4 - 2022493/11612160*Pi^8 - 6417/40*Pi^2*Zeta[5] + (623/48)*Pi^2*Zeta[3]^2 + (4599/320)*Pi^4*Zeta[3] + (11195/6)*Zeta[3]*Zeta[5] + (62755/12)*Zeta[3]^2 + (26335/8)*Pi^2 + 465751/6)*eps^5|>,
  "A71" -> <|"GGHIS" -> "A71", "figure" -> "figures/GGHIS_A71.svg", "props" -> 7, "planar" -> "non-planar",
     "description" -> "two-loop crossed vertex with bubble on line q-U1",
     "refs" -> "GHHS06 (4F3/5F4), LS11",
     "orders" -> {-5, 3},
     "closed" -> None,
     "series" -> (1/4)*eps^-5 + (1/2)*eps^-4 + (1 - 11/48*Pi^2)*eps^-3 + (-41/4*Zeta[3] - 11/24*Pi^2 + 2)*eps^-2 + (-41/2*Zeta[3] - 227/1152*Pi^4 - 11/12*Pi^2 + 4)*eps^-1 + (-1763/20*Zeta[5] - 41*Zeta[3] - 227/576*Pi^4 - 11/6*Pi^2 + 8 + (355/48)*Pi^2*Zeta[3])*eps^0 + (-1763/10*Zeta[5] - 82*Zeta[3] - 227/288*Pi^4 - 32759/483840*Pi^6 - 11/3*Pi^2 + 16 + (355/24)*Pi^2*Zeta[3] + (1649/8)*Zeta[3]^2)*eps^1 + (-16397/28*Zeta[7] - 1763/5*Zeta[5] - 164*Zeta[3] - 227/144*Pi^4 - 32759/241920*Pi^6 - 22/3*Pi^2 + 32 + (355/12)*Pi^2*Zeta[3] + (1649/4)*Zeta[3]^2 + (4971/80)*Pi^2*Zeta[5] + (43847/5760)*Pi^4*Zeta[3])*eps^2 + (-352*MZV[-6, -2] - 12251/96*Pi^2*Zeta[3]^2 - 16397/14*Zeta[7] - 3526/5*Zeta[5] - 328*Zeta[3] - 227/72*Pi^4 - 32759/120960*Pi^6 - 44/3*Pi^2 + 64 + (355/6)*Pi^2*Zeta[3] + (1649/2)*Zeta[3]^2 + (4971/40)*Pi^2*Zeta[5] + (43847/2880)*Pi^4*Zeta[3] + (23355197/116121600)*Pi^8 + (60603/20)*Zeta[3]*Zeta[5])*eps^3|>,
  "A72" -> <|"GGHIS" -> "A72", "figure" -> "figures/GGHIS_A72.svg", "props" -> 7, "planar" -> "non-planar",
     "description" -> "two-loop crossed vertex with bubble on line adjacent to p1",
     "refs" -> "GHHS06 (hypergeometric), LS11",
     "orders" -> {-3, 3},
     "closed" -> None,
     "series" -> ((1/12)*Pi^2)*eps^-3 + ((1/6)*Pi^2 + 2*Zeta[3])*eps^-2 + ((1/3)*Pi^2 + 4*Zeta[3] + (17/180)*Pi^4)*eps^-1 + (-9/4*Pi^2*Zeta[3] + (2/3)*Pi^2 + 8*Zeta[3] + 15*Zeta[5] + (17/90)*Pi^4)*eps^0 + (-75*Zeta[3]^2 - 9/2*Pi^2*Zeta[3] + (1199/362880)*Pi^6 + (4/3)*Pi^2 + 16*Zeta[3] + 30*Zeta[5] + (17/45)*Pi^4)*eps^1 + (-2959/720*Pi^4*Zeta[3] - 150*Zeta[3]^2 - 262/15*Pi^2*Zeta[5] - 1383/8*Zeta[7] - 9*Pi^2*Zeta[3] + (1199/181440)*Pi^6 + (8/3)*Pi^2 + 32*Zeta[3] + 60*Zeta[5] + (34/45)*Pi^4)*eps^2 + ((2912/9)*MZV[-6, -2] - 21377/75600*Pi^8 - 12493/15*Zeta[3]*Zeta[5] - 2959/360*Pi^4*Zeta[3] - 300*Zeta[3]^2 - 524/15*Pi^2*Zeta[5] - 1383/4*Zeta[7] - 18*Pi^2*Zeta[3] + (1199/90720)*Pi^6 + (16/3)*Pi^2 + 64*Zeta[3] + 120*Zeta[5] + (68/45)*Pi^4 + (883/24)*Pi^2*Zeta[3]^2)*eps^3|>,
  "A73" -> <|"GGHIS" -> "A73", "figure" -> "figures/GGHIS_A73.svg", "props" -> 7, "planar" -> "planar",
     "description" -> "genuine 7-line planar vertex",
     "refs" -> "HHM08, GGHIS10, LS11",
     "orders" -> {-1, 2},
     "closed" -> None,
     "series" -> (-10*Zeta[5] - 1/6*Pi^2*Zeta[3])*eps^-1 + (-119/2160*Pi^6 - 31/2*Zeta[3]^2)*eps^0 + (-1279/4*Zeta[7] - 23/3*Pi^2*Zeta[5] - 61/360*Pi^4*Zeta[3])*eps^1 + (-656/3*MZV[-6, -2] - 775*Zeta[3]*Zeta[5] - 48707/604800*Pi^8 + (149/24)*Pi^2*Zeta[3]^2)*eps^2|>,
  "A74" -> <|"GGHIS" -> "A74", "figure" -> "figures/GGHIS_A74.svg", "props" -> 7, "planar" -> "non-planar",
     "description" -> "genuine 7-line non-planar vertex",
     "refs" -> "HHM08, GGHIS10, LS11",
     "orders" -> {-2, 3},
     "closed" -> None,
     "series" -> (6*Zeta[3])*eps^-2 + ((11/90)*Pi^4 + 36*Zeta[3])*eps^-1 + (-7/2*Pi^2*Zeta[3] + 46*Zeta[5] + (11/15)*Pi^4 + 216*Zeta[3])*eps^0 + (-288*Zeta[3]^2 - 21*Pi^2*Zeta[3] - 109/1080*Pi^6 + 276*Zeta[5] + (22/5)*Pi^4 + 1296*Zeta[3])*eps^1 + (-1728*Zeta[3]^2 - 126*Pi^2*Zeta[3] - 8381/720*Pi^4*Zeta[3] - 1309*Zeta[7] - 109/180*Pi^6 - 11/2*Pi^2*Zeta[5] + 1656*Zeta[5] + (132/5)*Pi^4 + 7776*Zeta[3])*eps^2 + ((7712/3)*MZV[-6, -2] - 72323/43200*Pi^8 - 10368*Zeta[3]^2 - 756*Pi^2*Zeta[3] - 8381/120*Pi^4*Zeta[3] - 7854*Zeta[7] - 109/30*Pi^6 - 3668/5*Zeta[3]*Zeta[5] - 33*Pi^2*Zeta[5] + 160*Pi^2*Zeta[3]^2 + 9936*Zeta[5] + (792/5)*Pi^4 + 46656*Zeta[3])*eps^3|>,
  "A75" -> <|"GGHIS" -> "A75", "figure" -> "figures/GGHIS_A75.svg", "props" -> 7, "planar" -> "non-planar",
     "description" -> "genuine 7-line non-planar vertex (same graph as family-H 67/74, but different external legs)",
     "refs" -> "HHM08, GGHIS10, LS11",
     "orders" -> {0, 3},
     "closed" -> None,
     "series" -> (10*Zeta[5] + 2*Pi^2*Zeta[3])*eps^0 + (18*Zeta[3]^2 + 60*Zeta[5] + (11/162)*Pi^6 + 12*Pi^2*Zeta[3])*eps^1 + (-7/6*Pi^2*Zeta[5] + (13/30)*Pi^4*Zeta[3] + (307/2)*Zeta[7] + 108*Zeta[3]^2 + 360*Zeta[5] + (11/27)*Pi^6 + 72*Pi^2*Zeta[3])*eps^2 + (-992/9*MZV[-6, -2] - 185/2*Pi^2*Zeta[3]^2 - 1690/3*Zeta[3]*Zeta[5] - 7*Pi^2*Zeta[5] + (2431/680400)*Pi^8 + (13/5)*Pi^4*Zeta[3] + 921*Zeta[7] + 648*Zeta[3]^2 + 2160*Zeta[5] + (22/9)*Pi^6 + 432*Pi^2*Zeta[3])*eps^3|>,
  "A8" -> <|"GGHIS" -> "A81", "figure" -> "figures/GGHIS_A81.svg", "props" -> 8, "planar" -> "non-planar",
     "description" -> "genuine 8-line non-planar vertex",
     "refs" -> "HHM08, GGHIS10, LS11",
     "orders" -> {-2, 3},
     "closed" -> None,
     "series" -> (-8/3*Zeta[3])*eps^-2 + (-5/27*Pi^4 + 8*Zeta[3])*eps^-1 + (-352/3*Zeta[5] - 24*Zeta[3] + (5/9)*Pi^4 + (58/9)*Pi^2*Zeta[3])*eps^0 + (-58/3*Pi^2*Zeta[3] - 5/3*Pi^4 - 5261/34020*Pi^6 + 72*Zeta[3] + (340/3)*Zeta[3]^2 + 352*Zeta[5])*eps^1 + (-1056*Zeta[5] - 2443/3*Zeta[7] - 340*Zeta[3]^2 - 216*Zeta[3] + (256/9)*Pi^2*Zeta[5] + (5261/11340)*Pi^6 + 5*Pi^4 + 58*Pi^2*Zeta[3] + (5431/540)*Pi^4*Zeta[3])*eps^2 + (-9728/9*MZV[-6, -2] - 5431/180*Pi^4*Zeta[3] - 1775/9*Pi^2*Zeta[3]^2 - 174*Pi^2*Zeta[3] - 15*Pi^4 - 5261/3780*Pi^6 - 256/3*Pi^2*Zeta[5] + 648*Zeta[3] + 1020*Zeta[3]^2 + 2443*Zeta[7] + (9968/5)*Zeta[3]*Zeta[5] + 3168*Zeta[5] + (1625251/4082400)*Pi^8)*eps^3|>,
  "A91" -> <|"GGHIS" -> "A91", "figure" -> "figures/GGHIS_A91.svg", "props" -> 9, "planar" -> "planar",
     "description" -> "9-line planar ladder vertex",
     "refs" -> "HHKS09, LSS10, GGHIS10, LS11",
     "orders" -> {-5, 3},
     "closed" -> None,
     "series" -> (-1/18)*eps^-5 + (1/2)*eps^-4 + (-53/18 - 29/216*Pi^2)*eps^-3 + (-35/18*Zeta[3] + (149/216)*Pi^2 + 29/2)*eps^-2 + (-129/2 - 5473/25920*Pi^4 - 139/72*Pi^2 + (307/18)*Zeta[3])*eps^-1 + (-793/10*Zeta[5] - 1103/18*Zeta[3] - 871/216*Pi^2*Zeta[3] + (19/8)*Pi^2 + (3125/5184)*Pi^4 + 537/2)*eps^0 + (-2133/2 - 76801/186624*Pi^6 - 5521/36*Zeta[3]^2 - 2969/216*Pi^2*Zeta[3] - 4717/2880*Pi^4 + (97/8)*Pi^2 + (287/2)*Zeta[3] + (8251/30)*Zeta[5])*eps^1 + (-228799/126*Zeta[7] - 969/8*Pi^2 - 22487/30*Zeta[5] - 138403/25920*Pi^4*Zeta[3] - 799/4*Zeta[3]^2 - 195/2*Zeta[3] + (11987/1080)*Pi^2*Zeta[5] + (1333/320)*Pi^4 + (4286603/6531840)*Pi^6 + (5887/72)*Pi^2*Zeta[3] + 8181/2)*eps^2 + ((78272/81)*MZV[-6, -2] - 30537/2 - 325463531/313528320*Pi^8 - 2047/8*Pi^2*Zeta[3] - 295667/1080*Pi^2*Zeta[5] - 2685/2*Zeta[3] - 9163/960*Pi^4 - 38845/54*Zeta[3]*Zeta[5] - 125549/25920*Pi^4*Zeta[3] - 9379/16128*Pi^6 + (14613/10)*Zeta[5] + (70259/432)*Pi^2*Zeta[3]^2 + (73819/36)*Zeta[3]^2 + (542683/126)*Zeta[7] + (5589/8)*Pi^2)*eps^3|>,
  "A91n" -> <|"GGHIS" -> "A91 with numerator", "figure" -> "figures/GGHIS_A91.svg", "props" -> 9, "planar" -> "planar",
     "description" -> "UT replacement of A91 (numerator defined in 0902.3512)",
     "refs" -> "HHKS09, LS11",
     "orders" -> {-6, 2},
     "closed" -> None,
     "series" -> (-1/36)*eps^-6 + (0)*eps^-5 + (-7/144*Pi^2)*eps^-4 + (-55/36*Zeta[3])*eps^-3 + (-5329/51840*Pi^4)*eps^-2 + (-1171/432*Pi^2*Zeta[3] - 1199/60*Zeta[5])*eps^-1 + (-2606843/13063680*Pi^6 - 353/8*Zeta[3]^2)*eps^0 + (-180625/252*Zeta[7] - 46657/2160*Pi^2*Zeta[5] - 61831/51840*Pi^4*Zeta[3])*eps^1 + ((23680/81)*MZV[-6, -2] - 1527522487/3135283200*Pi^8 - 85381/108*Zeta[3]*Zeta[5] + (80579/864)*Pi^2*Zeta[3]^2)*eps^2|>,
  "A92" -> <|"GGHIS" -> "A92", "figure" -> "figures/GGHIS_A92.svg", "props" -> 9, "planar" -> "non-planar",
     "description" -> "9-line non-planar vertex",
     "refs" -> "HHKS09 (poles), LSS10a (finite), LS11",
     "orders" -> {-6, 2},
     "closed" -> None,
     "series" -> (2/9)*eps^-6 + (5/6)*eps^-5 + (-17/54*Pi^2 - 20/9)*eps^-4 + (-31/3*Zeta[3] - 181/216*Pi^2 + 50/9)*eps^-3 + (-119/432*Pi^4 - 347/18*Zeta[3] - 110/9 + (17/9)*Pi^2)*eps^-2 + (-2507/15*Zeta[5] - 19/6*Pi^2 - 163/960*Pi^4 + 170/9 + (514/9)*Zeta[3] + (341/36)*Pi^2*Zeta[3])*eps^-1 + (-2783/6*Zeta[5] - 195551/544320*Pi^6 - 1516/9*Zeta[3] - 1/2*Pi^2 + 130/9 + 29*Zeta[3]^2 + (943/1080)*Pi^4 + (737/24)*Pi^2*Zeta[3])*eps^0 + (-43329/14*Zeta[7] - 2416889/2177280*Pi^6 - 1357/18*Pi^2*Zeta[3] - 19169/36*Zeta[3]^2 - 8801/2160*Pi^4 - 2950/9 - 5963/4320*Pi^4*Zeta[3] + (83/2)*Pi^2 + (4444/9)*Zeta[3] + (2830/3)*Zeta[5] + (8183/60)*Pi^2*Zeta[5])*eps^1 + ((20752/9)*MZV[-6, -2] - 24178127/14515200*Pi^8 - 580805/42*Zeta[7] - 58517/1728*Pi^4*Zeta[3] - 569/2*Pi^2 - 6419/36*Pi^2*Zeta[3]^2 - 12916/9*Zeta[3] - 3112/3*Zeta[5] + (1214/9)*Zeta[3]^2 + (64733/54432)*Pi^6 + (7795/432)*Pi^4 + (1433/9)*Pi^2*Zeta[3] + 19090/9 + (146521/360)*Pi^2*Zeta[5] + (101288/15)*Zeta[3]*Zeta[5])*eps^2|>,
  "A92n" -> <|"GGHIS" -> "A92 with numerator", "figure" -> "figures/GGHIS_A92.svg", "props" -> 9, "planar" -> "non-planar",
     "description" -> "UT replacement of A92 (numerator defined in 0902.3512)",
     "refs" -> "HHKS09, LS11",
     "orders" -> {-6, 2},
     "closed" -> None,
     "series" -> (-2/9)*eps^-6 + (0)*eps^-5 + ((17/54)*Pi^2)*eps^-4 + ((31/3)*Zeta[3])*eps^-3 + ((119/432)*Pi^4)*eps^-2 + (-341/36*Pi^2*Zeta[3] + (2507/15)*Zeta[5])*eps^-1 + (-29*Zeta[3]^2 + (195551/544320)*Pi^6)*eps^0 + (-8183/60*Pi^2*Zeta[5] + (5963/4320)*Pi^4*Zeta[3] + (43329/14)*Zeta[7])*eps^1 + (-20752/9*MZV[-6, -2] - 101288/15*Zeta[3]*Zeta[5] + (6419/36)*Pi^2*Zeta[3]^2 + (24178127/14515200)*Pi^8)*eps^2|>,
  "A94" -> <|"GGHIS" -> "A94", "figure" -> "figures/GGHIS_A94.svg", "props" -> 9, "planar" -> "non-planar",
     "description" -> "9-line non-planar vertex",
     "refs" -> "HHKS09 (poles), LSS10a (finite), LS11",
     "orders" -> {-6, 2},
     "closed" -> None,
     "series" -> (1/9)*eps^-6 + (8/9)*eps^-5 + (-43/108*Pi^2 - 1)*eps^-4 + (-53/27*Pi^2 - 109/9*Zeta[3] - 14/9)*eps^-3 + (-608/9*Zeta[3] + (481/12960)*Pi^4 + 17 + (311/108)*Pi^2)*eps^-2 + (-84 - 3463/45*Zeta[5] - 85/108*Pi^4 - 11/18*Pi^2 + (949/9)*Zeta[3] + (2975/108)*Pi^2*Zeta[3])*eps^-1 + (-7868/15*Zeta[5] - 77/4*Pi^2 - 434/9*Zeta[3] + (2539/2592)*Pi^4 + 339 + (247613/466560)*Pi^6 + (3115/6)*Zeta[3]^2 + (299/3)*Pi^2*Zeta[3])*eps^0 + (-19499/108*Pi^2*Zeta[3] - 1242 - 589*Zeta[3] - 8564/63*Zeta[7] + (487/432)*Pi^4 + (30067/45)*Zeta[5] + (25567/30240)*Pi^6 + 112*Pi^2 + (38903/2592)*Pi^4*Zeta[3] + (113629/540)*Pi^2*Zeta[5] + (18512/9)*Zeta[3]^2)*eps^1 + (-76288/81*MZV[-6, -2] - 152299/216*Pi^2*Zeta[3]^2 - 26291/6*Zeta[3]^2 - 72537/14*Zeta[7] - 1887/4*Pi^2 - 9901847/3265920*Pi^6 - 491/32*Pi^4 + (7874/15)*Zeta[5] + (4019/18)*Pi^2*Zeta[3] + 4293 + 3756*Zeta[3] + (730841/135)*Zeta[3]*Zeta[5] + (9037/135)*Pi^4*Zeta[3] + (35728/45)*Pi^2*Zeta[5] + (30535087/31352832)*Pi^8)*eps^2|>,
  "A94n" -> <|"GGHIS" -> "A94 with numerator", "figure" -> "figures/GGHIS_A94.svg", "props" -> 9, "planar" -> "non-planar",
     "description" -> "UT replacement of A94 (numerator defined in 0902.3512)",
     "refs" -> "HHKS09, LS11",
     "orders" -> {-6, 2},
     "closed" -> None,
     "series" -> (1/9)*eps^-6 + (0)*eps^-5 + (-43/108*Pi^2)*eps^-4 + (-109/9*Zeta[3])*eps^-3 + ((481/12960)*Pi^4)*eps^-2 + (-3463/45*Zeta[5] + (2975/108)*Pi^2*Zeta[3])*eps^-1 + ((247613/466560)*Pi^6 + (3115/6)*Zeta[3]^2)*eps^0 + (-8564/63*Zeta[7] + (38903/2592)*Pi^4*Zeta[3] + (113629/540)*Pi^2*Zeta[5])*eps^1 + (-76288/81*MZV[-6, -2] - 152299/216*Pi^2*Zeta[3]^2 + (730841/135)*Zeta[3]*Zeta[5] + (30535087/31352832)*Pi^8)*eps^2|>,
  "B81" -> <|"GGHIS" -> "B81", "figure" -> "figures/GGHIS_B81.svg", "props" -> 8, "planar" -> "non-planar",
     "description" -> "non-planar three-loop propagator (NO, 'I17' of 1005.0362); 20 Zeta[5] at eps^0",
     "refs" -> "CKT80 (eps^0), Kazakov84 (eps^1), Bekavac05 (eps^2), LSS10 (eps^3..eps^5), GGHIS10",
     "orders" -> {0, 5},
     "closed" -> None,
     "series" -> 1/(1 - 2 eps) (20 Zeta[5] + eps (68 Zeta[3]^2 + 10 Pi^6/189) + eps^2 (34 Pi^4 Zeta[3]/15 - 5 Pi^2 Zeta[5] + 450 Zeta[7])
        + eps^3 (-9072/5 MZV[5, 3] - 2588 Zeta[3] Zeta[5] - 17 Pi^2 Zeta[3]^2 + 6487 Pi^8/10500)
        + eps^4 (-4897 Pi^6 Zeta[3]/630 - 6068 Zeta[3]^3/3 + 13063 Pi^4 Zeta[5]/120 - 225 Pi^2 Zeta[7]/2 + 88036 Zeta[9]/9)
        + eps^5 (2268/5 Pi^2 MZV[5, 3] + 42513 MZV[8, 2] - 145328 Zeta[3] Zeta[7] - 73394 Zeta[5]^2 + 647 Pi^2 Zeta[3] Zeta[5]
           - 11813 Pi^4 Zeta[3]^2/120 + 28138577 Pi^10/9355500))|>
|>;
FFMasterSeries[name_String] := FFMasters[name, "series"];
MZV /: N[MZV[-6, -2], ___] := -0.01476857327566665580645;
End[];
EndPackage[];
