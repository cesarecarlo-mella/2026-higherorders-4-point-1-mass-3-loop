(* Shared configuration for the per-family boundary pipeline (families/README.md). *)

$INT = "/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/";
$ALLINTS = $INT <> "all.ints/";
$SLC = $INT <> "NPL/DEQ_solutions/SLC/Solutions/";

(* names: all.ints uses the new labels, the SLC solution folders the legacy ones *)
inputfamall = {"A", "B1", "B2", "C", "D", "E1", "E2", "F1", "F2", "G1", "G2", "G3", "G4", "H", "I1", "I2", "I3", "B3", "B8", "A3", "B4", "A2"};
newlabels = {"A", "B1", "B2", "C", "D", "E1", "E2", "F1", "F2", "G1", "G2", "G3", "G4", "H", "I1", "I2", "I3", "N3", "L3", "M3", "P3", "R2"};
legacyL2newL = Thread[Rule[inputfamall, newlabels]];
newL2legacyL = Thread[Rule[newlabels, inputfamall]];
legacyName[f_String] := f /. newL2legacyL;

families = newlabels;
planarFamilies = {"A", "E1", "E2", "M3", "R2"};                         (* fully reducible to PL *)
rationalFamilies = {"B2", "C", "F2", "G2", "G3"};                        (* letters l[1..6] only *)
algebraicFamilies = Complement[families, planarFamilies, rationalFamilies]; (* sqrt / quadratic letters *)

(* inputs; all.ints is the reference (the SLC Atilde of D and L3 differ by +-I in the
   sqrt-normalized rows, bases/rootpos.m) *)
atildeFile[f_] := First[FileNames[{f <> "_atilde.m", f <> "_Atilde.m"}, $ALLINTS <> "a_tildes"], Missing[]];
lettersDictFile = $ALLINTS <> "a_tildes/letters.dict.m";
canonFile[f_] := $ALLINTS <> "bases/" <> f <> "_canon.m";
propsFile[f_] := $ALLINTS <> "bases/" <> f <> "_props.m";
bcFile[f_] := $ALLINTS <> "boundaries_all/" <> f <> "_bc.m";              (* BC[f, crossing, i], eps^0..6 *)
rootPos = Association[Get[$ALLINTS <> "bases/rootpos.m"]];
(* F1, G1, I2 have no shift map in all.ints: their Kira shifts (same propagators, checked) are
   taken from DEQ_solutions/Leading_Color, together with the PL -> BASISLC table found there *)
$LC = $INT <> "NPL/DEQ_solutions/Leading_Color/";
lcShiftDir = <|"F1" -> $LC <> "F1/NPL1_30389_F1_rotated_boundaries/toshift/",
   "G1" -> $LC <> "G1/NPL1_30629_G1_rotated_boundaries/toshift/", "I2" -> $LC <> "I2/NPL2_15275_I2_rotated_boundaries/toshift/"|>;
shiftMapFile[f_] := With[{p = $ALLINTS <> "result/" <> f <> "_shift.m"},
   Which[FileExistsQ[p], p, KeyExistsQ[lcShiftDir, f], lcShiftDir[f] <> "list.m", True, Missing[]]];
slcShiftedBasisFile[f_] := With[{p = $SLC <> legacyName[f] <> "/shiftedBasis/shiftedBasis.m"}, If[FileExistsQ[p], p, Missing[]]];

(* PL -> BASISLC reduction tables (PL[...] -> Sum BASISLC[i] coeff(s13, s23, d)) *)
plTableFiles := Select[Join[
   {$ALLINTS <> "result/kira_inputtoinvertPL.txt.m", $ALLINTS <> "result/kira_inputextra.txt.m",
    $SLC <> "kira_inputPLtoBASISLC_5.txt.m", $SLC <> "kira_inputtoPL_BASISLC.txt.m"},
   FileNames["kira_inputextra*.m", $SLC, 2], (# <> "kira_inputPLtoBASISLC.txt.m") & /@ Values[lcShiftDir]], FileExistsQ];
plHeads = {"PL", "PLx12", "PLx13", "PLx23", "PLx123", "PLx132"};
