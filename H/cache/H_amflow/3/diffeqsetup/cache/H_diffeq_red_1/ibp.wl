AppendTo[$Path,"/space/ttpmath/tools_src/blade/BLAddOns/BladeIBP"];
AppendTo[$Path,"/space/ttpmath/tools_src/finiteflow/mathlink"];
AppendTo[$LibraryPath,"/space/ttpmath/tools_src/finiteflow"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<BladeIBP.m;
<<FiniteFlow.m
<<topology.wl;
IntegralOrdering = 1;
CreateDirectory["results"];
Put[ZeroSectors[H],"results/zerosectors"];
If[NonZeroSectors[H]==={}, Print["this family is trivial."];Quit[]];

Module[{rank,maxdots,GetSectors,GetSeeds,cutposi,learn,maplearn,filters,prefix,famcount,selectedfiles,lirankmin},

rank=2;
maxdots=1;
GetSectors=Automatic;
cutposi=Position[{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},1]//Flatten;
learn={8 -> {1, 2}, 7 -> {1, 2}, 6 -> {1, 2}, 5 -> {1, 2}, 4 -> {1, 2}, _ -> {1, 2}};
maplearn=Map[{1,0}+#&,learn//Association]//Normal;
filters[learn_]:=Join[If[learn=!={},{If[-Plus@@Select[#[[2]],#<0&]>(Length@Select[#[[2]],#1>0&]/.learn)[[2]],False,True]&,
If[(Plus@@Select[#[[2]],#>0&]-Length@Select[#[[2]],#>0&])>(Length@Select[#[[2]],#1>0&]/.learn)[[1]],False,True]&},{}],
If[False,{If[(#[[2]][[cutposi]]=!=ConstantArray[1,Length[cutposi]]),False,True]&},{}]];
lirankmin[_]:=0;
GetSeeds["IBP"][sector_]:=GenSeeds[sector,{0,maxdots},0,rank,filters[learn]];
GetSeeds["LI"][sector_]:=GenSeeds[sector,{0,maxdots},lirankmin[sector],rank,filters[learn]];
GetSeeds["SR"][sector_]:=GenSeeds[sector,{0,0},0,rank,filters[learn]];
GetSeeds["Map"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank,filters[maplearn]];
GetSeeds["ExtMap"][sector_]:=GenSeeds[sector,{0,maxdots+1},0,rank];
FastGenIds[H,GetSeeds,"Directory"->"ibps","LaunchKernels"->50,"GetSectors"->GetSectors];
];

Module[{exints,usints,usZeroIntegrals},

exints = {BLInternalEX[1] -> BL[H, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[2] -> BL[H, {1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[3] -> BL[H, {1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[4] -> BL[H, {1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[5] -> BL[H, {1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[6] -> BL[H, {1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[53] -> BL[H, {1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[7] -> BL[H, {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1}], BLInternalEX[8] -> BL[H, {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1, 0, 0}], BLInternalEX[9] -> BL[H, {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, -1, 0}], BLInternalEX[10] -> BL[H, {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[11] -> BL[H, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[12] -> BL[H, {1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[13] -> BL[H, {1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[14] -> BL[H, {1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[15] -> BL[H, {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[16] -> BL[H, {1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[17] -> BL[H, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[18] -> BL[H, {1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[19] -> BL[H, {1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[20] -> BL[H, {1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[21] -> BL[H, {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[50] -> BL[H, {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[51] -> BL[H, {1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[52] -> BL[H, {1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[22] -> BL[H, {1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[23] -> BL[H, {1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[24] -> BL[H, {1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[25] -> BL[H, {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[26] -> BL[H, {1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[27] -> BL[H, {1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[28] -> BL[H, {1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[29] -> BL[H, {1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[30] -> BL[H, {1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[31] -> BL[H, {1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[32] -> BL[H, {1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[33] -> BL[H, {1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[34] -> BL[H, {1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[45] -> BL[H, {1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[46] -> BL[H, {1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[47] -> BL[H, {1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[48] -> BL[H, {1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[49] -> BL[H, {1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[35] -> BL[H, {0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1}], BLInternalEX[36] -> BL[H, {0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[37] -> BL[H, {1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[38] -> BL[H, {1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[39] -> BL[H, {1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[40] -> BL[H, {1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[41] -> BL[H, {1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[42] -> BL[H, {1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[43] -> BL[H, {1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}], BLInternalEX[44] -> BL[H, {0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}]};
Export["ibps/exints_def.m",exints];
Export["ibps/ids_H_exints.mx",(-#[[1]]+#[[2]])&/@exints,"MX"];
Export["ibps/ints_H_exints.mx",IntegralsIn[exints],"MX"];
usints = {};
Export["ibps/usints_def.m",usints];
Export["ibps/ids_H_usints.mx",(-#[[1]]+#[[2]])&/@usints,"MX"];
Export["ibps/ints_H_usints.mx",IntegralsIn[usints],"MX"];


If[!({eta, eps}==={}),SerializeFastIds[FileNames["ibps/ids_H_*.mx"], {s -> -1, t -> -3, pp4 -> -15}, {eta, eps}, "UserDefinedInts"->(First/@usints),"ExtraInts"->(First/@exints)]];
];
CloseKernels[];Pause[0.1];
Quit[];