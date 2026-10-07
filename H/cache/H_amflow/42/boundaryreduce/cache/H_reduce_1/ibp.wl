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
maxdots=0;
GetSectors=Automatic;
cutposi=Position[{0, 0},1]//Flatten;
learn={1 -> {0, 2}, 0 -> {0, 2}, -1 -> {0, 2}, _ -> {0, 2}};
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

exints = {};
Export["ibps/exints_def.m",exints];
Export["ibps/ids_H_exints.mx",(-#[[1]]+#[[2]])&/@exints,"MX"];
Export["ibps/ints_H_exints.mx",IntegralsIn[exints],"MX"];
usints = {};
Export["ibps/usints_def.m",usints];
Export["ibps/ids_H_usints.mx",(-#[[1]]+#[[2]])&/@usints,"MX"];
Export["ibps/ints_H_usints.mx",IntegralsIn[usints],"MX"];


If[!({eps}==={}),SerializeFastIds[FileNames["ibps/ids_H_*.mx"], {}, {eps}, "UserDefinedInts"->(First/@usints),"ExtraInts"->(First/@exints)]];
];
CloseKernels[];Pause[0.1];
Quit[];