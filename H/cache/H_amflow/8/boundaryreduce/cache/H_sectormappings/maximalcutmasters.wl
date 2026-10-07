AppendTo[$Path,"/space/ttpmath/tools_src/blade/BLAddOns/BladeIBP"];
AppendTo[$Path,"/space/ttpmath/tools_src/finiteflow/mathlink"];
AppendTo[$LibraryPath,"/space/ttpmath/tools_src/finiteflow"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<BladeIBP.m;
<<FiniteFlow.m;
<<topology.wl;


IntegralOrdering=1;
family=H;


Block[{mapdir},
mapdir = "results";
If[!DirectoryQ@mapdir,CreateDirectory[mapdir]];
Put[MappedSectors[family],FileNameJoin[{mapdir,"mappedsectors"}]];
Put[UniqueSectors[family],FileNameJoin[{mapdir,"uniquesectors"}]];
Put[ZeroSectors[family],FileNameJoin[{mapdir,"zerosectors"}]];
];


CloseKernels[];Pause[0.1];


Quit[];