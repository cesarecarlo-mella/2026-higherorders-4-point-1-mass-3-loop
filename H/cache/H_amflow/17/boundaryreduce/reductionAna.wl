If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
Get["/space/ttpmath/tools_src/blade/Blade.wl"];
Get["topology.wl"];


target = Get["target"];
preferred = Get["preferred"];


amfint[blint_]:=j[blint[[1]], Sequence@@blint[[2]]];
res = BLReduce[target, preferred]/.BLFamilyInf[family]["MastersRules"];
masters = Cases[res, _BL, Infinity]//DeleteDuplicates;
res = Coefficient[#, masters]&/@res;
res = Join[res, IdentityMatrix[Length[masters]]];
Put[{amfint/@masters, Thread[amfint/@Join[target, masters] -> res]}, "results/table"];


Quit[];