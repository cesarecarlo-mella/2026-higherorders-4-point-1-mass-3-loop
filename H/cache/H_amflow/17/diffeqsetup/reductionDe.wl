If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
Get["/space/ttpmath/tools_src/blade/Blade.wl"];
Get["topology.wl"];


preferred = Get["preferred"];


amfint[blint_]:=j[blint[[1]], Sequence@@blint[[2]]];
{masters, invs, diffeq} = BLDifferentialEquation[{eta}, preferred];
Put[{amfint/@masters, invs, diffeq}, "results/diffeq"];


Quit[];