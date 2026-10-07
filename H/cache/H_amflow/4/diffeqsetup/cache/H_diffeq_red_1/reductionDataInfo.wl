AppendTo[$Path,"/space/ttpmath/tools_src/blade/BLAddOns/BladeIBP"];
AppendTo[$Path,"/space/ttpmath/tools_src/finiteflow/mathlink"];
AppendTo[$LibraryPath,"/space/ttpmath/tools_src/finiteflow"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
If[!DirectoryQ["ibps"], Put[{{}, {}}, "results/table"]; Abort[]];

<<BladeIBP.m;
<<FiniteFlow.m
<<topology.wl;

IntegralOrdering = 1;


Module[
{resdir,usints,exints,all,target,eqsfiles,vars,eqs,ibplearn,mis,nonmis,indep,graph,in,ibps,depvarsorder,
time,pid,ps,eva,nonzero,sparseinfo, alleqs, fflearn, sol, rule,mem1,mem2},

resdir = "results";
If[!DirectoryQ@resdir,CreateDirectory[resdir]];

usints = Import["ibps/usints_def.m"];
exints = Import["ibps/exints_def.m"];
all = Join[First/@usints,GetAllInts[FileNames["ibps/ints_H_*.mx"]], First/@exints];
target = Get["target"];
If[target === 0,
target = Select[all, IntRank[#]<=Blade`Trimming`IBPRank && IntDots[#]<=Blade`Trimming`IBPDot &],
target = Select[EliminateZeroSectors[target], #=!=0&]];
If[target === {}, Print["target is empty."];Put[{{},{}}, "results/table"];Abort[]];

eqsfiles = FileNames["ibps/sids_*.json"];
vars = {eta, eps};

If[!(vars==={}),WriteSystemJSON[eqsfiles, all, target, vars, "FileName"->"system.json"]];

mem1=MemoryAvailable[];

FFNewGraph[graph,in,vars];
If[!(vars==={}),
FFAlgJSONSparseSolver[graph,ibps,{in},"system.json"],
eqsfiles = FileNames["ibps/ids_*.mx"];
eqs=Join@@((Import[#]/.{s -> -1, t -> -3, pp4 -> -15})&/@eqsfiles);
eqs=Thread[eqs==0];
eqs=DeleteCases[eqs,True];
FFAlgSparseSolver[graph,ibps,{in},{},eqs,all,"NeededVars"->target]];
FFSolverOnlyHomogeneous[graph,ibps];
FFGraphOutput[graph,ibps];
ibplearn = FFSparseSolverLearn[graph,all];
{nonmis, mis} = {"DepVars", "IndepVars"}/.ibplearn;
Print["number of integrals" -> Length/@{nonmis,mis}];
If[Length[mis]==0,Print["target is zero"];
Put[If[!(vars==={}),Range[Total[Import/@(StringReplace[#1,{"sids_"->"nids_",".json"->".m"}]&)/@eqsfiles]],Length[eqs]],FileNameJoin@{resdir, "indep"}];
Put[{nonmis,mis,MappedSectors[H],UniqueSectors[H],ConstantArray[{},Length[nonmis]]}, FileNameJoin@{resdir, "datainfo"}];
Abort[]];

FFSparseSolverMarkAndSweepEqs[graph,ibps];
mem2=MemoryAvailable[];
Print["memory used -> ",(mem1-mem2)/Power[1024,3]//N,"GB"];
FFSparseSolverDeleteUnneededEqs[graph,ibps];
Print["number of independent eqs."->FFSolverNIndepEqs[graph,ibps]];
indep = FFSolverIndepEqs[graph,ibps];
Put[indep,FileNameJoin@{resdir, "indep"}];

depvarsorder=FFNonZeroesGetInfo[graph,ibps];
Put[depvarsorder,FileNameJoin@{resdir, "depvarsorder"}];

time=AbsoluteTiming[FFGraphEvaluateMany[graph,RandomInteger[{1,FFPrimeNo[200]-1},{3,Length@vars}],"NThreads"->1];];
Print["average sample time for ibps (single core) ->"<>ToString[time[[1]]/3//N]];

pid=0;
ps = RandomInteger[{1, FFPrimeNo[pid]}, Length[vars]];
eva = FFGraphEvaluate[graph, ps, "PrimeNo"->pid];
eva = Join[Partition[eva, Length[mis]], IdentityMatrix[Length[mis]]];
nonzero[l0_]:=Position[l0, Except[0], Heads->False]//Flatten;
sparseinfo = nonzero/@eva;

Put[{FFPrimeNo[pid], eva}, FileNameJoin@{resdir, "data"}];
Put[{nonmis,mis,MappedSectors[H],UniqueSectors[H],sparseinfo}, FileNameJoin@{resdir, "datainfo"}];
Put[usints,FileNameJoin@{resdir, "usints"}];
Put[exints,FileNameJoin@{resdir,"exints"}];

FFDeleteGraph[graph];];
Quit[];