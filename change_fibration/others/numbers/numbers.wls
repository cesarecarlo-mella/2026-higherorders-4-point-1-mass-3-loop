#!/usr/bin/env wolframscript
(* ::Package:: *)

GetEnvironment["PATH"][[2]]
SetEnvironment["PATH"->GetEnvironment["PATH"][[2]]<>":/usr/local/bin"]
<<PolyLogTools`
tofix = Import["./tofix.m"];
todo = Import["./w5.m"];
numbers = Table[nothin, {i,1,Length[todo]}];
Do[numbers[[i]] = Ginsh[todo[[i]]- tofix[[i]], {delta -> 3/14 + I Log[2] 10^(-30),rho -> 5/13 + I Log[5] 10^(-30)}, PrecisionGoal-> 50]; Print[i]; If[Mod[i,20]===0, Export["./numbers.m", numbers]],{i,1,Length[todo]}]
Export["./numbers.m", numbers];
