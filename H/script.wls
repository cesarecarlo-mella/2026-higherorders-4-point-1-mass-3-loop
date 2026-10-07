SetDirectory["/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/NPL/DEQ_solutions/SLC/Solutions/I3"]


letters = {Log[-1 + x], Log[x], Log[-1 + y], Log[y], 
   Log[-(1 - x)^2 + y], Log[-1 + x + y], Log[x + y], Log[-x + x^2 + y],
    Log[x^2 + y + x y], Log[-1 + (1 + x) y], Log[1 + (-2 + x) y + y^2],
    Log[1 + x^2 + x (-2 + y)]};


sqroot = {Log[(x (1 - x - y) - I Sqrt[x (1 - x - y) y])/(
   x (1 - x - y) + I Sqrt[x (1 - x - y) y])], 
  Log[(x y - I Sqrt[x (1 - x - y) y])/(x y + I Sqrt[x (1 - x - y) y])]}
DSrules = {Log[(x (1 - x - y) - I Sqrt[x (1 - x - y) y])/(
     x (1 - x - y) + I Sqrt[x (1 - x - y) y])] -> (
    sqrt I y (-1 + 2 x + y))/(x y (1 - x - y) (-1 + x) (x + y)), 
   Log[(x y - I Sqrt[x (1 - x - y) y])/(
     x y + I Sqrt[x (1 - x - y) y])] -> (sqrt I y)/(
    x y (1 - x - y) (1 - x))};
DTrules = {Log[(x (1 - x - y) - I Sqrt[x (1 - x - y) y])/(
     x (1 - x - y) + I Sqrt[x (1 - x - y) y])] -> -((I x sqrt)/(
     x y (1 - x - y) (x + y))), 
   Log[(x y - I Sqrt[x (1 - x - y) y])/(
     x y + I Sqrt[x (1 - x - y) y])] -> (
    I x sqrt)/((1 - y) x y (1 - x - y))};

derLx = Thread[
   Rule[{l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], 
     l[11], l[12]}, D[letters, x]]];
derLy = Thread[
   Rule[{l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], 
     l[11], l[12]}, D[letters, y]]];
derLxroot =  {l[13] -> (sqrt I y (-1 + 2 x + y))/(
    x y (1 - x - y) (-1 + x) (x + y)), 
   l[14] -> (sqrt I y)/(x y (1 - x - y) (1 - x))};
derLyroot =  {l[13] -> -((I x sqrt)/(x y (1 - x - y) (x + y))), 
   l[14] -> (I x sqrt)/((1 - y) x y (1 - x - y))};


matrS = Import["./things/matS_apart.m"] /. {s -> x, t -> y};
matrT = Import["./things/matT_apart.m"] /. {s -> x, t -> y};


ansatz = Table[
   Sum[a[i, j][k] l[k], {k, 1, 14}], {i, 1, 333}, {j, 1, 333}];
counter = 0;
Do[
 tmp = Collect[
   Numerator[
    Together[(ansatz[[i, j]] /. derLx /. derLxroot) - 
      matrS[[i, j]]]], {x, y, sqrt}, FF];
 eqs = Union@Cases[tmp, _FF, \[Infinity]] /. 
   FF[a__] :> Thread[Equal[a, 0]];
 tmp = Collect[
   Numerator[
    Together[(ansatz[[i, j]] /. derLy /. derLyroot) - 
      matrT[[i, j]]]], {x, y, sqrt}, FF];
 eqs = Join[eqs, 
   Union@Cases[tmp, _FF, \[Infinity]] /. FF[a__] :> Thread[Equal[a, 0]]];
 sol = Flatten[Solve[eqs]];
 ansatz[[i, j]] = ansatz[[i, j]] /. sol;
 counter = counter + 1; 
 If[Mod[counter, 100] === 0, Print[i, "  ,   ", j]], {i, 1, 333}, {j, 
  1, 333}]

ansatz = Export["tmpansatz.m", ansatz];
