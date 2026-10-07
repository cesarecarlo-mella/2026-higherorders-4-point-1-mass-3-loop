(* Created with the Wolfram Language : www.wolfram.com *)
{G[-1, -1 + \[Delta]] -> G[0, \[Delta]], G[-1, -1, -1 + \[Delta]] -> 
  G[0, 0, \[Delta]], G[-1, 0, -1 + \[Delta]] -> 
  Pi^2/6 - I*Pi*G[0, \[Delta]] + G[0, 1, \[Delta]], 
 G[-1, 1, -1 + \[Delta]] -> Pi^2/12 + G[0, 2]*G[0, \[Delta]] - G[0, 0, 2] + 
   G[0, 2, \[Delta]], G[-1, -1, -1, -1 + \[Delta]] -> G[0, 0, 0, \[Delta]], 
 G[-1, -1, 0, -1 + \[Delta]] -> (Pi^2*G[0, \[Delta]])/6 - 
   I*Pi*G[0, 0, \[Delta]] + G[0, 0, 1, \[Delta]] + Zeta[3], 
 G[-1, 0, -1, -1 + \[Delta]] -> -1/6*(Pi^2*G[0, \[Delta]]) + 
   G[0, 1, 0, \[Delta]] - 2*Zeta[3], G[-1, 0, 0, -1 + \[Delta]] -> 
  (-1/6*I)*Pi^3 - (Pi^2*G[0, \[Delta]])/2 - I*Pi*G[0, 1, \[Delta]] + 
   G[0, 1, 1, \[Delta]] - Zeta[3], G[-1, 0, 1, -1 + \[Delta]] -> 
  (Pi^2*G[0, \[Delta]])/12 + G[0, 2]*G[0, 1, \[Delta]] + 
   G[0, 1, 2, \[Delta]] + (5*Zeta[3])/8, G[-1, 1, 0, -1 + \[Delta]] -> 
  (-1/12*I)*Pi^3 + Pi^2*(G[0, 2]/4 - G[0, \[Delta]]/12) + 
   I*Pi*(-(G[0, 2]*G[0, \[Delta]]) + G[0, 0, 2] - G[0, 2, \[Delta]]) + 
   G[0, 2, 1, \[Delta]] - (13*Zeta[3])/8, G[-1, -1, -1, -1, -1 + \[Delta]] -> 
  G[0, 0, 0, 0, \[Delta]], G[-1, -1, -1, 0, -1 + \[Delta]] -> 
  Pi^4/90 + (Pi^2*G[0, 0, \[Delta]])/6 - I*Pi*G[0, 0, 0, \[Delta]] + 
   G[0, 0, 0, 1, \[Delta]] + G[0, \[Delta]]*Zeta[3], 
 G[-1, -1, 0, -1, -1 + \[Delta]] -> -1/30*Pi^4 - (Pi^2*G[0, 0, \[Delta]])/6 + 
   G[0, 0, 1, 0, \[Delta]] - 2*G[0, \[Delta]]*Zeta[3], 
 G[-1, -1, 0, 0, -1 + \[Delta]] -> -1/360*Pi^4 - (I/6)*Pi^3*G[0, \[Delta]] - 
   (Pi^2*G[0, 0, \[Delta]])/2 - I*Pi*G[0, 0, 1, \[Delta]] + 
   G[0, 0, 1, 1, \[Delta]] - I*Pi*Zeta[3] - G[0, \[Delta]]*Zeta[3], 
 G[-1, -1, 0, 1, -1 + \[Delta]] -> -1/288*Pi^4 + 
   (Pi^2*G[0, 0, \[Delta]])/12 + G[0, 2]*G[0, 0, 1, \[Delta]] + 
   G[0, 0, 1, 2, \[Delta]] - (Pi^2*Log[2]^2)/24 + Log[2]^4/24 + 
   PolyLog[4, 1/2] + (5*G[0, \[Delta]]*Zeta[3])/8 + (7*Log[2]*Zeta[3])/8, 
 G[-1, 0, -1, -1, -1 + \[Delta]] -> Pi^4/30 + G[0, 1, 0, 0, \[Delta]] + 
   G[0, \[Delta]]*Zeta[3], G[-1, 0, -1, 0, -1 + \[Delta]] -> 
  (7*Pi^4)/360 + (I/6)*Pi^3*G[0, \[Delta]] + (Pi^2*G[0, 1, \[Delta]])/6 - 
   I*Pi*G[0, 1, 0, \[Delta]] + G[0, 1, 0, 1, \[Delta]] + (2*I)*Pi*Zeta[3] + 
   2*G[0, \[Delta]]*Zeta[3], G[-1, 0, 0, -1, -1 + \[Delta]] -> 
  -1/72*Pi^4 - (Pi^2*G[0, 1, \[Delta]])/6 + G[0, 1, 1, 0, \[Delta]] - 
   G[0, \[Delta]]*Zeta[3], G[-1, 0, 0, 0, -1 + \[Delta]] -> 
  (-13*Pi^4)/180 + (I/6)*Pi^3*G[0, \[Delta]] - (Pi^2*G[0, 1, \[Delta]])/2 - 
   I*Pi*G[0, 1, 1, \[Delta]] + G[0, 1, 1, 1, \[Delta]] + I*Pi*Zeta[3], 
 G[-1, 0, 0, 1, -1 + \[Delta]] -> (11*Pi^4)/360 + 
   (Pi^2*G[0, 1, \[Delta]])/12 + G[0, 2]*G[0, 1, 1, \[Delta]] + 
   G[0, 1, 1, 2, \[Delta]] + (Pi^2*Log[2]^2)/12 - Log[2]^4/12 - 
   2*PolyLog[4, 1/2] + (3*G[0, \[Delta]]*Zeta[3])/4 - (7*Log[2]*Zeta[3])/4, 
 G[-1, 0, 1, 0, -1 + \[Delta]] -> (-17*Pi^4)/1440 - 
   (I/12)*Pi^3*G[0, \[Delta]] - (Pi^2*G[0, 1, \[Delta]])/12 + 
   I*Pi*(-(G[0, 2]*G[0, 1, \[Delta]]) - G[0, 1, 2, \[Delta]]) + 
   G[0, 1, 2, 1, \[Delta]] - ((5*I)/8)*Pi*Zeta[3] - 
   (3*G[0, \[Delta]]*Zeta[3])/2, G[-1, 0, 1, 1, -1 + \[Delta]] -> 
  Pi^4/480 + G[0, 0, 2]*G[0, 1, \[Delta]] + G[0, 2]*G[0, 1, 2, \[Delta]] + 
   G[0, 1, 2, 2, \[Delta]] + (G[0, \[Delta]]*Zeta[3])/8, 
 G[-1, 1, 0, 0, -1 + \[Delta]] -> (-17*Pi^4)/360 + 
   I*Pi^3*(-1/4*G[0, 2] + G[0, \[Delta]]/12) + 
   Pi^2*(-1/2*(G[0, 2]*G[0, \[Delta]]) + G[0, 0, 2]/2 - 
     G[0, 2, \[Delta]]/2) - I*Pi*G[0, 2, 1, \[Delta]] + 
   G[0, 2, 1, 1, \[Delta]] - (Pi^2*Log[2]^2)/12 + Log[2]^4/12 + 
   2*PolyLog[4, 1/2] + ((13*I)/8)*Pi*Zeta[3] + 
   ((-7*G[0, 2])/4 + (3*G[0, \[Delta]])/4)*Zeta[3] + (7*Log[2]*Zeta[3])/4, 
 G[-1, 1, 0, 1, -1 + \[Delta]] -> 
  Pi^4/24 + Pi^2*((G[0, 2]*G[0, \[Delta]])/12 - G[0, 0, 2]/12 + 
     G[0, 2, \[Delta]]/12) + G[0, 2]*G[0, 2, 1, \[Delta]] + 
   G[0, 2, 1, 2, \[Delta]] + (Pi^2*Log[2]^2)/6 - Log[2]^4/6 - 
   4*PolyLog[4, 1/2] + ((7*G[0, 2])/8 - G[0, \[Delta]]/4)*Zeta[3] - 
   (7*Log[2]*Zeta[3])/2, G[-1, -1, -1, -1, -1, -1 + \[Delta]] -> 
  G[0, 0, 0, 0, 0, \[Delta]], G[-1, -1, -1, -1, 0, -1 + \[Delta]] -> 
  (Pi^4*G[0, \[Delta]])/90 + (Pi^2*G[0, 0, 0, \[Delta]])/6 - 
   I*Pi*G[0, 0, 0, 0, \[Delta]] + G[0, 0, 0, 0, 1, \[Delta]] + 
   G[0, 0, \[Delta]]*Zeta[3] + Zeta[5], 
 G[-1, -1, -1, 0, -1, -1 + \[Delta]] -> -1/30*(Pi^4*G[0, \[Delta]]) - 
   (Pi^2*G[0, 0, 0, \[Delta]])/6 + G[0, 0, 0, 1, 0, \[Delta]] - 
   2*G[0, 0, \[Delta]]*Zeta[3] - 4*Zeta[5], 
 G[-1, -1, -1, 0, 0, -1 + \[Delta]] -> (-1/90*I)*Pi^5 - 
   (Pi^4*G[0, \[Delta]])/360 - (I/6)*Pi^3*G[0, 0, \[Delta]] - 
   (Pi^2*G[0, 0, 0, \[Delta]])/2 - I*Pi*G[0, 0, 0, 1, \[Delta]] + 
   G[0, 0, 0, 1, 1, \[Delta]] + (Pi^2*Zeta[3])/6 - 
   I*Pi*G[0, \[Delta]]*Zeta[3] - G[0, 0, \[Delta]]*Zeta[3] - 2*Zeta[5], 
 G[-1, -1, 0, -1, -1, -1 + \[Delta]] -> (Pi^4*G[0, \[Delta]])/30 + 
   G[0, 0, 1, 0, 0, \[Delta]] + G[0, 0, \[Delta]]*Zeta[3] + 6*Zeta[5], 
 G[-1, -1, 0, -1, 0, -1 + \[Delta]] -> 
  (I/30)*Pi^5 + (7*Pi^4*G[0, \[Delta]])/360 + (I/6)*Pi^3*G[0, 0, \[Delta]] + 
   (Pi^2*G[0, 0, 1, \[Delta]])/6 - I*Pi*G[0, 0, 1, 0, \[Delta]] + 
   G[0, 0, 1, 0, 1, \[Delta]] - (Pi^2*Zeta[3])/3 + 
   (2*I)*Pi*G[0, \[Delta]]*Zeta[3] + 2*G[0, 0, \[Delta]]*Zeta[3] + 
   (11*Zeta[5])/2, G[-1, -1, 0, 0, -1, -1 + \[Delta]] -> 
  -1/72*(Pi^4*G[0, \[Delta]]) - (Pi^2*G[0, 0, 1, \[Delta]])/6 + 
   G[0, 0, 1, 1, 0, \[Delta]] - (Pi^2*Zeta[3])/6 - 
   G[0, 0, \[Delta]]*Zeta[3] + Zeta[5]/2, 
 G[-1, -1, 0, 0, 0, -1 + \[Delta]] -> (I/360)*Pi^5 - 
   (13*Pi^4*G[0, \[Delta]])/180 + (I/6)*Pi^3*G[0, 0, \[Delta]] - 
   (Pi^2*G[0, 0, 1, \[Delta]])/2 - I*Pi*G[0, 0, 1, 1, \[Delta]] + 
   G[0, 0, 1, 1, 1, \[Delta]] - (2*Pi^2*Zeta[3])/3 + 
   I*Pi*G[0, \[Delta]]*Zeta[3] + 2*Zeta[5], 
 G[-1, -1, 0, 1, 0, -1 + \[Delta]] -> (-17*Pi^4*G[0, \[Delta]])/1440 - 
   (I/12)*Pi^3*G[0, 0, \[Delta]] - (Pi^2*G[0, 0, 1, \[Delta]])/12 + 
   I*Pi*(-(G[0, 2]*G[0, 0, 1, \[Delta]]) - G[0, 0, 1, 2, \[Delta]]) + 
   G[0, 0, 1, 2, 1, \[Delta]] + (Pi^2*Zeta[3])/32 - 
   ((5*I)/8)*Pi*G[0, \[Delta]]*Zeta[3] - (3*G[0, 0, \[Delta]]*Zeta[3])/2 + 
   I*(Pi^5/288 + (Pi^3*Log[2]^2)/24 - (Pi*Log[2]^4)/24 - Pi*PolyLog[4, 1/2] - 
     (7*Pi*Log[2]*Zeta[3])/8) - (41*Zeta[5])/32, 
 G[-1, 0, -1, -1, -1, -1 + \[Delta]] -> -1/90*(Pi^4*G[0, \[Delta]]) + 
   G[0, 1, 0, 0, 0, \[Delta]] - 4*Zeta[5], 
 G[-1, 0, -1, -1, 0, -1 + \[Delta]] -> (-1/30*I)*Pi^5 - 
   (Pi^4*G[0, \[Delta]])/72 + (Pi^2*G[0, 1, 0, \[Delta]])/6 - 
   I*Pi*G[0, 1, 0, 0, \[Delta]] + G[0, 1, 0, 0, 1, \[Delta]] + 
   (Pi^2*Zeta[3])/6 - I*Pi*G[0, \[Delta]]*Zeta[3] + 
   G[0, 1, \[Delta]]*Zeta[3] - (9*Zeta[5])/2, 
 G[-1, 0, -1, 0, -1, -1 + \[Delta]] -> (Pi^4*G[0, \[Delta]])/120 - 
   (Pi^2*G[0, 1, 0, \[Delta]])/6 + G[0, 1, 0, 1, 0, \[Delta]] + 
   (Pi^2*Zeta[3])/3 - 2*G[0, 1, \[Delta]]*Zeta[3] - 2*Zeta[5], 
 G[-1, 0, -1, 0, 0, -1 + \[Delta]] -> ((-7*I)/360)*Pi^5 + 
   (Pi^4*G[0, \[Delta]])/20 - (I/6)*Pi^3*G[0, 1, \[Delta]] - 
   (Pi^2*G[0, 1, 0, \[Delta]])/2 - I*Pi*G[0, 1, 0, 1, \[Delta]] + 
   G[0, 1, 0, 1, 1, \[Delta]] + (4*Pi^2*Zeta[3])/3 - 
   (2*I)*Pi*G[0, \[Delta]]*Zeta[3] - G[0, 1, \[Delta]]*Zeta[3] - 
   (11*Zeta[5])/2, G[-1, 0, 0, -1, -1, -1 + \[Delta]] -> 
  (Pi^4*G[0, \[Delta]])/360 + G[0, 1, 1, 0, 0, \[Delta]] + 
   G[0, 1, \[Delta]]*Zeta[3] + Zeta[5]/2, 
 G[-1, 0, 0, -1, 0, -1 + \[Delta]] -> (I/72)*Pi^5 + 
   (Pi^4*G[0, \[Delta]])/30 + (I/6)*Pi^3*G[0, 1, \[Delta]] + 
   (Pi^2*G[0, 1, 1, \[Delta]])/6 - I*Pi*G[0, 1, 1, 0, \[Delta]] + 
   G[0, 1, 1, 0, 1, \[Delta]] - (Pi^2*Zeta[3])/6 + 
   I*Pi*G[0, \[Delta]]*Zeta[3] + 2*G[0, 1, \[Delta]]*Zeta[3] + (9*Zeta[5])/2, 
 G[-1, 0, 0, 0, -1, -1 + \[Delta]] -> -1/90*(Pi^4*G[0, \[Delta]]) - 
   (Pi^2*G[0, 1, 1, \[Delta]])/6 + G[0, 1, 1, 1, 0, \[Delta]] + 
   (Pi^2*Zeta[3])/6 - G[0, 1, \[Delta]]*Zeta[3] - 3*Zeta[5], 
 G[-1, 0, 0, 0, 0, -1 + \[Delta]] -> (I/60)*Pi^5 + (Pi^4*G[0, \[Delta]])/24 + 
   (I/6)*Pi^3*G[0, 1, \[Delta]] - (Pi^2*G[0, 1, 1, \[Delta]])/2 - 
   I*Pi*G[0, 1, 1, 1, \[Delta]] + G[0, 1, 1, 1, 1, \[Delta]] + 
   (Pi^2*Zeta[3])/2 - Zeta[5], G[-1, 0, 0, 1, 0, -1 + \[Delta]] -> 
  (-7*Pi^4*G[0, \[Delta]])/240 - (I/12)*Pi^3*G[0, 1, \[Delta]] - 
   (Pi^2*G[0, 1, 1, \[Delta]])/12 + I*Pi*(-(G[0, 2]*G[0, 1, 1, \[Delta]]) - 
     G[0, 1, 1, 2, \[Delta]]) + G[0, 1, 1, 2, 1, \[Delta]] + 
   (13*Pi^2*Zeta[3])/48 - ((3*I)/4)*Pi*G[0, \[Delta]]*Zeta[3] - 
   (3*G[0, 1, \[Delta]]*Zeta[3])/2 + I*((-11*Pi^5)/360 - (Pi^3*Log[2]^2)/12 + 
     (Pi*Log[2]^4)/12 + 2*Pi*PolyLog[4, 1/2] + (7*Pi*Log[2]*Zeta[3])/4) - 
   (83*Zeta[5])/16, G[-1, 0, 1, 0, -1, -1 + \[Delta]] -> 
  -1/480*(Pi^4*G[0, \[Delta]]) + Pi^2*(-1/6*(G[0, 2]*G[0, 1, \[Delta]]) - 
     G[0, 1, 2, \[Delta]]/6) + G[0, 1, 2, 1, 0, \[Delta]] - 
   (5*Pi^2*Zeta[3])/48 + (5*G[0, 1, \[Delta]]*Zeta[3])/8 + (29*Zeta[5])/32, 
 G[-1, 0, 1, 0, 0, -1 + \[Delta]] -> ((17*I)/1440)*Pi^5 - 
   (Pi^4*G[0, \[Delta]])/80 + (I/12)*Pi^3*G[0, 1, \[Delta]] + 
   Pi^2*(-1/2*(G[0, 2]*G[0, 1, \[Delta]]) - G[0, 1, 2, \[Delta]]/2) - 
   I*Pi*G[0, 1, 2, 1, \[Delta]] + G[0, 1, 2, 1, 1, \[Delta]] - 
   (9*Pi^2*Zeta[3])/16 + ((3*I)/2)*Pi*G[0, \[Delta]]*Zeta[3] + 
   (3*G[0, 1, \[Delta]]*Zeta[3])/4 + (67*Zeta[5])/16, 
 G[-1, 0, 1, 1, 0, -1 + \[Delta]] -> (-1/480*I)*Pi^5 - 
   (Pi^4*G[0, \[Delta]])/288 + Pi^2*(-1/12*(G[0, 2]*G[0, 1, \[Delta]]) - 
     G[0, 1, 2, \[Delta]]/12) + I*Pi*(-(G[0, 0, 2]*G[0, 1, \[Delta]]) - 
     G[0, 2]*G[0, 1, 2, \[Delta]] - G[0, 1, 2, 2, \[Delta]]) + 
   G[0, 1, 2, 2, 1, \[Delta]] + (Pi^2*Zeta[3])/48 - 
   (I/8)*Pi*G[0, \[Delta]]*Zeta[3] + (G[0, 1, \[Delta]]*Zeta[3])/8 - 
   (5*Zeta[5])/8, G[-1, -1, -1, -1, 0, -1 + \[Delta]] -> 
  (Pi^4*G[0, \[Delta]])/90 + (Pi^2*G[0, 0, 0, \[Delta]])/6 - 
   I*Pi*G[0, 0, 0, 0, \[Delta]] + G[0, 0, 0, 0, 1, \[Delta]] + 
   G[0, 0, \[Delta]]*Zeta[3] + Zeta[5], G[-1, -1, -1, 0, 0, -1 + \[Delta]] -> 
  (-1/90*I)*Pi^5 - (Pi^4*G[0, \[Delta]])/360 - (I/6)*Pi^3*G[0, 0, \[Delta]] - 
   (Pi^2*G[0, 0, 0, \[Delta]])/2 - I*Pi*G[0, 0, 0, 1, \[Delta]] + 
   G[0, 0, 0, 1, 1, \[Delta]] + (Pi^2*Zeta[3])/6 - 
   I*Pi*G[0, \[Delta]]*Zeta[3] - G[0, 0, \[Delta]]*Zeta[3] - 2*Zeta[5], 
 G[-1, -1, 0, -1, 0, -1 + \[Delta]] -> 
  (I/30)*Pi^5 + (7*Pi^4*G[0, \[Delta]])/360 + (I/6)*Pi^3*G[0, 0, \[Delta]] + 
   (Pi^2*G[0, 0, 1, \[Delta]])/6 - I*Pi*G[0, 0, 1, 0, \[Delta]] + 
   G[0, 0, 1, 0, 1, \[Delta]] - (Pi^2*Zeta[3])/3 + 
   (2*I)*Pi*G[0, \[Delta]]*Zeta[3] + 2*G[0, 0, \[Delta]]*Zeta[3] + 
   (11*Zeta[5])/2, G[-1, -1, 0, 0, 0, -1 + \[Delta]] -> 
  (I/360)*Pi^5 - (13*Pi^4*G[0, \[Delta]])/180 + 
   (I/6)*Pi^3*G[0, 0, \[Delta]] - (Pi^2*G[0, 0, 1, \[Delta]])/2 - 
   I*Pi*G[0, 0, 1, 1, \[Delta]] + G[0, 0, 1, 1, 1, \[Delta]] - 
   (2*Pi^2*Zeta[3])/3 + I*Pi*G[0, \[Delta]]*Zeta[3] + 2*Zeta[5], 
 G[-1, -1, 0, 1, 0, -1 + \[Delta]] -> (-17*Pi^4*G[0, \[Delta]])/1440 - 
   (I/12)*Pi^3*G[0, 0, \[Delta]] - (Pi^2*G[0, 0, 1, \[Delta]])/12 + 
   I*Pi*(-(G[0, 2]*G[0, 0, 1, \[Delta]]) - G[0, 0, 1, 2, \[Delta]]) + 
   G[0, 0, 1, 2, 1, \[Delta]] + (Pi^2*Zeta[3])/32 - 
   ((5*I)/8)*Pi*G[0, \[Delta]]*Zeta[3] - (3*G[0, 0, \[Delta]]*Zeta[3])/2 + 
   I*(Pi^5/288 + (Pi^3*Log[2]^2)/24 - (Pi*Log[2]^4)/24 - Pi*PolyLog[4, 1/2] - 
     (7*Pi*Log[2]*Zeta[3])/8) - (41*Zeta[5])/32, 
 G[-1, 0, -1, -1, 0, -1 + \[Delta]] -> (-1/30*I)*Pi^5 - 
   (Pi^4*G[0, \[Delta]])/72 + (Pi^2*G[0, 1, 0, \[Delta]])/6 - 
   I*Pi*G[0, 1, 0, 0, \[Delta]] + G[0, 1, 0, 0, 1, \[Delta]] + 
   (Pi^2*Zeta[3])/6 - I*Pi*G[0, \[Delta]]*Zeta[3] + 
   G[0, 1, \[Delta]]*Zeta[3] - (9*Zeta[5])/2, 
 G[-1, 0, -1, 0, 0, -1 + \[Delta]] -> ((-7*I)/360)*Pi^5 + 
   (Pi^4*G[0, \[Delta]])/20 - (I/6)*Pi^3*G[0, 1, \[Delta]] - 
   (Pi^2*G[0, 1, 0, \[Delta]])/2 - I*Pi*G[0, 1, 0, 1, \[Delta]] + 
   G[0, 1, 0, 1, 1, \[Delta]] + (4*Pi^2*Zeta[3])/3 - 
   (2*I)*Pi*G[0, \[Delta]]*Zeta[3] - G[0, 1, \[Delta]]*Zeta[3] - 
   (11*Zeta[5])/2, G[-1, 0, 0, -1, 0, -1 + \[Delta]] -> 
  (I/72)*Pi^5 + (Pi^4*G[0, \[Delta]])/30 + (I/6)*Pi^3*G[0, 1, \[Delta]] + 
   (Pi^2*G[0, 1, 1, \[Delta]])/6 - I*Pi*G[0, 1, 1, 0, \[Delta]] + 
   G[0, 1, 1, 0, 1, \[Delta]] - (Pi^2*Zeta[3])/6 + 
   I*Pi*G[0, \[Delta]]*Zeta[3] + 2*G[0, 1, \[Delta]]*Zeta[3] + (9*Zeta[5])/2, 
 G[-1, 0, 0, 0, 0, -1 + \[Delta]] -> (I/60)*Pi^5 + (Pi^4*G[0, \[Delta]])/24 + 
   (I/6)*Pi^3*G[0, 1, \[Delta]] - (Pi^2*G[0, 1, 1, \[Delta]])/2 - 
   I*Pi*G[0, 1, 1, 1, \[Delta]] + G[0, 1, 1, 1, 1, \[Delta]] + 
   (Pi^2*Zeta[3])/2 - Zeta[5], G[-1, 0, 0, 1, 0, -1 + \[Delta]] -> 
  (-7*Pi^4*G[0, \[Delta]])/240 - (I/12)*Pi^3*G[0, 1, \[Delta]] - 
   (Pi^2*G[0, 1, 1, \[Delta]])/12 + I*Pi*(-(G[0, 2]*G[0, 1, 1, \[Delta]]) - 
     G[0, 1, 1, 2, \[Delta]]) + G[0, 1, 1, 2, 1, \[Delta]] + 
   (13*Pi^2*Zeta[3])/48 - ((3*I)/4)*Pi*G[0, \[Delta]]*Zeta[3] - 
   (3*G[0, 1, \[Delta]]*Zeta[3])/2 + I*((-11*Pi^5)/360 - (Pi^3*Log[2]^2)/12 + 
     (Pi*Log[2]^4)/12 + 2*Pi*PolyLog[4, 1/2] + (7*Pi*Log[2]*Zeta[3])/4) - 
   (83*Zeta[5])/16, G[-1, 0, 1, 0, 0, -1 + \[Delta]] -> 
  ((17*I)/1440)*Pi^5 - (Pi^4*G[0, \[Delta]])/80 + 
   (I/12)*Pi^3*G[0, 1, \[Delta]] + Pi^2*(-1/2*(G[0, 2]*G[0, 1, \[Delta]]) - 
     G[0, 1, 2, \[Delta]]/2) - I*Pi*G[0, 1, 2, 1, \[Delta]] + 
   G[0, 1, 2, 1, 1, \[Delta]] - (9*Pi^2*Zeta[3])/16 + 
   ((3*I)/2)*Pi*G[0, \[Delta]]*Zeta[3] + (3*G[0, 1, \[Delta]]*Zeta[3])/4 + 
   (67*Zeta[5])/16, G[-1, 0, 1, 1, 0, -1 + \[Delta]] -> 
  (-1/480*I)*Pi^5 - (Pi^4*G[0, \[Delta]])/288 + 
   Pi^2*(-1/12*(G[0, 2]*G[0, 1, \[Delta]]) - G[0, 1, 2, \[Delta]]/12) + 
   I*Pi*(-(G[0, 0, 2]*G[0, 1, \[Delta]]) - G[0, 2]*G[0, 1, 2, \[Delta]] - 
     G[0, 1, 2, 2, \[Delta]]) + G[0, 1, 2, 2, 1, \[Delta]] + 
   (Pi^2*Zeta[3])/48 - (I/8)*Pi*G[0, \[Delta]]*Zeta[3] + 
   (G[0, 1, \[Delta]]*Zeta[3])/8 - (5*Zeta[5])/8, 
 G[0, -1, -1, -1, 0, -1 + \[Delta]] -> 
  (I/90)*Pi^5 + (Pi^4*G[1, \[Delta]])/90 + (Pi^2*G[1, 0, 0, \[Delta]])/6 - 
   I*Pi*G[1, 0, 0, 0, \[Delta]] + G[1, 0, 0, 0, 1, \[Delta]] - 
   (Pi^2*Zeta[3])/6 + G[1, 0, \[Delta]]*Zeta[3] + 3*Zeta[5], 
 G[0, -1, -1, 0, 0, -1 + \[Delta]] -> (I/72)*Pi^5 - 
   (Pi^4*G[1, \[Delta]])/360 - (I/6)*Pi^3*G[1, 0, \[Delta]] - 
   (Pi^2*G[1, 0, 0, \[Delta]])/2 - I*Pi*G[1, 0, 0, 1, \[Delta]] + 
   G[1, 0, 0, 1, 1, \[Delta]] - (Pi^2*Zeta[3])/3 - 
   I*Pi*G[1, \[Delta]]*Zeta[3] - G[1, 0, \[Delta]]*Zeta[3] - Zeta[5]/2, 
 G[0, -1, 0, -1, 0, -1 + \[Delta]] -> (-1/120*I)*Pi^5 + 
   (7*Pi^4*G[1, \[Delta]])/360 + (I/6)*Pi^3*G[1, 0, \[Delta]] + 
   (Pi^2*G[1, 0, 1, \[Delta]])/6 - I*Pi*G[1, 0, 1, 0, \[Delta]] + 
   G[1, 0, 1, 0, 1, \[Delta]] - (Pi^2*Zeta[3])/3 + 
   (2*I)*Pi*G[1, \[Delta]]*Zeta[3] + 2*G[1, 0, \[Delta]]*Zeta[3] + 2*Zeta[5], 
 G[0, -1, 0, 0, 0, -1 + \[Delta]] -> (I/180)*Pi^5 - 
   (13*Pi^4*G[1, \[Delta]])/180 + (I/6)*Pi^3*G[1, 0, \[Delta]] - 
   (Pi^2*G[1, 0, 1, \[Delta]])/2 - I*Pi*G[1, 0, 1, 1, \[Delta]] + 
   G[1, 0, 1, 1, 1, \[Delta]] - Pi^2*Zeta[3] + I*Pi*G[1, \[Delta]]*Zeta[3] + 
   4*Zeta[5], G[0, -1, 0, 1, 0, -1 + \[Delta]] -> 
  (-17*Pi^4*G[1, \[Delta]])/1440 - (I/12)*Pi^3*G[1, 0, \[Delta]] - 
   (Pi^2*G[1, 0, 1, \[Delta]])/12 + I*Pi*(-(G[0, 2]*G[1, 0, 1, \[Delta]]) - 
     G[1, 0, 1, 2, \[Delta]]) + G[1, 0, 1, 2, 1, \[Delta]] - 
   (Pi^2*Zeta[3])/24 - ((5*I)/8)*Pi*G[1, \[Delta]]*Zeta[3] - 
   (3*G[1, 0, \[Delta]]*Zeta[3])/2 + I*((71*Pi^5)/1440 + (Pi^3*Log[2]^2)/6 - 
     (Pi*Log[2]^4)/6 - 4*Pi*PolyLog[4, 1/2] - (7*Pi*Log[2]*Zeta[3])/2) + 
   2*Zeta[5], G[0, 0, -1, -1, 0, -1 + \[Delta]] -> 
  (-1/360*I)*Pi^5 - (Pi^4*G[1, \[Delta]])/72 + (Pi^2*G[1, 1, 0, \[Delta]])/
    6 - I*Pi*G[1, 1, 0, 0, \[Delta]] + G[1, 1, 0, 0, 1, \[Delta]] - 
   I*Pi*G[1, \[Delta]]*Zeta[3] + G[1, 1, \[Delta]]*Zeta[3] - Zeta[5]/2, 
 G[0, 0, -1, 0, 0, -1 + \[Delta]] -> (-1/30*I)*Pi^5 + 
   (Pi^4*G[1, \[Delta]])/20 - (I/6)*Pi^3*G[1, 1, \[Delta]] - 
   (Pi^2*G[1, 1, 0, \[Delta]])/2 - I*Pi*G[1, 1, 0, 1, \[Delta]] + 
   G[1, 1, 0, 1, 1, \[Delta]] + (Pi^2*Zeta[3])/2 - 
   (2*I)*Pi*G[1, \[Delta]]*Zeta[3] - G[1, 1, \[Delta]]*Zeta[3] - 6*Zeta[5], 
 G[0, 0, 0, -1, 0, -1 + \[Delta]] -> (I/90)*Pi^5 + (Pi^4*G[1, \[Delta]])/30 + 
   (I/6)*Pi^3*G[1, 1, \[Delta]] + (Pi^2*G[1, 1, 1, \[Delta]])/6 - 
   I*Pi*G[1, 1, 1, 0, \[Delta]] + G[1, 1, 1, 0, 1, \[Delta]] + 
   I*Pi*G[1, \[Delta]]*Zeta[3] + 2*G[1, 1, \[Delta]]*Zeta[3] + 4*Zeta[5], 
 G[0, 0, 0, 0, 0, -1 + \[Delta]] -> (-1/120*I)*Pi^5 + 
   (Pi^4*G[1, \[Delta]])/24 + (I/6)*Pi^3*G[1, 1, \[Delta]] - 
   (Pi^2*G[1, 1, 1, \[Delta]])/2 - I*Pi*G[1, 1, 1, 1, \[Delta]] + 
   G[1, 1, 1, 1, 1, \[Delta]], G[0, 0, 0, 1, 0, -1 + \[Delta]] -> 
  ((-7*I)/720)*Pi^5 - (7*Pi^4*G[1, \[Delta]])/240 - 
   (I/12)*Pi^3*G[1, 1, \[Delta]] - (Pi^2*G[1, 1, 1, \[Delta]])/12 + 
   I*Pi*(-(G[0, 2]*G[1, 1, 1, \[Delta]]) - G[1, 1, 1, 2, \[Delta]]) + 
   G[1, 1, 1, 2, 1, \[Delta]] - ((3*I)/4)*Pi*G[1, \[Delta]]*Zeta[3] - 
   (3*G[1, 1, \[Delta]]*Zeta[3])/2 - (15*Zeta[5])/4, 
 G[0, 0, 1, 0, 0, -1 + \[Delta]] -> ((7*I)/240)*Pi^5 - 
   (Pi^4*G[1, \[Delta]])/80 + (I/12)*Pi^3*G[1, 1, \[Delta]] + 
   Pi^2*(-1/2*(G[0, 2]*G[1, 1, \[Delta]]) - G[1, 1, 2, \[Delta]]/2) - 
   I*Pi*G[1, 1, 2, 1, \[Delta]] + G[1, 1, 2, 1, 1, \[Delta]] - 
   (3*Pi^2*Zeta[3])/8 + ((3*I)/2)*Pi*G[1, \[Delta]]*Zeta[3] + 
   (3*G[1, 1, \[Delta]]*Zeta[3])/4 + (45*Zeta[5])/8, 
 G[0, 0, 1, 1, 0, -1 + \[Delta]] -> -1/288*(Pi^4*G[1, \[Delta]]) + 
   Pi^2*(-1/12*(G[0, 2]*G[1, 1, \[Delta]]) - G[1, 1, 2, \[Delta]]/12) + 
   I*Pi*(-(G[0, 0, 2]*G[1, 1, \[Delta]]) - G[0, 2]*G[1, 1, 2, \[Delta]] - 
     G[1, 1, 2, 2, \[Delta]]) + G[1, 1, 2, 2, 1, \[Delta]] - 
   (7*Pi^2*Zeta[3])/48 - (I/8)*Pi*G[1, \[Delta]]*Zeta[3] + 
   (G[1, 1, \[Delta]]*Zeta[3])/8 + I*(Pi^5/48 + (Pi^3*Log[2]^2)/12 - 
     (Pi*Log[2]^4)/12 - 2*Pi*PolyLog[4, 1/2] - (7*Pi*Log[2]*Zeta[3])/4) + 
   (23*Zeta[5])/16, G[0, 1, 0, -1, 0, -1 + \[Delta]] -> 
  (I/480)*Pi^5 + (17*Pi^4*G[1, \[Delta]])/480 + 
   I*Pi^3*((G[0, 2]*G[1, \[Delta]])/6 + G[1, 2, \[Delta]]/6) + 
   (Pi^2*G[1, 2, 1, \[Delta]])/6 - I*Pi*G[1, 2, 1, 0, \[Delta]] + 
   G[1, 2, 1, 0, 1, \[Delta]] - (Pi^2*Zeta[3])/8 - 
   ((5*I)/8)*Pi*G[1, \[Delta]]*Zeta[3] + 
   (2*G[0, 2]*G[1, \[Delta]] + 2*G[1, 2, \[Delta]])*Zeta[3] - 
   4*G[1, \[Delta]]*(-1/24*(Pi^2*Log[2]^2) + Log[2]^4/24 + PolyLog[4, 1/2] + 
     (7*Log[2]*Zeta[3])/8) + 2*Zeta[5], G[0, 1, 0, 0, 0, -1 + \[Delta]] -> 
  ((-11*I)/720)*Pi^5 + (23*Pi^4*G[1, \[Delta]])/720 + 
   I*Pi^3*((G[0, 2]*G[1, \[Delta]])/6 + G[1, 2, \[Delta]]/6) - 
   (Pi^2*G[1, 2, 1, \[Delta]])/2 - I*Pi*G[1, 2, 1, 1, \[Delta]] + 
   G[1, 2, 1, 1, 1, \[Delta]] + (3*Pi^2*Zeta[3])/4 - 
   ((3*I)/4)*Pi*G[1, \[Delta]]*Zeta[3] - (15*Zeta[5])/4, 
 G[0, 1, 0, 1, 0, -1 + \[Delta]] -> (-11*Pi^4*G[1, \[Delta]])/288 + 
   I*Pi^3*(-1/12*(G[0, 2]*G[1, \[Delta]]) - G[1, 2, \[Delta]]/12) - 
   (Pi^2*G[1, 2, 1, \[Delta]])/12 + I*Pi*(-(G[0, 2]*G[1, 2, 1, \[Delta]]) - 
     G[1, 2, 1, 2, \[Delta]]) + G[1, 2, 1, 2, 1, \[Delta]] + 
   (11*Pi^2*Zeta[3])/24 + (I/4)*Pi*G[1, \[Delta]]*Zeta[3] + 
   ((-3*G[0, 2]*G[1, \[Delta]])/2 - (3*G[1, 2, \[Delta]])/2)*Zeta[3] + 
   4*G[1, \[Delta]]*(-1/24*(Pi^2*Log[2]^2) + Log[2]^4/24 + PolyLog[4, 1/2] + 
     (7*Log[2]*Zeta[3])/8) + I*((-13*Pi^5)/288 - (Pi^3*Log[2]^2)/6 + 
     (Pi*Log[2]^4)/6 + 4*Pi*PolyLog[4, 1/2] + (7*Pi*Log[2]*Zeta[3])/2) - 
   (23*Zeta[5])/4, G[0, 1, 1, 0, 0, -1 + \[Delta]] -> 
  (I/288)*Pi^5 + (Pi^4*G[1, \[Delta]])/48 + 
   I*Pi^3*((G[0, 2]*G[1, \[Delta]])/12 + G[1, 2, \[Delta]]/12) + 
   Pi^2*(-1/2*(G[1, \[Delta]]*G[0, 0, 2]) - (G[0, 2]*G[1, 2, \[Delta]])/2 - 
     G[1, 2, 2, \[Delta]]/2) - I*Pi*G[1, 2, 2, 1, \[Delta]] + 
   G[1, 2, 2, 1, 1, \[Delta]] - (7*Pi^2*Zeta[3])/48 - 
   (I/8)*Pi*G[1, \[Delta]]*Zeta[3] + ((3*G[0, 2]*G[1, \[Delta]])/4 + 
     (3*G[1, 2, \[Delta]])/4)*Zeta[3] - 2*G[1, \[Delta]]*
    (-1/24*(Pi^2*Log[2]^2) + Log[2]^4/24 + PolyLog[4, 1/2] + 
     (7*Log[2]*Zeta[3])/8) + (23*Zeta[5])/16, 
 G[0, 1, 1, 1, 0, -1 + \[Delta]] -> -1/90*(Pi^4*G[1, \[Delta]]) + 
   Pi^2*(-1/12*(G[1, \[Delta]]*G[0, 0, 2]) - (G[0, 2]*G[1, 2, \[Delta]])/12 - 
     G[1, 2, 2, \[Delta]]/12) + I*Pi*(-(G[0, 0, 2]*G[1, 2, \[Delta]]) - 
     G[1, \[Delta]]*G[0, 0, 0, 2] - G[0, 2]*G[1, 2, 2, \[Delta]] - 
     G[1, 2, 2, 2, \[Delta]]) + G[1, 2, 2, 2, 1, \[Delta]] + 
   (5*Pi^2*Zeta[3])/32 + ((G[0, 2]*G[1, \[Delta]])/8 + G[1, 2, \[Delta]]/8)*
    Zeta[3] + G[1, \[Delta]]*(-1/24*(Pi^2*Log[2]^2) + Log[2]^4/24 + 
     PolyLog[4, 1/2] + (7*Log[2]*Zeta[3])/8) + 
   I*(-1/90*Pi^5 - (Pi^3*Log[2]^2)/24 + (Pi*Log[2]^4)/24 + 
     Pi*PolyLog[4, 1/2] + (7*Pi*Log[2]*Zeta[3])/8) - (59*Zeta[5])/32, 
 G[1, 0, 0, 1, 0, -1 + \[Delta]] -> (I/288)*Pi^5 + 
   Pi^4*((-7*G[0, 2])/240 - (7*G[2, \[Delta]])/240) - 
   (I/12)*Pi^3*G[2, 1, \[Delta]] - (Pi^2*G[2, 1, 1, \[Delta]])/12 + 
   I*Pi*(-(G[0, 2]*G[2, 1, 1, \[Delta]]) - G[2, 1, 1, 2, \[Delta]]) + 
   G[2, 1, 1, 2, 1, \[Delta]] - (Pi^2*Zeta[3])/16 + 
   I*Pi*((-3*G[0, 2])/4 - (3*G[2, \[Delta]])/4)*Zeta[3] - 
   (3*G[2, 1, \[Delta]]*Zeta[3])/2 + (51*Zeta[5])/32, 
 G[1, 0, 1, 0, 0, -1 + \[Delta]] -> 
  Pi^4*(-1/80*G[0, 2] - G[2, \[Delta]]/80) + (I/12)*Pi^3*G[2, 1, \[Delta]] + 
   Pi^2*(-1/2*(G[0, 2]*G[2, 1, \[Delta]]) - G[2, 1, 2, \[Delta]]/2) - 
   I*Pi*G[2, 1, 2, 1, \[Delta]] + G[2, 1, 2, 1, 1, \[Delta]] - 
   (Pi^2*Zeta[3])/24 + I*Pi*((3*G[0, 2])/2 + (3*G[2, \[Delta]])/2)*Zeta[3] + 
   (3*G[2, 1, \[Delta]]*Zeta[3])/4 + I*((11*Pi^5)/288 + (Pi^3*Log[2]^2)/6 - 
     (Pi*Log[2]^4)/6 - 4*Pi*PolyLog[4, 1/2] - (7*Pi*Log[2]*Zeta[3])/2) + 
   (41*Zeta[5])/32, G[1, 0, 1, 1, 0, -1 + \[Delta]] -> 
  Pi^4*(-1/288*G[0, 2] - G[2, \[Delta]]/288) + 
   Pi^2*(-1/12*(G[0, 2]*G[2, 1, \[Delta]]) - G[2, 1, 2, \[Delta]]/12) + 
   I*Pi*(-(G[0, 0, 2]*G[2, 1, \[Delta]]) - G[0, 2]*G[2, 1, 2, \[Delta]] - 
     G[2, 1, 2, 2, \[Delta]]) + G[2, 1, 2, 2, 1, \[Delta]] - 
   (11*Pi^2*Zeta[3])/48 + I*Pi*(-1/8*G[0, 2] - G[2, \[Delta]]/8)*Zeta[3] + 
   (G[2, 1, \[Delta]]*Zeta[3])/8 + I*(Pi^5/30 + (Pi^3*Log[2]^2)/8 - 
     (Pi*Log[2]^4)/8 - 3*Pi*PolyLog[4, 1/2] - (21*Pi*Log[2]*Zeta[3])/8) + 
   (177*Zeta[5])/64, G[-1, -1, -1, -1, -1, 0, -1 + \[Delta]] -> 
  Pi^6/945 + (Pi^4*G[0, 0, \[Delta]])/90 + (Pi^2*G[0, 0, 0, 0, \[Delta]])/6 - 
   I*Pi*G[0, 0, 0, 0, 0, \[Delta]] + G[0, 0, 0, 0, 0, 1, \[Delta]] + 
   G[0, 0, 0, \[Delta]]*Zeta[3] + G[0, \[Delta]]*Zeta[5], 
 G[-1, -1, -1, -1, 0, 0, -1 + \[Delta]] -> 
  -1/2*(Pi^2*G[0, 0, 0, 0, \[Delta]]) - I*Pi*G[0, 0, 0, 0, 1, \[Delta]] + 
   G[0, 0, 0, 0, 1, 1, \[Delta]] + G[0, 0, 0, \[Delta]]*
    ((-1/6*I)*Pi^3 - Zeta[3]) + G[0, 0, \[Delta]]*
    (-1/360*Pi^4 - I*Pi*Zeta[3]) + (-Pi^6 + 630*Zeta[3]^2)/1260 + 
   G[0, \[Delta]]*((-1/90*I)*Pi^5 + (Pi^2*Zeta[3])/6 - 2*Zeta[5]) - 
   I*Pi*Zeta[5], G[-1, -1, -1, 0, -1, 0, -1 + \[Delta]] -> 
  (Pi^2*G[0, 0, 0, 1, \[Delta]])/6 - I*Pi*G[0, 0, 0, 1, 0, \[Delta]] + 
   G[0, 0, 0, 1, 0, 1, \[Delta]] + G[0, 0, 0, \[Delta]]*
    ((I/6)*Pi^3 + 2*Zeta[3]) + G[0, 0, \[Delta]]*
    ((7*Pi^4)/360 + (2*I)*Pi*Zeta[3]) + (37*Pi^6 - 11340*Zeta[3]^2)/11340 + 
   (4*I)*Pi*Zeta[5] + G[0, \[Delta]]*((I/30)*Pi^5 - (Pi^2*Zeta[3])/3 + 
     (11*Zeta[5])/2), G[-1, -1, -1, 0, 0, 0, -1 + \[Delta]] -> 
  (I/6)*Pi^3*G[0, 0, 0, \[Delta]] - (Pi^2*G[0, 0, 0, 1, \[Delta]])/2 - 
   I*Pi*G[0, 0, 0, 1, 1, \[Delta]] + G[0, 0, 0, 1, 1, 1, \[Delta]] + 
   G[0, 0, \[Delta]]*((-13*Pi^4)/180 + I*Pi*Zeta[3]) + 
   (-61*Pi^6 - 15120*Zeta[3]^2)/15120 + G[0, \[Delta]]*
    ((I/360)*Pi^5 - (2*Pi^2*Zeta[3])/3 + 2*Zeta[5]) + 
   (I/6)*(-(Pi^3*Zeta[3]) + 12*Pi*Zeta[5]), 
 G[-1, -1, -1, 0, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 0, 0, 1, \[Delta]] - 
   I*Pi*G[0, 0, 0, 1, 2, \[Delta]] + G[0, 0, 0, 1, 2, 1, \[Delta]] + 
   G[0, 0, 0, \[Delta]]*((-1/12*I)*Pi^3 - (3*Zeta[3])/2) + 
   G[0, 0, \[Delta]]*((-17*Pi^4)/1440 - ((5*I)/8)*Pi*Zeta[3]) + 
   (-1102*Pi^6 - 3780*Pi^4*Log[2]^2 + 3780*Pi^2*Log[2]^4 + 
     90720*Pi^2*PolyLog[4, 1/2] + 79380*Pi^2*Log[2]*Zeta[3] - 
     144585*Zeta[3]^2)/362880 + G[0, \[Delta]]*
    ((I/288)*Pi^5 + (I/24)*Pi^3*Log[2]^2 - (I/24)*Pi*Log[2]^4 - 
     I*Pi*PolyLog[4, 1/2] + (Pi^2*Zeta[3])/32 - ((7*I)/8)*Pi*Log[2]*Zeta[3] - 
     (41*Zeta[5])/32) + (I/1440)*(-40*Pi^3*Log[2]^3 + 48*Pi*Log[2]^5 + 
     1440*Pi*Log[2]*PolyLog[4, 1/2] + 1440*Pi*PolyLog[5, 1/2] - 
     105*Pi^3*Zeta[3] + 630*Pi*Log[2]^2*Zeta[3] - 1215*Pi*Zeta[5]), 
 G[-1, -1, 0, -1, -1, 0, -1 + \[Delta]] -> (Pi^2*G[0, 0, 1, 0, \[Delta]])/6 - 
   I*Pi*G[0, 0, 1, 0, 0, \[Delta]] + G[0, 0, 1, 0, 0, 1, \[Delta]] + 
   G[0, 0, 1, \[Delta]]*Zeta[3] + G[0, 0, \[Delta]]*
    (-1/72*Pi^4 - I*Pi*Zeta[3]) + (-19*Pi^6 + 1890*Zeta[3]^2)/3780 + 
   G[0, \[Delta]]*((-1/30*I)*Pi^5 + (Pi^2*Zeta[3])/6 - (9*Zeta[5])/2) - 
   (6*I)*Pi*Zeta[5], G[-1, -1, 0, -1, 0, 0, -1 + \[Delta]] -> 
  -1/2*(Pi^2*G[0, 0, 1, 0, \[Delta]]) - I*Pi*G[0, 0, 1, 0, 1, \[Delta]] + 
   G[0, 0, 1, 0, 1, 1, \[Delta]] + G[0, 0, 1, \[Delta]]*
    ((-1/6*I)*Pi^3 - Zeta[3]) + G[0, 0, \[Delta]]*
    (Pi^4/20 - (2*I)*Pi*Zeta[3]) + (79*Pi^6 + 12960*Zeta[3]^2)/6480 + 
   G[0, \[Delta]]*(((-7*I)/360)*Pi^5 + (4*Pi^2*Zeta[3])/3 - (11*Zeta[5])/2) + 
   (I/6)*(2*Pi^3*Zeta[3] - 33*Pi*Zeta[5]), 
 G[-1, -1, 0, 0, -1, 0, -1 + \[Delta]] -> (Pi^2*G[0, 0, 1, 1, \[Delta]])/6 - 
   I*Pi*G[0, 0, 1, 1, 0, \[Delta]] + G[0, 0, 1, 1, 0, 1, \[Delta]] + 
   G[0, 0, 1, \[Delta]]*((I/6)*Pi^3 + 2*Zeta[3]) + 
   G[0, 0, \[Delta]]*(Pi^4/30 + I*Pi*Zeta[3]) + (17*Pi^6 + 4536*Zeta[3]^2)/
    9072 + G[0, \[Delta]]*((I/72)*Pi^5 - (Pi^2*Zeta[3])/6 + (9*Zeta[5])/2) + 
   (I/6)*(Pi^3*Zeta[3] - 3*Pi*Zeta[5]), 
 G[-1, -1, 0, 0, 0, 0, -1 + \[Delta]] -> (Pi^4*G[0, 0, \[Delta]])/24 + 
   (I/6)*Pi^3*G[0, 0, 1, \[Delta]] - (Pi^2*G[0, 0, 1, 1, \[Delta]])/2 - 
   I*Pi*G[0, 0, 1, 1, 1, \[Delta]] + G[0, 0, 1, 1, 1, 1, \[Delta]] + 
   (Pi^6 + 840*Zeta[3]^2)/1680 + G[0, \[Delta]]*
    ((I/60)*Pi^5 + (Pi^2*Zeta[3])/2 - Zeta[5]) + 
   (I/3)*(Pi^3*Zeta[3] - 6*Pi*Zeta[5]), 
 G[-1, -1, 0, 0, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 0, 1, 1, \[Delta]] - 
   I*Pi*G[0, 0, 1, 1, 2, \[Delta]] + G[0, 0, 1, 1, 2, 1, \[Delta]] + 
   G[0, 0, 1, \[Delta]]*((-1/12*I)*Pi^3 - (3*Zeta[3])/2) + 
   G[0, 0, \[Delta]]*((-7*Pi^4)/240 - ((3*I)/4)*Pi*Zeta[3]) + 
   (331*Pi^6 + 3780*Pi^4*Log[2]^2 - 3780*Pi^2*Log[2]^4 - 
     90720*Pi^2*PolyLog[4, 1/2] - 79380*Pi^2*Log[2]*Zeta[3] + 
     189945*Zeta[3]^2)/181440 + G[0, \[Delta]]*(((-11*I)/360)*Pi^5 - 
     (I/12)*Pi^3*Log[2]^2 + (I/12)*Pi*Log[2]^4 + (2*I)*Pi*PolyLog[4, 1/2] + 
     (13*Pi^2*Zeta[3])/48 + ((7*I)/4)*Pi*Log[2]*Zeta[3] - (83*Zeta[5])/16) + 
   (I/1440)*(80*Pi^3*Log[2]^3 - 96*Pi*Log[2]^5 - 
     2880*Pi*Log[2]*PolyLog[4, 1/2] - 2880*Pi*PolyLog[5, 1/2] + 
     90*Pi^3*Zeta[3] - 1260*Pi*Log[2]^2*Zeta[3] + 675*Pi*Zeta[5]), 
 G[-1, -1, 0, 1, 0, 0, -1 + \[Delta]] -> 
  -1/2*(Pi^2*G[0, 0, 1, 2, \[Delta]]) - I*Pi*G[0, 0, 1, 2, 1, \[Delta]] + 
   G[0, 0, 1, 2, 1, 1, \[Delta]] + G[0, 0, 1, \[Delta]]*
    ((I/12)*Pi^3 - (Pi^2*G[0, 2])/2 + (3*Zeta[3])/4) + 
   G[0, 0, \[Delta]]*(-1/80*Pi^4 + ((3*I)/2)*Pi*Zeta[3]) + 
   (718*Pi^6 + 3780*Pi^4*Log[2]^2 - 3780*Pi^2*Log[2]^4 - 
     90720*Pi^2*PolyLog[4, 1/2] - 79380*Pi^2*Log[2]*Zeta[3] - 
     133245*Zeta[3]^2)/181440 + G[0, \[Delta]]*(((17*I)/1440)*Pi^5 - 
     (9*Pi^2*Zeta[3])/16 + (67*Zeta[5])/16) + 
   (I/32)*(-(Pi^3*Zeta[3]) + 41*Pi*Zeta[5]), 
 G[-1, -1, 0, 1, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 0, 1, 2, \[Delta]] - 
   I*Pi*G[0, 0, 1, 2, 2, \[Delta]] + G[0, 0, 1, 2, 2, 1, \[Delta]] + 
   G[0, 0, 1, \[Delta]]*(-1/12*(Pi^2*G[0, 2]) - I*Pi*G[0, 0, 2] + 
     Zeta[3]/8) + G[0, 0, \[Delta]]*(-1/288*Pi^4 - (I/8)*Pi*Zeta[3]) + 
   (12*Pi^6 + 140*Pi^4*Log[2]^2 - 140*Pi^2*Log[2]^4 - 
     3360*Pi^2*PolyLog[4, 1/2] - 2940*Pi^2*Log[2]*Zeta[3] + 4725*Zeta[3]^2)/
    40320 + G[0, \[Delta]]*((-1/480*I)*Pi^5 + (Pi^2*Zeta[3])/48 - 
     (5*Zeta[5])/8) + (I/1440)*(40*Pi^3*Log[2]^3 - 48*Pi*Log[2]^5 - 
     1440*Pi*Log[2]*PolyLog[4, 1/2] - 1440*Pi*PolyLog[5, 1/2] - 
     630*Pi*Log[2]^2*Zeta[3] + 1125*Pi*Zeta[5]), 
 G[-1, 0, -1, -1, -1, 0, -1 + \[Delta]] -> (Pi^4*G[0, 1, \[Delta]])/90 + 
   (Pi^2*G[0, 1, 0, 0, \[Delta]])/6 - I*Pi*G[0, 1, 0, 0, 0, \[Delta]] + 
   G[0, 1, 0, 0, 0, 1, \[Delta]] + G[0, 1, 0, \[Delta]]*Zeta[3] + 
   (59*Pi^6 - 11340*Zeta[3]^2)/11340 + (4*I)*Pi*Zeta[5] + 
   G[0, \[Delta]]*((I/90)*Pi^5 - (Pi^2*Zeta[3])/6 + 3*Zeta[5]), 
 G[-1, 0, -1, -1, 0, 0, -1 + \[Delta]] -> 
  -1/2*(Pi^2*G[0, 1, 0, 0, \[Delta]]) - I*Pi*G[0, 1, 0, 0, 1, \[Delta]] + 
   G[0, 1, 0, 0, 1, 1, \[Delta]] + G[0, 1, 0, \[Delta]]*
    ((-1/6*I)*Pi^3 - Zeta[3]) + G[0, 1, \[Delta]]*
    (-1/360*Pi^4 - I*Pi*Zeta[3]) + (-671*Pi^6 + 22680*Zeta[3]^2)/45360 + 
   G[0, \[Delta]]*((I/72)*Pi^5 - (Pi^2*Zeta[3])/3 - Zeta[5]/2) + 
   (I/6)*(-(Pi^3*Zeta[3]) + 27*Pi*Zeta[5]), 
 G[-1, 0, -1, 0, -1, 0, -1 + \[Delta]] -> (Pi^2*G[0, 1, 0, 1, \[Delta]])/6 - 
   I*Pi*G[0, 1, 0, 1, 0, \[Delta]] + G[0, 1, 0, 1, 0, 1, \[Delta]] + 
   G[0, 1, 0, \[Delta]]*((I/6)*Pi^3 + 2*Zeta[3]) + 
   G[0, 1, \[Delta]]*((7*Pi^4)/360 + (2*I)*Pi*Zeta[3]) + 
   (31*Pi^6 - 60480*Zeta[3]^2)/15120 + G[0, \[Delta]]*
    ((-1/120*I)*Pi^5 - (Pi^2*Zeta[3])/3 + 2*Zeta[5]) + 
   (I/3)*(-(Pi^3*Zeta[3]) + 6*Pi*Zeta[5]), 
 G[-1, 0, -1, 0, 0, 0, -1 + \[Delta]] -> (I/6)*Pi^3*G[0, 1, 0, \[Delta]] - 
   (Pi^2*G[0, 1, 0, 1, \[Delta]])/2 - I*Pi*G[0, 1, 0, 1, 1, \[Delta]] + 
   G[0, 1, 0, 1, 1, 1, \[Delta]] + G[0, 1, \[Delta]]*
    ((-13*Pi^4)/180 + I*Pi*Zeta[3]) + (-293*Pi^6 - 45360*Zeta[3]^2)/45360 + 
   G[0, \[Delta]]*((I/180)*Pi^5 - Pi^2*Zeta[3] + 4*Zeta[5]) + 
   (I/6)*(-4*Pi^3*Zeta[3] + 33*Pi*Zeta[5]), 
 G[-1, 0, -1, 0, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 1, 0, 1, \[Delta]] - 
   I*Pi*G[0, 1, 0, 1, 2, \[Delta]] + G[0, 1, 0, 1, 2, 1, \[Delta]] + 
   G[0, 1, 0, \[Delta]]*((-1/12*I)*Pi^3 - (3*Zeta[3])/2) + 
   G[0, 1, \[Delta]]*((-17*Pi^4)/1440 - ((5*I)/8)*Pi*Zeta[3]) + 
   (-191*Pi^6 - 1008*Pi^4*Log[2]^2 + 1008*Pi^2*Log[2]^4 + 
     24192*Pi^2*PolyLog[4, 1/2] + 21168*Pi^2*Log[2]*Zeta[3] + 3780*Zeta[3]^2)/
    24192 + G[0, \[Delta]]*(((71*I)/1440)*Pi^5 + (I/6)*Pi^3*Log[2]^2 - 
     (I/6)*Pi*Log[2]^4 - (4*I)*Pi*PolyLog[4, 1/2] - (Pi^2*Zeta[3])/24 - 
     ((7*I)/2)*Pi*Log[2]*Zeta[3] + 2*Zeta[5]) + 
   (I/2880)*(-320*Pi^3*Log[2]^3 + 384*Pi*Log[2]^5 + 
     11520*Pi*Log[2]*PolyLog[4, 1/2] + 11520*Pi*PolyLog[5, 1/2] - 
     360*Pi^3*Zeta[3] + 5040*Pi*Log[2]^2*Zeta[3] - 4005*Pi*Zeta[5]), 
 G[-1, 0, 0, -1, -1, 0, -1 + \[Delta]] -> (-13*Pi^6)/15120 + 
   (Pi^2*G[0, 1, 1, 0, \[Delta]])/6 - I*Pi*G[0, 1, 1, 0, 0, \[Delta]] + 
   G[0, 1, 1, 0, 0, 1, \[Delta]] + G[0, 1, 1, \[Delta]]*Zeta[3] + 
   G[0, 1, \[Delta]]*(-1/72*Pi^4 - I*Pi*Zeta[3]) + 
   G[0, \[Delta]]*((-1/360*I)*Pi^5 - Zeta[5]/2) - (I/2)*Pi*Zeta[5], 
 G[-1, 0, 0, -1, 0, 0, -1 + \[Delta]] -> 
  -1/2*(Pi^2*G[0, 1, 1, 0, \[Delta]]) - I*Pi*G[0, 1, 1, 0, 1, \[Delta]] + 
   G[0, 1, 1, 0, 1, 1, \[Delta]] + G[0, 1, 1, \[Delta]]*
    ((-1/6*I)*Pi^3 - Zeta[3]) + G[0, 1, \[Delta]]*
    (Pi^4/20 - (2*I)*Pi*Zeta[3]) + (29*Pi^6 + 7560*Zeta[3]^2)/15120 + 
   G[0, \[Delta]]*((-1/30*I)*Pi^5 + (Pi^2*Zeta[3])/2 - 6*Zeta[5]) + 
   (I/6)*(Pi^3*Zeta[3] - 27*Pi*Zeta[5]), 
 G[-1, 0, 0, 0, -1, 0, -1 + \[Delta]] -> (Pi^2*G[0, 1, 1, 1, \[Delta]])/6 - 
   I*Pi*G[0, 1, 1, 1, 0, \[Delta]] + G[0, 1, 1, 1, 0, 1, \[Delta]] + 
   G[0, 1, 1, \[Delta]]*((I/6)*Pi^3 + 2*Zeta[3]) + 
   G[0, 1, \[Delta]]*(Pi^4/30 + I*Pi*Zeta[3]) + (59*Pi^6 - 11340*Zeta[3]^2)/
    11340 + G[0, \[Delta]]*((I/90)*Pi^5 + 4*Zeta[5]) + 
   (I/6)*(-(Pi^3*Zeta[3]) + 18*Pi*Zeta[5]), 
 G[-1, 0, 0, 0, 0, 0, -1 + \[Delta]] -> (37*Pi^6)/15120 - 
   (I/120)*Pi^5*G[0, \[Delta]] + (Pi^4*G[0, 1, \[Delta]])/24 + 
   (I/6)*Pi^3*G[0, 1, 1, \[Delta]] - (Pi^2*G[0, 1, 1, 1, \[Delta]])/2 - 
   I*Pi*G[0, 1, 1, 1, 1, \[Delta]] + G[0, 1, 1, 1, 1, 1, \[Delta]] + 
   (I/6)*(-(Pi^3*Zeta[3]) + 6*Pi*Zeta[5]), 
 G[-1, 0, 0, 0, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 1, 1, 1, \[Delta]] - 
   I*Pi*G[0, 1, 1, 1, 2, \[Delta]] + G[0, 1, 1, 1, 2, 1, \[Delta]] + 
   G[0, 1, 1, \[Delta]]*((-1/12*I)*Pi^3 - (3*Zeta[3])/2) + 
   G[0, 1, \[Delta]]*((-7*Pi^4)/240 - ((3*I)/4)*Pi*Zeta[3]) + 
   (-811*Pi^6 + 136080*Zeta[3]^2)/181440 + 
   G[0, \[Delta]]*(((-7*I)/720)*Pi^5 - (15*Zeta[5])/4) + 
   (I/96)*(8*Pi^3*Zeta[3] - 177*Pi*Zeta[5]), 
 G[-1, 0, 0, 1, 0, 0, -1 + \[Delta]] -> -1/2*(Pi^2*G[0, 1, 1, 2, \[Delta]]) - 
   I*Pi*G[0, 1, 1, 2, 1, \[Delta]] + G[0, 1, 1, 2, 1, 1, \[Delta]] + 
   G[0, 1, 1, \[Delta]]*((I/12)*Pi^3 - (Pi^2*G[0, 2])/2 + (3*Zeta[3])/4) + 
   G[0, 1, \[Delta]]*(-1/80*Pi^4 + ((3*I)/2)*Pi*Zeta[3]) + 
   (-299*Pi^6 - 1260*Pi^4*Log[2]^2 + 1260*Pi^2*Log[2]^4 + 
     30240*Pi^2*PolyLog[4, 1/2] + 26460*Pi^2*Log[2]*Zeta[3] - 
     31185*Zeta[3]^2)/30240 + G[0, \[Delta]]*(((7*I)/240)*Pi^5 - 
     (3*Pi^2*Zeta[3])/8 + (45*Zeta[5])/8) + 
   (I/48)*(-13*Pi^3*Zeta[3] + 249*Pi*Zeta[5]), 
 G[-1, 0, 0, 1, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 1, 1, 2, \[Delta]] - 
   I*Pi*G[0, 1, 1, 2, 2, \[Delta]] + G[0, 1, 1, 2, 2, 1, \[Delta]] + 
   G[0, 1, 1, \[Delta]]*(-1/12*(Pi^2*G[0, 2]) - I*Pi*G[0, 0, 2] + 
     Zeta[3]/8) + G[0, 1, \[Delta]]*(-1/288*Pi^4 - (I/8)*Pi*Zeta[3]) + 
   (-53*Pi^6 - 280*Pi^4*Log[2]^2 + 280*Pi^2*Log[2]^4 + 
     6720*Pi^2*PolyLog[4, 1/2] + 5880*Pi^2*Log[2]*Zeta[3] - 15435*Zeta[3]^2)/
    10080 + G[0, \[Delta]]*((I/48)*Pi^5 + (I/12)*Pi^3*Log[2]^2 - 
     (I/12)*Pi*Log[2]^4 - (2*I)*Pi*PolyLog[4, 1/2] - (7*Pi^2*Zeta[3])/48 - 
     ((7*I)/4)*Pi*Log[2]*Zeta[3] + (23*Zeta[5])/16) + 
   (I/720)*(-80*Pi^3*Log[2]^3 + 96*Pi*Log[2]^5 + 
     2880*Pi*Log[2]*PolyLog[4, 1/2] + 2880*Pi*PolyLog[5, 1/2] - 
     105*Pi^3*Zeta[3] + 1260*Pi*Log[2]^2*Zeta[3] - 1755*Pi*Zeta[5]), 
 G[-1, 0, 1, 0, -1, 0, -1 + \[Delta]] -> (Pi^2*G[0, 1, 2, 1, \[Delta]])/6 - 
   I*Pi*G[0, 1, 2, 1, 0, \[Delta]] + G[0, 1, 2, 1, 0, 1, \[Delta]] + 
   G[0, 1, 2, \[Delta]]*((I/6)*Pi^3 + 2*Zeta[3]) + 
   G[0, 1, \[Delta]]*((17*Pi^4)/480 + (I/6)*Pi^3*G[0, 2] + 
     (Pi^2*Log[2]^2)/6 - Log[2]^4/6 - 4*PolyLog[4, 1/2] - 
     ((5*I)/8)*Pi*Zeta[3] + 2*G[0, 2]*Zeta[3] - (7*Log[2]*Zeta[3])/2) + 
   (17*Pi^6 + 4725*Zeta[3]^2)/30240 + G[0, \[Delta]]*
    ((I/480)*Pi^5 - (Pi^2*Zeta[3])/8 + 2*Zeta[5]) + 
   (I/96)*(10*Pi^3*Zeta[3] - 87*Pi*Zeta[5]), 
 G[-1, 0, 1, 0, 0, 0, -1 + \[Delta]] -> (I/6)*Pi^3*G[0, 1, 2, \[Delta]] - 
   (Pi^2*G[0, 1, 2, 1, \[Delta]])/2 - I*Pi*G[0, 1, 2, 1, 1, \[Delta]] + 
   G[0, 1, 2, 1, 1, 1, \[Delta]] + G[0, 1, \[Delta]]*
    ((23*Pi^4)/720 + (I/6)*Pi^3*G[0, 2] - ((3*I)/4)*Pi*Zeta[3]) + 
   (583*Pi^6 + 136080*Zeta[3]^2)/181440 + G[0, \[Delta]]*
    (((-11*I)/720)*Pi^5 + (3*Pi^2*Zeta[3])/4 - (15*Zeta[5])/4) + 
   (I/48)*(17*Pi^3*Zeta[3] - 201*Pi*Zeta[5]), 
 G[-1, 0, 1, 0, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 1, 2, 1, \[Delta]] - 
   I*Pi*G[0, 1, 2, 1, 2, \[Delta]] + G[0, 1, 2, 1, 2, 1, \[Delta]] + 
   G[0, 1, 2, \[Delta]]*((-1/12*I)*Pi^3 - (3*Zeta[3])/2) + 
   G[0, 1, \[Delta]]*((-11*Pi^4)/288 - (I/12)*Pi^3*G[0, 2] - 
     (Pi^2*Log[2]^2)/6 + Log[2]^4/6 + 4*PolyLog[4, 1/2] + (I/4)*Pi*Zeta[3] - 
     (3*G[0, 2]*Zeta[3])/2 + (7*Log[2]*Zeta[3])/2) + 
   (49*Pi^6 + 360*Pi^4*Log[2]^2 - 360*Pi^2*Log[2]^4 - 
     8640*Pi^2*PolyLog[4, 1/2] - 7560*Pi^2*Log[2]*Zeta[3] + 25920*Zeta[3]^2)/
    8640 + G[0, \[Delta]]*(((-13*I)/288)*Pi^5 - (I/6)*Pi^3*Log[2]^2 + 
     (I/6)*Pi*Log[2]^4 + (4*I)*Pi*PolyLog[4, 1/2] + (11*Pi^2*Zeta[3])/24 + 
     ((7*I)/2)*Pi*Log[2]*Zeta[3] - (23*Zeta[5])/4) + 
   (I/1440)*(160*Pi^3*Log[2]^3 - 192*Pi*Log[2]^5 - 
     5760*Pi*Log[2]*PolyLog[4, 1/2] - 5760*Pi*PolyLog[5, 1/2] + 
     345*Pi^3*Zeta[3] - 2520*Pi*Log[2]^2*Zeta[3] + 1485*Pi*Zeta[5]), 
 G[-1, 0, 1, 1, 0, 0, -1 + \[Delta]] -> -1/2*(Pi^2*G[0, 1, 2, 2, \[Delta]]) - 
   I*Pi*G[0, 1, 2, 2, 1, \[Delta]] + G[0, 1, 2, 2, 1, 1, \[Delta]] + 
   G[0, 1, 2, \[Delta]]*((I/12)*Pi^3 - (Pi^2*G[0, 2])/2 + (3*Zeta[3])/4) + 
   G[0, 1, \[Delta]]*(Pi^4/48 + (I/12)*Pi^3*G[0, 2] - (Pi^2*G[0, 0, 2])/2 + 
     (Pi^2*Log[2]^2)/12 - Log[2]^4/12 - 2*PolyLog[4, 1/2] - 
     (I/8)*Pi*Zeta[3] + (3*G[0, 2]*Zeta[3])/4 - (7*Log[2]*Zeta[3])/4) + 
   (2*Pi^6 - 8505*Zeta[3]^2)/22680 + G[0, \[Delta]]*
    ((I/288)*Pi^5 - (7*Pi^2*Zeta[3])/48 + (23*Zeta[5])/16) + 
   (I/48)*(-(Pi^3*Zeta[3]) + 30*Pi*Zeta[5]), 
 G[-1, 0, 1, 1, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 1, 2, 2, \[Delta]] - 
   I*Pi*G[0, 1, 2, 2, 2, \[Delta]] + G[0, 1, 2, 2, 2, 1, \[Delta]] + 
   G[0, 1, 2, \[Delta]]*(-1/12*(Pi^2*G[0, 2]) - I*Pi*G[0, 0, 2] + 
     Zeta[3]/8) + G[0, 1, \[Delta]]*(-1/90*Pi^4 - (Pi^2*G[0, 0, 2])/12 - 
     I*Pi*G[0, 0, 0, 2] - (Pi^2*Log[2]^2)/24 + Log[2]^4/24 + 
     PolyLog[4, 1/2] + (G[0, 2]*Zeta[3])/8 + (7*Log[2]*Zeta[3])/8) + 
   (394*Pi^6 + 3780*Pi^4*Log[2]^2 - 3780*Pi^2*Log[2]^4 - 
     90720*Pi^2*PolyLog[4, 1/2] - 79380*Pi^2*Log[2]*Zeta[3] + 
     371385*Zeta[3]^2)/362880 + G[0, \[Delta]]*((-1/90*I)*Pi^5 - 
     (I/24)*Pi^3*Log[2]^2 + (I/24)*Pi*Log[2]^4 + I*Pi*PolyLog[4, 1/2] + 
     (5*Pi^2*Zeta[3])/32 + ((7*I)/8)*Pi*Log[2]*Zeta[3] - (59*Zeta[5])/32) + 
   (I/1440)*(40*Pi^3*Log[2]^3 - 48*Pi*Log[2]^5 - 
     1440*Pi*Log[2]*PolyLog[4, 1/2] - 1440*Pi*PolyLog[5, 1/2] + 
     105*Pi^3*Zeta[3] - 630*Pi*Log[2]^2*Zeta[3] + 180*Pi*Zeta[5]), 
 G[-1, 1, 0, 0, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 2, 1, 1, \[Delta]] - 
   I*Pi*G[0, 2, 1, 1, 2, \[Delta]] + G[0, 2, 1, 1, 2, 1, \[Delta]] + 
   G[0, 2, 1, \[Delta]]*((-1/12*I)*Pi^3 - (3*Zeta[3])/2) + 
   G[0, 2, \[Delta]]*((-7*Pi^4)/240 - ((3*I)/4)*Pi*Zeta[3]) + 
   (-1270080*HPLs6 - 871*Pi^6 + 7812*Pi^4*Log[2]^2 - 2520*Pi^2*Log[2]^4 - 
     60480*Pi^2*PolyLog[4, 1/2] + 68040*Pi^2*Log[2]*Zeta[3] + 
     697410*Zeta[3]^2)/362880 + G[0, \[Delta]]*
    ((I/288)*Pi^5 - (7*Pi^4*G[0, 2])/240 - (Pi^2*Zeta[3])/16 - 
     ((3*I)/4)*Pi*G[0, 2]*Zeta[3] + (51*Zeta[5])/32) + 
   (I/1440)*(-49*Pi^5*Log[2] + 40*Pi^3*Log[2]^3 - 72*Pi*Log[2]^5 - 
     2880*Pi*Log[2]*PolyLog[4, 1/2] - 5760*Pi*PolyLog[5, 1/2] + 
     105*Pi^3*Zeta[3] + 540*Pi*Log[2]^2*Zeta[3] + 4860*Pi*Zeta[5]), 
 G[-1, 1, 0, 1, 0, 0, -1 + \[Delta]] -> -1/2*(Pi^2*G[0, 2, 1, 2, \[Delta]]) - 
   I*Pi*G[0, 2, 1, 2, 1, \[Delta]] + G[0, 2, 1, 2, 1, 1, \[Delta]] + 
   G[0, 2, 1, \[Delta]]*((I/12)*Pi^3 - (Pi^2*G[0, 2])/2 + (3*Zeta[3])/4) + 
   G[0, 2, \[Delta]]*(-1/80*Pi^4 + ((3*I)/2)*Pi*Zeta[3]) + 
   (272160*HPLs6 - 4013*Pi^6 - 19026*Pi^4*Log[2]^2 + 20160*Pi^2*Log[2]^4 + 
     483840*Pi^2*PolyLog[4, 1/2] + 328860*Pi^2*Log[2]*Zeta[3] - 
     552825*Zeta[3]^2)/181440 + G[0, \[Delta]]*(((11*I)/288)*Pi^5 - 
     (Pi^4*G[0, 2])/80 + (I/6)*Pi^3*Log[2]^2 - (I/6)*Pi*Log[2]^4 - 
     (4*I)*Pi*PolyLog[4, 1/2] - (Pi^2*Zeta[3])/24 + 
     ((3*I)/2)*Pi*G[0, 2]*Zeta[3] - ((7*I)/2)*Pi*Log[2]*Zeta[3] + 
     (41*Zeta[5])/32) + (I/1440)*(-38*Pi^5*Log[2] - 80*Pi^3*Log[2]^3 + 
     48*Pi*Log[2]^5 - 5760*Pi*PolyLog[5, 1/2] - 495*Pi^3*Zeta[3] + 
     1440*Pi*Log[2]^2*Zeta[3] + 11610*Pi*Zeta[5]), 
 G[-1, 1, 0, 1, 1, 0, -1 + \[Delta]] -> 
  (-1/12*Pi^2 - I*Pi*G[0, 2])*G[0, 2, 1, 2, \[Delta]] - 
   I*Pi*G[0, 2, 1, 2, 2, \[Delta]] + G[0, 2, 1, 2, 2, 1, \[Delta]] + 
   G[0, 2, 1, \[Delta]]*(-1/12*(Pi^2*G[0, 2]) - I*Pi*G[0, 0, 2] + 
     Zeta[3]/8) + G[0, 2, \[Delta]]*(-1/288*Pi^4 - (I/8)*Pi*Zeta[3]) + 
   (-211680*HPLs6 - 737*Pi^6 - 3570*Pi^4*Log[2]^2 + 3780*Pi^2*Log[2]^4 + 
     90720*Pi^2*PolyLog[4, 1/2] + 109620*Pi^2*Log[2]*Zeta[3] - 
     210735*Zeta[3]^2)/120960 + G[0, \[Delta]]*
    ((I/30)*Pi^5 - (Pi^4*G[0, 2])/288 + (I/8)*Pi^3*Log[2]^2 - 
     (I/8)*Pi*Log[2]^4 - (3*I)*Pi*PolyLog[4, 1/2] - (11*Pi^2*Zeta[3])/48 - 
     (I/8)*Pi*G[0, 2]*Zeta[3] - ((21*I)/8)*Pi*Log[2]*Zeta[3] + 
     (177*Zeta[5])/64) + (I/480)*(-17*Pi^5*Log[2] - 60*Pi^3*Log[2]^3 + 
     60*Pi*Log[2]^5 + 1440*Pi*Log[2]*PolyLog[4, 1/2] - 110*Pi^3*Zeta[3] + 
     1290*Pi*Log[2]^2*Zeta[3] + 1305*Pi*Zeta[5])}
