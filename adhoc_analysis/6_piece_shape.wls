#!/usr/bin/env wolframscript

(* rows are the allowed motions *)

rows = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    {-1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, -1},
    {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
    {-1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, -1}
};

For[rid = 1, rid <=Length[rows], rid++,
    Print[rid];
    Print[rows[[rid]]];
];

rowsRank = MatrixRank[rows, Tolerance -> 10^-8];
Print["rows rank: ", rowsRank];

mean = Mean[rows];
Print[mean];

relativeMotion = Map[ # - mean& , rows];
Map[Print, relativeMotion];

inner = relativeMotion[[1]].relativeMotion[[2]];
Print[inner];
norm = Norm[relativeMotion[[1]]];
Print[norm];
angle = ArcCos[ inner/norm^2 ] * 180.0 / Pi;
Print[angle];

innerFromApex = rows[[1]].rows[[2]];
Print[innerFromApex];
normFromApex = Norm[rows[[1]]];
Print[normFromApex];
angleFromApex = ArcCos[ innerFromApex/normFromApex^2 ] * 180.0 / Pi;
Print[angleFromApex];

(* 
rows rank = 4

The shape of the 4D region is a 4D pyramid. The 3D base is a regular tetrahedron. So, the 4D shape is like a 5-cell, extended from one of the cells. However, this 5-cell isn't the regular 5-cell. The apex in the different 4-th dimension is not symmetrical with other vertices.

from mean:
-4
2*Sqrt[3]
109.47122063449069

from apex:
2
3*Sqrt[2]
83.62062979155719


 *)