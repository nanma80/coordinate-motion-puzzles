#!/usr/bin/env wolframscript

(* 

 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

PairsToGlue = 0;
(* PairsToGlue = 1; *)
(* PairsToGlue = 3; *)

findVertices = False;
(* findVertices = True; *)

computeDimensions = False;
(* computeDimensions = True; *)

If[forceUniformExpansion && PairsToGlue > 0,
    Print["Cannot both force uniform expansion and glue any pair of pieces at the same time. Will not glue any."];
    PairsToGlue = 0;
];

If[findVertices && (forceUniformExpansion || PairsToGlue > 0 ),
    Print["Finding vertices of the feasible set. We are setting uniform expansion and gluing to false"];
    PairsToGlue = 0;
    forceUniformExpansion = False;
];

If[computeDimensions && (forceUniformExpansion || PairsToGlue > 0 ),
    Print["Computing dimensions. We are setting uniform expansion and gluing to false"];
    PairsToGlue = 0;
    forceUniformExpansion = False;
    findVertices = False;
];

phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;
findVerticesBatchSize = 100;
getAngle[v1_, v2_] := ArcCos[v1.v2/Norm[v1]/Norm[v2]]/Pi*180.;

shape = "RhombicIcosahedron";

outputFolder = "..\\output";
dataFolder = "data";
vertexFileName = shape <> "_vertices.json";
vertexFileNameFull = FileNameJoin[{outputFolder, dataFolder, vertexFileName}];
Print[vertexFileNameFull];

feasibleSetVertices = Import[vertexFileNameFull, "RawJSON"];
Print["Load this many vertices of the feasible set: ", feasibleSetVertices//Length];

allGroupings = Map[#[["grouping"]]&, feasibleSetVertices];
Print["All groupings:"];
Map[Print, Tally[allGroupings]];

motionNorms = {};
flatMotions = {};

For[vid = 1, vid <= Length[feasibleSetVertices], vid++,
    vertex = feasibleSetVertices[[vid]];
    pieceAndMotion = vertex[["pieceAndMotion"]];
    motionVectors = Map[Last, pieceAndMotion];
    motion = Flatten[ motionVectors ];

    AppendTo[flatMotions, motion];
];

avgMotion = Total[flatMotions]/Length[flatMotions];
(* Print[avgMotion]; *)

getCentralAngleById[id1_, id2_] := Round[ getAngle[flatMotions[[id1]] - avgMotion, flatMotions[[id2]] - avgMotion], 0.001];

(* all vertices form a decagon in a plane *)
Print[getCentralAngleById[3, 1]];
Print[getCentralAngleById[1, 2]];
Print[getCentralAngleById[2, 4]];
Print[getCentralAngleById[4, 6]];
Print[getCentralAngleById[6, 8]];
Print[getCentralAngleById[8, 10]];
Print[getCentralAngleById[10, 9]];
Print[getCentralAngleById[9, 7]];
Print[getCentralAngleById[7, 5]];
Print[getCentralAngleById[5, 3]];

Print[];
Print["Two adjacent vertices"];
Print[feasibleSetVertices[[1]]];
Print[feasibleSetVertices[[2]]];

(* 

The two adjacent vertices are related in this way.

Each motion represented by a vertex has two triplets. Viewed from the north pole, the two triplets are different by 3/10 * 360, or 3/5 pi. Between these two triplets, one is closer to the north pole and the other is closer to the south pole.

Two adjacent vertices share a common triplet, but the other triplets are mirrored with respect to the common triplet, viewed from the north pole.

If one vertex has the triplets at angles 0 and 3/5 pi, an adjacent vertex has the triplets at angles 0 and 3/5 pi. Another adjacent vertex has the triplets at angles 3/5 pi and 6/5 pi. 

If we trace the triplets by adjacent vertices, we should draw a decagram (star polygon) {10/3}.

 *)