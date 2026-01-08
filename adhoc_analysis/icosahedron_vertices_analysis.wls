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

shape = "Icosahedron";

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

(* 131 *)
focusedGrouping = {5, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1};

(* 60 *)
(* focusedGrouping = {5, 5, 2, 2, 2, 1, 1, 1, 1}; *)

(* focusedGrouping = {5, 5, 2, 2, 1, 1, 1, 1, 1, 1}; *)



focusedVertices = {};
motionNorms = {};
countInFocusedGroup = 0;

For[vid = 1, vid <= Length[feasibleSetVertices], vid++,
    vertex = feasibleSetVertices[[vid]];
    If[vertex[["grouping"]] == focusedGrouping,
        countInFocusedGroup++;
        pieceAndMotion = vertex[["pieceAndMotion"]];
        motionVectors = Map[Last, pieceAndMotion];
        motion = Flatten[ motionVectors ];
        normMotion = Norm[motion];
        If[normMotion < 4, Continue[]];
        Print[countInFocusedGroup];

        Print[{countInFocusedGroup, Norm[motion]}];
        scaledMotionVectors = Map[#*Norm[#]&, motionVectors];
        Print[Total[scaledMotionVectors]];
        motionVectorTally = Tally[ Round[motionVectors, 10000 epsilon] ];
        
        Print[Select[motionVectorTally, #[[2]] == 5&]];
        (* Print[Select[motionVectorTally, #[[2]] == 1&]]; *)
        (* Print[motionVectorTally]; *)


        AppendTo[motionNorms, Norm[motion]];
        AppendTo[ focusedVertices,  motion];
    ];
];

Print[Length[focusedVertices]];

(* sameTestThreshold = epsilon * 100;
Print["sameTestThreshold: ", sameTestThreshold]; *)
(* dedupedVertices = Union[focusedVertices, SameTest -> (EuclideanDistance[#1, #2] < sameTestThreshold &) ];
Print["after dedup"];
Print[Length[dedupedVertices]];
 *)
(* Print[Tally[Round[motionNorms, 1000 epsilon]]]; *)

(* 
Three variants for grouping signature {5, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1}

{{4.062289000000001, 60}, {3.9085030000000005, 60}, {3.8221850000000006, 11}}


 *)
