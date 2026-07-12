#!/usr/bin/env wolframscript

(* 

 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

PairsToGlue = 0;
(* PairsToGlue = 1; *)
(* PairsToGlue = 2; *)

findVertices = False;
(* findVertices = True; *)

computeDimensions = False;
computeDimensions = True;

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
findVerticesBatchSize = 300;
(* findVerticesBatchSize = 400000; *)

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)

(* shape = "Octahedron"; *)
shape = "Icosahedron";

(* shape = "RhombicDodecahedron"; *)
(* shape = "RhombicTriacontahedron"; *)

(* shape = "Cuboctahedron"; *)
(* shape = "Icosidodecahedron"; *)

(* shape = "SquareAntiprism"; *)
(* shape = "RhombicIcosahedron"; *)

outputFolder = FileNameJoin[{"..", "output"}];
dataFolder = "data";
vertexFileName = shape <> "_diagonal_constraint_vertices.json";
vertexFileNameFull = FileNameJoin[{outputFolder, dataFolder, vertexFileName}];

Print["The shape is: ", shape];
isExtreme[v_, verts_] := Module[
  {others, n, vars, res, minval},
  others = Select[verts, EuclideanDistance[v, #] > 10 epsilon&];
  n = Length[others];
  If[n == 0, Return[True]];
  vars = Table[Unique["w"], {n}];
  res = Quiet @ NMinimize[
    {
      Norm[vars.others - v]^2,
      Total[vars] == 1.0 && And @@ Thread[vars >= 0]
    },
    vars
  ];
  (* Print[res]; *)
  minval = res[[1]];
  minval > 10 epsilon
];

vertices = PolyhedronData[shape, "VertexCoordinates"];
faces = PolyhedronData[shape, "FaceIndices"];
edges = PolyhedronData[shape, "EdgeIndices"];
(* Print["vertices: ", vertices];
Print["faces: ", faces];
Print["edges: ", edges]; *)
faceEdges = Table[
   Sort /@ Table[
     {v[[i]], v[[Mod[i, Length[v]] + 1]]},
     {i, Length[v]}
   ],
   {v, faces}
];
edgeToFaces = Table[
   edge = Sort[edges[[i]]];
   matchingFaces = First /@ Position[faceEdges, edge];
   {edge, matchingFaces},
   {i, Length[edges]}
];

(*
    each item in edgeToFaces describes an edge, like:
    {{{1, 2}, {13, 14}}
    It means the edge is between vertices with index 1 to vertext index 2. The edge is shared between faces 13 and 14
    It should work for convex polyhedra. No guarantee for 
*)

faceCenters = Mean /@ (vertices[[#]] & /@ faces);
faceCenters = N[faceCenters];
faceCentersFlat = Flatten[faceCenters];

(* For[fid = 1, fid <= Length[faceCenters], fid++,
    normalizationFactor = Norm[faceCenters[[1]]];
    Print[fid, ": ", faceCenters[[fid]]/normalizationFactor];
];
 *)
faceCount = Length[faces];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];

(* each face becomes a piece *)
pieceCoordinates = Array[x, {faceCount, 3}];
pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {};
offsetCondition = pieceCoordinatesFlat.faceCentersFlat == Norm[faceCentersFlat]^2;
If[!computeDimensions && !findVertices,
    AppendTo[conditions, offsetCondition];
];

If[forceUniformExpansion,
    conditions = {};
    For[pid = 1, pid <= faceCount, pid++,
        AppendTo[conditions, pieceCoordinates[[pid]][[1]] == faceCenters[[pid]][[1]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[2]] == faceCenters[[pid]][[2]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[3]] == faceCenters[[pid]][[3]] ];
    ];
];

(* center of mass constraint *)
For[dim = 1, dim <= 3, dim++,
    AppendTo[conditions, Total[Table[pieceCoordinates[[k]][[dim]], {k, 1, faceCount }]] == 0];
];

Print["edgeToFaces: ", edgeToFaces];

For[edgeId = 1, edgeId <= Length[edgeToFaces], edgeId++,
    edgeIndex = edgeToFaces[[edgeId]];
    edgeVertex1 = edgeIndex[[1]][[1]];
    edgeVertex2 = edgeIndex[[1]][[2]];
    face1 = edgeIndex[[2]][[1]];
    face2 = edgeIndex[[2]][[2]];
    edgeVector = N[vertices[[edgeVertex1]] - vertices[[edgeVertex2]]];
    normalVector = Cross[edgeVector, N[vertices[[edgeVertex1]]]];
    (* there's no guarantee normalVector will point to face1 or face2. We quickly check here to ensure it's pointing to face1 *)
    If[normalVector.( faceCenters[[face1]] - faceCenters[[face2]] ) < 0, normalVector = - normalVector];

    noOverlapCondition = ( pieceCoordinates[[face1]] - pieceCoordinates[[face2]] ).normalVector >= - epsilon * 0.00001 ;
    crossProduct = Cross[pieceCoordinates[[face1]] - pieceCoordinates[[face2]], normalVector];

    AppendTo[conditions, noOverlapCondition];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, - crossProduct[[dim]] >= -epsilon * 0.00001 ];
        AppendTo[conditions, crossProduct[[dim]] >= -epsilon * 0.00001 ];
    ];
];

(* for testing *)
(* facesToIgnore = {2, 4, 6, 7};
For[fid = 1, fid <= Length[facesToIgnore], fid++,
    pid = facesToIgnore[[fid]];
    Print["pid: ", pid];
    AppendTo[conditions, pieceCoordinates[[pid]][[1]] == 0 ];
    AppendTo[conditions, pieceCoordinates[[pid]][[2]] == 0 ];
    AppendTo[conditions, pieceCoordinates[[pid]][[3]] == 0 ];
]; *)

(* diagonalFacePairs contains only face indices. There's no "edge" *)
(* dim = 2 *)
(* diagonalFacePairs = { {2, 4}, {2, 5}, {4, 17}, {5, 17}, {7, 9}, {7, 10}, {12, 10}, {12, 9} }; *)

(* dim = 1 *)
(* diagonalFacePairs = { {2, 4}, {2, 5}, {4, 17}, {5, 17}, {7, 9}, {7, 10}, {12, 10}, {12, 9}, {4, 13}, {4, 16}, {10, 13}, {10, 16} }; *)

(* dim = 1, 2 diagonals per face *)
(* diagonalFacePairs = {
    {2, 4}, {2, 5}, {4, 17}, {5, 17}, {7, 9}, {7, 10}, {12, 10}, {12, 9}, 
    {4, 13}, {4, 16}, {10, 13}, {10, 16}, {5, 11}, {5, 18}, {9, 11}, {9, 18},
    {11, 2}, {11, 12}, {13, 2}, {13, 12}, {16, 7}, {16, 17}, {18, 7}, {18, 17}
}; *)

(* dim = 1, 1 diagonal per face *)
diagonalFacePairs = {
    {2, 4}, {5, 17}, {7, 10}, {12, 9}, 
    {4, 13}, {10, 16}, {5, 18}, {9, 11},
    {11, 12}, {13, 2}, {16, 7}, {18, 17}
};



Print["diagonalFacePairs: ", diagonalFacePairs];
rotatedConstraintAngle = Pi/2;

For[pairId = 1, pairId <= Length[diagonalFacePairs], pairId++,
    pair = diagonalFacePairs[[pairId]];
    face1 = pair[[1]];
    face2 = pair[[2]];
    vector12 = N[faceCenters[[face1]] - faceCenters[[face2]]];
    middlePoint = N[faceCenters[[face1]] + faceCenters[[face2]]];
    vector12rotated = RotationMatrix[rotatedConstraintAngle, middlePoint].vector12;

    noOverlapConditionRotated = ( pieceCoordinates[[face1]] - pieceCoordinates[[face2]] ).vector12rotated >= - epsilon * 0.00001 ;
    AppendTo[conditions, noOverlapConditionRotated];
];


For[gid = 1, gid <= PairsToGlue, gid++,
    edgeIndex = edgeToFaces[[gid]];
    gluedPair = edgeIndex[[2]];
    Print["Gluing a pair of pieces: "];
    Print[gluedPair];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];

Print["condition count:"];
Print[Length[conditions]];

computeDimensionsModule[] := Module[
    {},
    tightnessThreshold = 0.1;
    rows = {};

    For[cid = 1, cid <= Length[conditions], cid++,
        If[Mod[cid, 20] == 0,
            Print[{DateString["ISODateTime"], "Checking constraints", cid, Length[conditions]}];
        ];
        condition = conditions[[cid]];
        (* Print[condition]; *)
        If[condition[[0]] == GreaterEqual,
            newCondition = condition;
            newCondition[[2]] = tightnessThreshold;
            (* Print[{"modified condition", newCondition}]; *)
            addedConditions = Append[conditions, newCondition];
            combinedConditions = Apply[ And, addedConditions];
            solutions = FindInstance[
                combinedConditions,
                pieceCoordinatesFlat,
                Reals
            ];
            (* Print["solution count:"] *)
            (* Print[Length[solutions]]; *)
            If[Length[solutions] == 0,
                (* Print["Tight inequality!"]; *)
                row = Coefficient[condition[[1]], #] & /@ pieceCoordinatesFlat;
                AppendTo[rows, row];
            ];
        ];

        If[condition[[0]] == Equal,
            (* Print["Equality constraint!"]; *)
            row = Coefficient[condition[[1]], #] & /@ pieceCoordinatesFlat;
            AppendTo[rows, row];
        ];
    ];

    rowsRank = MatrixRank[rows, Tolerance -> 10^-6];
    Print["rows rank: ", rowsRank];

    coneDimension = Length[rows[[1]] ] - rowsRank;
    Print["coneDimension: ", coneDimension];


    sv = SingularValueList[N[rows], Tolerance -> 0];
    tol = 10^-10 * Max[sv];
    numericalRank = Count[sv, _?(# > tol &)];

    Print[sv];
    Print["numerical Rank: ", numericalRank];
    numericalConeDimension = Length[rows[[1]] ] - numericalRank;
    Print["numerical Dimension: ", numericalConeDimension];
    nullSpace = NullSpace[rows];
    Print["nullSpace dimension: ", Length[nullSpace]];
    nullSpace
];

If[computeDimensions,
    computeDimensionsModule[],

    If[findVertices,
        (* find vertices *)
        Print["To find vertices, first compute dimensions and a set of base vectors"];
        nullSpace = computeDimensionsModule[];
        AppendTo[conditions, offsetCondition];

        Print["Starting to find vertices of the feasible set"];
        discoveredVertices = {};

        For[directionId = 1, directionId <= findVerticesBatchSize, directionId++,
            If[Mod[directionId, 10000] == 0,
                Print[{DateString["ISODateTime"], "Find vertex using random directions", directionId, findVerticesBatchSize}];
            ];
            randomDirectionInNullSpace = RandomVariate[NormalDistribution[], Length[nullSpace]];
            randomDirection = randomDirectionInNullSpace.nullSpace;
            combinedConditions = Apply[ And, conditions];

            maxSolution = NMaximize[
                {randomDirection.pieceCoordinatesFlat, combinedConditions},
                pieceCoordinatesFlat
            ];
            If[Head[maxSolution] === NMaximize,
                Print["Find vertices: maxSolution is not evaluated. Skipping"];
                Print[maxSolution];
                Continue[];
            ];
            If[Length[maxSolution] == 0,
                Print["Find vertices: maxSolution is empty. Cannot find vertices. Skipping"];
                Continue[];
            ];
            If[maxSolution[[1]] == -Infinity || maxSolution[[1]] == Infinity,
                Print["Find vertices: feasible set is empty or the value is infinity. Cannot find vertices. Skipping"];
                Print[maxSolution];
                Continue[];
            ];
            preciseValues = pieceCoordinates/.maxSolution[[2]];
            If[NumberQ[ preciseValues[[1]][[1]] ] == False,
                Print["Find vertices: preciseValues is not a number. Skipping"];
                Print[preciseValues];
                Print[maxSolution];
                Continue[];
            ];
            values = Round[preciseValues, 10000 epsilon];
            pieceAndMotion = {};
            For[pid = 1, pid <= faceCount, pid++,
                AppendTo[pieceAndMotion, {pid, Round[ faceCenters[[pid]], epsilon ], values[[pid]]}];
            ];
            pieceAndMotion = SortBy[pieceAndMotion, Last];
            motionTally = Tally[Map[Last, pieceAndMotion]];
            tallyCounts = ReverseSort[ Map[#[[2]]&, motionTally] ];

            AppendTo[discoveredVertices, { preciseValues, tallyCounts } ];
        ];
        Print["before dedup"];
        Print[Length[discoveredVertices]];

        (* dedup *)
        discoveredVertices = Union[discoveredVertices, SameTest -> (EuclideanDistance[Flatten[#1[[1]]], Flatten[#2[[1]]] ] < epsilon * 1000&) ];
        Print["after dedup"];
        Print[Length[discoveredVertices]];
        For[vid = 1, vid <= Length[discoveredVertices], vid++,
            Print["groups: ", discoveredVertices[[vid]][[2]] ];
        ];

        allDiscoveredVerticesValues = Map[Flatten[First[#]]&, discoveredVertices];
        feasibleSetVertices = {};
        For[vid = 1, vid <= Length[discoveredVertices], vid++,
            If[Mod[vid, 10] == 0,
                Print[{DateString["ISODateTime"], "Checking vertex", vid, Length[discoveredVertices]}];
            ];
            vert = discoveredVertices[[vid]];
            
            vertexIsTrue = isExtreme[Flatten[vert[[1]]], allDiscoveredVerticesValues ];
            If[vertexIsTrue,
                Print["True vertex:"];
                Print[vert[[2]]];
                Print[vert[[1]]];
                vertexToSave = <||>;
                vertexToSave["grouping"] = vert[[2]];
                vertexToSave["pieceAndMotion"] = Transpose[ {faceCenters, vert[[1]]} ];
                AppendTo[feasibleSetVertices, vertexToSave],
                Print["False vertex. Ignore"];
                Print[vert[[2]]];
                Print[vert[[1]]];
            ];
        ];
        Print["after true vertex check. sort by groupings"];
        feasibleSetVertices = SortBy[feasibleSetVertices, #[["grouping"]]&];
        Print["after sorting. all groupings"];
        For[vid = 1, vid <= Length[feasibleSetVertices], vid++,
            Print[ feasibleSetVertices[[vid]][["grouping"]] ];
        ];
        groupings = Map[#[["grouping"]]&, feasibleSetVertices ];
        groupingsTally = Tally[groupings];
        Print["tally of groupings"];
        Map[Print, groupingsTally];

        Print["true vertex count"];
        Print[Length[feasibleSetVertices]];
        If[FileExistsQ[vertexFileNameFull], DeleteFile[vertexFileNameFull]];
        Export[vertexFileNameFull, feasibleSetVertices, "RawJSON", Compact -> True];
        ,

        (* find instance *)
        combinedConditions = Apply[ And, conditions];

        Print["Starting FindInstance"];
        solutions = FindInstance[
            combinedConditions,
            pieceCoordinatesFlat,
            Reals
        ];

        (* Print[solutions]; *)
        Print["Found these many solutions:"]
        Print[solutions//Length];

        For[sid = 1, sid <= Length[solutions], sid++,
            Print["solution #", sid];
            solution = solutions[[sid]];
            values = pieceCoordinates/.solution;
            values = Round[values, epsilon];
            (* Print[values]; *)
            (* values = Round[values, 1]; *)
            Print["piece centers and motions"];
            pieceAndMotion = {};
            For[pid = 1, pid <= faceCount, pid++,
                AppendTo[pieceAndMotion, {pid, Round[ faceCenters[[pid]], epsilon ], values[[pid]]}];
            ];
            pieceAndMotion = SortBy[pieceAndMotion, Last];
            Map[Print, pieceAndMotion];
            (* Print[pieceAndMotion]; *)

            motionTally = Tally[Map[Last, pieceAndMotion]];
            tallyCounts = Map[#[[2]]&, motionTally];
            Print["Groups"];
            Print[ tallyCounts ];
            Print["Group count"]
            Print[ tallyCounts//Length ];
        ];
    ];
];

Print["The shape is: ", shape];




(* 
results:

Original dimension, no diagonal: dim = 9

dim = 2:
diagonalFacePairs = { {2, 4}, {2, 5}, {4, 17}, {5, 17}, {7, 9}, {7, 10}, {12, 10}, {12, 9} };

dim = 1:
diagonalFacePairs = { {2, 4}, {2, 5}, {4, 17}, {5, 17}, {7, 9}, {7, 10}, {12, 10}, {12, 9}, {4, 13}, {4, 16}, {10, 13}, {10, 16} };

dim = 1:
2 diagonals per face
diagonalFacePairs = {
    {2, 4}, {2, 5}, {4, 17}, {5, 17}, {7, 9}, {7, 10}, {12, 10}, {12, 9}, 
    {4, 13}, {4, 16}, {10, 13}, {10, 16}, {5, 11}, {5, 18}, {9, 11}, {9, 18},
    {11, 2}, {11, 12}, {13, 2}, {13, 12}, {16, 7}, {16, 17}, {18, 7}, {18, 17}
};


(* dim = 1, 1 diagonal per face *)
diagonalFacePairs = {
    {2, 4}, {5, 17}, {7, 10}, {12, 9}, 
    {4, 13}, {10, 16}, {5, 18}, {9, 11},
    {11, 12}, {13, 2}, {16, 7}, {18, 17}
};


 *)