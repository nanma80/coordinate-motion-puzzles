#!/usr/bin/env wolframscript

(* 
derived from the generic script perp_constraints. 
It's an older version of it, using distances to calculate neighbor pairs
Focus on Icosahedron
Can we add extra constraints to make the dimension < 9 but > 1?
 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

glueFirstPair = False;
(* glueFirstPair = True; *)

computeDimensions = False;
computeDimensions = True;

If[forceUniformExpansion && glueFirstPair,
    Print["Cannot both force uniform expansion and glue first pair at the same time. Will not glue first pair."];
    glueFirstPair = False;
];

If[computeDimensions && (forceUniformExpansion || glueFirstPair ),
    Print["Computing dimensions. We are setting uniform expansion and gluing to false"];
    glueFirstPair = False;
    forceUniformExpansion = False;
];

phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;

shape = "Icosahedron";

Print["shape is: ", shape];

vertices = PolyhedronData[shape, "VertexCoordinates"];
faces = PolyhedronData[shape, "FaceIndices"];

faceCenters = Mean /@ (vertices[[#]] & /@ faces);
faceCenters = N[faceCenters];
faceCentersFlat = Flatten[faceCenters];

faceCount = Length[faceCenters];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];
For[fid = 1, fid <= Length[faceCenters], fid++,
    Print[fid, " ", faceCenters[[fid]]];
];

minFaceCenterDistance = 1000.0;

(* each face becomes a piece *)
pieceCoordinates = Array[x, {faceCount, 3}];
pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {};

If[!computeDimensions,
    (* AppendTo[conditions, pieceCoordinates[[1]][[3]] == 1.0 ]; *)
    AppendTo[conditions, pieceCoordinatesFlat.faceCentersFlat == 10.0];
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

For[fid = 2, fid <= faceCount, fid++,
    distance = N[ EuclideanDistance[faceCenters[[1]], faceCenters[[fid]]] ];
    minFaceCenterDistance = Min[minFaceCenterDistance, distance];
];

secondMinFaceCenterDistance = 1000.0;
For[fid = 2, fid <= faceCount, fid++,
    distance = N[ EuclideanDistance[faceCenters[[1]], faceCenters[[fid]]] ];
    If[distance > minFaceCenterDistance * 1.01,
        secondMinFaceCenterDistance = Min[secondMinFaceCenterDistance, distance];
    ];
];

neighborFacePairs = {};
secondNeighborFacePairs = {};

For[f1id = 1, f1id <= faceCount, f1id++,
    For[f2id = f1id + 1, f2id <= faceCount, f2id++,
        distance = N[ EuclideanDistance[faceCenters[[f1id]], faceCenters[[f2id]]] ];
        If[distance < minFaceCenterDistance + epsilon,
            (* Print["neighboring faces: ", f1id, " - ", f2id]; *)
            AppendTo[neighborFacePairs, {f1id, f2id}];
        ];

        If[distance < secondMinFaceCenterDistance + epsilon && distance > secondMinFaceCenterDistance - epsilon,
            AppendTo[secondNeighborFacePairs, {f1id, f2id}];
        ];
    ];
];

Print["neighborFacePairs"];
Print[Sort[neighborFacePairs]];

Print["secondNeighborFacePairs"];
Print[Sort[secondNeighborFacePairs]];

For[pairId = 1, pairId <= Length[neighborFacePairs], pairId++,
    pair = neighborFacePairs[[pairId]];
    piece1Id = pair[[1]];
    piece2Id = pair[[2]];
    vector12 = N[faceCenters[[piece1Id]] - faceCenters[[piece2Id]]];

    noOverlapCondition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).vector12 >= - epsilon * 0.00001 ;
    crossProduct = Cross[pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]], vector12];

    AppendTo[conditions, noOverlapCondition];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, - crossProduct[[dim]] >= -epsilon * 0.00001 ];
        AppendTo[conditions, crossProduct[[dim]] >= -epsilon * 0.00001 ];
    ];
];

For[pairId = 1, pairId <= Length[secondNeighborFacePairs], pairId++,
    pair = secondNeighborFacePairs[[pairId]];
    piece1Id = pair[[1]];
    piece2Id = pair[[2]];
    vector12 = N[faceCenters[[piece1Id]] - faceCenters[[piece2Id]]];
    edgeCenter = N[faceCenters[[piece1Id]] + faceCenters[[piece2Id]]];
    edgeCenter = edgeCenter/Norm[edgeCenter]*Norm[faceCenters[[1]]];
    distancesFromEdgeCenterToFaceCenters = Map[ EuclideanDistance[edgeCenter, #]&, faceCenters ];
    closestFaceId = Select[Range[Length[faceCenters]], distancesFromEdgeCenterToFaceCenters[[#]] < Min[distancesFromEdgeCenterToFaceCenters] + epsilon&][[1]];
    newFaceNorm = edgeCenter - faceCenters[[closestFaceId]];
    tripleProduct = Cross[vector12, newFaceNorm].edgeCenter;

    (* AppendTo[conditions, Sign[tripleProduct] * ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).newFaceNorm >= - epsilon * 0.00001]; *)
    AppendTo[conditions, ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).edgeCenter >= - epsilon * 0.00001];
];

gluedGroups = {
    (* 8 way CM *)
(*     {1, 2}, {1, 3}, {1, 4}, {1, 5}, {6, 7}, {6, 8}, {6, 9}, {6, 10}, 
    {11, 19}, {12, 20} *)

    (* 5-fold symmetry, piece 1 going up strictly, 5+5+1+1+1...+1*)
    (* {1, 2}, {1, 3}, {1, 4}, {1, 5}, {6, 7}, {6, 8}, {6, 9}, {6, 10} *)

    (* 5-fold symmetry, sum of pieces 1~5 going up. 10 pairs *)
    (* {1, 11}, {2, 12}, {3, 13}, {4, 14}, {5, 15}, {6, 16}, {7, 17}, {8, 18}, {9, 19}, {10, 20} *)

    (* 3-fold symmetry *)
    (* {1, 11}, {3, 4}, {12, 20} *)

    (* 3-fold symmetry, centered around face 12 *)
    (* {1,2,3,4,5}, {19,9,11,18,8}, {20,10,13,16,6} *)

    (* explore *)
    (* {12, 19, 20}, {1, 2, 3} *)

    (* test, no solution *)
    (* {1, 2} *)
};

gluedPairs = Join @@ (Thread[{#[[1]], Rest[#]}] & /@ gluedGroups);


If[glueFirstPair,
    Print["Gluing groups of pieces: "];
    Print[gluedGroups];

    Print["Gluing pairs of pieces: "];
    Print[gluedPairs];
    For[pid = 1, pid <= Length[gluedPairs], pid++,
        gluedPair = gluedPairs[[pid]];
        For[dim = 1, dim <= 3, dim++,
            AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
        ];
    ];
];

Print["condition count:"];
Print[Length[conditions]];
(* Print[conditions]; *)

If[computeDimensions,
    tightnessThreshold = 0.1;
    rows = {};

    For[cid = 1, cid <= Length[conditions], cid++,
        If[Mod[cid, 20] == 0,
            Print[{"Condition", cid, Length[conditions]}];
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
    Print["numerical Dimension: ", numericalConeDimension],

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
        Print[ tallyCounts//Length ];
    ];
];

Print["shape is: ", shape];




(* 

If we use first neighbors + second neighbors to construct cross product constraints

first neighbors: 30 pairs, corresponding to 30 edges of the icosahedron
second neighbors: 60 pairs

dimension = 1. Uniform expansion

piece centers and motions
{12, {-0.742344243, 0., 0.141775135}, {-1.236067977, 0., 0.236067977}}
{20, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-1., -0.7265425280000001, -0.236067977}}
{19, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-1., 0.7265425280000001, -0.236067977}}
{2, {-0.458793973, 0., 0.6005691080000001}, {-0.7639320230000001, 0., 1.}}
{10, {-0.371172121, -0.26967233100000004, -0.6005691080000001}, {-0.6180339890000001, -0.449027977, -1.}}
{9, {-0.371172121, 0.26967233100000004, -0.6005691080000001}, {-0.6180339890000001, 0.449027977, -1.}}
{13, {-0.22939698700000002, -0.70601133, 0.141775135}, {-0.38196601100000005, -1.175570505, 0.236067977}}
{11, {-0.22939698700000002, 0.70601133, 0.141775135}, {-0.38196601100000005, 1.175570505, 0.236067977}}
{3, {-0.141775135, -0.43633899800000003, 0.6005691080000001}, {-0.236067977, -0.7265425280000001, 1.}}
{1, {-0.141775135, 0.43633899800000003, 0.6005691080000001}, {-0.236067977, 0.7265425280000001, 1.}}
{6, {0.141775135, -0.43633899800000003, -0.6005691080000001}, {0.236067977, -0.7265425280000001, -1.}}
{8, {0.141775135, 0.43633899800000003, -0.6005691080000001}, {0.236067977, 0.7265425280000001, -1.}}
{16, {0.22939698700000002, -0.70601133, -0.141775135}, {0.38196601100000005, -1.175570505, -0.236067977}}
{18, {0.22939698700000002, 0.70601133, -0.141775135}, {0.38196601100000005, 1.175570505, -0.236067977}}
{4, {0.371172121, -0.26967233100000004, 0.6005691080000001}, {0.6180339890000001, -0.449027977, 1.}}
{5, {0.371172121, 0.26967233100000004, 0.6005691080000001}, {0.6180339890000001, 0.449027977, 1.}}
{7, {0.458793973, 0., -0.6005691080000001}, {0.7639320230000001, 0., -1.}}
{14, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {1., -0.7265425280000001, 0.236067977}}
{15, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {1., 0.7265425280000001, 0.236067977}}
{17, {0.742344243, 0., -0.141775135}, {1.236067977, 0., -0.236067977}}
Groups
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
20


 *)