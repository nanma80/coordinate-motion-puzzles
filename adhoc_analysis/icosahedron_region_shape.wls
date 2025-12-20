#!/usr/bin/env wolframscript

(* 
derived from the generic script perp_constraints. Focus on Icosahedron
What's the shape of the 9D region?
9D: cone with 8D base.
What's the shape of the 8D base?

noOverlapCondition creates the boundaries of the 8D region. There are 30 of them.
crossProduct condition restrict dimensions
If I flip one noOverlapCondition (pairId == 1), < -0.1 (up to -0.5)
there's still a solution. So every condition is active. 30 facets. I don't know any 8D shape with 30 facets.

 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

glueFirstPair = False;
(* glueFirstPair = True; *)

computeDimensions = False;
(* computeDimensions = True; *)

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

(* shape = "RhombicDodecahedron"; *)
(* shape = "RhombicTriacontahedron"; *)

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)

(* shape = "Octahedron"; *)
shape = "Icosahedron";

(* shape = "TetrakisHexahedron"; *)


Print["shape is: ", shape];

vertices = PolyhedronData[shape, "VertexCoordinates"];
faces = PolyhedronData[shape, "FaceIndices"];

faceCenters = Mean /@ (vertices[[#]] & /@ faces);
faceCenters = N[faceCenters];
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
    AppendTo[conditions, pieceCoordinates[[1]][[3]] == 1.0 ];
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

neighborFacePairs = {};

For[f1id = 1, f1id <= faceCount, f1id++,
    For[f2id = f1id + 1, f2id <= faceCount, f2id++,
        distance = N[ EuclideanDistance[faceCenters[[f1id]], faceCenters[[f2id]]] ];
        If[distance < minFaceCenterDistance + epsilon,
            (* Print["neighboring faces: ", f1id, " - ", f2id]; *)
            AppendTo[neighborFacePairs, {f1id, f2id}];
        ];
    ];
];

Print["neighborFacePairs"];
Print[Sort[neighborFacePairs]];
Print["Count: ", Length[neighborFacePairs]];

For[pairId = 1, pairId <= Length[neighborFacePairs], pairId++,
    pair = neighborFacePairs[[pairId]];
    piece1Id = pair[[1]];
    piece2Id = pair[[2]];
    vector12 = N[faceCenters[[piece1Id]] - faceCenters[[piece2Id]]];

    noOverlapCondition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).vector12 >= - epsilon * 0.00001 ;
    oppositeOfNoOverlapCondition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).vector12 < - 0.1 ;
    crossProduct = Cross[pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]], vector12];

    If[pairId == 23,
        AppendTo[conditions, oppositeOfNoOverlapCondition],
        AppendTo[conditions, noOverlapCondition]
    ];
    (* AppendTo[conditions, noOverlapCondition]; *)
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, - crossProduct[[dim]] >= -epsilon * 0.00001 ];
        AppendTo[conditions, crossProduct[[dim]] >= -epsilon * 0.00001 ];
    ];
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
    {12, 19, 20}, {1, 2, 3}
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
=====
Gluing pairs of pieces:
{{1, 2}, {1, 3}, {1, 4}, {1, 5}, {6, 7}, {6, 8}, {6, 9}, {6, 10}, {11, 19}, {12, 20}}

piece centers and motions
{12, {-0.742344243, 0., 0.141775135}, {-7.472135955000001, -3.0776835370000004, 3.2360679770000003}}
{20, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-7.472135955000001, -3.0776835370000004, 3.2360679770000003}}
{11, {-0.22939698700000002, 0.70601133, 0.141775135}, {-4.236067977, 6.881909602, -3.2360679770000003}}
{19, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-4.236067977, 6.881909602, -3.2360679770000003}}
{13, {-0.22939698700000002, -0.70601133, 0.141775135}, {-2.2360679770000003, -6.881909602, 7.236067977}}
{1, {-0.141775135, 0.43633899800000003, 0.6005691080000001}, {-1., -3.0776835370000004, 13.708203932}}
{2, {-0.458793973, 0., 0.6005691080000001}, {-1., -3.0776835370000004, 13.708203932}}
{3, {-0.141775135, -0.43633899800000003, 0.6005691080000001}, {-1., -3.0776835370000004, 13.708203932}}
{4, {0.371172121, -0.26967233100000004, 0.6005691080000001}, {-1., -3.0776835370000004, 13.708203932}}
{5, {0.371172121, 0.26967233100000004, 0.6005691080000001}, {-1., -3.0776835370000004, 13.708203932}}
{6, {0.141775135, -0.43633899800000003, -0.6005691080000001}, {1., 3.0776835370000004, -13.708203932}}
{7, {0.458793973, 0., -0.6005691080000001}, {1., 3.0776835370000004, -13.708203932}}
{8, {0.141775135, 0.43633899800000003, -0.6005691080000001}, {1., 3.0776835370000004, -13.708203932}}
{9, {-0.371172121, 0.26967233100000004, -0.6005691080000001}, {1., 3.0776835370000004, -13.708203932}}
{10, {-0.371172121, -0.26967233100000004, -0.6005691080000001}, {1., 3.0776835370000004, -13.708203932}}
{18, {0.22939698700000002, 0.70601133, -0.141775135}, {2.2360679770000003, 6.881909602, -7.236067977}}
{14, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {4.236067977, -6.881909602, 3.2360679770000003}}
{16, {0.22939698700000002, -0.70601133, -0.141775135}, {4.236067977, -6.881909602, 3.2360679770000003}}
{15, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {7.472135955000001, 3.0776835370000004, -3.2360679770000003}}
{17, {0.742344243, 0., -0.141775135}, {7.472135955000001, 3.0776835370000004, -3.2360679770000003}}
Groups
{2, 2, 1, 5, 5, 1, 2, 2}
8

This motion doesn't have much symmetry. The top and bottom pieces form caps but they don't move straight up and down.
=====
5-fold sym motion, 5+5 +1 +1... +1
{1, 2}, {1, 3}, {1, 4}, {1, 5}, {6, 7}, {6, 8}, {6, 9}, {6, 10}
Also, piece 1 going up, strictly

piece centers and motions
{12, {-0.742344243, 0., 0.141775135}, {-0.552786405, 0., 0.105572809}}
{20, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-0.447213595, -0.32491969600000004, -0.105572809}}
{19, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-0.447213595, 0.32491969600000004, -0.105572809}}
{13, {-0.22939698700000002, -0.70601133, 0.141775135}, {-0.17082039300000001, -0.525731112, 0.105572809}}
{11, {-0.22939698700000002, 0.70601133, 0.141775135}, {-0.17082039300000001, 0.525731112, 0.105572809}}
{6, {0.141775135, -0.43633899800000003, -0.6005691080000001}, {0., 0., -1.}}
{7, {0.458793973, 0., -0.6005691080000001}, {0., 0., -1.}}
{8, {0.141775135, 0.43633899800000003, -0.6005691080000001}, {0., 0., -1.}}
{9, {-0.371172121, 0.26967233100000004, -0.6005691080000001}, {0., 0., -1.}}
{10, {-0.371172121, -0.26967233100000004, -0.6005691080000001}, {0., 0., -1.}}
{1, {-0.141775135, 0.43633899800000003, 0.6005691080000001}, {0., 0., 1.}}
{2, {-0.458793973, 0., 0.6005691080000001}, {0., 0., 1.}}
{3, {-0.141775135, -0.43633899800000003, 0.6005691080000001}, {0., 0., 1.}}
{4, {0.371172121, -0.26967233100000004, 0.6005691080000001}, {0., 0., 1.}}
{5, {0.371172121, 0.26967233100000004, 0.6005691080000001}, {0., 0., 1.}}
{16, {0.22939698700000002, -0.70601133, -0.141775135}, {0.17082039300000001, -0.525731112, -0.105572809}}
{18, {0.22939698700000002, 0.70601133, -0.141775135}, {0.17082039300000001, 0.525731112, -0.105572809}}
{14, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {0.447213595, -0.32491969600000004, 0.105572809}}
{15, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {0.447213595, 0.32491969600000004, 0.105572809}}
{17, {0.742344243, 0., -0.141775135}, {0.552786405, 0., -0.105572809}}
Groups
{1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1}
12

The middle layer pieces move outwards by the norm direction of the faces. The top and bottom pieces form two caps and move up and down.


=======
10 pairs:
{1, 11}, {2, 12}, {3, 13}, {4, 14}, {5, 15}, {6, 16}, {7, 17}, {8, 18}, {9, 19}, {10, 20}

piece centers and motions
{2, {-0.458793973, 0., 0.6005691080000001}, {-5.236067977, 0., 1.}}
{12, {-0.742344243, 0., 0.141775135}, {-5.236067977, 0., 1.}}
{10, {-0.371172121, -0.26967233100000004, -0.6005691080000001}, {-4.236067977, -3.0776835370000004, -1.}}
{20, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-4.236067977, -3.0776835370000004, -1.}}
{9, {-0.371172121, 0.26967233100000004, -0.6005691080000001}, {-4.236067977, 3.0776835370000004, -1.}}
{19, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-4.236067977, 3.0776835370000004, -1.}}
{3, {-0.141775135, -0.43633899800000003, 0.6005691080000001}, {-1.6180339890000002, -4.97979657, 1.}}
{13, {-0.22939698700000002, -0.70601133, 0.141775135}, {-1.6180339890000002, -4.97979657, 1.}}
{1, {-0.141775135, 0.43633899800000003, 0.6005691080000001}, {-1.6180339890000002, 4.97979657, 1.}}
{11, {-0.22939698700000002, 0.70601133, 0.141775135}, {-1.6180339890000002, 4.97979657, 1.}}
{6, {0.141775135, -0.43633899800000003, -0.6005691080000001}, {1.6180339890000002, -4.97979657, -1.}}
{16, {0.22939698700000002, -0.70601133, -0.141775135}, {1.6180339890000002, -4.97979657, -1.}}
{8, {0.141775135, 0.43633899800000003, -0.6005691080000001}, {1.6180339890000002, 4.97979657, -1.}}
{18, {0.22939698700000002, 0.70601133, -0.141775135}, {1.6180339890000002, 4.97979657, -1.}}
{4, {0.371172121, -0.26967233100000004, 0.6005691080000001}, {4.236067977, -3.0776835370000004, 1.}}
{14, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {4.236067977, -3.0776835370000004, 1.}}
{5, {0.371172121, 0.26967233100000004, 0.6005691080000001}, {4.236067977, 3.0776835370000004, 1.}}
{15, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {4.236067977, 3.0776835370000004, 1.}}
{7, {0.458793973, 0., -0.6005691080000001}, {5.236067977, 0., -1.}}
{17, {0.742344243, 0., -0.141775135}, {5.236067977, 0., -1.}}
Groups
{2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
10

The middle layer pieces move outwards by the norm direction of the faces. The top and bottom pieces follow the middle layer pieces they glue to.

Remark: unfortunately, the second way to pair 20 faces into 10 pairs cannot be used to disassemble

====
3-fold symmetry, 3 groups of 5s around one face
{{1, 11}, {3, 4}, {12, 20}}

center:
{2, {-0.458793973, 0., 0.6005691080000001}, {-0.7639320230000001, 0., 1.}}

three groups of fives
{9, {-0.371172121, 0.26967233100000004, -0.6005691080000001}, {-2.198212717, 0., -1.320714913}}
{10, {-0.371172121, -0.26967233100000004, -0.6005691080000001}, {-2.198212717, 0., -1.320714913}}
{12, {-0.742344243, 0., 0.141775135}, {-2.198212717, 0., -1.320714913}}
{19, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-2.198212717, 0., -1.320714913}}
{20, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-2.198212717, 0., -1.320714913}}

{3, {-0.141775135, -0.43633899800000003, 0.6005691080000001}, {0.8396425430000001, -2.2071310410000002, 1.}}
{4, {0.371172121, -0.26967233100000004, 0.6005691080000001}, {0.8396425430000001, -2.2071310410000002, 1.}}
{13, {-0.22939698700000002, -0.70601133, 0.141775135}, {0.8396425430000001, -2.2071310410000002, 1.}}
{14, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {0.8396425430000001, -2.2071310410000002, 1.}}
{16, {0.22939698700000002, -0.70601133, -0.141775135}, {0.8396425430000001, -2.2071310410000002, 1.}}

{1, {-0.141775135, 0.43633899800000003, 0.6005691080000001}, {0.8396425430000001, 2.2071310410000002, 1.}}
{5, {0.371172121, 0.26967233100000004, 0.6005691080000001}, {0.8396425430000001, 2.2071310410000002, 1.}}
{11, {-0.22939698700000002, 0.70601133, 0.141775135}, {0.8396425430000001, 2.2071310410000002, 1.}}
{15, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {0.8396425430000001, 2.2071310410000002, 1.}}
{18, {0.22939698700000002, 0.70601133, -0.141775135}, {0.8396425430000001, 2.2071310410000002, 1.}}

{17, {0.742344243, 0., -0.141775135}, {1.5567828910000001, 0., -0.43428069500000005}}

{6, {0.141775135, -0.43633899800000003, -0.6005691080000001}, {0.396425434, -0.84304904, -1.320714913}}

{8, {0.141775135, 0.43633899800000003, -0.6005691080000001}, {0.396425434, 0.84304904, -1.320714913}}

anti center
{7, {0.458793973, 0., -0.6005691080000001}, {1.008936415, 0., -1.320714913}}

Groups
{5, 1, 1, 1, 5, 5, 1, 1}
8

 *)