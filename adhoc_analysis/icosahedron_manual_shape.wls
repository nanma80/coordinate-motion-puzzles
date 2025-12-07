#!/usr/bin/env wolframscript

(* 
derived from the generic script perp_constraints. Focus on Icosahedron
 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

glueFirstPair = False;
glueFirstPair = True;

computeDimensions = False;
computeDimensions = True;

gluedGroups = {
    (* 10 pairs, dimension 1 *)
    (* {10, 12}, {18, 20}, {2, 6}, {14, 17}, {5, 8}, {9, 11}, {1, 3}, {4, 7}, {13, 16}, {15, 19} *)

    (* upper and bottom caps, dimension 3 *)
    (* {12, 18, 6, 8, 14}, {3, 7, 9, 13, 15} *)

    (* upper and bottom caps and 1 pairs, dimension 2 *)
    (* {12, 18, 6, 8, 14}, {3, 7, 9, 13, 15}, {1, 2} *)

    (* upper and bottom caps and 4 pairs. 2 pairs force 2 extra pairs, 5 & 16 by themselves. 8 way CM *)
    (* {12, 18, 6, 8, 14}, {3, 7, 9, 13, 15}, {1, 2}, {4, 10}, {11, 17}, {19, 20} *)

    (* 3 groups of 5 + 5 individual. 8 way CM *)
    {12, 18, 6, 8, 14}, {10, 16, 13, 7, 4}, {15, 19, 17, 11, 9}
};


If[forceUniformExpansion && glueFirstPair,
    Print["Cannot both force uniform expansion and glue first pair at the same time. Will not glue first pair."];
    glueFirstPair = False;
];

If[computeDimensions && (forceUniformExpansion || glueFirstPair ),
    Print["Computing dimensions. We are setting uniform expansion and gluing to false"];
    (* glueFirstPair = False; *)
    forceUniformExpansion = False;
];

phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;
axisNorthPole = {0, 3+5 phi, 2+3 phi};

(* 

First 2 vectors: two pieces. The first piece goes over (outside of) the second
Third vector: norm of the contacting surface of the two pieces, going from second to first piece

Coordinates are from a vZome file

central surface, small right triangle, normal vector:
Cross[ {0,-1-2phi,3+5phi} - {0,0,3+5phi}, {2+3phi,0,3+5phi} - {0,0,3+5phi}]
same direction as {0,0,1}

pyramid normal vector:
Cross[{2+3phi,0,3+5phi} - {-2-3phi,0,3+5phi}, {0,0,0} - {-2-3phi,0,3+5phi} ]
same direction as {0, 1, 0}

central divider normal vector:
Cross[ {0, 1+2 phi, 3+5 phi}-{0,0,3+5phi}, {0,-1-phi,4+7phi}-{0,0,3+5phi} ]
same direction as {1,0,0}

We are ignoring the outer faces of the legs because I don't think they are active.
 *)

upperFaceCenter = {0, 1+2 phi, 3+5 phi};
lowerFaceCenter = {0, -1-2 phi, 3+5 phi};
lowerLeftFaceCenter = {-2-3phi, -2-3phi, 2+3phi};

faceTemplates = {
    {upperFaceCenter, lowerFaceCenter, {0,0,1} } (* central surface, small right triangle *)
    , { upperFaceCenter, lowerFaceCenter, {0, 1, 0} } (* pyramid to origin *)
    , { upperFaceCenter, lowerFaceCenter, {1, 0, 0} } (* central divider *)
    , { upperFaceCenter, lowerLeftFaceCenter, {1, 0, 0} } (* central divider, with farther pieces *)
}//N;

faces = {};
rotations = { IdentityMatrix[3] };
generators = { RotationMatrix[2 Pi/5, {phi, 1, 0}], RotationMatrix[Pi, {0, 0, 1}], RotationMatrix[Pi, {1, 0, 0}] };

For[round = 1, round <= 3, round++,
    For[gid = 1, gid <= Length[generators], gid ++,
        generator = generators[[gid]];
        newRotations = {};
        For[rid = 1, rid <= Length[rotations], rid ++,
            rotation = rotations[[rid]];
            newRotation = generator.rotation;
            AppendTo[newRotations, newRotation];
        ];
        rotations = Union[rotations, newRotations, SameTest -> (Max[Abs[N[#1] - N[#2]]] < epsilon &)];
        Print[Length[rotations]];
    ];
];

For[tid = 1, tid <= Length[faceTemplates], tid ++,
    faceTemplate = faceTemplates[[tid]];
    For[rid = 1, rid <= Length[rotations], rid ++,
        rotation = rotations[[rid]];
        newFace = Map[#.rotation&, faceTemplate];
        newFace = FullSimplify[newFace];
        AppendTo[faces, newFace];
    ];    
];

pieces = Map[First, faces];
roundedPieces = Round[pieces, epsilon];
roundedPieces = Sort[Union[roundedPieces]];
uniquePieces = {};
For[pid = 1, pid <= Length[roundedPieces], pid++,
    uniquePiece = First[Select[pieces, EuclideanDistance[#, roundedPieces[[pid]]] < 10 * epsilon & ]];
    AppendTo[uniquePieces, uniquePiece];
];

Print["Rounded piece count"];
Print[Length[roundedPieces]];

Print["piece count"];
Print[Length[pieces]];

Print["unique piece count"];
Print[Length[uniquePieces]];

pieceLookup = <||>;

Print["Rounded pieces"];
(* Print[pieces]; *)
For[pid = 1, pid <= Length[roundedPieces], pid ++,
    piece = roundedPieces[[pid]];
    Print[{pid, piece, piece.axisNorthPole } ];
    pieceLookup[piece] = pid;
];

minFaceCenterDistance = 10000.0;
For[fid = 2, fid <= Length[uniquePieces], fid++,
    distance = N[ EuclideanDistance[uniquePieces[[1]], uniquePieces[[fid]]] ];
    minFaceCenterDistance = Min[minFaceCenterDistance, distance];
];

neighborFacePairs = {};

For[f1id = 1, f1id <= Length[uniquePieces], f1id++,
    For[f2id = f1id + 1, f2id <= Length[uniquePieces], f2id++,
        distance = N[ EuclideanDistance[uniquePieces[[f1id]], uniquePieces[[f2id]]] ];
        If[distance < minFaceCenterDistance * 1.01,
            Print["neighboring faces: ", f1id, " - ", f2id];
            AppendTo[neighborFacePairs, {f1id, f2id}];
        ];
    ];
];


(* Exit[]; *)

(* Print[pieceLookup]; *)

pieceCoordinates = Array[x, {Length[roundedPieces], 3}];

pieceCoordinatesFlat = Flatten[pieceCoordinates];

faceCount = Length[uniquePieces];
pieceCount = faceCount; (* for convenience *)

(* faces = Union[faces]; *)

Print["face (constraint) count"];
Print[Length[faces]];

conditions = {};
If[!computeDimensions,
    conditions = {
        pieceCoordinates[[1]][[1]] == -1.0
    };
];

If[forceUniformExpansion,
    conditions = {};
    For[pid = 1, pid <= Length[uniquePieces], pid++,
        AppendTo[conditions, pieceCoordinates[[pid]][[1]] == uniquePieces[[pid]][[1]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[2]] == uniquePieces[[pid]][[2]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[3]] == uniquePieces[[pid]][[3]] ];
    ];
];

(* hack: force piece #10 to move along its direction *)
(* conditions = {};
pid = 10;
AppendTo[conditions, pieceCoordinates[[pid]][[1]] == uniquePieces[[pid]][[1]] ];
AppendTo[conditions, pieceCoordinates[[pid]][[2]] == uniquePieces[[pid]][[2]] ];
AppendTo[conditions, pieceCoordinates[[pid]][[3]] == uniquePieces[[pid]][[3]] ]; *)


gluedPairs = Join @@ (Thread[{#[[1]], Rest[#]}] & /@ gluedGroups);

If[glueFirstPair,
    For[pid = 1, pid <= Length[gluedPairs], pid++,
        gluedPair = gluedPairs[[pid]];
        distance = EuclideanDistance[ uniquePieces[[gluedPair[[1]]]], uniquePieces[[gluedPair[[2]]]] ];
        gluedPairCenter = Mean[{uniquePieces[[gluedPair[[1]]]], uniquePieces[[gluedPair[[2]]]]}];
        inner = axisNorthPole.gluedPairCenter;
        Print["gluedPair: ", gluedPair, " distance: ", distance, " inner: ", inner];
        For[dim = 1, dim <= 3, dim++,
            AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] - pieceCoordinates[[gluedPair[[2]]]][[dim]] == 0 ];
        ];
    ];
];

For[fid = 1, fid <= Length[faces], fid++,
    face = faces[[fid]];
    piece1Id = pieceLookup[Round[face[[1]], epsilon]];
    piece2Id = pieceLookup[Round[face[[2]], epsilon]];
    condition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).face[[3]] >= - epsilon * 0.00001 ;
    AppendTo[conditions, condition];
];

For[dim = 1, dim <= 3, dim++,
    AppendTo[conditions, Total[Table[pieceCoordinates[[k]][[dim]], {k, 1, Length[roundedPieces]}]] == 0];
];

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

    (* not computing dimensions *)
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
        Print[sid];
        solution = solutions[[sid]];
        values = pieceCoordinates/.solution;
        values = Round[values, 0.000001];
        Print[values];
        (* values = Round[values, 1]; *)
        Print["piece id, piece, and motion"];
        pieceAndMotion = {};
        For[pid = 1, pid <= Length[roundedPieces], pid++,
            AppendTo[pieceAndMotion, {pid, roundedPieces[[pid]], values[[pid]]}];
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


(* 
with three constraints, up to central divider
compute dimension: 9. Already as same as the prep analysis.

with four constraints, up to central divider with farther pieces, the dimension = 9.

Rounded pieces
{1, {-11.090169944000001, 0., -4.236067977}, -4.236067977}
{2, {-11.090169944000001, 0., 4.236067977}, 4.236067977}, front left
{3, {-6.854101966, -6.854101966, -6.854101966}, -6.854101966}
{4, {-6.854101966, -6.854101966, 6.854101966}, 6.854101966}
{5, {-6.854101966, 6.854101966, -6.854101966}, -6.854101966}: lower up left back
{6, {-6.854101966, 6.854101966, 6.854101966}, 6.854101966}: up left front
{7, {-4.236067977, -11.090169944000001, 0.}, 0.}
{8, {-4.236067977, 11.090169944000001, 0.}, 0.}: up left back
{9, {0., -4.236067977, -11.090169944000001}, -11.090169944000001}, lower back center
{10, {0., -4.236067977, 11.090169944000001}, 11.090169944000001}, lower front center
{11, {0., 4.236067977, -11.090169944000001}, -11.090169944000001}, upper back center
{12, {0., 4.236067977, 11.090169944000001}, 11.090169944000001}, upper front center
{13, {4.236067977, -11.090169944000001, 0.}, 0.}
{14, {4.236067977, 11.090169944000001, 0.}, 0.}: up right back
{15, {6.854101966, -6.854101966, -6.854101966}, -6.854101966}
{16, {6.854101966, -6.854101966, 6.854101966}, 6.854101966}
{17, {6.854101966, 6.854101966, -6.854101966}, -6.854101966}: lower up right back
{18, {6.854101966, 6.854101966, 6.854101966}, 6.854101966}: up right front
{19, {11.090169944000001, 0., -4.236067977}, -4.236067977}
{20, {11.090169944000001, 0., 4.236067977}, 4.236067977}, front right


upper and lower caps, not very symmetrical:
piece id, piece, and motion
{1, {-11.090169944000001, 0., -4.236067977}, {-6.854101999999999, 2.6180339999999998, 1.}}
{2, {-11.090169944000001, 0., 4.236067977}, {-6.854101999999999, 2.6180339999999998, 1.}}
{5, {-6.854101966, 6.854101966, -6.854101966}, {-3.6180339999999998, 7.854101999999999, -1.}}
{4, {-6.854101966, -6.854101966, 6.854101966}, {-1.618034, -5.854102, 4.2360679999999995}}
{10, {0., -4.236067977, 11.090169944000001}, {-1.618034, -5.854102, 4.2360679999999995}}
{6, {-6.854101966, 6.854101966, 6.854101966}, {-1.618034, 11.090169999999999, 4.2360679999999995}}
{8, {-4.236067977, 11.090169944000001, 0.}, {-1.618034, 11.090169999999999, 4.2360679999999995}}
{12, {0., 4.236067977, 11.090169944000001}, {-1.618034, 11.090169999999999, 4.2360679999999995}}
{14, {4.236067977, 11.090169944000001, 0.}, {-1.618034, 11.090169999999999, 4.2360679999999995}}
{18, {6.854101966, 6.854101966, 6.854101966}, {-1.618034, 11.090169999999999, 4.2360679999999995}}
{3, {-6.854101966, -6.854101966, -6.854101966}, {1.618034, -11.090169999999999, -4.2360679999999995}}
{7, {-4.236067977, -11.090169944000001, 0.}, {1.618034, -11.090169999999999, -4.2360679999999995}}
{9, {0., -4.236067977, -11.090169944000001}, {1.618034, -11.090169999999999, -4.2360679999999995}}
{13, {4.236067977, -11.090169944000001, 0.}, {1.618034, -11.090169999999999, -4.2360679999999995}}
{15, {6.854101966, -6.854101966, -6.854101966}, {1.618034, -11.090169999999999, -4.2360679999999995}}
{11, {0., 4.236067977, -11.090169944000001}, {1.618034, 5.854102, -4.2360679999999995}}
{17, {6.854101966, 6.854101966, -6.854101966}, {1.618034, 5.854102, -4.2360679999999995}}
{16, {6.854101966, -6.854101966, 6.854101966}, {3.6180339999999998, -7.854101999999999, 1.}}
{19, {11.090169944000001, 0., -4.236067977}, {6.854101999999999, -2.6180339999999998, -1.}}
{20, {11.090169944000001, 0., 4.236067977}, {6.854101999999999, -2.6180339999999998, -1.}}
Groups
{2, 1, 2, 5, 5, 2, 1, 2}
8

upper and lower caps, also force piece #10 to move along its direction to force symmetry
piece id, piece, and motion
{1, {-11.090169944000001, 0., -4.236067977}, {-11.090169999999999, 0., -4.2360679999999995}}
{2, {-11.090169944000001, 0., 4.236067977}, {-11.090169999999999, 0., 4.2360679999999995}}
{4, {-6.854101966, -6.854101966, 6.854101966}, {-6.854101999999999, -6.854101999999999, 6.854101999999999}}
{5, {-6.854101966, 6.854101966, -6.854101966}, {-6.854101999999999, 6.854101999999999, -6.854101999999999}}
{3, {-6.854101966, -6.854101966, -6.854101966}, {0., -17.944271999999998, -11.090169999999999}}
{7, {-4.236067977, -11.090169944000001, 0.}, {0., -17.944271999999998, -11.090169999999999}}
{9, {0., -4.236067977, -11.090169944000001}, {0., -17.944271999999998, -11.090169999999999}}
{13, {4.236067977, -11.090169944000001, 0.}, {0., -17.944271999999998, -11.090169999999999}}
{15, {6.854101966, -6.854101966, -6.854101966}, {0., -17.944271999999998, -11.090169999999999}}
{10, {0., -4.236067977, 11.090169944000001}, {0., -4.2360679999999995, 11.090169999999999}}
{11, {0., 4.236067977, -11.090169944000001}, {0., 4.2360679999999995, -11.090169999999999}}
{6, {-6.854101966, 6.854101966, 6.854101966}, {0., 17.944271999999998, 11.090169999999999}}
{8, {-4.236067977, 11.090169944000001, 0.}, {0., 17.944271999999998, 11.090169999999999}}
{12, {0., 4.236067977, 11.090169944000001}, {0., 17.944271999999998, 11.090169999999999}}
{14, {4.236067977, 11.090169944000001, 0.}, {0., 17.944271999999998, 11.090169999999999}}
{18, {6.854101966, 6.854101966, 6.854101966}, {0., 17.944271999999998, 11.090169999999999}}
{16, {6.854101966, -6.854101966, 6.854101966}, {6.854101999999999, -6.854101999999999, 6.854101999999999}}
{17, {6.854101966, 6.854101966, -6.854101966}, {6.854101999999999, 6.854101999999999, -6.854101999999999}}
{19, {11.090169944000001, 0., -4.236067977}, {11.090169999999999, 0., -4.2360679999999995}}
{20, {11.090169944000001, 0., 4.236067977}, {11.090169999999999, 0., 4.2360679999999995}}
Groups
{1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1}
12

Working 10 pairs:
piece id, piece, and motion
{1, {-11.090169944000001, 0., -4.236067977}, {-1., 0., -0.38196599999999997}}
{3, {-6.854101966, -6.854101966, -6.854101966}, {-1., 0., -0.38196599999999997}}
{2, {-11.090169944000001, 0., 4.236067977}, {-1., 0., 0.38196599999999997}}
{6, {-6.854101966, 6.854101966, 6.854101966}, {-1., 0., 0.38196599999999997}}
{4, {-6.854101966, -6.854101966, 6.854101966}, {-0.618034, -0.618034, 0.618034}}
{7, {-4.236067977, -11.090169944000001, 0.}, {-0.618034, -0.618034, 0.618034}}
{5, {-6.854101966, 6.854101966, -6.854101966}, {-0.618034, 0.618034, -0.618034}}
{8, {-4.236067977, 11.090169944000001, 0.}, {-0.618034, 0.618034, -0.618034}}
{10, {0., -4.236067977, 11.090169944000001}, {0., -0.38196599999999997, 1.}}
{12, {0., 4.236067977, 11.090169944000001}, {0., -0.38196599999999997, 1.}}
{9, {0., -4.236067977, -11.090169944000001}, {0., 0.38196599999999997, -1.}}
{11, {0., 4.236067977, -11.090169944000001}, {0., 0.38196599999999997, -1.}}
{13, {4.236067977, -11.090169944000001, 0.}, {0.618034, -0.618034, 0.618034}}
{16, {6.854101966, -6.854101966, 6.854101966}, {0.618034, -0.618034, 0.618034}}
{14, {4.236067977, 11.090169944000001, 0.}, {0.618034, 0.618034, -0.618034}}
{17, {6.854101966, 6.854101966, -6.854101966}, {0.618034, 0.618034, -0.618034}}
{15, {6.854101966, -6.854101966, -6.854101966}, {1., 0., -0.38196599999999997}}
{19, {11.090169944000001, 0., -4.236067977}, {1., 0., -0.38196599999999997}}
{18, {6.854101966, 6.854101966, 6.854101966}, {1., 0., 0.38196599999999997}}
{20, {11.090169944000001, 0., 4.236067977}, {1., 0., 0.38196599999999997}}
Groups
{2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
10

dimension:
no gluing: 9

(* 10 pairs *)
{10, 12}, {18, 20}, {2, 6}, {14, 17}, {5, 8}, {9, 11}, {1, 3}, {4, 7}, {13, 16}, {15, 19}
dimension: 1

(* upper and bottom caps *)
{12, 18, 6, 8, 14}, {3, 7, 9, 13, 15}
dimension: 3

(* upper and bottom caps and 1 pairs *)
{12, 18, 6, 8, 14}, {3, 7, 9, 13, 15}, {1, 2}
dimension: 2

(* upper and bottom caps and 2 pairs *)
{12, 18, 6, 8, 14}, {3, 7, 9, 13, 15}, {1, 2}, {4, 10}
dimension: 1

(* upper and bottom caps and 4 pairs. 2 pairs force 2 extra pairs, 5 & 16 by themselves *)
{12, 18, 6, 8, 14}, {3, 7, 9, 13, 15}, {1, 2}, {4, 10}, {11, 17}, {19, 20}
dimension: 1

(* 3 groups of 5 + 5 individual *)
{12, 18, 6, 8, 14}, {10, 16, 13, 7, 4}, {15, 19, 17, 11, 9}
dimension: 1

 *)