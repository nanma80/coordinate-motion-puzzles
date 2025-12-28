#!/usr/bin/env wolframscript

(* 
Focus on antiprisms and see what motions are allowed
Dim = 3 for all. Why?
 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

PairsToGlue = 0;
PairsToGlue = 1;

computeDimensions = False;
(* computeDimensions = True; *)

If[forceUniformExpansion && PairsToGlue > 0,
    Print["Cannot both force uniform expansion and glue any pair of pieces at the same time. Will not glue any."];
    PairsToGlue = 0;
];

If[computeDimensions && (forceUniformExpansion || PairsToGlue > 0 ),
    Print["Computing dimensions. We are setting uniform expansion and gluing to false"];
    PairsToGlue = 0;
    forceUniformExpansion = False;
];

phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;

(* shape = "SquareAntiprism"; *)
(* shape = "PentagonalAntiprism"; *)
(* shape = "HexagonalAntiprism"; *)
(* shape = "HeptagonalAntiprism"; *)
(* shape = "OctagonalAntiprism"; *)
(* shape = "NonagonalAntiprism"; *)
shape = "DecagonalAntiprism";

gluedGroupsForAllShapes = <||>;

gluedGroupsForAllShapes["SquareAntiprism"] = {
    (* {2, 1, 1, 1, 1, 2, 1, 1}: 8 groups *)
    (* {3, 7}, {5, 8} *)

    (* {1, 1, 1, 1, 2, 1, 2, 1}: 8 groups *)
    {3, 7}, {6, 9}
};

gluedGroupsForAllShapes["PentagonalAntiprism"] = {
    (* {1, 2, 2, 1, 1, 2, 2, 1}: 8 groups *)
    {3, 8}, {4, 9}, {6, 10}, {7, 11}
};

gluedGroupsForAllShapes["HexagonalAntiprism"] = {
    (* {2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}: 12 groups *)
    {3, 9}, {6, 11}

    (* {1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1}: 12 groups *)
    (* {3, 9}, {7, 12} *)
};

gluedGroupsForAllShapes["HeptagonalAntiprism"] = {
    (* {1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1}: 12 groups *)
    {3, 10}, {4, 11}, {7, 13}, {8, 14}
};

gluedGroupsForAllShapes["OctagonalAntiprism"] = {
    (* {2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}: 16 groups *)
    {3, 11}, {7, 14}

    (* {1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}: 16 groups *)
    (* {3, 11}, {8, 15} *)
};

gluedGroupsForAllShapes["NonagonalAntiprism"] = {
    (* {1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1}: 16 groups *)
    {3, 12}, {4, 13}, {8, 16}, {9, 17}
};

gluedGroupsForAllShapes["DecagonalAntiprism"] = {
    (* {2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}: 20 groups *)
    {3, 13}, {8, 17}
};

gluedGroups = gluedGroupsForAllShapes[shape];

gluedPairs = Join @@ (Thread[{#[[1]], Rest[#]}] & /@ gluedGroups);

Print["The shape is: ", shape];

vertices = PolyhedronData[shape, "VertexCoordinates"];
faces = PolyhedronData[shape, "FaceIndices"];
edges = PolyhedronData[shape, "EdgeIndices"];
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
Print["Edges defined by indices of neighboring faces"];
Print[Sort[Map[Last, edgeToFaces]]];

faceCenters = Mean /@ (vertices[[#]] & /@ faces);
faceCenters = N[faceCenters];

faceCount = Length[faces];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];

For[fid = 1, fid <= Length[faceCenters], fid++,
    firstTwoDims = faceCenters[[fid]][[1]] + I faceCenters[[fid]][[2]];
    arg = If[Abs[firstTwoDims] > epsilon, Arg[ firstTwoDims ]/Pi*180.0, None];
    Print[{fid, arg, faceCenters[[fid]]}];
];

(* each face becomes a piece *)
pieceCoordinates = Array[x, {faceCount, 3}];
pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {};

If[!computeDimensions,
    If[Abs[faceCenters[[1]][[1]]] > epsilon, 
        AppendTo[conditions, pieceCoordinates[[1]][[1]] == 1.0 * Sign[ faceCenters[[1]][[1]] ] ],
        If[Abs[faceCenters[[1]][[2]]] > epsilon, 
            AppendTo[conditions, pieceCoordinates[[1]][[2]] == 1.0 * Sign[ faceCenters[[1]][[2]] ] ],
            AppendTo[conditions, pieceCoordinates[[1]][[3]] == 1.0 * Sign[ faceCenters[[1]][[3]] ] ]
        ];
    ];
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

If[PairsToGlue > 0,
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

If[computeDimensions,
    tightnessThreshold = 0.1;
    rows = {};

    For[cid = 1, cid <= Length[conditions], cid++,
        If[Mod[cid, 20] == 0,
            Print[{DateString["ISODateTime"], "Condition", cid, Length[conditions]}];
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
        Print["Group count"]
        Print[ tallyCounts//Length ];
    ];

];

Print["The shape is: ", shape];




(* 
results:

SquareAntiprism. Two maximal ways to glue. Equivalent. Side pieces: 2 doubles, 4 singles.
2 doubles separated by 1 single on one side, 3 singles on the other side

{3, 7}, {5, 8}
{up 45, down 90}, {up -135, down 180}
piece centers and motions
{5, {-0.40236892700000004, -0.40236892700000004, 0.140149403}, {-0.738878573, -0.246292858, 0.171572875}}
{8, {-0.569035594, 0., -0.140149403}, {-0.738878573, -0.246292858, 0.171572875}}
{4, {-0.40236892700000004, 0.40236892700000004, 0.140149403}, {-0.45032854200000005, 0.45032854200000005, 0.656854249}}
{2, {0., 0., 0.42044820800000005}, {-0.246292858, 0.246292858, 1.}}
{9, {0., -0.569035594, -0.140149403}, {0.246292858, -0.654364226, -0.5147186260000001}}
{1, {0., 0., -0.42044820800000005}, {0.246292858, -0.246292858, -1.}}
{3, {0.40236892700000004, 0.40236892700000004, 0.140149403}, {0.246292858, 0.738878573, 0.171572875}}
{7, {0., 0.569035594, -0.140149403}, {0.246292858, 0.738878573, 0.171572875}}
{6, {0.40236892700000004, -0.40236892700000004, 0.140149403}, {0.534842889, -0.534842889, -0.31370849900000003}}
{10, {0.569035594, 0., -0.140149403}, {0.654364226, -0.246292858, -0.5147186260000001}}
Groups
{2, 1, 1, 1, 1, 2, 1, 1}
Group count
8

{3, 7}, {6, 9}
{up 45, down 90}, {up -45, down -90}
Symmetric to {3, 7}, {5, 8}


PentagonalAntiprism. Just like the icosahedron puzzle, with two quintuples on top and bottom.
Side pieces: 4 doubles, 2 singles. Arranged like: double, double, single, double, double, single

{3, 8}, {4, 9}, {6, 10}, {7, 11}
piece centers and motions
{5, {-0.742344243, 0., 0.141775135}, {-0.5278640450000001, 0., 0.5278640450000001}}
{6, {-0.22939698700000002, -0.70601133, 0.141775135}, {-0.38196601100000005, -0.449027977, 0.236067977}}
{10, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-0.38196601100000005, -0.449027977, 0.236067977}}
{4, {-0.22939698700000002, 0.70601133, 0.141775135}, {-0.38196601100000005, 0.449027977, 0.236067977}}
{9, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-0.38196601100000005, 0.449027977, 0.236067977}}
{2, {0., 0., 0.425325404}, {-0.236067977, 0., 1.}}
{1, {0., 0., -0.425325404}, {0.236067977, 0., -1.}}
{7, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {0.38196601100000005, -0.449027977, -0.236067977}}
{11, {0.22939698700000002, -0.70601133, -0.141775135}, {0.38196601100000005, -0.449027977, -0.236067977}}
{3, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {0.38196601100000005, 0.449027977, -0.236067977}}
{8, {0.22939698700000002, 0.70601133, -0.141775135}, {0.38196601100000005, 0.449027977, -0.236067977}}
{12, {0.742344243, 0., -0.141775135}, {0.5278640450000001, 0., -0.5278640450000001}}
Groups
{1, 2, 2, 1, 1, 2, 2, 1}
Group count
8

HexagonalAntiprism. Side pieces: 2 doubles, 8 singles.
2 doubles separated by 3 singles on one side, 5 singles on the other side
{3, 9}, {6, 11}

piece centers and motions
{6, {-0.7886751350000001, -0.45534180100000005, 0.142599946}, {-0.463266057, -0.11462862100000001, 0.07179677000000001}}
{11, {-0.9106836030000001, 0., -0.142599946}, {-0.463266057, -0.11462862100000001, 0.07179677000000001}}
{5, {-0.7886751350000001, 0.45534180100000005, 0.142599946}, {-0.35686740100000003, 0.28245657100000005, 0.32050807600000003}}
{10, {-0.45534180100000005, 0.7886751350000001, -0.142599946}, {-0.234008815, 0.40531515700000004, 0.215390309}}
{12, {-0.45534180100000005, -0.7886751350000001, -0.142599946}, {-0.127610158, -0.45028452, -0.215390309}}
{7, {0., -0.9106836030000001, 0.142599946}, {-0.066180865, -0.46674445000000003, -0.176914536}}
{2, {0., 0., 0.427799839}, {-0.066180865, 0.11462862100000001, 1.}}
{4, {0., 0.9106836030000001, 0.142599946}, {-0.066180865, 0.45028452, 0.32050807600000003}}
{1, {0., 0., -0.427799839}, {0.066180865, -0.11462862100000001, -1.}}
{13, {0.45534180100000005, -0.7886751350000001, -0.142599946}, {0.22450567000000002, -0.388855227, -0.358983849}}
{3, {0.7886751350000001, 0.45534180100000005, 0.142599946}, {0.330904327, 0.34388586400000004, 0.07179677000000001}}
{9, {0.45534180100000005, 0.7886751350000001, -0.142599946}, {0.330904327, 0.34388586400000004, 0.07179677000000001}}
{8, {0.7886751350000001, -0.45534180100000005, 0.142599946}, {0.43730298300000003, -0.176057914, -0.176914536}}
{14, {0.9106836030000001, 0., -0.142599946}, {0.45376291300000005, -0.11462862100000001, -0.215390309}}
Groups
{2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}
Group count
12

HeptagonalAntiprism. Side pieces: 4 doubles, 6 singles.
2 doubles, 3 singles, 2 doubles, 3 singles
{3, 10}, {4, 11}, {7, 13}, {8, 14}

piece centers and motions
{6, {-1.076301277, 0., 0.143078866}, {-0.38708935600000005, 0.044302149000000006, 0.28620826400000005}}
{12, {-0.969713943, 0.466989623, -0.143078866}, {-0.357865765, 0.172339069, 0.20775094300000002}}
{7, {-0.6710628700000001, -0.8414862230000001, 0.143078866}, {-0.32142454200000004, -0.243394199, 0.10991626400000001}}
{13, {-0.969713943, -0.466989623, -0.143078866}, {-0.32142454200000004, -0.243394199, 0.10991626400000001}}
{5, {-0.6710628700000001, 0.8414862230000001, 0.143078866}, {-0.275983081, 0.275016707, 0.28620826400000005}}
{2, {0., 0., 0.429236598}, {-0.09199436000000001, 0.044302149000000006, 1.}}
{4, {0.23949956500000003, 1.049316157, 0.143078866}, {-0.010111676, 0.403053627, 0.10991626400000001}}
{11, {-0.23949956500000003, 1.049316157, -0.143078866}, {-0.010111676, 0.403053627, 0.10991626400000001}}
{8, {0.23949956500000003, -1.049316157, 0.143078866}, {0.010111676, -0.403053627, -0.10991626400000001}}
{14, {-0.23949956500000003, -1.049316157, -0.143078866}, {0.010111676, -0.403053627, -0.10991626400000001}}
{1, {0., 0., -0.429236598}, {0.09199436000000001, -0.044302149000000006, -1.}}
{15, {0.6710628700000001, -0.8414862230000001, -0.143078866}, {0.275983081, -0.275016707, -0.28620826400000005}}
{3, {0.969713943, 0.466989623, 0.143078866}, {0.32142454200000004, 0.243394199, -0.10991626400000001}}
{10, {0.6710628700000001, 0.8414862230000001, -0.143078866}, {0.32142454200000004, 0.243394199, -0.10991626400000001}}
{9, {0.969713943, -0.466989623, 0.143078866}, {0.357865765, -0.172339069, -0.20775094300000002}}
{16, {1.076301277, 0., -0.143078866}, {0.38708935600000005, -0.044302149000000006, -0.28620826400000005}}
Groups
{1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1}
Group count
12

OctagonalAntiprism. Side pieces: 2 doubles, 12 singles.
1 double, 5 singles, 1 double, 7 singles
{3, 11}, {7, 14}

piece centers and motions
{7, {-1.1458497600000002, -0.47462651100000003, 0.143382595}, {-0.34332009700000005, -0.065486101, 0.039566130000000005}}
{14, {-1.240258842, 0., -0.143382595}, {-0.34332009700000005, -0.065486101, 0.039566130000000005}}
{6, {-1.1458497600000002, 0.47462651100000003, 0.143382595}, {-0.29518229500000004, 0.17651897200000002, 0.18578348}}
{13, {-0.876995438, 0.876995438, -0.143382595}, {-0.263017644, 0.224656774, 0.151476045}}
{15, {-0.876995438, -0.876995438, -0.143382595}, {-0.194940512, -0.287551844, -0.11869839}}
{8, {-0.47462651100000003, -1.1458497600000002, 0.143382595}, {-0.17803676400000001, -0.298846567, -0.10665122}}
{5, {-0.47462651100000003, 1.1458497600000002, 0.143382595}, {-0.129898961, 0.313603835, 0.24634869}}
{2, {0., 0., 0.430147785}, {-0.027125231000000003, 0.065486101, 1.}}
{16, {0., -1.240258842, -0.143382595}, {0.027125231000000003, -0.339655825, -0.230608305}}
{1, {0., 0., -0.430147785}, {0.027125231000000003, -0.065486101, -1.}}
{12, {0., 1.240258842, -0.143382595}, {0.027125231000000003, 0.34483788800000004, 0.151476045}}
{4, {0.47462651100000003, 1.1458497600000002, 0.143382595}, {0.08390764, 0.33354316500000003, 0.18578348}}
{9, {0.47462651100000003, -1.1458497600000002, 0.143382595}, {0.132045442, -0.318785898, -0.16721643}}
{17, {0.876995438, -0.876995438, -0.143382595}, {0.220992502, -0.259353372, -0.230608305}}
{3, {1.1458497600000002, 0.47462651100000003, 0.143382595}, {0.289069635, 0.196458303, 0.039566130000000005}}
{11, {0.876995438, 0.876995438, -0.143382595}, {0.289069635, 0.196458303, 0.039566130000000005}}
{10, {1.1458497600000002, -0.47462651100000003, 0.143382595}, {0.337207437, -0.08542543100000001, -0.10665122}}
{18, {1.240258842, 0., -0.143382595}, {0.341173616, -0.065486101, -0.11869839}}
Groups
{2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}
Group count
16

NonagonalAntiprism. Side pieces: 4 doubles, 10 singles.
2 doubles, 5 singles, 2 doubles, 5 singles
{3, 12}, {4, 13}, {8, 16}, {9, 17}

piece centers and motions
{7, {-1.4031265400000001, 0., 0.14358768100000002}, {-0.30244048900000003, 0.037246564, 0.177051744}}
{15, {-1.318507656, 0.47989754, -0.14358768100000002}, {-0.293867129, 0.085868499, 0.147955904}}
{8, {-1.074857289, -0.9019123550000001, 0.14358768100000002}, {-0.269181126, -0.151376656, 0.06417777200000001}}
{16, {-1.318507656, -0.47989754, -0.14358768100000002}, {-0.269181126, -0.151376656, 0.06417777200000001}}
{6, {-1.074857289, 0.9019123550000001, 0.14358768100000002}, {-0.23135997400000002, 0.19413406900000002, 0.22162915700000002}}
{14, {-0.7015632700000001, 1.215143228, -0.14358768100000002}, {-0.135593455, 0.274491719, 0.147955904}}
{9, {-0.243650367, -1.3818098950000002, 0.14358768100000002}, {-0.102334093, -0.29137794, -0.06417777200000001}}
{17, {-0.7015632700000001, -1.215143228, -0.14358768100000002}, {-0.102334093, -0.29137794, -0.06417777200000001}}
{5, {-0.243650367, 1.3818098950000002, 0.14358768100000002}, {-0.089198944, 0.29137794, 0.177051744}}
{2, {0., 0., 0.43076304400000004}, {-0.044388726, 0.037246564, 1.}}
{1, {0., 0., -0.43076304400000004}, {0.044388726, -0.037246564, -1.}}
{18, {0.243650367, -1.3818098950000002, -0.14358768100000002}, {0.089198944, -0.29137794, -0.177051744}}
{4, {0.7015632700000001, 1.215143228, 0.14358768100000002}, {0.102334093, 0.29137794, 0.06417777200000001}}
{13, {0.243650367, 1.3818098950000002, -0.14358768100000002}, {0.102334093, 0.29137794, 0.06417777200000001}}
{10, {0.7015632700000001, -1.215143228, 0.14358768100000002}, {0.135593455, -0.274491719, -0.147955904}}
{19, {1.074857289, -0.9019123550000001, -0.14358768100000002}, {0.23135997400000002, -0.19413406900000002, -0.22162915700000002}}
{3, {1.318507656, 0.47989754, 0.14358768100000002}, {0.269181126, 0.151376656, -0.06417777200000001}}
{12, {1.074857289, 0.9019123550000001, -0.14358768100000002}, {0.269181126, 0.151376656, -0.06417777200000001}}
{11, {1.318507656, -0.47989754, 0.14358768100000002}, {0.293867129, -0.085868499, -0.147955904}}
{20, {1.4031265400000001, 0., -0.14358768100000002}, {0.30244048900000003, -0.037246564, -0.177051744}}
Groups
{1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1}
Group count
16

DecagonalAntiprism. Side pieces: 2 doubles, 16 singles.
1 double, 7 singles, 1 double, 9 singles
{3, 13}, {8, 17}

piece centers and motions
{8, {-1.488630917, -0.483685505, 0.143732834}, {-0.273524569, -0.042208714, 0.025085631}}
{17, {-1.565239175, 0., -0.143732834}, {-0.273524569, -0.042208714, 0.025085631}}
{7, {-1.488630917, 0.483685505, 0.143732834}, {-0.24809254400000003, 0.118362775, 0.120517042}}
{16, {-1.2663050930000002, 0.920024504, -0.143732834}, {-0.23706958600000003, 0.139996548, 0.106264438}}
{18, {-1.2663050930000002, -0.920024504, -0.143732834}, {-0.195919704, -0.19451683700000003, -0.075256893}}
{9, {-0.920024504, -1.2663050930000002, 0.143732834}, {-0.190003787, -0.20043275300000002, -0.07034578000000001}}
{6, {-0.920024504, 1.2663050930000002, 0.143732834}, {-0.148853906, 0.22821222800000002, 0.17949689800000002}}
{15, {-0.483685505, 1.488630917, -0.143732834}, {-0.08476146300000001, 0.260868959, 0.13727198300000001}}
{19, {-0.483685505, -1.488630917, -0.143732834}, {-0.059329437000000006, -0.26701466, -0.1564357}}
{10, {0., -1.565239175, 0.143732834}, {-0.013714442, -0.27423936600000004, -0.12932563600000002}}
{2, {0., 0., 0.43119850200000004}, {-0.013714442, 0.042208714, 1.}}
{5, {0., 1.565239175, 0.143732834}, {-0.013714442, 0.272121701, 0.17949689800000002}}
{1, {0., 0., -0.43119850200000004}, {0.013714442, -0.042208714, -1.}}
{20, {0.483685505, -1.488630917, -0.143732834}, {0.08407339200000001, -0.258751294, -0.187443245}}
{14, {0.483685505, 1.488630917, -0.143732834}, {0.10950541700000001, 0.252605593, 0.106264438}}
{4, {0.920024504, 1.2663050930000002, 0.143732834}, {0.13113919000000002, 0.24158263500000002, 0.120517042}}
{11, {0.920024504, -1.2663050930000002, 0.143732834}, {0.17228907200000002, -0.21380316000000002, -0.12932563600000002}}
{21, {1.2663050930000002, -0.920024504, -0.143732834}, {0.204945803, -0.18114643000000002, -0.1564357}}
{3, {1.488630917, 0.483685505, 0.143732834}, {0.246095684, 0.126626141, 0.025085631}}
{13, {1.2663050930000002, 0.920024504, -0.143732834}, {0.246095684, 0.126626141, 0.025085631}}
{12, {1.488630917, -0.483685505, 0.143732834}, {0.27152771000000003, -0.05047208, -0.07034578000000001}}
{22, {1.565239175, 0., -0.143732834}, {0.27283649800000004, -0.042208714, -0.075256893}}
Groups
{2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1}
Group count
20


It's still not clear why such motions + all symmetrical ones, span a 3D region.
Why do we have the pattern of singles/doubles

 *)