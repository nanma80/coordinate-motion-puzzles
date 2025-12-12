#!/usr/bin/env wolframscript

(* 
Why does SmallRhombicuboctahedron have dimension 6 rather than 3?
 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

(* PairsToGlue is not used in this script. We glue based on gluedGroups *)
PairsToGlue = 0;

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

shape = "SmallRhombicuboctahedron";

gluedGroups = {
    (* bottom cap, dim = 2 *)
    (* {1, 15, 16, 17, 18, 23, 24, 25, 26} *)

    (* top and bottom caps, dim = 1 *)
    (* {1, 15, 16, 17, 18, 23, 24, 25, 26}, {2, 11, 12, 13, 14, 19, 20, 21, 22} *)

    (* bottom cap + 4 pairs *)
    {1, 15, 16, 17, 18, 23, 24, 25, 26}, {7, 19}, {8, 20}, {9, 21}, {10, 22}
};

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

(*
    each item in edgeToFaces describes an edge, like:
    {{{1, 2}, {13, 14}}
    It means the edge is between vertices with index 1 to vertext index 2. The edge is shared between faces 13 and 14
    It should work for convex polyhedra. No guarantee for 
*)

faceCenters = Mean /@ (vertices[[#]] & /@ faces);
faceCenters = N[faceCenters];

faceCount = Length[faces];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];

minFaceCenterDistance = 1000.0;

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

For[pid = 1, pid <= Length[gluedPairs], pid++,
    gluedPair = gluedPairs[[pid]];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] - pieceCoordinates[[gluedPair[[2]]]][[dim]] == 0 ];
    ];
];

(* For[gid = 1, gid <= PairsToGlue, gid++,
    edgeIndex = edgeToFaces[[gid]];
    gluedPair = edgeIndex[[2]];
    Print["Gluing a pair of pieces: "];
    Print[gluedPair];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];
 *)
Print["condition count:"];
Print[Length[conditions]];

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
        Print["Group count"]
        Print[ tallyCounts//Length ];
    ];

];

Print["The shape is: ", shape];




(* 

piece centers and motions
{4, {-1.207106781, 0., 0.}, {-0.41421356200000004, 0., 0.}}
{9, {-0.853553391, -0.853553391, 0.}, {-0.29289321900000004, -0.29289321900000004, 0.}}
{10, {-0.853553391, 0.853553391, 0.}, {-0.29289321900000004, 0.29289321900000004, 0.}}
{5, {0., -1.207106781, 0.}, {0., -0.41421356200000004, 0.}}

{1, {0., 0., -1.207106781}, {0., 0., -1.}}
{15, {0., -0.853553391, -0.853553391}, {0., 0., -1.}}
{16, {0., 0.853553391, -0.853553391}, {0., 0., -1.}}
{17, {0.853553391, 0., -0.853553391}, {0., 0., -1.}}
{18, {-0.853553391, 0., -0.853553391}, {0., 0., -1.}}
{23, {0.73570226, -0.73570226, -0.73570226}, {0., 0., -1.}}
{24, {0.73570226, 0.73570226, -0.73570226}, {0., 0., -1.}}
{25, {-0.73570226, -0.73570226, -0.73570226}, {0., 0., -1.}}
{26, {-0.73570226, 0.73570226, -0.73570226}, {0., 0., -1.}}

{2, {0., 0., 1.207106781}, {0., 0., 1.}}
{11, {0., -0.853553391, 0.853553391}, {0., 0., 1.}}
{12, {0., 0.853553391, 0.853553391}, {0., 0., 1.}}
{13, {0.853553391, 0., 0.853553391}, {0., 0., 1.}}
{14, {-0.853553391, 0., 0.853553391}, {0., 0., 1.}}
{19, {0.73570226, -0.73570226, 0.73570226}, {0., 0., 1.}}
{20, {0.73570226, 0.73570226, 0.73570226}, {0., 0., 1.}}
{21, {-0.73570226, -0.73570226, 0.73570226}, {0., 0., 1.}}
{22, {-0.73570226, 0.73570226, 0.73570226}, {0., 0., 1.}}

{6, {0., 1.207106781, 0.}, {0., 0.41421356200000004, 0.}}
{7, {0.853553391, -0.853553391, 0.}, {0.29289321900000004, -0.29289321900000004, 0.}}
{8, {0.853553391, 0.853553391, 0.}, {0.29289321900000004, 0.29289321900000004, 0.}}
{3, {1.207106781, 0., 0.}, {0.41421356200000004, 0., 0.}}
Groups
{1, 1, 1, 1, 9, 9, 1, 1, 1, 1}
Group count
10

if we glue
    {1, 15, 16, 17, 18, 23, 24, 25, 26}
dimension = 2

Alternative motion, with only one cap
piece centers and motions
bottom cap (9), the whole bottom layer
{1, {0., 0., -1.207106781}, {0., 0., -1.}}
{15, {0., -0.853553391, -0.853553391}, {0., 0., -1.}}
{16, {0., 0.853553391, -0.853553391}, {0., 0., -1.}}
{17, {0.853553391, 0., -0.853553391}, {0., 0., -1.}}
{18, {-0.853553391, 0., -0.853553391}, {0., 0., -1.}}
{23, {0.73570226, -0.73570226, -0.73570226}, {0., 0., -1.}}
{24, {0.73570226, 0.73570226, -0.73570226}, {0., 0., -1.}}
{25, {-0.73570226, -0.73570226, -0.73570226}, {0., 0., -1.}}
{26, {-0.73570226, 0.73570226, -0.73570226}, {0., 0., -1.}}

4 pairs: corner on top + edge below the corner
{9, {-0.853553391, -0.853553391, 0.}, {-0.43303287700000004, -0.43303287700000004, 0.47846672}}
{21, {-0.73570226, -0.73570226, 0.73570226}, {-0.43303287700000004, -0.43303287700000004, 0.47846672}}

{10, {-0.853553391, 0.853553391, 0.}, {-0.43303287700000004, 0.43303287700000004, 0.47846672}}
{22, {-0.73570226, 0.73570226, 0.73570226}, {-0.43303287700000004, 0.43303287700000004, 0.47846672}}

{7, {0.853553391, -0.853553391, 0.}, {0.43303287700000004, -0.43303287700000004, 0.47846672}}
{19, {0.73570226, -0.73570226, 0.73570226}, {0.43303287700000004, -0.43303287700000004, 0.47846672}}

{8, {0.853553391, 0.853553391, 0.}, {0.43303287700000004, 0.43303287700000004, 0.47846672}}
{20, {0.73570226, 0.73570226, 0.73570226}, {0.43303287700000004, 0.43303287700000004, 0.47846672}}

singles:
{4, {-1.207106781, 0., 0.}, {-0.612400967, 0., 0.47846672}}
{14, {-0.853553391, 0., 0.853553391}, {-0.5598652700000001, 0., 0.605299113}}
{5, {0., -1.207106781, 0.}, {0., -0.612400967, 0.47846672}}
{11, {0., -0.853553391, 0.853553391}, {0., -0.5598652700000001, 0.605299113}}
{2, {0., 0., 1.207106781}, {0., 0., 0.837202901}}
{12, {0., 0.853553391, 0.853553391}, {0., 0.5598652700000001, 0.605299113}}
{6, {0., 1.207106781, 0.}, {0., 0.612400967, 0.47846672}}
{13, {0.853553391, 0., 0.853553391}, {0.5598652700000001, 0., 0.605299113}}
{3, {1.207106781, 0., 0.}, {0.612400967, 0., 0.47846672}}

Groups
{1, 1, 2, 2, 1, 1, 9, 1, 1, 1, 2, 2, 1, 1}
Group count
14

 *)