#!/usr/bin/env wolframscript

(* 

 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

glueFirstPair = False;
glueFirstPair = True;

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

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)

(* shape = "Octahedron"; *)
(* shape = "Icosahedron"; *)

(* shape = "RhombicDodecahedron"; *)
(* shape = "RhombicTriacontahedron"; *)

(* shape = "Cuboctahedron"; *)
(* shape = "Icosidodecahedron"; *)

shape = "DeltoidalIcositetrahedron";


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

If[glueFirstPair,
    edgeIndex = edgeToFaces[[1]];
    gluedPair = edgeIndex[[2]];
    Print["Gluing a pair of pieces: "];
    Print[gluedPair];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];

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
            AppendTo[pieceAndMotion, { Round[ faceCenters[[pid]], epsilon ], values[[pid]]}];
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

For rhombic triacontahedron and rhombic dodecahedron, we see interesting non-uniform expansion solutions like the physical puzzle. 

RhombicDodecahedron: dimension = 2
piece centers and motions
{{-0.816496581, 0., 0.}, {-1., 0., -0.707106781}}
{{-0.40824829, -0.40824829, -0.5773502690000001}, {-1., 0., -0.707106781}}
{{-0.40824829, 0.40824829, -0.5773502690000001}, {-1., 0., -0.707106781}}
{{-0.40824829, -0.40824829, 0.5773502690000001}, {0., -1., 0.707106781}}
{{0., -0.816496581, 0.}, {0., -1., 0.707106781}}
{{0.40824829, -0.40824829, 0.5773502690000001}, {0., -1., 0.707106781}}
{{-0.40824829, 0.40824829, 0.5773502690000001}, {0., 1., 0.707106781}}
{{0., 0.816496581, 0.}, {0., 1., 0.707106781}}
{{0.40824829, 0.40824829, 0.5773502690000001}, {0., 1., 0.707106781}}
{{0.40824829, -0.40824829, -0.5773502690000001}, {1., 0., -0.707106781}}
{{0.40824829, 0.40824829, -0.5773502690000001}, {1., 0., -0.707106781}}
{{0.816496581, 0., 0.}, {1., 0., -0.707106781}}
Groups
{3, 3, 3, 3}
Group count
4

RhombicTriacontahedron: dimension = 4
piece centers and motions
{{-1.309016994, -0.425325404, 0.}, {-0.7847047290000001, -0.49259709100000004, 0.139676038}}
{{-1.309016994, 0.425325404, 0.}, {-0.7847047290000001, 0.695558255, 0.139676038}}
{{-0.947213595, 0.68819096, 0.723606798}, {-0.7847047290000001, 0.695558255, 0.139676038}}
{{-0.809016994, 1.113516364, 0.}, {-0.7847047290000001, 0.695558255, 0.139676038}}
{{-0.947213595, -0.68819096, 0.723606798}, {-0.591677192, -0.632839806, 0.525731112}}
{{-1.170820393, 0., -0.723606798}, {-0.591677192, 0.101480582, -0.871029268}}
{{-0.5854101970000001, -0.425325404, -1.170820393}, {-0.591677192, 0.101480582, -0.871029268}}
{{-0.5854101970000001, 0.425325404, -1.170820393}, {-0.591677192, 0.101480582, -0.871029268}}
{{-0.809016994, -1.113516364, 0.}, {-0.517947233, -0.859757285, 0.139676038}}
{{-0.723606798, 0., 1.170820393}, {-0.47237961300000003, -0.265679612, 0.764326269}}
{{-0.361803399, 1.113516364, -0.723606798}, {-0.39864965500000005, 0.695558255, -0.484974194}}
{{-0.361803399, -1.113516364, -0.723606798}, {-0.279352076, -0.859757285, -0.24637903600000002}}
{{-0.22360679800000002, -0.68819096, 1.170820393}, {-0.20562211800000002, -0.632839806, 0.764326269}}
{{0., -1.37638192, 0.}, {-0.086324539, -1., 0.139676038}}
{{0., 1.37638192, 0.}, {-0.086324539, 0.922475734, 0.139676038}}
{{0.361803399, -1.113516364, 0.723606798}, {0.10670299800000001, -0.859757285, 0.525731112}}
{{0.22360679800000002, 0.68819096, -1.170820393}, {0.10670299800000001, 0.328398061, -0.871029268}}
{{0.5854101970000001, -0.425325404, 1.170820393}, {0.226000577, -0.49259709100000004, 0.764326269}}
{{-0.22360679800000002, 0.68819096, 1.170820393}, {0.226000577, 0.695558255, 0.764326269}}
{{0.361803399, 1.113516364, 0.723606798}, {0.226000577, 0.695558255, 0.764326269}}
{{0.5854101970000001, 0.425325404, 1.170820393}, {0.226000577, 0.695558255, 0.764326269}}
{{0.809016994, -1.113516364, 0.}, {0.345298156, -0.859757285, 0.139676038}}
{{0.22360679800000002, -0.68819096, -1.170820393}, {0.538325693, -0.265679612, -0.871029268}}
{{0.723606798, 0., -1.170820393}, {0.538325693, -0.265679612, -0.871029268}}
{{0.947213595, -0.68819096, -0.723606798}, {0.538325693, -0.265679612, -0.871029268}}
{{0.809016994, 1.113516364, 0.}, {0.612055651, 0.695558255, 0.139676038}}
{{0.947213595, 0.68819096, -0.723606798}, {0.7313532300000001, 0.328398061, -0.484974194}}
{{1.170820393, 0., 0.723606798}, {1.043678345, 0.101480582, 0.139676038}}
{{1.309016994, -0.425325404, 0.}, {1.043678345, 0.101480582, 0.139676038}}
{{1.309016994, 0.425325404, 0.}, {1.043678345, 0.101480582, 0.139676038}}
Groups
{1, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 3, 1, 1, 3}
Group count
20

Open question: are there more motions beyond the ones described by gluing pieces? Are the gluing pieces motion the vertices?

For Cube, Tetrahedron, Dodecahedron, after we glue any pair, no solution is found.
For these 3 shapes, dimension = 1. The only allowed motions are uniform expansions.

For Octahedron, the 8 pieces can split in two halves along equator. It's a boring solution. 
Dimension = 3: splitting in the x, y, z directions.
piece centers and motions
{{-0.23570226000000002, -0.23570226000000002, -0.23570226000000002}, {-1., 0., 0.}}
{{-0.23570226000000002, -0.23570226000000002, 0.23570226000000002}, {-1., 0., 0.}}
{{-0.23570226000000002, 0.23570226000000002, -0.23570226000000002}, {-1., 0., 0.}}
{{-0.23570226000000002, 0.23570226000000002, 0.23570226000000002}, {-1., 0., 0.}}
{{0.23570226000000002, -0.23570226000000002, -0.23570226000000002}, {1., 0., 0.}}
{{0.23570226000000002, -0.23570226000000002, 0.23570226000000002}, {1., 0., 0.}}
{{0.23570226000000002, 0.23570226000000002, -0.23570226000000002}, {1., 0., 0.}}
{{0.23570226000000002, 0.23570226000000002, 0.23570226000000002}, {1., 0., 0.}}
Groups
{4, 4}
Group count
2

For Icosahedron, more cases are allowed. dimension = 9
piece centers and motions
{{-0.742344243, 0., 0.141775135}, {-1., 0.7265425280000001, 1.145898034}}
{{-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-1., 0.7265425280000001, 1.145898034}}
{{-0.458793973, 0., 0.6005691080000001}, {-1., 0.7265425280000001, 1.145898034}}
{{-0.22939698700000002, 0.70601133, 0.141775135}, {-1., 0.7265425280000001, 1.145898034}}
{{-0.141775135, 0.43633899800000003, 0.6005691080000001}, {-1., 0.7265425280000001, 1.145898034}}
{{-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-0.47213595500000005, -0.8980559530000001, 0.090169944}}
{{-0.22939698700000002, -0.70601133, 0.141775135}, {-0.47213595500000005, -0.8980559530000001, 0.090169944}}
{{-0.141775135, -0.43633899800000003, 0.6005691080000001}, {-0.27050983100000003, -0.277514551, 1.145898034}}
{{0.371172121, -0.26967233100000004, 0.6005691080000001}, {0.180339887, -0.13102459800000002, 1.145898034}}
{{0.371172121, 0.26967233100000004, 0.6005691080000001}, {0.180339887, 0.34302685, 1.145898034}}
{{-0.371172121, -0.26967233100000004, -0.6005691080000001}, {0.38196601100000005, -0.277514551, -1.6180339890000002}}
{{-0.371172121, 0.26967233100000004, -0.6005691080000001}, {0.38196601100000005, -0.277514551, -1.6180339890000002}}
{{0.141775135, -0.43633899800000003, -0.6005691080000001}, {0.38196601100000005, -0.277514551, -1.6180339890000002}}
{{0.141775135, 0.43633899800000003, -0.6005691080000001}, {0.38196601100000005, -0.277514551, -1.6180339890000002}}
{{0.458793973, 0., -0.6005691080000001}, {0.38196601100000005, -0.277514551, -1.6180339890000002}}
{{0.22939698700000002, -0.70601133, -0.141775135}, {0.5835921350000001, -0.8980559530000001, -0.562305899}}
{{0.22939698700000002, 0.70601133, -0.141775135}, {0.7082039320000001, 0.7265425280000001, 0.090169944}}
{{0.6005691080000001, 0.43633899800000003, 0.141775135}, {0.7082039320000001, 0.7265425280000001, 0.090169944}}
{{0.6005691080000001, -0.43633899800000003, 0.141775135}, {0.9098300560000001, -0.6610302290000001, -0.31308230400000003}}
{{0.742344243, 0., -0.141775135}, {1.034441854, -0.277514551, -0.562305899}}
Groups
{5, 2, 1, 1, 1, 5, 1, 2, 1, 1}
Group count
10

Cuboctahedron: dimension = 4
piece centers and motions
{{-0.666666667, 0., -0.471404521}, {-1.414213562, 0., 0.}}
{{-0.666666667, 0., 0.471404521}, {-1.414213562, 0., 0.}}
{{-0.5, -0.5, 0.}, {-1.414213562, 0., 0.}}
{{-0.5, 0.5, 0.}, {-1.414213562, 0., 0.}}
{{0., -0.666666667, -0.471404521}, {0., 0., -1.}}
{{0., 0., -0.707106781}, {0., 0., -1.}}
{{0., 0.666666667, -0.471404521}, {0., 0., -1.}}
{{0., -0.666666667, 0.471404521}, {0., 0., 1.}}
{{0., 0., 0.707106781}, {0., 0., 1.}}
{{0., 0.666666667, 0.471404521}, {0., 0., 1.}}
{{0.5, -0.5, 0.}, {1.414213562, 0., 0.}}
{{0.5, 0.5, 0.}, {1.414213562, 0., 0.}}
{{0.666666667, 0., -0.471404521}, {1.414213562, 0., 0.}}
{{0.666666667, 0., 0.471404521}, {1.414213562, 0., 0.}}
Groups
{4, 3, 3, 4}
Group count
4


Icosidodecahedron: dimension = 6
piece centers and motions
{{-1.484688486, 0., -0.283550269}, {-1.236067977, 0., -0.236067977}}
{{-1.231073415, 0., 0.6155367070000001}, {-1.236067977, 0., -0.236067977}}
{{-1.2011382160000001, 0.8726779960000001, 0.283550269}, {-1.236067977, 0., -0.236067977}}
{{-0.995959314, 0.723606798, -0.6155367070000001}, {-1.236067977, 0., -0.236067977}}
{{-0.742344243, -0.5393446630000001, 1.2011382160000001}, {-1.236067977, 0., -0.236067977}}
{{-0.283550269, 0.8726779960000001, -1.2011382160000001}, {-1.236067977, 0., -0.236067977}}
{{-1.2011382160000001, -0.8726779960000001, 0.283550269}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{-0.995959314, -0.723606798, -0.6155367070000001}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{-0.917587947, 0., -1.2011382160000001}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{-0.458793973, -1.412022659, -0.283550269}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{-0.380422607, -1.170820393, 0.6155367070000001}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{-0.283550269, -0.8726779960000001, -1.2011382160000001}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{0., 0., -1.37638192}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{0.380422607, -1.170820393, -0.6155367070000001}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{0.458793973, -1.412022659, 0.283550269}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{0.742344243, -0.5393446630000001, -1.2011382160000001}, {-0.7639320230000001, -1.4530850560000002, -1.}}
{{-0.742344243, 0.5393446630000001, 1.2011382160000001}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{-0.458793973, 1.412022659, -0.283550269}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{-0.380422607, 1.170820393, 0.6155367070000001}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0., 0., 1.37638192}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0.283550269, 0.8726779960000001, 1.2011382160000001}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0.380422607, 1.170820393, -0.6155367070000001}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0.458793973, 1.412022659, 0.283550269}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0.917587947, 0., 1.2011382160000001}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0.995959314, 0.723606798, 0.6155367070000001}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{1.2011382160000001, 0.8726779960000001, -0.283550269}, {0.7639320230000001, 1.4530850560000002, 1.}}
{{0.283550269, -0.8726779960000001, 1.2011382160000001}, {1.236067977, 0., 0.236067977}}
{{0.742344243, 0.5393446630000001, -1.2011382160000001}, {1.236067977, 0., 0.236067977}}
{{0.995959314, -0.723606798, 0.6155367070000001}, {1.236067977, 0., 0.236067977}}
{{1.2011382160000001, -0.8726779960000001, -0.283550269}, {1.236067977, 0., 0.236067977}}
{{1.231073415, 0., -0.6155367070000001}, {1.236067977, 0., 0.236067977}}
{{1.484688486, 0., 0.283550269}, {1.236067977, 0., 0.236067977}}
Groups
{6, 10, 10, 6}
Group count
4

DeltoidalIcositetrahedron (dual of Rhombicuboctahedron): dimension = 4
piece centers and motions
{{-1.250520274, -0.5295892790000001, -0.5295892790000001}, {-1., 0., 0.}}
{{-1.250520274, -0.5295892790000001, 0.5295892790000001}, {-1., 0., 0.}}
{{-1.250520274, 0.5295892790000001, -0.5295892790000001}, {-1., 0., 0.}}
{{-1.250520274, 0.5295892790000001, 0.5295892790000001}, {-1., 0., 0.}}
{{-0.5295892790000001, -1.250520274, -0.5295892790000001}, {-1., 0., 0.}}
{{-0.5295892790000001, -1.250520274, 0.5295892790000001}, {-1., 0., 0.}}
{{-0.5295892790000001, -0.5295892790000001, -1.250520274}, {-1., 0., 0.}}
{{-0.5295892790000001, -0.5295892790000001, 1.250520274}, {-1., 0., 0.}}
{{-0.5295892790000001, 0.5295892790000001, -1.250520274}, {-1., 0., 0.}}
{{-0.5295892790000001, 0.5295892790000001, 1.250520274}, {-1., 0., 0.}}
{{-0.5295892790000001, 1.250520274, -0.5295892790000001}, {-1., 0., 0.}}
{{-0.5295892790000001, 1.250520274, 0.5295892790000001}, {-1., 0., 0.}}
{{0.5295892790000001, -1.250520274, -0.5295892790000001}, {1., 0., 0.}}
{{0.5295892790000001, -1.250520274, 0.5295892790000001}, {1., 0., 0.}}
{{0.5295892790000001, -0.5295892790000001, -1.250520274}, {1., 0., 0.}}
{{0.5295892790000001, -0.5295892790000001, 1.250520274}, {1., 0., 0.}}
{{0.5295892790000001, 0.5295892790000001, -1.250520274}, {1., 0., 0.}}
{{0.5295892790000001, 0.5295892790000001, 1.250520274}, {1., 0., 0.}}
{{0.5295892790000001, 1.250520274, -0.5295892790000001}, {1., 0., 0.}}
{{0.5295892790000001, 1.250520274, 0.5295892790000001}, {1., 0., 0.}}
{{1.250520274, -0.5295892790000001, -0.5295892790000001}, {1., 0., 0.}}
{{1.250520274, -0.5295892790000001, 0.5295892790000001}, {1., 0., 0.}}
{{1.250520274, 0.5295892790000001, -0.5295892790000001}, {1., 0., 0.}}
{{1.250520274, 0.5295892790000001, 0.5295892790000001}, {1., 0., 0.}}
Groups
{12, 12}
Group count
2

some catalan solids:
DeltoidalIcositetrahedron (dual of Rhombicuboctahedron): dimension = 4
DeltoidalHexecontahedron (dual of rhombicosidodecahedron): dimension = 10
TetrakisHexahedron: dimension = 11
PentagonalIcositetrahedron (dual of snub cube): dim = 1


some archimedian solids:

TruncatedTetrahedron: dim = 1
TruncatedCube: dim = 1
TruncatedOctahedron: dim = 1
TruncatedIcosahedron: dim = 1
TruncatedDodecahedron: dim = 1

Cuboctahedron: dim = 4
Icosidodecahedron: dim = 6

SnubCube: dim = 15

 *)