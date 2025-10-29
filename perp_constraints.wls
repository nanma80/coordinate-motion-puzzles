#!/usr/bin/env wolframscript

(* 

 *)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

glueFirstPair = False;
glueFirstPair = True;

computeDimensions = False;
computeDimensions = True;

If[forceUniformExpansion && glueFirstPair,
    Print["Cannot both force uniform expansion and glue first pair at the same time. Will not glue first pair."];
    glueFirstPair = False;
];

If[computeDimensions && (forceUniformExpansion || glueFirstPair ),
    Print["Computing dimensions. We are setting uniform expansion and glueing to false"];
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


Print["shape is: ", shape];

vertices = PolyhedronData[shape, "VertexCoordinates"];
faces = PolyhedronData[shape, "FaceIndices"];

faceCenters = Mean /@ (vertices[[#]] & /@ faces);
faceCenters = N[faceCenters];
faceCount = Length[faceCenters];
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

If[glueFirstPair,
    gluedPair = neighborFacePairs[[1]];
    Print["Glueing a pair of pieces: "];
    Print[gluedPair];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
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
(*     Print[];
    Print["rows"];
    Print[rows];
    Print[];
    Print["rows row count: ", Length[rows]];
    Print["rows col count: ", Length[rows[[1]] ]]; *)

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
        Print["piece and motion"];
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
        Print[ tallyCounts//Length ];
    ];

];




(* 
results:

For rhombic triacontahedron and rhombic dodecahedron, we see interesting non-uniform expansion solutions like the physical puzzle. 
RhombicDodecahedron: dimension = 2
RhombicTriacontahedron: dimension = 4
Question: are there more motions beyond the ones described by gluing pieces? Are the gluing pieces motion the vertices?

For Cube, Tetrahedron, Dodecahedron, after we glue any pair, no solution is found. It suggests only uniform expansion is allowed.
For these 3 shapes, dimension = 1

For Octahedron, the 8 pieces can split in two halves along equator. It's a boring solution. 
Dimension = 3: splitting in the x, y, z directions.

piece and motion
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

For Icosahedron, more cases are allowed. dimension = 9
piece and motion
{{-0.742344243, 0., 0.141775135}, {-1., 0.5156624860000001, 0.42479043600000005}}
{{-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-1., 0.5156624860000001, 0.42479043600000005}}
{{-0.458793973, 0., 0.6005691080000001}, {-1., 0.5156624860000001, 0.42479043600000005}}
{{-0.22939698700000002, 0.70601133, 0.141775135}, {-1., 0.5156624860000001, 0.42479043600000005}}
{{-0.141775135, 0.43633899800000003, 0.6005691080000001}, {-1., 0.5156624860000001, 0.42479043600000005}}
{{-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-0.668173698, -0.5055938600000001, -0.238862168}}
{{-0.371172121, -0.26967233100000004, -0.6005691080000001}, {-0.54422152, -0.415537331, -0.486766525}}
{{-0.371172121, 0.26967233100000004, -0.6005691080000001}, {-0.54422152, 0.184520037, -0.486766525}}
{{-0.22939698700000002, -0.70601133, 0.141775135}, {-0.355499061, -0.732765281, 0.}}
{{-0.141775135, -0.43633899800000003, 0.6005691080000001}, {-0.274371307, -0.48307972800000004, 0.42479043600000005}}
{{0.141775135, -0.43633899800000003, -0.6005691080000001}, {0.27716549700000004, -0.682422151, -0.486766525}}
{{0.22939698700000002, -0.70601133, -0.141775135}, {0.293522972, -0.732765281, -0.40111767600000003}}
{{0.141775135, 0.43633899800000003, -0.6005691080000001}, {0.474930145, 0.5156624860000001, -0.486766525}}
{{0.22939698700000002, 0.70601133, -0.141775135}, {0.474930145, 0.5156624860000001, -0.486766525}}
{{0.458793973, 0., -0.6005691080000001}, {0.811277542, 0.052720010000000005, -0.486766525}}
{{0.371172121, -0.26967233100000004, 0.6005691080000001}, {0.8997205810000001, -0.10159414800000001, 0.42479043600000005}}
{{0.371172121, 0.26967233100000004, 0.6005691080000001}, {0.8997205810000001, -0.10159414800000001, 0.42479043600000005}}
{{0.6005691080000001, -0.43633899800000003, 0.141775135}, {1.030988045, -0.19696554300000002, 0.16225550800000002}}
{{0.6005691080000001, 0.43633899800000003, 0.141775135}, {1.1121157990000001, 0.052720010000000005, 0.}}
{{0.742344243, 0., -0.141775135}, {1.1121157990000001, 0.052720010000000005, 0.}}
Groups
{5, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2}
13


 *)