#!/usr/bin/env wolframscript

(* 
Two types of constraints:
cross product (equality); nooverlap condition (inequalities)
if both constraints are rotated, if it's not 0 or Pi, always no solution. Nothing to see 
if only the nooverlap condition is rotated, there are solutions, dimension is not changed. 

interesting case, when rotateAngle = Pi/2;
{8, 1, 1, 1, 1, 8}
Group count
6

*)

forceUniformExpansion = False;
(* forceUniformExpansion = True; *)

PairsToGlue = 0;
(* PairsToGlue = 1; *)
PairsToGlue = 18;

computeDimensions = False;
(* computeDimensions = True; *)

rotateAngle = 0;
rotateAngle = Pi/2;

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

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)

(* shape = "Octahedron"; *)
shape = "Icosahedron";

(* shape = "RhombicDodecahedron"; *)
(* shape = "RhombicTriacontahedron"; *)

(* shape = "Cuboctahedron"; *)
(* shape = "Icosidodecahedron"; *)

(* shape = "PentagonalHexecontahedron"; *)
(* shape = "SmallRhombicosidodecahedron"; *)


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
    edgeCenter = N[vertices[[edgeVertex1]] + vertices[[edgeVertex2]]]/2;
    normalVector = Cross[edgeVector, N[vertices[[edgeVertex1]]]];
    (* there's no guarantee normalVector will point to face1 or face2. We quickly check here to ensure it's pointing to face1 *)
    If[normalVector.( faceCenters[[face1]] - faceCenters[[face2]] ) < 0, normalVector = - normalVector];
    (* twist the normal vector *)
    rotatedNormalVector = RotationMatrix[rotateAngle, edgeCenter].normalVector;

    noOverlapCondition = ( pieceCoordinates[[face1]] - pieceCoordinates[[face2]] ).rotatedNormalVector >= - epsilon * 0.00001 ;
    crossProduct = Cross[pieceCoordinates[[face1]] - pieceCoordinates[[face2]], normalVector];

    AppendTo[conditions, noOverlapCondition];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, - crossProduct[[dim]] >= -epsilon * 0.00001 ];
        AppendTo[conditions, crossProduct[[dim]] >= -epsilon * 0.00001 ];
    ];
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

Rotated angle = Pi/2,
Apply rotation only to no overlap condition:

piece centers and motions
{{-0.458793973, 0., 0.6005691080000001}, {-1., 0., -1.6180339890000002}}
{{-0.141775135, -0.43633899800000003, 0.6005691080000001}, {-1., 0., -1.6180339890000002}}
{{-0.141775135, 0.43633899800000003, 0.6005691080000001}, {-1., 0., -1.6180339890000002}}
{{0.371172121, -0.26967233100000004, 0.6005691080000001}, {-1., 0., -1.6180339890000002}}
{{0.371172121, 0.26967233100000004, 0.6005691080000001}, {-1., 0., -1.6180339890000002}}
{{0.6005691080000001, -0.43633899800000003, 0.141775135}, {-1., 0., -1.6180339890000002}}
{{0.6005691080000001, 0.43633899800000003, 0.141775135}, {-1., 0., -1.6180339890000002}}
{{0.742344243, 0., -0.141775135}, {-1., 0., -1.6180339890000002}}
{{-0.22939698700000002, 0.70601133, 0.141775135}, {-0.6180339890000001, -1.175570505, 0.38196601100000005}}
{{-0.22939698700000002, -0.70601133, 0.141775135}, {-0.6180339890000001, 1.175570505, 0.38196601100000005}}
{{0.22939698700000002, 0.70601133, -0.141775135}, {0.6180339890000001, -1.175570505, -0.38196601100000005}}
{{0.22939698700000002, -0.70601133, -0.141775135}, {0.6180339890000001, 1.175570505, -0.38196601100000005}}
{{-0.742344243, 0., 0.141775135}, {1., 0., 1.6180339890000002}}
{{-0.6005691080000001, -0.43633899800000003, -0.141775135}, {1., 0., 1.6180339890000002}}
{{-0.6005691080000001, 0.43633899800000003, -0.141775135}, {1., 0., 1.6180339890000002}}
{{-0.371172121, -0.26967233100000004, -0.6005691080000001}, {1., 0., 1.6180339890000002}}
{{-0.371172121, 0.26967233100000004, -0.6005691080000001}, {1., 0., 1.6180339890000002}}
{{0.141775135, -0.43633899800000003, -0.6005691080000001}, {1., 0., 1.6180339890000002}}
{{0.141775135, 0.43633899800000003, -0.6005691080000001}, {1., 0., 1.6180339890000002}}
{{0.458793973, 0., -0.6005691080000001}, {1., 0., 1.6180339890000002}}
Groups
{8, 1, 1, 1, 1, 8}
Group count
6

 *)