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
findVerticesBatchSize = 1000;

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)

(* shape = "Octahedron"; *)
shape = "Icosahedron";

(* shape = "RhombicDodecahedron"; *)
(* shape = "RhombicTriacontahedron"; *)

(* shape = "Cuboctahedron"; *)
(* shape = "Icosidodecahedron"; *)

(* shape = "NonagonalAntiprism"; *)

outputFolder = "output";
dataFolder = "data";
vertexFileName = shape <> "_vertices.wl";
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

faceCount = Length[faces];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];

(* each face becomes a piece *)
pieceCoordinates = Array[x, {faceCount, 3}];
pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {};

If[!computeDimensions,
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

For[gid = 1, gid <= PairsToGlue, gid++,
    edgeIndex = edgeToFaces[[gid]];
    gluedPair = edgeIndex[[2]];
    If[ StringContainsQ[shape, "Antiprism"] && (gluedPair[[1]] < 3 || gluedPair[[2]] < 3), 
        Print["Skipping a pair of pieces for antiprism: "];
        Print[gluedPair];
        Continue[];
    ];
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
    Print["numerical Dimension: ", numericalConeDimension],

    If[findVertices,
        (* find vertices *)
        Print["Starting to find vertices of the feasible set"];
        discoveredVertices = {};

        For[directionId = 1, directionId <= findVerticesBatchSize, directionId++,
            randomDirection = RandomVariate[NormalDistribution[], Length[pieceCoordinatesFlat]];
            combinedConditions = Apply[ And, conditions];

            maxSolution = NMaximize[
                {randomDirection.pieceCoordinatesFlat, combinedConditions},
                pieceCoordinatesFlat
            ];
            If[maxSolution[[1]] == -Infinity,
                (* Print["Find vertices: feasible set is empty. Cannot find vertices. Skipping"]; *)
                Continue[];
            ];
            preciseValues = pieceCoordinates/.maxSolution[[2]];
            values = Round[preciseValues, 1000 epsilon];
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
            If[Mod[vid, 20] == 0,
                Print[{DateString["ISODateTime"], "Checking vertex", vid, Length[discoveredVertices]}];
            ];
            vert = discoveredVertices[[vid]];
            
            vertexIsTrue = isExtreme[Flatten[vert[[1]]], allDiscoveredVerticesValues ];
            If[vertexIsTrue,
                Print["True vertex:"];
                Print[vert[[2]]];
                Print[vert[[1]]];
                AppendTo[feasibleSetVertices, vert];
            ];
        ];
        Print["after true vertex check"];
        For[vid = 1, vid <= Length[feasibleSetVertices], vid++,
            Print[ feasibleSetVertices[[vid]][[2]] ];
        ];
        Print["true vertex count"];
        Print[Length[feasibleSetVertices]];
        Save[vertexFileNameFull, feasibleSetVertices];
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

RhombicTriacontahedron: dimension = 4. The region has 12 vertices. All symmetricaly. One of them is as follows.
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

Cuboctahedron: dimension = 4, split in half
piece centers and motions
{{-0.666666667, 0., -0.471404521}, {-1.414213562, 0., 1.}}
{{-0.666666667, 0., 0.471404521}, {-1.414213562, 0., 1.}}
{{-0.5, -0.5, 0.}, {-1.414213562, 0., 1.}}
{{-0.5, 0.5, 0.}, {-1.414213562, 0., 1.}}
{{0., -0.666666667, 0.471404521}, {-1.414213562, 0., 1.}}
{{0., 0., 0.707106781}, {-1.414213562, 0., 1.}}
{{0., 0.666666667, 0.471404521}, {-1.414213562, 0., 1.}}
{{0., -0.666666667, -0.471404521}, {1.414213562, 0., -1.}}
{{0., 0., -0.707106781}, {1.414213562, 0., -1.}}
{{0., 0.666666667, -0.471404521}, {1.414213562, 0., -1.}}
{{0.5, -0.5, 0.}, {1.414213562, 0., -1.}}
{{0.5, 0.5, 0.}, {1.414213562, 0., -1.}}
{{0.666666667, 0., -0.471404521}, {1.414213562, 0., -1.}}
{{0.666666667, 0., 0.471404521}, {1.414213562, 0., -1.}}
Groups
{7, 7}
Group count
2

Icosidodecahedron: dimension = 6, split in half
piece centers and motions
{{-1.484688486, 0., -0.283550269}, {-1.6180339890000002, -1.175570505, -1.}}
{{-1.231073415, 0., 0.6155367070000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-1.2011382160000001, -0.8726779960000001, 0.283550269}, {-1.6180339890000002, -1.175570505, -1.}}
{{-1.2011382160000001, 0.8726779960000001, 0.283550269}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.995959314, -0.723606798, -0.6155367070000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.995959314, 0.723606798, -0.6155367070000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.917587947, 0., -1.2011382160000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.742344243, -0.5393446630000001, 1.2011382160000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.458793973, -1.412022659, -0.283550269}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.380422607, -1.170820393, 0.6155367070000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.283550269, -0.8726779960000001, -1.2011382160000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.283550269, 0.8726779960000001, -1.2011382160000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{0., 0., -1.37638192}, {-1.6180339890000002, -1.175570505, -1.}}
{{0.380422607, -1.170820393, -0.6155367070000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{0.458793973, -1.412022659, 0.283550269}, {-1.6180339890000002, -1.175570505, -1.}}
{{0.742344243, -0.5393446630000001, -1.2011382160000001}, {-1.6180339890000002, -1.175570505, -1.}}
{{-0.742344243, 0.5393446630000001, 1.2011382160000001}, {1.6180339890000002, 1.175570505, 1.}}
{{-0.458793973, 1.412022659, -0.283550269}, {1.6180339890000002, 1.175570505, 1.}}
{{-0.380422607, 1.170820393, 0.6155367070000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0., 0., 1.37638192}, {1.6180339890000002, 1.175570505, 1.}}
{{0.283550269, -0.8726779960000001, 1.2011382160000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0.283550269, 0.8726779960000001, 1.2011382160000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0.380422607, 1.170820393, -0.6155367070000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0.458793973, 1.412022659, 0.283550269}, {1.6180339890000002, 1.175570505, 1.}}
{{0.742344243, 0.5393446630000001, -1.2011382160000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0.917587947, 0., 1.2011382160000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0.995959314, -0.723606798, 0.6155367070000001}, {1.6180339890000002, 1.175570505, 1.}}
{{0.995959314, 0.723606798, 0.6155367070000001}, {1.6180339890000002, 1.175570505, 1.}}
{{1.2011382160000001, -0.8726779960000001, -0.283550269}, {1.6180339890000002, 1.175570505, 1.}}
{{1.2011382160000001, 0.8726779960000001, -0.283550269}, {1.6180339890000002, 1.175570505, 1.}}
{{1.231073415, 0., -0.6155367070000001}, {1.6180339890000002, 1.175570505, 1.}}
{{1.484688486, 0., 0.283550269}, {1.6180339890000002, 1.175570505, 1.}}
Groups
{16, 16}
Group count
2

DeltoidalIcositetrahedron (dual of Rhombicuboctahedron): dimension = 4, split in half
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

SmallRhombicuboctahedron (aka Rhombicuboctahedron): dim = 6
piece centers and motions
{{-1.207106781, 0., 0.}, {-0.41421356200000004, 0., 0.}}
{{-0.853553391, -0.853553391, 0.}, {-0.29289321900000004, -0.29289321900000004, 0.}}
{{-0.853553391, 0.853553391, 0.}, {-0.29289321900000004, 0.29289321900000004, 0.}}
{{0., -1.207106781, 0.}, {0., -0.41421356200000004, 0.}}
{{-0.853553391, 0., -0.853553391}, {0., 0., -1.}}
{{-0.73570226, -0.73570226, -0.73570226}, {0., 0., -1.}}
{{-0.73570226, 0.73570226, -0.73570226}, {0., 0., -1.}}
{{0., -0.853553391, -0.853553391}, {0., 0., -1.}}
{{0., 0., -1.207106781}, {0., 0., -1.}}
{{0., 0.853553391, -0.853553391}, {0., 0., -1.}}
{{0.73570226, -0.73570226, -0.73570226}, {0., 0., -1.}}
{{0.73570226, 0.73570226, -0.73570226}, {0., 0., -1.}}
{{0.853553391, 0., -0.853553391}, {0., 0., -1.}}
{{-0.853553391, 0., 0.853553391}, {0., 0., 1.}}
{{-0.73570226, -0.73570226, 0.73570226}, {0., 0., 1.}}
{{-0.73570226, 0.73570226, 0.73570226}, {0., 0., 1.}}
{{0., -0.853553391, 0.853553391}, {0., 0., 1.}}
{{0., 0., 1.207106781}, {0., 0., 1.}}
{{0., 0.853553391, 0.853553391}, {0., 0., 1.}}
{{0.73570226, -0.73570226, 0.73570226}, {0., 0., 1.}}
{{0.73570226, 0.73570226, 0.73570226}, {0., 0., 1.}}
{{0.853553391, 0., 0.853553391}, {0., 0., 1.}}
{{0., 1.207106781, 0.}, {0., 0.41421356200000004, 0.}}
{{0.853553391, -0.853553391, 0.}, {0.29289321900000004, -0.29289321900000004, 0.}}
{{0.853553391, 0.853553391, 0.}, {0.29289321900000004, 0.29289321900000004, 0.}}
{{1.207106781, 0., 0.}, {0.41421356200000004, 0., 0.}}
Groups
{1, 1, 1, 1, 9, 9, 1, 1, 1, 1}
Group count
10

SquareAntiprism: (this motion may not represent an extreme case)
piece centers and motions
{8, {-0.569035594, 0., -0.140149403}, {-0.6368607310000001, 0., -0.656854249}}
{5, {-0.40236892700000004, -0.40236892700000004, 0.140149403}, {-0.34831070000000003, -0.6966213990000001, -0.171572875}}
{9, {0., -0.569035594, -0.140149403}, {-0.34831070000000003, -0.6966213990000001, -0.171572875}}
{1, {0., 0., -0.42044820800000005}, {-0.34831070000000003, 0., -1.}}
{4, {-0.40236892700000004, 0.40236892700000004, 0.140149403}, {-0.34831070000000003, 0.6966213990000001, -0.171572875}}
{7, {0., 0.569035594, -0.140149403}, {-0.34831070000000003, 0.6966213990000001, -0.171572875}}
{2, {0., 0., 0.42044820800000005}, {0.34831070000000003, 0., 1.}}
{6, {0.40236892700000004, -0.40236892700000004, 0.140149403}, {0.6368607310000001, -0.288550032, 0.5147186260000001}}
{3, {0.40236892700000004, 0.40236892700000004, 0.140149403}, {0.6368607310000001, 0.288550032, 0.5147186260000001}}
{10, {0.569035594, 0., -0.140149403}, {0.7563820680000001, 0., 0.31370849900000003}}
Groups
{1, 2, 1, 2, 1, 1, 1, 1}
Group count
8

PentagonalAntiprism: (this motion may not represent an extreme case)
piece centers and motions
{5, {-0.742344243, 0., 0.141775135}, {-0.5450849720000001, -0.22451398800000003, 0.236067977}}
{10, {-0.6005691080000001, -0.43633899800000003, -0.141775135}, {-0.5450849720000001, -0.22451398800000003, 0.236067977}}
{4, {-0.22939698700000002, 0.70601133, 0.141775135}, {-0.30901699400000004, 0.5020285400000001, -0.236067977}}
{9, {-0.6005691080000001, 0.43633899800000003, -0.141775135}, {-0.30901699400000004, 0.5020285400000001, -0.236067977}}
{6, {-0.22939698700000002, -0.70601133, 0.141775135}, {-0.163118961, -0.5020285400000001, 0.5278640450000001}}
{2, {0., 0., 0.425325404}, {-0.072949017, -0.22451398800000003, 1.}}
{1, {0., 0., -0.425325404}, {0.072949017, 0.22451398800000003, -1.}}
{8, {0.22939698700000002, 0.70601133, -0.141775135}, {0.163118961, 0.5020285400000001, -0.5278640450000001}}
{7, {0.6005691080000001, -0.43633899800000003, 0.141775135}, {0.30901699400000004, -0.5020285400000001, 0.236067977}}
{11, {0.22939698700000002, -0.70601133, -0.141775135}, {0.30901699400000004, -0.5020285400000001, 0.236067977}}
{3, {0.6005691080000001, 0.43633899800000003, 0.141775135}, {0.5450849720000001, 0.22451398800000003, -0.236067977}}
{12, {0.742344243, 0., -0.141775135}, {0.5450849720000001, 0.22451398800000003, -0.236067977}}
Groups
{2, 2, 1, 1, 1, 1, 2, 2}
Group count
8

HexagonalAntiprism: (this motion may not represent an extreme case)
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

Platonic solids:
Tetrahedron: dim = 1
Cube: dim = 1
Octahedron: dim = 3
Dodecahedron: dim = 1
Icosahedron: dim = 9

Kepler-Poinsot (nonconvex):

SmallStellatedDodecahedron: dim = 1 
GreatDodecahedron: dim = 1
GreatStellatedDodecahedron: dim = 1 
GreatIcosahedron: dim = 9 (same grouping as icosahedron? maybe)

All Catalan solids:

RhombicDodecahedron: dim = 2, Splinter 12
RhombicTriacontahedron: dim = 4, Splinter 30

(Triakis/Tetrakis/Pentakis xxx: dual of truncated xxx)
TriakisTetrahedron: dim = 5
TetrakisHexahedron: dim = 11
TriakisOctahedron: dim = 11
TriakisIcosahedron (633 constraints): dim = 29
PentakisDodecahedron (633 constraints): dim = 29

DeltoidalIcositetrahedron (dual of Rhombicuboctahedron): dim = 4, split in half
DeltoidalHexecontahedron (dual of rhombicosidodecahedron): dim = 10
PentagonalIcositetrahedron (dual of snub cube): dim = 1
PentagonalHexecontahedron (dual of snub dodecahedron): dim = 1

DisdyakisDodecahedron (dual of truncated cuboctahedron aka great rhombicuboctahedron, 507): dim = 23
DisdyakisTriacontahedron (dual of truncated icosidodecahedron aka great rhombicosidodecahedron, 1263): dim = 59

All Archimedean solids:
command: Map[Print, PolyhedronData["Archimedean"]]

TruncatedTetrahedron: dim = 1
TruncatedCube: dim = 1
TruncatedOctahedron: dim = 1
TruncatedIcosahedron: dim = 1
TruncatedDodecahedron: dim = 1

Cuboctahedron: dim = 4, split in half
Icosidodecahedron: dim = 6, split in half

SmallRhombicuboctahedron (aka Rhombicuboctahedron): dim = 6, see motion above
SmallRhombicosidodecahedron (aka Rhombicosidodecahedron): dim = 12
GreatRhombicuboctahedron (aka truncated cuboctahedron): dim = 1
GreatRhombicosidodecahedron (aka truncated icosidodecahedron): dim = 1

SnubCube: dim = 15
SnubDodecahedron (1053 constraints, dim should be large, slow, ongoing): dim = 33

Solids with dihedral symmetry:
Useful command:
Table[ PolyhedronData[{"Antiprism", k}, "AlternateStandardNames"], {k, 3, 10}]

TriangularPrism: dim = 1
TriangularDipyramid: dim = 2, split up/down vs radially 3-way
TriangularAntiprism (octahedron): dim = 3
TriangularTrapezohedron: name not found for PolyhedronData, cube, likely dim = 1

SquarePrism (cube): dim = 1
SquareDipyramid (octahedron): dim = 3
SquareAntiprism: dim = 3
SquareTrapezohedron: bad name, likely dim = 1

PentagonalPrism: dim = 1
PentagonalDipyramid: dim = 4
PentagonalAntiprism: dim = 3
PentagonalTrapezohedron: name not found for PolyhedronData, likely dim = 1

HexagonalPrism: dim = 1
HexagonalDipyramid: name not found for PolyhedronData, likely dim = 5
HexagonalAntiprism: dim = 3
HexagonalTrapezohedron: name not found for PolyhedronData, likely dim = 1

HeptagonalPrism: dim = 1
HeptagonalDipyramid: name not found for PolyhedronData, likely dim = 6
HeptagonalAntiprism: dim = 3
HeptagonalTrapezohedron: name not found for PolyhedronData, likely dim = 1

OctagonalAntiprism: dim = 3

NonagonalAntiprism: dim = 3

DecagonalAntiprism: dim = 3

Q: why are the dimensions of all antiprisms 3? There must be a common reason

SnubSquareAntiprism: dim = 11

 *)