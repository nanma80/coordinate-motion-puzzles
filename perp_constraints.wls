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
findVerticesBatchSize = 30000;

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)

(* shape = "Octahedron"; *)
(* shape = "Icosahedron"; *)

(* shape = "RhombicDodecahedron"; *)
(* shape = "RhombicTriacontahedron"; *)

(* shape = "Cuboctahedron"; *)
(* shape = "Icosidodecahedron"; *)

shape = "TetrahedronTenCompound";

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
    AppendTo[conditions, pieceCoordinatesFlat.faceCentersFlat == Norm[faceCentersFlat]^2];
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

See results.txt

 *)