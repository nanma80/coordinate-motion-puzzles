#!/usr/bin/env wolframscript

(* 

This script analyzes an abstract puzzle defined in a Bridges 2026 paper by George Bell and Nan Ma. It has 3 modes:
- if computeDimensions = False and findExtremeRays = False, it computes an arbitrary solution. In this mode, if we set PairsToGlue = 1, it takes a pair of pieces and glue them together
- if computeDimensions = True, it computes the dimension of the convex cone that contains the the valid motions
- if computeDimensions = False and findExtremeRays = True, it computes the extreme rays and outputs the different types of groupings.

To run this script, install https://www.wolfram.com/wolframscript/
Alternatively, you can copy the script into Mathematica.

*)

PairsToGlue = 0;
(* PairsToGlue = 1; *)

computeDimensions = False;
(* computeDimensions = True; *)

findExtremeRays = False;
findExtremeRays = True;

If[findExtremeRays && PairsToGlue > 0,
    Print["Finding vertices of the feasible set. We are setting gluing to false"];
    PairsToGlue = 0;
];

If[computeDimensions,
    Print["Computing dimensions. We are setting findExtremeRays and gluing to false"];
    PairsToGlue = 0;
    findExtremeRays = False;
];

(* shape = "Tetrahedron"; *)
(* shape = "Cube"; *)
(* shape = "Dodecahedron"; *)
(* shape = "Octahedron"; *)
(* shape = "Icosahedron"; *)

shape = "RhombicDodecahedron";
(* shape = "RhombicTriacontahedron"; *)


phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;
findExtremeRaysBatchSize = If[shape == "Icosahedron", 1000, 300];


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
offsetCondition = pieceCoordinatesFlat.faceCentersFlat == Norm[faceCentersFlat]^2;
If[!computeDimensions && !findExtremeRays,
    AppendTo[conditions, offsetCondition];
];

(* center of mass constraint. Equation (2) in the paper *)
For[dim = 1, dim <= 3, dim++,
    AppendTo[conditions, Total[Table[pieceCoordinates[[k]][[dim]], {k, 1, faceCount }]] == 0];
];

For[edgeId = 1, edgeId <= Length[edgeToFaces], edgeId++,
    edgeIndex = edgeToFaces[[edgeId]];
    face1 = edgeIndex[[2]][[1]];
    face2 = edgeIndex[[2]][[2]];
    
    vectorBetweenPieceCenters = N[faceCenters[[face1]] - faceCenters[[face2]]];

    (* Pieces cannot overlap. Inequalities (3) in the paper *)
    noOverlapCondition = ( pieceCoordinates[[face1]] - pieceCoordinates[[face2]] ).vectorBetweenPieceCenters >= - epsilon * 0.00001 ;
    AppendTo[conditions, noOverlapCondition];

    (* Pieces move along the vector between piece centers. Equation (1) in the paper *)
    crossProduct = Cross[pieceCoordinates[[face1]] - pieceCoordinates[[face2]], vectorBetweenPieceCenters];

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

Print["condition count: ", Length[conditions]];

computeDimensionsModule[] := Module[
    {},
    tightnessThreshold = 0.1;
    rows = {};

    For[cid = 1, cid <= Length[conditions], cid++,
        condition = conditions[[cid]];
        If[condition[[0]] == GreaterEqual,
            newCondition = condition;
            newCondition[[2]] = tightnessThreshold;
            addedConditions = Append[conditions, newCondition];
            combinedConditions = Apply[ And, addedConditions];
            solutions = FindInstance[
                combinedConditions,
                pieceCoordinatesFlat,
                Reals
            ];
            If[Length[solutions] == 0,
                row = Coefficient[condition[[1]], #] & /@ pieceCoordinatesFlat;
                AppendTo[rows, row];
            ];
        ];

        If[condition[[0]] == Equal,
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

    Print["numerical Rank: ", numericalRank];
    numericalConeDimension = Length[rows[[1]] ] - numericalRank;
    Print["numerical Dimension: ", numericalConeDimension];
    nullSpace = NullSpace[rows];
    Print["nullSpace dimension: ", Length[nullSpace]];
    nullSpace
];

computeExtremeRays[] := Module[
    {},
    Print["To find extreme rays, first compute dimensions and a set of base vectors"];
    nullSpace = computeDimensionsModule[];
    AppendTo[conditions, offsetCondition];

    Print["Starting to find extreme rays of the feasible set"];
    discoveredVertices = {};

    For[directionId = 1, directionId <= findExtremeRaysBatchSize, directionId++,
        If[Mod[directionId, 10000] == 0,
            Print[{DateString["ISODateTime"], "Find vertex using random directions", directionId, findExtremeRaysBatchSize}];
        ];
        randomDirectionInNullSpace = RandomVariate[NormalDistribution[], Length[nullSpace]];
        randomDirection = randomDirectionInNullSpace.nullSpace;
        combinedConditions = Apply[ And, conditions];

        maxSolution = NMaximize[
            {randomDirection.pieceCoordinatesFlat, combinedConditions},
            pieceCoordinatesFlat
        ];
        If[Head[maxSolution] === NMaximize,
            Print["find extreme rays: maxSolution is not evaluated. Skipping"];
            Print[maxSolution];
            Continue[];
        ];
        If[Length[maxSolution] == 0,
            Print["find extreme rays: maxSolution is empty. Cannot find extreme rays. Skipping"];
            Continue[];
        ];
        If[maxSolution[[1]] == -Infinity || maxSolution[[1]] == Infinity,
            Print["find extreme rays: feasible set is empty or the value is infinity. Cannot find extreme rays. Skipping"];
            Print[maxSolution];
            Continue[];
        ];
        preciseValues = pieceCoordinates/.maxSolution[[2]];
        If[NumberQ[ preciseValues[[1]][[1]] ] == False,
            Print["find extreme rays: preciseValues is not a number. Skipping"];
            Print[preciseValues];
            Print[maxSolution];
            Continue[];
        ];
        values = Round[preciseValues, 10000 epsilon];
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
        If[Mod[vid, 10] == 0,
            Print[{DateString["ISODateTime"], "Checking extreme rays", vid, Length[discoveredVertices]}];
        ];
        vert = discoveredVertices[[vid]];
        
        vertexIsTrue = isExtreme[Flatten[vert[[1]]], allDiscoveredVerticesValues ];
        If[vertexIsTrue,
            Print["True extreme ray:"];
            Print[vert[[2]]];
            Print[vert[[1]]];
            vertexToSave = <||>;
            vertexToSave["grouping"] = vert[[2]];
            vertexToSave["pieceAndMotion"] = Transpose[ {faceCenters, vert[[1]]} ];
            AppendTo[feasibleSetVertices, vertexToSave],
            Print["False extreme ray. Ignore"];
            Print[vert[[2]]];
            Print[vert[[1]]];
        ];
    ];
    Print["after true extreme ray check. sort by groupings"];
    feasibleSetVertices = SortBy[feasibleSetVertices, #[["grouping"]]&];
    Print["after sorting. all groupings"];
    For[vid = 1, vid <= Length[feasibleSetVertices], vid++,
        Print[ feasibleSetVertices[[vid]][["grouping"]] ];
    ];
    groupings = Map[#[["grouping"]]&, feasibleSetVertices ];
    groupingsTally = Tally[groupings];
    Print["tally of groupings"];
    Map[Print, groupingsTally];

    Print["true extreme ray count"];
    Print[Length[feasibleSetVertices]];
];

findOneSolution[] := Module[
    {},
    combinedConditions = Apply[ And, conditions];

    Print["Starting FindInstance"];
    solutions = FindInstance[
        combinedConditions,
        pieceCoordinatesFlat,
        Reals
    ];

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

        motionTally = Tally[Map[Last, pieceAndMotion]];
        tallyCounts = Map[#[[2]]&, motionTally];
        Print["Groups"];
        Print[ tallyCounts ];
        Print["Group count"]
        Print[ tallyCounts//Length ];
    ];
];

If[computeDimensions,
    computeDimensionsModule[],

    If[findExtremeRays,
        computeExtremeRays[],
        findOneSolution[]
    ];
];
