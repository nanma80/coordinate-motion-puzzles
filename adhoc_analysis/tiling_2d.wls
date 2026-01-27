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
findVerticesBatchSize = 300;
(* findVerticesBatchSize = 400000; *)
targetPieceCount = 131;

(* shape = "HexagonalTiling"; *)
(* shape = "SquareTiling"; *)
(* shape = "TriangularTiling"; *)

(* shape = "TruncatedHexagonalTiling"; *)
(* shape = "TruncatedSquareTiling"; *)
(* shape = "SmallRhombitrihexagonalTiling"; *)
(* shape = "GreatRhombitrihexagonalTiling"; *)
(* shape = "SnubHexagonalTiling"; *)
(* shape = "SnubSquareTiling"; *)
(* shape = "TrihexagonalTiling"; *)
(* shape = "ElongatedTriangularTiling"; *)

shape = "RhombilleTiling";
(* shape = "CairoPentagonalTiling"; *)
(* shape = "FloretPentagonalTiling"; *)

outputFolder = FileNameJoin[{"..", "output"}];
dataFolder = "data";
vertexFileName = shape <> "_vertices.json";
vertexFileNameFull = FileNameJoin[{outputFolder, dataFolder, vertexFileName}];

Print["shape: ", shape];
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
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001;
sameCenter[set1_, set2_] := ApproxSamePoint[Total[N[set1]], Total[N[set2]]];
sameEdge[e1_, e2_] := sameCenter[e1, e2] && ((ApproxSamePoint[e1[[1]], e2[[1]]] && ApproxSamePoint[e1[[2]], e2[[2]]] ) || (ApproxSamePoint[e1[[1]], e2[[2]]] && ApproxSamePoint[e1[[2]], e2[[1]]] ));

tiling = Entity["PeriodicTiling", shape];

cell = tiling["PrimitiveUnit"];
primitivePolygons = Cases[cell, Polygon[p_] :> p, All];

lat = tiling["TranslationVectors"][];

pieceCount = 0;

For[polygonRadius = 1, pieceCount < targetPieceCount, polygonRadius++,
    Print["Building faces. Radius = ", polygonRadius, " pieceCount so far: ", pieceCount];
    nmax = Ceiling[polygonRadius/Min[Norm /@ lat]];
    polygonsTable = Table[
        Map[# + i lat[[1]] + j lat[[2]] &, primitivePolygons, {2}],
        {i, -nmax, nmax}, {j, -nmax, nmax}
    ];

    flatPolygons = FullSimplify[Flatten[polygonsTable, 2]];
    polygons = Select[flatPolygons, Norm[Mean[#]] < polygonRadius &];
    pieceCount = Length[polygons];
];

faceCenters = Map[Mean, polygons];
faceCenters = Map[{N[#[[1]]], N[#[[2]]], 0.0}&, faceCenters];
faceCentersFlat = Flatten[faceCenters];

faceCount = Length[faceCenters];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];
Print[];

edgesByFace = Map[getEdgesFromFace, N[polygons]];

neighborFacePairs = {};
edgePerpDirections = {};

For[f1id = 1, f1id <= faceCount, f1id++,
    For[f2id = f1id + 1, f2id <= faceCount, f2id++,
        edgesFace1 = edgesByFace[[f1id]];
        edgesFace2 = edgesByFace[[f2id]];
        neighborFlag = False;
        For[e1id = 1, e1id <= Length[edgesFace1] && neighborFlag == False, e1id++,
            For[e2id = 1, e2id <= Length[edgesFace2] && neighborFlag == False, e2id++,
                If[sameEdge[edgesFace1[[e1id]], edgesFace2[[e2id]]],
                    neighborFlag = True;
                    perpDirection = (edgesFace1[[e1id]][[1]] - edgesFace1[[e1id]][[2]] ).RotationMatrix[Pi/2];
                    AppendTo[neighborFacePairs, {f1id, f2id}];
                    AppendTo[edgePerpDirections, perpDirection];
                ];
            ];
        ];
    ];
];

Print["edge count: ", Length[neighborFacePairs]];

r = 0.1;
ballsPlot = Sphere[#, r] & /@ faceCenters;

rc = 0.04;
edgesPlot = Cylinder[
   {faceCenters[[#[[1]]]], faceCenters[[#[[2]]]]},
   rc
 ] & /@ neighborFacePairs;


plot = Graphics3D[
  {
    (* balls *)
    Directive[LightBlue, Specularity[White, 20]],
    ballsPlot,

    (* connections *)
    Directive[GrayLevel[0.3]],
    edgesPlot
  },
  Boxed -> False,
  ViewPoint -> {0, 0, polygonRadius},
  Lighting -> "Neutral"
];

Export["../output/test_plot/"<> shape <> "_" <> ToString[pieceCount] <> ".png", plot];

(* Exit[]; *)

(* each face becomes a piece *)
pieceCoordinates = Array[x, {faceCount, 3}];
pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {};
offsetCondition = pieceCoordinatesFlat.faceCentersFlat == Norm[faceCentersFlat]^2;
If[!computeDimensions && !findVertices,
    AppendTo[conditions, offsetCondition];
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

For[pairId = 1, pairId <= Length[neighborFacePairs], pairId++,
    pair = neighborFacePairs[[pairId]];
    face1 = pair[[1]];
    face2 = pair[[2]];
    vector12 = N[faceCenters[[face1]] - faceCenters[[face2]]];
    edgePerpDirection = edgePerpDirections[[pairId]];
    normalVector = {edgePerpDirection[[1]], edgePerpDirection[[2]], 0.0 };

    (* there's no guarantee normalVector will point to face1 or face2. We quickly check here to ensure it's pointing to face1 *)
    If[normalVector.vector12 < 0, normalVector = - normalVector];

    noOverlapCondition = ( pieceCoordinates[[face1]] - pieceCoordinates[[face2]] ).normalVector >= - epsilon * 0.00001 ;
    crossProduct = Cross[pieceCoordinates[[face1]] - pieceCoordinates[[face2]], normalVector];

    AppendTo[conditions, noOverlapCondition];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, - crossProduct[[dim]] >= -epsilon * 0.00001 ];
        AppendTo[conditions, crossProduct[[dim]] >= -epsilon * 0.00001 ];
    ];
];

For[gid = 1, gid <= PairsToGlue, gid++,
    gluedPair = neighborFacePairs[[gid]];
    Print["Gluing a pair of pieces: "];
    Print[gluedPair];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];

Print["condition count:"];
Print[Length[conditions]];

computeDimensionsModule[] := Module[
    {},
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
    Print["numerical Dimension: ", numericalConeDimension];
    nullSpace = NullSpace[rows];
    Print["nullSpace dimension: ", Length[nullSpace]];
    nullSpace
];

If[computeDimensions,
    computeDimensionsModule[],

    If[findVertices,
        (* find vertices *)
        Print["To find vertices, first compute dimensions and a set of base vectors"];
        nullSpace = computeDimensionsModule[];
        AppendTo[conditions, offsetCondition];

        Print["Starting to find vertices of the feasible set"];
        discoveredVertices = {};

        For[directionId = 1, directionId <= findVerticesBatchSize, directionId++,
            If[Mod[directionId, 10000] == 0,
                Print[{DateString["ISODateTime"], "Find vertex using random directions", directionId, findVerticesBatchSize}];
            ];
            randomDirectionInNullSpace = RandomVariate[NormalDistribution[], Length[nullSpace]];
            randomDirection = randomDirectionInNullSpace.nullSpace;
            combinedConditions = Apply[ And, conditions];

            maxSolution = NMaximize[
                {randomDirection.pieceCoordinatesFlat, combinedConditions},
                pieceCoordinatesFlat
            ];
            If[Head[maxSolution] === NMaximize,
                Print["Find vertices: maxSolution is not evaluated. Skipping"];
                Print[maxSolution];
                Continue[];
            ];
            If[Length[maxSolution] == 0,
                Print["Find vertices: maxSolution is empty. Cannot find vertices. Skipping"];
                Continue[];
            ];
            If[maxSolution[[1]] == -Infinity || maxSolution[[1]] == Infinity,
                Print["Find vertices: feasible set is empty or the value is infinity. Cannot find vertices. Skipping"];
                Print[maxSolution];
                Continue[];
            ];
            preciseValues = pieceCoordinates/.maxSolution[[2]];
            If[NumberQ[ preciseValues[[1]][[1]] ] == False,
                Print["Find vertices: preciseValues is not a number. Skipping"];
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
                Print[{DateString["ISODateTime"], "Checking vertex", vid, Length[discoveredVertices]}];
            ];
            vert = discoveredVertices[[vid]];
            
            vertexIsTrue = isExtreme[Flatten[vert[[1]]], allDiscoveredVerticesValues ];
            If[vertexIsTrue,
                Print["True vertex:"];
                Print[vert[[2]]];
                Print[vert[[1]]];
                vertexToSave = <||>;
                vertexToSave["grouping"] = vert[[2]];
                vertexToSave["pieceAndMotion"] = Transpose[ {faceCenters, vert[[1]]} ];
                AppendTo[feasibleSetVertices, vertexToSave],
                Print["False vertex. Ignore"];
                Print[vert[[2]]];
                Print[vert[[1]]];
            ];
        ];
        Print["after true vertex check. sort by groupings"];
        feasibleSetVertices = SortBy[feasibleSetVertices, #[["grouping"]]&];
        Print["after sorting. all groupings"];
        For[vid = 1, vid <= Length[feasibleSetVertices], vid++,
            Print[ feasibleSetVertices[[vid]][["grouping"]] ];
        ];
        groupings = Map[#[["grouping"]]&, feasibleSetVertices ];
        groupingsTally = Tally[groupings];
        Print["tally of groupings"];
        Map[Print, groupingsTally];

        Print["true vertex count"];
        Print[Length[feasibleSetVertices]];
        If[FileExistsQ[vertexFileNameFull], DeleteFile[vertexFileNameFull]];
        Export[vertexFileNameFull, feasibleSetVertices, "RawJSON", Compact -> True];
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

Print["shape: ", shape];
Print["piece (face) count: ", pieceCount];



(* 
If dim increases when the number of pieces go up, is this due to longer boundaries?
If dim decreases when the number of pieces go up, it means the dim is over estimated?

HexagonalTiling:
    dim = 1 for piece count = 12
    dim = 1 for piece count = 30
    dim = 1 Should be correct, because all vertices are degree 3

SquareTiling:
    dim = 6 for piece count = 12
    dim = 10 for piece count = 32
    dim -> infinity when piece count -> infinity

TriangularTiling:
    dim = 16 for piece count = 24
    dim = 38 for piece count = 76
    dim -> infinity when piece count -> infinity

RhombilleTiling:
    dim = 7 for piece count = 14
    dim = 7 for piece count = 20
    dim = 14 for piece count = 50
    dim = 14 for piece count = 60
    dim = 20 for piece count = 104
    dim = 24 for piece count = 130
    dim = 26 for piece count = 188
    likely dim -> infinity when piece count -> infinity

TruncatedHexagonalTiling:
    dim = 1 for piece count = 19
    dim = 13 for piece count = 31
    dim = 1 for piece count = 55
    dim = 1 Should be correct, because all vertices are degree 3

TruncatedSquareTiling:
    dim = 1 for piece count = 21

SmallRhombitrihexagonalTiling:
    dim = 6 for piece count = 13
    dim = 7 for piece count = 25
    dim = 13 for piece count = 43
    dim = 19 for piece count = 79
    dim = 31 for piece count = 109
    dim = 37 for piece count = 145
    likely dim -> infinity when piece count -> infinity. probably increasing

GreatRhombitrihexagonalTiling
    dim = 1 for piece count = 25
    Should be 1, because all vertices are degree 3

SnubHexagonalTiling
    dim = 18 for piece count 37
    very likely dim -> infinity when piece count -> infinity

SnubSquareTiling
    dim = 27 for piece count 48
    very likely dim -> infinity when piece count -> infinity

TrihexagonalTiling
    dim = 12 for piece count 37
    dim -> infinity when piece count -> infinity

ElongatedTriangularTiling
    dim = 14 for piece count 21
    dim -> infinity when piece count -> infinity

CairoPentagonalTiling
    dim = 2 for piece count 14
    dim = 2 for piece count 32
    dim = 2 for piece count 52
    dim = 2 for piece count 84
    dim = 2 when piece count -> infinity

FloretPentagonalTiling
    dim = 1 for piece count 24
    dim = 1 Should be correct, degree-3 vertices are connected


full list supported by wolfram:

{SquareTiling, HexagonalTiling, TriangularTiling, TrihexagonalTiling, SmallRhombitrihexagonalTiling, TruncatedSquareTiling, SnubSquareTiling, TruncatedHexagonalTiling, ElongatedTriangularTiling, GreatRhombitrihexagonalTiling, SnubHexagonalTiling, RhombilleTiling, DeltoidalTrihexagonalTiling, TetrakisSquareTiling, CairoPentagonalTiling, TriakisTriangularTiling, PrismaticPentagonalTiling, BisectedHexagonalTiling, FloretPentagonalTiling, DellaRobbiaWeaveTiling, PortugalTiling, StackBondTiling, HerringboneTiling, BasketweaveTiling, PersianHexagonalWeaveTiling, HopscotchTiling, StretcherBondTiling, PinwheelTiling, BrickworkSquareTiling, ChickenwireTiling, CorridorTiling, CorridorHorizontalTiling, BrickweaveTiling, TrellisTiling, HeeschIsohedralTiling, ArbitraryTriangleTiling, ArbitraryQuadrilateralTiling, PPentominoTiling, ChevronTiling, ShingleTiling, ZigzagTiling, KiteTiling, FalseCubicTiling, TrihexAndHexTiling, GlideReflectionTiling, PentagonType1Tiling, PentagonType2Tiling, PentagonType3Tiling, PentagonType4Tiling, PentagonType5Tiling, PentagonType6Tiling, PentagonType7Tiling, PentagonType8Tiling, PentagonType9Tiling, PentagonType10Tiling, PentagonType11Tiling, PentagonType12Tiling, PentagonType13Tiling, PentagonType14Tiling, PentagonType15Tiling}


 *)