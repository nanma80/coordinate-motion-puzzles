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
targetPieceCount = 10;

shape = "RhombilleTiling";
(* shape = "HexagonalTiling"; *)
(* shape = "SquareTiling"; *)
(* shape = "TriangularTiling"; *)
(* shape = "TruncatedHexagonalTiling"; *)
(* shape = "TruncatedSquareTiling"; *)

(* not reliable below this line: edges are incomplete *)
(* shape = "SmallRhombitrihexagonalTiling"; *)
(* shape = "GreatRhombitrihexagonalTiling"; *)
(* shape = "SnubHexagonalTiling"; *)
(* shape = "SnubSquareTiling"; *)
(* shape = "TrihexagonalTiling"; *)
(* shape = "ElongatedTriangularTiling"; *)


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

tiling = Entity["PeriodicTiling", shape];

cell = tiling["PrimitiveUnit"];

centers0 = FullSimplify[ Mean /@ (Cases[cell, Polygon[p_] :> p, All]) ];
Print["centers0: ", centers0];
lat = tiling["TranslationVectors"][];
(* poly = First @ Cases[cell, _Polygon, All];
polys = Cases[cell, _Polygon, All];


 *)
faceCenters = RegionCentroid /@ translatedPolys;
adjacentQ[p_, q_] := RegionDimension[RegionIntersection[p, q]] == 1;

neighborFacePairs =
  Flatten[
    Table[
      If[adjacentQ[translatedPolys[[i]], translatedPolys[[j]]],
        {i, j},
        Nothing
      ],
      {i, Length[translatedPolys]}, {j, i + 1, Length[translatedPolys]}
    ],
    1
  ];


pieceCount = 0;

For[polygonRadius = 1, pieceCount < targetPieceCount, polygonRadius++,
    Print["Building faces. Radius = ", polygonRadius, " pieceCount so far: ", pieceCount];
    nmax = Ceiling[polygonRadius/Min[Norm /@ lat]];
    Print["nmax: ", nmax];
    faceCenters =
      Select[
        Flatten[
          Table[
            # + i lat[[1]] + j lat[[2]] & /@ centers0,
            {i, -nmax, nmax}, {j, -nmax, nmax}
          ],
          2
        ],
        Norm[#] <= polygonRadius &
      ];
    Print["faceCenters: ", faceCenters];
    pieceCount = Length[faceCenters];
];

faceCenters = Map[{N[#[[1]]], N[#[[2]]], 0.0}&, faceCenters];
faceCentersFlat = Flatten[faceCenters];
Print["piece (face) centers"];
Print[faceCenters];

faceCount = Length[faceCenters];
pieceCount = faceCount; (* for convenience *)
Print["piece (face) count: ", pieceCount];

closestIndex = First[Ordering[Norm /@ faceCenters, 1]];
Print["closestIndex: ", closestIndex];

minFaceCenterDistance = 1000000;

For[fid = 1, fid <= faceCount, fid++,
    If[fid == closestIndex, Continue[]];
    distance = N[ EuclideanDistance[faceCenters[[closestIndex]], faceCenters[[fid]]] ];
    minFaceCenterDistance = Min[minFaceCenterDistance, distance];
];
Print["minFaceCenterDistance: ", minFaceCenterDistance];

neighborFacePairs = {};

For[f1id = 1, f1id <= faceCount, f1id++,
    For[f2id = f1id + 1, f2id <= faceCount, f2id++,
        distance = N[ EuclideanDistance[faceCenters[[f1id]], faceCenters[[f2id]]] ];
        If[distance < minFaceCenterDistance + epsilon,
            AppendTo[neighborFacePairs, {f1id, f2id}];
        ];
    ];
];


Print["edges by face ids:"];
Print[neighborFacePairs];
Print["edge count: ", Length[neighborFacePairs]];

edgesToCloestIndex = Select[neighborFacePairs, #[[1]] == closestIndex || #[[2]] == closestIndex & ];
Print["Edgest to the closest face: ", edgesToCloestIndex, " edge count to the closest face: ", Length[edgesToCloestIndex]];

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

Export["test_plot/"<> shape <> ".png", plot];

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

    noOverlapCondition = ( pieceCoordinates[[face1]] - pieceCoordinates[[face2]] ).vector12 >= - epsilon * 0.00001 ;
    crossProduct = Cross[pieceCoordinates[[face1]] - pieceCoordinates[[face2]], vector12];

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

Print["The shape is: ", shape];




(* 

HexagonalTiling:
    dim = 1 for target piece count = 10, piece count = 12
    dim = 1 for target piece count = 20, piece count = 30
    dim = 1 unchanged when piece count -> infinity

SquareTiling:
    dim = 6 for targetPieceCount = 10; piece count = 12
    dim = 10 for targetPieceCount = 20; piece count = 32
    dim -> infinity when piece count -> infinity

TriangularTiling:
    dim = 16 for targetPieceCount = 10; piece count = 24
    dim = 38 for targetPieceCount = 30; piece count = 76
    dim -> infinity when piece count -> infinity

RhombilleTiling:
    dim = 7 for targetPieceCount = 10; piece count = 16
    dim = 7 for targetPieceCount = 20; piece count = 20
    dim = 14 for targetPieceCount = 30; piece count = 52
    dim -> infinity when piece count -> infinity

Not reliable below this line


TruncatedHexagonalTiling: dim = 1
TruncatedSquareTiling: 

SmallRhombitrihexagonalTiling:

GreatRhombitrihexagonalTiling

shape = "SnubHexagonalTiling";

shape = "ElongatedTriangularTiling";

shape = "SnubSquareTiling";
shape = "TrihexagonalTiling";


full list supported by wolfram:

{SquareTiling, HexagonalTiling, TriangularTiling, TrihexagonalTiling, SmallRhombitrihexagonalTiling, TruncatedSquareTiling, SnubSquareTiling, TruncatedHexagonalTiling, ElongatedTriangularTiling, GreatRhombitrihexagonalTiling, SnubHexagonalTiling, RhombilleTiling, DeltoidalTrihexagonalTiling, TetrakisSquareTiling, CairoPentagonalTiling, TriakisTriangularTiling, PrismaticPentagonalTiling, BisectedHexagonalTiling, FloretPentagonalTiling, DellaRobbiaWeaveTiling, PortugalTiling, StackBondTiling, HerringboneTiling, BasketweaveTiling, PersianHexagonalWeaveTiling, HopscotchTiling, StretcherBondTiling, PinwheelTiling, BrickworkSquareTiling, ChickenwireTiling, CorridorTiling, CorridorHorizontalTiling, BrickweaveTiling, TrellisTiling, HeeschIsohedralTiling, ArbitraryTriangleTiling, ArbitraryQuadrilateralTiling, PPentominoTiling, ChevronTiling, ShingleTiling, ZigzagTiling, KiteTiling, FalseCubicTiling, TrihexAndHexTiling, GlideReflectionTiling, PentagonType1Tiling, PentagonType2Tiling, PentagonType3Tiling, PentagonType4Tiling, PentagonType5Tiling, PentagonType6Tiling, PentagonType7Tiling, PentagonType8Tiling, PentagonType9Tiling, PentagonType10Tiling, PentagonType11Tiling, PentagonType12Tiling, PentagonType13Tiling, PentagonType14Tiling, PentagonType15Tiling}


 *)