#!/usr/bin/env wolframscript

(* 
v3, start from scratch

If we only keep the first 3 constraints, the dimension is the same: 2

 *)
faceTemplates = {
    {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}} (* central surface, trapezoid *)
    ,{{1, 0, 1}, {0, 1, 1}, {1, -1, 0}} (* pyramid*)
    ,{{1, 0, 1}, {0, 1, 1}, {-1, -1, 1}} (* central small surface, trapezoid, diff from v2 *)
    (* ,{{1, 0, 1}, {-1, 1, 0}, {-1, -1, 1}} (* front-up and back-right, extended small surface *) *)
    (* ,{{1, 0, 1}, {0, 1, 1}, {0, -1, 1}} (* stopper of cavity, should be inactive *) *)
    (* ,{{1, 0, 1}, {0, 1, 1}, {1, 1, 0}} (* extruded face, optional, doesn't change behavior *) *)
};
faces = {};
rotations = { IdentityMatrix[3] };
generators = { RotationMatrix[Pi/2, {0, 0, 1}], RotationMatrix[Pi/2, {1, 0, 0}] };

For[round = 1, round <= 4, round++,
    For[gid = 1, gid <= Length[generators], gid ++,
        generator = generators[[gid]];
        newRotations = {};
        For[rid = 1, rid <= Length[rotations], rid ++,
            rotation = rotations[[rid]];
            newRotation = generator.rotation;
            AppendTo[newRotations, newRotation];
        ];
        rotations = Union[rotations, newRotations];
        Print[Length[rotations]];
    ];
];

For[tid = 1, tid <= Length[faceTemplates], tid ++,
    faceTemplate = faceTemplates[[tid]];
    For[rid = 1, rid <= Length[rotations], rid ++,
        rotation = rotations[[rid]];
        newFace = Map[#.rotation&, faceTemplate];
        AppendTo[faces, newFace];
    ];    
];


pieces = Map[First, faces];
pieces = Sort[Union[pieces]];

Print["piece count"];
Print[Length[pieces]];

(* faces = Union[faces]; *)

Print["face count"];
Print[Length[faces]];

pieceLookup = <||>;

Print["pieces"];
(* Print[pieces]; *)
For[pid = 1, pid <= Length[pieces], pid ++,
    piece = pieces[[pid]];
    Print[{pid, piece}];
    pieceLookup[piece] = pid;
];

(* Exit[]; *)

(* Print[pieceLookup]; *)

pieceCoordinates = Array[x, {Length[pieces], 3}];

pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {
};

For[fid = 1, fid <= Length[faces], fid++,
    face = faces[[fid]];
    piece1Id = pieceLookup[face[[1]]];
    piece2Id = pieceLookup[face[[2]]];
    condition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).face[[3]] >= 0 ;
    AppendTo[conditions, condition];
];

For[dim = 1, dim <= 3, dim++,
    AppendTo[conditions, Total[Table[pieceCoordinates[[k]][[dim]], {k, 1, Length[pieces]}]] == 0];
];


tightnessThreshold = 0.1;
rows = {};

For[cid = 1, cid <= Length[conditions], cid++,
    Print[{"Condition", cid, Length[conditions]}];
    condition = conditions[[cid]];
    Print[condition];
    If[condition[[0]] == GreaterEqual,
        newCondition = condition;
        newCondition[[2]] = tightnessThreshold;
        Print[{"modified condition", newCondition}];
        addedConditions = Append[conditions, newCondition];
        combinedConditions = Apply[ And, addedConditions];
        solutions = FindInstance[
            combinedConditions,
            pieceCoordinatesFlat,
            Reals
        ];
        Print["solution count:"]
        Print[Length[solutions]];
        If[Length[solutions] == 0,
            Print["Tight inequality!"];
            row = Coefficient[condition[[1]], #] & /@ pieceCoordinatesFlat;
            AppendTo[rows, row];
        ];
    ];

    If[condition[[0]] == Equal,
        Print["Equality constraint!"];
        row = Coefficient[condition[[1]], #] & /@ pieceCoordinatesFlat;
        AppendTo[rows, row];
    ];
];
Print[];
Print["rows row count: ", Length[rows]];
Print["rows col count: ", Length[rows[[1]] ]];

rowsRank = MatrixRank[rows, Tolerance -> 10^-8];
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



(* 
As same as v2

rows row count: 75
rows col count: 36
rows rank: 34
coneDimension: 2
{5.685468786330354, 5.685468786330354, 5.685468786330354, 5.522046047070366, 5.5220460470703605, 5.186288078271688, 5.186288078271688, 5.186288078271687, 5.081009134649531, 5.081009134649531, 5.081009134649531, 4.8267754611579825, 4.826775461157981, 4.47213595499958, 4.4357011485968325, 4.43570114859683, 4.435701148596829, 3.741657386773942, 3.7416573867739413, 3.741657386773939, 3.6309883666551763, 3.6309883666551754, 3.630988366655174, 3.464101615137755, 3.4641016151377544, 3.4641016151377535, 2.8651781971803687, 2.8651781971803683, 2.8139544151232725, 2.8139544151232725, 2.813954415123271, 2.0453229998237523, 2.0453229998237514, 2.045322999823751, 1.4085897194468968*^-15, 5.400908741023951*^-16}
numerical Rank: 34
numerical Dimension: 2


 *)
