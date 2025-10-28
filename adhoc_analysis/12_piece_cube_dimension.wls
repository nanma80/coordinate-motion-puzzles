#!/usr/bin/env wolframscript

(* 
(+x, +y, +z): Front, Right, Up

faceUFoverUR = {1, 1, 2};

First 2 vectors: two pieces. The first piece goes over (outside of) the second
Third vector: norm of the contacting surface of the two pieces, going from second to first piece

UF and RB are over UR. UF and RB has a contact surface.
,{{1, 0, 1}, {-1, 1, 0}, {0, 1, -1}}
After adding this constraint, 5+1+1+5 no longer valid
 *)
faceTemplates = {
    {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}}
    ,{{1, 0, 1}, {0, 1, 1}, {1, -1, 0}}
    ,{{1, 0, 1}, {0, 1, 1}, {0, -1, 1}}
    ,{{1, 0, 1}, {-1, 1, 0}, {0, -1, 1}}
    (* ,{{1, 0, 1}, {1, 0, -1}, {0, -1, 1}} *)
    (* ,{{1, 0, 1}, {0, -1, 1}, {1, 1, 0}} *)
    (* ,{{1, 0, 1}, {-1, 1, 0}, {0, -1, 1}} *)
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

rows row count: 27
rows col count: 36
rows rank: 26
coneDimension: 10
{5.291502622129181, 5.2915026221291805, 5.231569255668225, 5.231569255668224, 5.231569255668224, 4.472135954999581, 4.47213595499958, 4.472135954999577, 3.7416573867739418, 3.741657386773941, 3.74165738677394, 3.4641016151377544, 3.4641016151377544, 3.4641016151377535, 2.613125929752753, 2.613125929752753, 2.6131259297527523, 2.4494897427831783, 2.449489742783177, 1.9999999999999998, 1.621938076236889, 1.6219380762368887, 1.6219380762368882, 1.0823922002923942, 1.0823922002923942, 1.0823922002923938, 5.740845412558077*^-16}
numerical Rank: 26
numerical Dimension: 10

 *)
