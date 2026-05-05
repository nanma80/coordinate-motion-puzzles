#!/usr/bin/env wolframscript

(* 
(+x, +y, +z): Front, Right, Up

faceUFoverUR = {1, 1, 2};

First 2 vectors: two pieces. The first piece goes over (outside of) the second
Third vector: norm of the contacting surface of the two pieces, going from second to first piece

UF and RB are over UR. UF and RB has a contact surface.
,{{1, 0, 1}, {-1, 1, 0}, {0, 1, -1}}
After adding this constraint, 5+1+1+5 no longer valid

{1, 0, 1}, {0, 1, 1}, 

phi = (Sqrt[5] + 1)/2;

norm1 = Cross[{19-29/2 phi, 65-75/2 phi, -38+29phi} - {-8/5-6/5phi, 16/5+12/5phi, 16/5+12/5phi}, {0, 8/3+2phi, 8/3+2phi} - {-8/5-6/5phi, 16/5+12/5phi, 16/5+12/5phi}]
norm normalize to {1, 1, 2}. Consistent. Brown

{1, 0, 1}, {0, 1, 1},
norm2 = Cross[{0, 8/3+2phi, 8/3+2phi}- {8/5+6/5phi, 16/5+12/5phi, 16/5+12/5phi}, {8/3+2phi, 8/3+2phi, 8/3+2phi} - {8/5+6/5phi, 16/5+12/5phi, 16/5+12/5phi}]
norm2 ~ {0, -1, 1}, center green


{1, 0, 1}, {0, 1, 1},
norm3 = Cross[{122/3-27phi, 130/3-25phi, 16/3+4phi}-{19-29/2 phi, 65-75/2 phi, -38+29phi}, {84-52phi,0,-38+29phi}-{19-29/2 phi, 65-75/2 phi, -38+29phi}]
~ {1, 1, 0}, outside green

dim = 2 !
There's no slack compared to splinter 12

 *)



phi = (Sqrt[5] + 1)/2;

faceTemplates = {
    (* new analysis *)
    {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}}
    ,{{1, 0, 1}, {0, 1, 1}, {0, -1, 1}}
    ,{{1, 0, 1}, {0, 1, 1}, {1, 1, 0}}

    (* ,{{1, 0, 1}, {0, 1, 1}, {1, -1, 0}} *)
    (* ,{{1, 0, 1}, {-1, 1, 0}, {0, -1, 1}} *)
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

(* specified constraints *)
For[fid = 1, fid <= Length[faces], fid++,
    face = faces[[fid]];
    piece1Id = pieceLookup[face[[1]]];
    piece2Id = pieceLookup[face[[2]]];
    condition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).face[[3]] >= 0 ;
    AppendTo[conditions, condition];
];

(* center of gravity is at the origin *)
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

rows row count: 51
rows col count: 36
rows rank: 34
coneDimension: 2
{5.515632792514014, 5.515632792514013, 5.492315091161141, 5.492315091161141, 5.492315091161141, 4.828427124746191, 4.82842712474619, 4.828427124746189, 4.1486265861038865, 4.148626586103886, 4.148626586103885, 3.4641016151377553, 3.4641016151377553, 3.464101615137755, 3.464101615137755, 3.464101615137755, 3.464101615137755, 3.4641016151377544, 3.4641016151377544, 2.8284271247461903, 2.4154657810458753, 2.4154657810458753, 2.4154657810458744, 2.000000000000001, 1.9999999999999998, 1.9999999999999993, 1.6699992362489338, 1.6699992362489333, 1.6699992362489333, 1.256103060319512, 1.2561030603195118, 0.8284271247461905, 0.8284271247461902, 0.8284271247461897, 1.1263065357759656*^-15, 8.449237438836995*^-16}
numerical Rank: 34
numerical Dimension: 2

 *)