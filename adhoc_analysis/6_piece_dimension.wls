#!/usr/bin/env wolframscript

(* 
(+x, +y, +z): Right, Back, Up

faceFrontUpOutside = {0, -1, 1};
faceFrontUpInside = {0, -1, -1};

First 2 vectors: two pieces. The first piece goes over (outside of) the second
Third vector: norm of the contacting surface of the two pieces, going from second to first piece

 *)
faceTemplates = {
    {{0, -1, 0}, {0, 0, 1}, {0, -1, 1}},
    {{0, -1, 0}, {0, 0, 1}, {0, -1, -1}}
};
(* faceTemplate = {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}}; *)
faces = {};
rotations = { IdentityMatrix[3] };
generators = { RotationMatrix[Pi, {0, 0, 1}], RotationMatrix[2 Pi/3, {1, 1, 1}] };

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

rows row count: 15
rows col count: 18
rows rank: 14
coneDimension: 4
{2.4494897427831788, 2.4494897427831783, 2.4494897427831783, 2.4494897427831783, 2.4494897427831783, 2.449489742783178, 2.449489742783178, 2.4494897427831774, 2., 1.9999999999999998, 1.9999999999999991, 1.4142135623730951, 1.4142135623730947, 1.4142135623730947, 9.540047058232763*^-17}
numerical Rank: 14
numerical Dimension: 4

 *)