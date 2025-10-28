#!/usr/bin/env wolframscript

(* 
(+x, +y, +z): Front, Right, Up

faceUFoverUR = {1, 1, 2};

First 2 vectors: two pieces. The first piece goes over (outside of) the second
Third vector: norm of the contacting surface of the two pieces, going from second to first piece
 *)
faceTemplates = {
    {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}},
    {{1, 0, 1}, {1, 0, -1}, {0, -1, 1}}
};
(* faceTemplate = {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}}; *)
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

gluedPairs = {
    {3, 6},
    {8, 11},
    {11, 12},
    {1, 5},
    {1, 2}
};

conditions = {
    (* pieceCoordinates[[12]][[1]] == 1.0, *)
    pieceCoordinates[[1]][[3]] == -1.0,
    pieceCoordinates[[1]][[2]] == -1.0,
    pieceCoordinates[[1]][[1]] == -1.0
};

For[pid = 1, pid <= Length[gluedPairs], pid++,
    gluedPair = gluedPairs[[pid]];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];

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

combinedConditions = Apply[ And, conditions];

solutions = FindInstance[
    combinedConditions,
    pieceCoordinatesFlat,
    Reals
];

(* Print[solutions]; *)
Print["Found these many solutions:"]
Print[solutions//Length];

For[sid = 1, sid <= Length[solutions], sid++,
    Print[sid];
    solution = solutions[[sid]];
    values = pieceCoordinates/.solution;
    values = Round[values, 0.000001];
    Print[values];
    (* values = Round[values, 1]; *)
    Print["piece and motion"];
    pieceAndMotion = {};
    For[pid = 1, pid <= Length[pieces], pid++,
        AppendTo[pieceAndMotion, { pieces[[pid]], values[[pid]]}];
    ];
    pieceAndMotion = SortBy[pieceAndMotion, Last];
    Map[Print, pieceAndMotion];

    motionTally = Tally[Map[Last, pieceAndMotion]];
    tallyCounts = Map[#[[2]]&, motionTally];
    Print["Groups"];
    Print[ tallyCounts ];
    Print[ tallyCounts//Length ];
];


(* 

Grouping 1,2,3,5,6

piece and motion
{{-1, -1, 0}, {-3, -3, 0}}
{{-1, 0, -1}, {-3, -3, 0}}
{{-1, 0, 1}, {-3, -3, 0}}
{{0, -1, -1}, {-3, -3, 0}}
{{0, -1, 1}, {-3, -3, 0}}

{{-1, 1, 0}, {-1, 1, 0}}

{{1, -1, 0}, {1, -1, 0}}

{{0, 1, -1}, {3, 3, 0}}
{{0, 1, 1}, {3, 3, 0}}
{{1, 0, -1}, {3, 3, 0}}
{{1, 0, 1}, {3, 3, 0}}
{{1, 1, 0}, {3, 3, 0}}
Groups
{5, 1, 1, 5}
4


Grouping 1,2,3,4
piece and motion
{{-1, -1, 0}, {-3, 0, 0}}
{{-1, 0, -1}, {-3, 0, 0}}
{{-1, 0, 1}, {-3, 0, 0}}
{{-1, 1, 0}, {-3, 0, 0}}
{{0, -1, -1}, {0, -1, -1}}
{{0, -1, 1}, {0, -1, 1}}
{{0, 1, -1}, {0, 1, -1}}
{{0, 1, 1}, {0, 1, 1}}
{{1, -1, 0}, {3, 0, 0}}
{{1, 0, -1}, {3, 0, 0}}
{{1, 0, 1}, {3, 0, 0}}
{{1, 1, 0}, {3, 0, 0}}
Groups
{4, 1, 1, 1, 1, 4}
6

Grouping 1, 2, 5 (3-fold face), {8, 11, 12}, and {3, 4}
piece and motion
{{-1, -1, 0}, {-1., -1., -1.}}
{{-1, 0, -1}, {-1., -1., -1.}}
{{0, -1, -1}, {-1., -1., -1.}}
{{-1, 0, 1}, {-1., 0.411765, 0.411765}}
{{-1, 1, 0}, {-1., 0.411765, 0.411765}}
{{0, -1, 1}, {0.411765, -1., 0.411765}}
{{1, -1, 0}, {0.411765, -1., 0.411765}}
{{0, 1, -1}, {0.411765, 0.411765, -1.}}
{{1, 0, -1}, {0.411765, 0.411765, -1.}}
{{0, 1, 1}, {1.117647, 1.117647, 1.117647}}
{{1, 0, 1}, {1.117647, 1.117647, 1.117647}}
{{1, 1, 0}, {1.117647, 1.117647, 1.117647}}
Groups
{3, 2, 2, 2, 3}
5

Grouping 1, 2, 5 (3-fold face), {8, 11, 12}, and {3, 6}
piece and motion
{{-1, -1, 0}, {-1., -1., -1.}}
{{-1, 0, -1}, {-1., -1., -1.}}
{{0, -1, -1}, {-1., -1., -1.}}
{{-1, 0, 1}, {-0.368421, -0.368421, 0.894737}}
{{0, -1, 1}, {-0.368421, -0.368421, 0.894737}}
{{-1, 1, 0}, {-0.368421, 0.894737, -0.368421}}
{{0, 1, -1}, {-0.368421, 0.894737, -0.368421}}
{{1, -1, 0}, {0.894737, -0.368421, -0.368421}}
{{1, 0, -1}, {0.894737, -0.368421, -0.368421}}
{{0, 1, 1}, {0.894737, 0.894737, 0.894737}}
{{1, 0, 1}, {0.894737, 0.894737, 0.894737}}
{{1, 1, 0}, {0.894737, 0.894737, 0.894737}}
Groups
{3, 2, 2, 2, 3}

Grouping 1,2,5, and 8,11,12, and forcing these two groups to go opposite in same speed
piece and motion
{{-1, -1, 0}, {-1., -1., -1.}}
{{-1, 0, -1}, {-1., -1., -1.}}
{{0, -1, -1}, {-1., -1., -1.}}
{{-1, 0, 1}, {-0.666667, 0., 0.666667}}
{{-1, 1, 0}, {-0.666667, 0.666667, 0.}}
{{0, -1, 1}, {0., -0.666667, 0.666667}}
{{0, 1, -1}, {0., 0.666667, -0.666667}}
{{1, -1, 0}, {0.666667, -0.666667, 0.}}
{{1, 0, -1}, {0.666667, 0., -0.666667}}
{{0, 1, 1}, {1., 1., 1.}}
{{1, 0, 1}, {1., 1., 1.}}
{{1, 1, 0}, {1., 1., 1.}}
Groups
{3, 1, 1, 1, 1, 1, 1, 3}

 *)