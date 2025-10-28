#!/usr/bin/env wolframscript

(* 
(+x, +y, +z): Front, Right, Up

faceUFoverUR = {1, 1, 2};

First 2 vectors: two pieces. The first piece goes over (outside of) the second
Third vector: norm of the contacting surface of the two pieces, going from second to first piece

New model started from scratch.

constraint 1:
central surface
vzome: 
phi = (Sqrt[5] + 1)/2;
p1={0, 4+20/3 phi, 4+20/3 phi}
p2={0,0,6+10phi}
p3={-(6+10phi),6+10phi,6+10phi}

norm = Cross[p2-p1, p3-p1] ~ {-1,-1,-2}

constraint 2: pyramid, as same as before

constraint 3: central surface, small trapezoid
(-,-,+)
p1={0, 4+20/3 phi, 4+20/3 phi}
p2={3+5phi,3+5phi,6+10phi}
p3={6+10phi,6+10phi,6+10phi}*2/3 + {-(6+10phi),0,6+10phi}*1/3
norm = Cross[p2-p1, p3-p1]

 ~ {-1,-1,-2}

If we only keep the first 3 constraints, the dimension is the same: 2

 *)
faceTemplates = {
    {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}} (* central surface, trapezoid *)
    ,{{1, 0, 1}, {0, 1, 1}, {1, -1, 0}} (* pyramid *)
    ,{{1, 0, 1}, {0, 1, 1}, {-1, -1, 1}} (* central small surface, trapezoid, diff from v2, split *)
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

gluedPairs = {
    {1, 2}
};

conditions = {
    (* pieceCoordinates[[12]][[1]] == 1.0, *)
    (* pieceCoordinates[[1]][[3]] == -1.0, *)
    (* pieceCoordinates[[1]][[2]] == -1.0, *)
    pieceCoordinates[[1]][[1]] == -1.0
};

(* uniform explosion *)
If[False,
    For[pid = 1, pid <= Length[pieces], pid++,
        AppendTo[conditions, pieceCoordinates[[pid]][[1]] == pieces[[pid]][[1]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[2]] == pieces[[pid]][[2]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[3]] == pieces[[pid]][[3]] ];
    ];
];

(* glued pairs move together *)
For[pid = 1, pid <= Length[gluedPairs], pid++,
    gluedPair = gluedPairs[[pid]];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];

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

Uniform expansion: works

Also, 4 groups of 3:
{{-1, -1, 0}, {-1., -1., -1.}}
{{-1, 0, -1}, {-1., -1., -1.}}
{{0, -1, -1}, {-1., -1., -1.}}
{{-1, 0, 1}, {-1., 1., 1.}}
{{-1, 1, 0}, {-1., 1., 1.}}
{{0, 1, 1}, {-1., 1., 1.}}
{{0, -1, 1}, {1., -1., 1.}}
{{1, -1, 0}, {1., -1., 1.}}
{{1, 0, 1}, {1., -1., 1.}}
{{0, 1, -1}, {1., 1., -1.}}
{{1, 0, -1}, {1., 1., -1.}}
{{1, 1, 0}, {1., 1., -1.}}
Groups
{3, 3, 3, 3}
4


 *)