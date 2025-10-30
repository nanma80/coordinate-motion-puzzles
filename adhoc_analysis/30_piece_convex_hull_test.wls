#!/usr/bin/env wolframscript

(* 
(+x, +y, +z): Front, Right, Up

icosahedral vertices:
{phi, 1, 0}, {1, 0, phi}, {0, phi, 1}

edges:
Edge between {phi, 1, 0} and {phi, -1, 0}: {2 phi, 0, 0}
Edge between {phi, 1, 0} and {1, 0, phi}: {phi + 1, 1, phi}
sum of two edges:
{3 phi + 1, 1, phi}
Another edge over: {phi + 1, -1, -phi}

Rhombus in front: in plane: x = phi
+/- {phi, 1, 0}, +/- {phi, 0, phi - 1}

pyramid surface determined by:
{phi, 1, 0}, {phi, 0, phi - 1}, {0, 0, 0}

Cross[{phi, 0, phi - 1}, {phi, 1, 0}]

For v3 let's go with cuts with both blue struts
 *)


phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;
faceTemplates = {
    {{phi + 1, 1, phi}, {2 phi, 0, 0}, {3 phi + 1, 1, phi} } (* central surface *)
    , {{phi + 1, 1, phi}, {2 phi, 0, 0}, Cross[{phi, 0, phi - 1}, {phi, 1, 0}] } (* pyramid *)
    (* , {{phi + 1, 1, phi}, {2 phi, 0, 0}, {0, 1, 0} } (* old split. not used in v3*) *)
    , {{phi + 1, 1, phi}, {2 phi, 0, 0}, {0, phi, -1} } (* split, computed by cross product, as same as the line above *)
    , {{phi + 1, 1, phi}, {2 phi, 0, 0}, {1+phi, -1-2phi, 2+phi} } (* extruded surface, is this needed? Adding it doesn't change dimensionality *)
}//N;

faces = {};
rotations = { IdentityMatrix[3] };
generators = { RotationMatrix[2 Pi/5, {phi, 1, 0}], RotationMatrix[Pi, {0, 0, 1}], RotationMatrix[Pi, {1, 0, 0}] };

For[round = 1, round <= 3, round++,
    For[gid = 1, gid <= Length[generators], gid ++,
        generator = generators[[gid]];
        newRotations = {};
        For[rid = 1, rid <= Length[rotations], rid ++,
            rotation = rotations[[rid]];
            newRotation = generator.rotation;
            AppendTo[newRotations, newRotation];
        ];
        rotations = Union[rotations, newRotations, SameTest -> (Max[Abs[N[#1] - N[#2]]] < epsilon &)];
        Print[Length[rotations]];
    ];
];

For[tid = 1, tid <= Length[faceTemplates], tid ++,
    faceTemplate = faceTemplates[[tid]];
    For[rid = 1, rid <= Length[rotations], rid ++,
        rotation = rotations[[rid]];
        newFace = Map[#.rotation&, faceTemplate];
        newFace = FullSimplify[newFace];
        AppendTo[faces, newFace];
    ];    
];

pieces = Map[First, faces];
roundedPieces = Round[pieces, epsilon];
roundedPieces = Sort[Union[roundedPieces]];
uniquePieces = {};
For[pid = 1, pid <= Length[roundedPieces], pid++,
    uniquePiece = First[Select[pieces, EuclideanDistance[#, roundedPieces[[pid]]] < 10 * epsilon & ]];
    AppendTo[uniquePieces, uniquePiece];
];

Print["Rounded piece count"];
Print[Length[roundedPieces]];

Print["piece count"];
Print[Length[pieces]];

Print["unique piece count"];
Print[Length[uniquePieces]];

(* faces = Union[faces]; *)

Print["face (constraint) count"];
Print[Length[faces]];

(* Exit[]; *)

pieceLookup = <||>;

Print["Rounded pieces"];
(* Print[pieces]; *)
For[pid = 1, pid <= Length[roundedPieces], pid ++,
    piece = roundedPieces[[pid]];
    Print[{pid, piece, piece.{0, 0, 1}} ];
    pieceLookup[piece] = pid;
];

(* Exit[]; *)

(* Print[pieceLookup]; *)

pieceCoordinates = Array[x, {Length[roundedPieces], 3}];

pieceCoordinatesFlat = Flatten[pieceCoordinates];

gluedPairs = {
};

allMotions = {{-4 - Sqrt[5], 0, (1 + Sqrt[5])/2, (-3*(1 + Sqrt[5]))/2, -Sqrt[5], -2, -4 - Sqrt[5], 0, (1 + Sqrt[5])/2, (-3*(1 + Sqrt[5]))/2, Sqrt[5], -2, -4 - Sqrt[5], 0, (1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, 1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], (-3 + Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], 1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], (-3 + Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, 1, 0, (-9 - Sqrt[5])/2, 1, 0, (1 + 5*Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], (7 - Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], (7 - Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, 3, 1 + Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, 3, 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], (7 + Sqrt[5])/2, (-5 + Sqrt[5])/2, (1 + 3*Sqrt[5])/2, 1 + Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], (7 + Sqrt[5])/2, (5 - Sqrt[5])/2, (1 + 3*Sqrt[5])/2, 1 + 2*Sqrt[5], 0, (1 + Sqrt[5])/2}, {-4 - Sqrt[5], 0, (-1 - Sqrt[5])/2, -4 - Sqrt[5], 0, (-1 - Sqrt[5])/2, (-3*(1 + Sqrt[5]))/2, -Sqrt[5], 2, -4 - Sqrt[5], 0, (-1 - Sqrt[5])/2, (-3*(1 + Sqrt[5]))/2, Sqrt[5], 2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], 1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], (-3 + Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], 1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, 1, 0, (-1 - 5*Sqrt[5])/2, 1, 0, (9 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (7 - Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], (7 - Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, -3, 1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, -3, 1 + Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], (7 + Sqrt[5])/2, (-5 + Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, 1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], (7 + Sqrt[5])/2, (5 - Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, 1 + Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 + 2*Sqrt[5], 0, (-1 - Sqrt[5])/2}, {-1 - 2*Sqrt[5], 0, (-1 - Sqrt[5])/2, (-7 - Sqrt[5])/2, (-5 + Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], (-7 - Sqrt[5])/2, (5 - Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, -1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, -3, -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, -3, -1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], (-7 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], (-7 + Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, -1, 0, (-1 - 5*Sqrt[5])/2, -1, 0, (9 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 - Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], -1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], (3 - Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], -1 + Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, 4 + Sqrt[5], 0, (-1 - Sqrt[5])/2, (3*(1 + Sqrt[5]))/2, -Sqrt[5], 2, 4 + Sqrt[5], 0, (-1 - Sqrt[5])/2, (3*(1 + Sqrt[5]))/2, Sqrt[5], 2, 4 + Sqrt[5], 0, (-1 - Sqrt[5])/2}, {-1 - 2*Sqrt[5], 0, (1 + Sqrt[5])/2, -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], (-7 - Sqrt[5])/2, (-5 + Sqrt[5])/2, (1 + 3*Sqrt[5])/2, -1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], (-7 - Sqrt[5])/2, (5 - Sqrt[5])/2, (1 + 3*Sqrt[5])/2, -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, 3, -1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, 3, -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], (-7 + Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], (-7 + Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, -1, 0, (-9 - Sqrt[5])/2, -1, 0, (1 + 5*Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, -1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], (3 - Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], -1 + Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], (3 - Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3*(1 + Sqrt[5]))/2, -Sqrt[5], -2, 4 + Sqrt[5], 0, (1 + Sqrt[5])/2, (3*(1 + Sqrt[5]))/2, Sqrt[5], -2, 4 + Sqrt[5], 0, (1 + Sqrt[5])/2, 4 + Sqrt[5], 0, (1 + Sqrt[5])/2}, {(-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - Sqrt[5])/2, 3, 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], -Sqrt[5], 3 + Sqrt[5], (-3 + Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], (7 - Sqrt[5])/2, -Sqrt[5], -2, (-3*(1 + Sqrt[5]))/2, (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], 0, (1 + Sqrt[5])/2, -4 - Sqrt[5], (-5 + Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (7 + Sqrt[5])/2, 0, (-9 - Sqrt[5])/2, 1, 0, (1 + Sqrt[5])/2, -4 - Sqrt[5], 0, (1 + Sqrt[5])/2, 1 + 2*Sqrt[5], 0, (1 + 5*Sqrt[5])/2, 1, Sqrt[5], -2, (-3*(1 + Sqrt[5]))/2, (5 + Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], 0, (1 + Sqrt[5])/2, -4 - Sqrt[5], (5 - Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (7 + Sqrt[5])/2, (5 + Sqrt[5])/2, -2 - Sqrt[5], 1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], Sqrt[5], 3 + Sqrt[5], (-3 + Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (7 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + Sqrt[5])/2, 3, 1 + Sqrt[5], (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2}, {(-9 - Sqrt[5])/2, -1, 0, -2 - Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 + Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, -2, (3*(1 + Sqrt[5]))/2, -Sqrt[5], -2, (3*(1 + Sqrt[5]))/2, Sqrt[5], -2 - Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 + Sqrt[5])/2, -1 - 2*Sqrt[5], 0, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 + Sqrt[5])/2, 4 + Sqrt[5], 0, 3, -1 - Sqrt[5], (-5 - Sqrt[5])/2, 3, -1 - Sqrt[5], (5 + Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, (-5 + Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, (5 - Sqrt[5])/2, (1 + Sqrt[5])/2, 4 + Sqrt[5], 0, (1 + Sqrt[5])/2, 4 + Sqrt[5], 0, 3 + Sqrt[5], (-7 + Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], (-7 + Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (3 - Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], 2/(3 + Sqrt[5]), Sqrt[5], (1 + 5*Sqrt[5])/2, -1, 0}, {(-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - Sqrt[5])/2, 3, -1 - Sqrt[5], (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 + Sqrt[5], -Sqrt[5], 3 + Sqrt[5], (-7 + Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], (3 - Sqrt[5])/2, (-5 - Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], -Sqrt[5], -2, (3*(1 + Sqrt[5]))/2, (-5 + Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, 0, (1 + Sqrt[5])/2, 4 + Sqrt[5], 0, (-9 - Sqrt[5])/2, -1, 0, (1 + Sqrt[5])/2, -1 - 2*Sqrt[5], 0, (1 + Sqrt[5])/2, 4 + Sqrt[5], 0, (1 + 5*Sqrt[5])/2, -1, (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], Sqrt[5], -2, (3*(1 + Sqrt[5]))/2, (5 - Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, 0, (1 + Sqrt[5])/2, 4 + Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 + Sqrt[5], Sqrt[5], 3 + Sqrt[5], (-7 + Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (3 - Sqrt[5])/2, (5 + Sqrt[5])/2, -2 - Sqrt[5], -1 - Sqrt[5], (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + Sqrt[5])/2, 3, -1 - Sqrt[5], (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (3 + Sqrt[5])/2}, {(-9 - Sqrt[5])/2, 1, 0, -2 - Sqrt[5], 1 - Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 - Sqrt[5], (5 + Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, -2, (-3*(1 + Sqrt[5]))/2, -Sqrt[5], -2, (-3*(1 + Sqrt[5]))/2, Sqrt[5], -2 - Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, -2 - Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, (1 + Sqrt[5])/2, -4 - Sqrt[5], 0, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (1 + Sqrt[5])/2, 1 + 2*Sqrt[5], 0, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, 3, 1 + Sqrt[5], (-5 - Sqrt[5])/2, 3, 1 + Sqrt[5], (5 + Sqrt[5])/2, (1 + Sqrt[5])/2, -4 - Sqrt[5], 0, (1 + Sqrt[5])/2, -4 - Sqrt[5], 0, (1 + 3*Sqrt[5])/2, (7 + Sqrt[5])/2, (-5 + Sqrt[5])/2, (1 + 3*Sqrt[5])/2, (7 + Sqrt[5])/2, (5 - Sqrt[5])/2, 3 + Sqrt[5], (-3 + Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], (-3 + Sqrt[5])/2, Sqrt[5], 3 + Sqrt[5], (7 - Sqrt[5])/2, -Sqrt[5], 3 + Sqrt[5], (7 - Sqrt[5])/2, Sqrt[5], (1 + 5*Sqrt[5])/2, 1, 0}, {(-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - Sqrt[5])/2, -3, 1 + Sqrt[5], (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], -Sqrt[5], -3 - Sqrt[5], (-3 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], (7 - Sqrt[5])/2, (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], 0, (-1 - Sqrt[5])/2, -4 - Sqrt[5], (-5 + Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, (7 + Sqrt[5])/2, -Sqrt[5], 2, (-3*(1 + Sqrt[5]))/2, (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], 0, (-1 - 5*Sqrt[5])/2, 1, 0, (-1 - Sqrt[5])/2, -4 - Sqrt[5], 0, (-1 - Sqrt[5])/2, 1 + 2*Sqrt[5], 0, (9 + Sqrt[5])/2, 1, 0, (-1 - Sqrt[5])/2, -4 - Sqrt[5], (5 - Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, (7 + Sqrt[5])/2, Sqrt[5], 2, (-3*(1 + Sqrt[5]))/2, (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], Sqrt[5], -3 - Sqrt[5], (-3 + Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], (7 - Sqrt[5])/2, (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + Sqrt[5])/2, -3, 1 + Sqrt[5], (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2}, {(-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - Sqrt[5])/2, -3, -1 - Sqrt[5], (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], (-7 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], (3 - Sqrt[5])/2, (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 + Sqrt[5], (-5 + Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, 0, (-1 - Sqrt[5])/2, 4 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], -Sqrt[5], 2, (3*(1 + Sqrt[5]))/2, 0, (-1 - 5*Sqrt[5])/2, -1, 0, (-1 - Sqrt[5])/2, -1 - 2*Sqrt[5], 0, (-1 - Sqrt[5])/2, 4 + Sqrt[5], 0, (9 + Sqrt[5])/2, -1, (5 - Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, 0, (-1 - Sqrt[5])/2, 4 + Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], Sqrt[5], 2, (3*(1 + Sqrt[5]))/2, Sqrt[5], -3 - Sqrt[5], (-7 + Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], (3 - Sqrt[5])/2, (5 + Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], -1 + Sqrt[5], (5 + Sqrt[5])/2, -3, -1 - Sqrt[5], (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2}, {(-1 - 5*Sqrt[5])/2, -1, 0, -3 - Sqrt[5], (-7 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], (-7 + Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], (3 - Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], 2/(3 + Sqrt[5]), Sqrt[5], (-1 - 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, (-5 + Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, (-7 - Sqrt[5])/2, (5 - Sqrt[5])/2, (-1 - Sqrt[5])/2, 4 + Sqrt[5], 0, (-1 - Sqrt[5])/2, 4 + Sqrt[5], 0, -3, -1 - Sqrt[5], (-5 - Sqrt[5])/2, -3, -1 - Sqrt[5], (5 + Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 - Sqrt[5])/2, -1 - 2*Sqrt[5], 0, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 - Sqrt[5])/2, 4 + Sqrt[5], 0, 2 + Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (5 + 3*Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, 2, (3*(1 + Sqrt[5]))/2, -Sqrt[5], 2, (3*(1 + Sqrt[5]))/2, Sqrt[5], 2 + Sqrt[5], -1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], -1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], -1 + Sqrt[5], (5 + Sqrt[5])/2, (9 + Sqrt[5])/2, -1, 0}, {(-1 - 5*Sqrt[5])/2, 1, 0, -3 - Sqrt[5], (-3 + Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], (-3 + Sqrt[5])/2, Sqrt[5], -3 - Sqrt[5], (7 - Sqrt[5])/2, -Sqrt[5], -3 - Sqrt[5], (7 - Sqrt[5])/2, Sqrt[5], (-1 - Sqrt[5])/2, -4 - Sqrt[5], 0, (-1 - Sqrt[5])/2, -4 - Sqrt[5], 0, (-1 - 3*Sqrt[5])/2, (7 + Sqrt[5])/2, (-5 + Sqrt[5])/2, (-1 - 3*Sqrt[5])/2, (7 + Sqrt[5])/2, (5 - Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, -3, 1 + Sqrt[5], (-5 - Sqrt[5])/2, -3, 1 + Sqrt[5], (5 + Sqrt[5])/2, (-1 - Sqrt[5])/2, -4 - Sqrt[5], 0, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, (-1 - Sqrt[5])/2, 1 + 2*Sqrt[5], 0, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (-5 - 3*Sqrt[5])/2, (-1 + Sqrt[5])/2, (-3 - Sqrt[5])/2, (5 + 3*Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, 2, (-3*(1 + Sqrt[5]))/2, -Sqrt[5], 2, (-3*(1 + Sqrt[5]))/2, Sqrt[5], 2 + Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 - Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 - Sqrt[5], (5 + Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (-5 - Sqrt[5])/2, 2 + Sqrt[5], 1 + Sqrt[5], (5 + Sqrt[5])/2, (9 + Sqrt[5])/2, 1, 0}};

conditions = {
    (* pieceCoordinates[[30]][[1]] == 1.0, *)
    (* pieceCoordinates[[1]][[3]] == -1.0, *)
    (* pieceCoordinates[[1]][[2]] == -1.0/phi, *)
    (* pieceCoordinates[[1]][[1]] == -1.0 *)
};

(* Exit[]; *)
(* uniform explosion *)
If[False,
    For[pid = 1, pid <= Length[uniquePieces], pid++,
        AppendTo[conditions, pieceCoordinates[[pid]][[1]] == uniquePieces[[pid]][[1]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[2]] == uniquePieces[[pid]][[2]] ];
        AppendTo[conditions, pieceCoordinates[[pid]][[3]] == uniquePieces[[pid]][[3]] ];
    ];
];

For[pid = 1, pid <= Length[gluedPairs], pid++,
    gluedPair = gluedPairs[[pid]];
    For[dim = 1, dim <= 3, dim++,
        AppendTo[conditions, pieceCoordinates[[gluedPair[[1]]]][[dim]] == pieceCoordinates[[gluedPair[[2]]]][[dim]] ];
    ];
];

For[fid = 1, fid <= Length[faces], fid++,
    face = faces[[fid]];
    piece1Id = pieceLookup[Round[face[[1]], epsilon]];
    piece2Id = pieceLookup[Round[face[[2]], epsilon]];
    condition = ( pieceCoordinates[[piece1Id]] - pieceCoordinates[[piece2Id]] ).face[[3]] >= - epsilon * 0.00001 ;
    AppendTo[conditions, condition];
];

For[dim = 1, dim <= 3, dim++,
    AppendTo[conditions, Total[Table[pieceCoordinates[[k]][[dim]], {k, 1, Length[roundedPieces]}]] == 0];
];

weightVariables = Array[weight, Length[allMotions] ];
allVariables = Join[pieceCoordinatesFlat, weightVariables];

For[mid = 1, mid <= Length[allMotions[[1]] ], mid++,
    combinedMotionComponent = Total[Table[N[allMotions[[k]][[mid]]] * weightVariables[[k]], {k, 1, Length[allMotions] }]];
    AppendTo[conditions, pieceCoordinatesFlat[[mid]] - combinedMotionComponent >= -epsilon * 0.00001 ];
    AppendTo[conditions, -(pieceCoordinatesFlat[[mid]] - combinedMotionComponent) >= -epsilon * 0.00001 ];
];

AppendTo[conditions, Total[weightVariables] == 1.0 ];

For[wid = 1, wid <= Length[weightVariables], wid++,
    If[wid != 1 && wid != 2 && wid != 12 && wid != 3,
        AppendTo[conditions, weightVariables[[wid]] >= 0 ];
        AppendTo[conditions, weightVariables[[wid]] <= epsilon * 0.00001 ],
        If[wid == 1 || wid == 2 || wid == 12,
            AppendTo[conditions, weightVariables[[wid]] >= 0 ],
            (* wid == 3, opposite *)
            (* AppendTo[conditions, weightVariables[[wid]] >= 0 ]; *)
            AppendTo[conditions, weightVariables[[wid]] < - epsilon ];
        ];
    ];
];

combinedConditions = Apply[ And, conditions];

Print["Starting FindInstance"];
solutions = FindInstance[
    combinedConditions,
    allVariables,
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

    weightsValues = weightVariables/.solution;
    Print["weights:"];
    Print[weightsValues];

    Print["flat motion"];
    Print[pieceCoordinatesFlat/.solution];

    (* values = Round[values, 1]; *)
    Print["piece and motion"];
    pieceAndMotion = {};
    For[pid = 1, pid <= Length[roundedPieces], pid++,
        AppendTo[pieceAndMotion, { roundedPieces[[pid]], values[[pid]]}];
    ];
    pieceAndMotion = SortBy[pieceAndMotion, Last];
    Map[Print, pieceAndMotion];
    (* Print[pieceAndMotion]; *)

    motionTally = Tally[Map[Last, pieceAndMotion]];
    tallyCounts = Map[#[[2]]&, motionTally];
    Print["Groups"];
    Print[ tallyCounts ];
    Print[ tallyCounts//Length ];
];

(* 

weights:
{-1.0000000000000124, 0.9999999999999898, 0., 0.9999999999999662, 1.0000000000000002*^-14, 0., 1.0000000000000002*^-14, 6.120122866394084*^-15, 0., 1.0000000000000002*^-14, 1.0000000000000002*^-14, 1.0000000000000002*^-14}

flat motion
{-5.472135954999579, 2.9023437232074247*^-15, -1.618033988749936, -4.618033988749927, -1.3819660112500642, -3.854101966249656, -3.236067977499782, -3.6180339887499136, 4.236067977499776, -4.618033988749928, 1.38196601125007, -3.854101966249654, -3.236067977499779, 3.618033988749912, 4.2360679774997765, -3.236067977499786, -3.6180339887499136, -2.9999999999999507, -3.2360679774997867, -3.6180339887499113, 4.236067977499776, -3.236067977499783, 3.618033988749912, -2.999999999999947, -3.236067977499784, 3.6180339887499113, 4.236067977499774, -2.3819660112500687, -2.2360679774997707, -5.236067977499801, -3.236067977499787, -3.6180339887499082, 4.236067977499777, -2.3819660112500727, 2.236067977499774, -5.236067977499791, -3.2360679774997845, 3.6180339887499104, 4.236067977499775, 2.6180339887498483, -5.854101966249683, 0.6180339887499084, -0.9999999999999767, -1.0743221089275384*^-15, -6.0901699437494745, -0.9999999999999859, -8.354043586480537*^-16, 5.618033988749886, 2.618033988749842, 5.854101966249678, 0.6180339887499057, 0.38196601125011487, -2.2360679774997756, -5.2360679774998, 1.2360679774998091, -3.6180339887499073, 4.236067977499777, 0.38196601125012436, 2.2360679774997827, -5.236067977499793, 1.236067977499813, 3.6180339887499025, 4.236067977499772, 2.6180339887498962, -5.85410196624967, 0.6180339887498737, 2.6180339887498962, -5.854101966249668, 0.6180339887499371, 2.6180339887499033, 5.854101966249657, 0.6180339887498663, 2.618033988749902, 5.854101966249661, 0.6180339887499294, 6.236067977499777, 1.4600842229138146*^-15, -1.618033988749912, 4.854101966249661, -2.236067977499815, 2.0000000000000195, 6.236067977499779, 1.925039820564095*^-15, -1.6180339887499153, 4.85410196624966, 2.2360679774998236, 2.000000000000017, 6.236067977499776, 1.3172040745721043*^-15, -1.618033988749917}


allMotions: 1, 2, 12 are on the same small facet. 3, 4, 6 on the opposite facet

Table[ {k, ArcCos[ allMotions[[2]].allMotions[[k]]/Norm[allMotions[[1]]]^2]/Pi*180//N}, {k, 1,12}]

face1212 = (allMotions[[1]]+allMotions[[2]]+allMotions[[12]])/3.;

Table[ {k, ArcCos[ allMotions[[k]].face1212 /Norm[allMotions[[1]]]/Norm[face1212] ] /Pi*180//N}, {k, 1,12}]

if we require positive weights on 1, 2, 12, and negative weight on 3, no solution. So, no solution lies beyond the face of 1,2,12, proving that we are not missing any solution beyond the solution of five triplets.


 *)