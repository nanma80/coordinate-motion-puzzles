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
== {-0.618034, 1., 1.61803}
== {-phi+1, 1, phi}

v2: On each pentagonal face, extrude a bit orthogonally from the face.

vzome after scaling down a few times.
icosa vertices:
{1+2phi, 1+phi, 0}

norm vector of the new constraint:
{1+3phi, 0, 1}-{2phi, 1+2phi, -1-phi}
={1+phi, -1-2phi, 2+phi}
This is proportional to:
Cross[ {7+12phi, 2+3phi, -3-5phi} - {6+10phi, 4+6phi, 0}, {8+64/5 phi, 24/5+8 phi, 0} - {6+10phi, 4+6phi, 0}]

Pattern from 12 piece v3:
The third constraint = Cross[ pyramid, central surface ]

intuition:
pyr x central = split
extended(extrude) x central = central x split

split = {0, 6 + 2 Sqrt[5], -2 - 2 Sqrt[5]} ~ {0, phi, -1} == {0., 1.61803, -1.}

vzome file is phi^2 scaled up.

For v3 let's go with cuts with both blue struts
 *)



phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;
faceTemplates = {
    {{phi + 1, 1, phi}, {2 phi, 0, 0}, {3 phi + 1, 1, phi} }  (* central surface *)
    , {{phi + 1, 1, phi}, {2 phi, 0, 0}, Cross[{phi, 0, phi - 1}, {phi, 1, 0}] } (* pyramid *)
    (* , {{phi + 1, 1, phi}, {2 phi, 0, 0}, Cross[Cross[{phi, 0, phi - 1}, {phi, 1, 0}], {3 phi + 1, 1, phi}] } (* split, computed by cross product *) *)
    , {{phi + 1, 1, phi}, {2 phi, 0, 0}, {0, phi, -1} } (* split, computed by cross product, as same as the line above *)
    , {{phi + 1, 1, phi}, {2 phi, 0, 0}, {1+phi, -1-2phi, 2+phi} } (* extruded surface, is this needed? Adding it doesn't change dimensionality *)
    (* , {{0, 2 phi, 0}, {2 phi, 0, 0}, Cross[{5 + 8 phi, 6 + 9 phi, -3 - 5 phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi}, {31/5 + 48/5 phi, 31/5 + 48/5 phi, -31/5 - 48/5 phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi}] }  (* end of wings *) *)
    (* , {{0, 2 phi, 0}, {2 phi, 0, 0}, -Cross[{5+8 phi, 6 + 9phi, -3-5 phi}-{57/10 + 46/5 phi, 33/5+21/2 phi, -7/2-57/10 phi}, {28/5+44/5 phi, 4 + 8 phi, -16/5-28/5 phi} -{57/10 + 46/5 phi, 33/5+21/2 phi, -7/2-57/10 phi} ] }  (* end of wings 1/2, blue *) *)

    (* , {{0, 2 phi, 0}, {2 phi, 0, 0}, Cross[{28/5+44/5 phi, 4 + 8 phi, -16/5-28/5 phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi},  {31/5 + 48/5 phi, 31/5 + 48/5 phi, -31/5 - 48/5 phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi} ] }  (* end of wings 2/2, orange *) *)
    , {{0, 2 phi, 0}, {2 phi, 0, 0}, -Cross[{5+8 phi, 6 + 9phi, -3-5 phi}-{57/10 + 46/5 phi, 33/5+21/2 phi, -7/2-57/10 phi}, {28/5+44/5 phi, 4 + 8 phi, -16/5-28/5 phi} -{57/10 + 46/5 phi, 33/5+21/2 phi, -7/2-57/10 phi} ] }  (* end of wings all blue style, 1/2 *)

    , {{0, 2 phi, 0}, {2 phi, 0, 0}, Cross[{6+8 phi, 4+6 phi, -2-6 phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi},  {31/5 + 48/5 phi, 31/5 + 48/5 phi, -31/5 - 48/5 phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi} ] }  (* end of wings all blue style, 2/2 *)

    (* , {{0, 2 phi, 0}, {2 phi, 0, 0}, Cross[{32/5 +52/5 phi, 36/5+12 phi, -4-32/5 phi } - {5+8 phi, 6+9 phi, -3 -5 phi}, {24/5+42/5phi, 26/5+38/5phi,-18/5-24/5phi} - {5+8 phi, 6+9 phi, -3 -5 phi} ] }  (* end of wings deep cut, 1/3 inner *) *)

    (* , {{0, 2 phi, 0}, {2 phi, 0, 0}, Cross[{44/5+72/5phi, 16/5+28/5phi, -28/5-44/5phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi},  {32/5+56/5phi,32/5+56/5phi,-32/5-56/5phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi} ] }  (* end of wings deep cut, 2/3 outer *) *)

    (* , {{0, 2 phi, 0}, {2 phi, 0, 0}, Cross[{ 24/5+42/5phi, 26/5+38/5phi,-18/5-24/5phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi},{44/5+72/5phi, 16/5+28/5phi, -28/5-44/5phi} - {4 + 8 phi, 4 + 8 phi, -4 - 8 phi} ] }  (* end of wings deep cut, 3/3 buttom *) *)

    (* , {{phi + 1, 1, phi}, {2 phi, 0, 0}, {0, 1, 0} } (* old central split *) *)
    (* , {{phi, phi+1, 1}, {phi + 1, -1, -phi}, {0, 1, 0} } (* smallest surface at the end *) *)
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
    Print[{pid, piece, piece.{1,1,1}} ];
    pieceLookup[piece] = pid;
];

(* Exit[]; *)

(* Print[pieceLookup]; *)

pieceCoordinates = Array[x, {Length[roundedPieces], 3}];

pieceCoordinatesFlat = Flatten[pieceCoordinates];

conditions = {
};

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
Print["rows"];
Print[rows];
Print[];
Print["rows row count: ", Length[rows]];
Print["rows col count: ", Length[rows[[1]] ]];

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

(* 

rows row count: 123
rows col count: 90
rows rank: 86
coneDimension: 4
{28.51372206018814, 28.513722060188126, 28.513722060188126, 27.95403204309489, 27.954032043094877, 27.95403204309487, 27.954032043094863, 23.707842330143624, 23.70784233014358, 23.707842330143553, 23.707842330143553, 22.67159038809749, 22.67159038809748, 22.671590388097478, 22.671590388097474, 22.671590388097457, 22.392626217580315, 22.39262621758029, 22.39262621758029, 21.01276138494568, 21.01276138494568, 21.01276138494567, 18.30672416732552, 18.306724167325513, 18.3067241673255, 18.3067241673255, 18.306724167325495, 15.962178834491652, 15.962178834491636, 15.962178834491633, 15.928452503904007, 15.928452503904007, 15.928452503903992, 14.417596710647079, 14.417596710647071, 14.417596710647071, 14.41759671064707, 14.417596710647066, 14.330481127008762, 14.33048112700876, 14.33048112700876, 14.330481127008753, 14.330481127008753, 13.908880114639864, 13.908880114639862, 13.908880114639858, 13.908880114639857, 13.908880114639851, 13.45579365369639, 13.455793653696386, 13.455793653696384, 13.220545949703467, 13.220545949703464, 13.220545949703462, 13.220545949703457, 12.084437368129127, 12.08443736812912, 12.084437368129118, 12.084437368129116, 9.60321112460192, 9.603211124601918, 9.603211124601918, 9.603211124601916, 9.603211124601915, 9.139428273658298, 9.139428273658293, 9.139428273658291, 9.139428273658286, 8.45470522932122, 8.454705229321213, 8.45470522932121, 6.387577774962743, 6.387577774962742, 6.3875777749627405, 6.38757777496274, 6.387577774962738, 6.2461931317654384, 6.246193131765438, 6.246193131765438, 5.687060849489612, 5.687060849489612, 5.68706084948961, 5.687060849489608, 5.477225575051661, 5.477225575051661, 5.477225575051659, 9.25651436198871*^-15, 6.027482380995763*^-15, 5.397319876117236*^-15, 2.3866043456715638*^-15}
numerical Rank: 86
numerical Dimension: 4
 *)
