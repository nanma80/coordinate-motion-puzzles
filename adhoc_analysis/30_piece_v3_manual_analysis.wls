#!/usr/bin/env wolframscript

(* 
Analyze why 12 different mode only share 4 dimensions

 *)


phi = (Sqrt[5] + 1)/2;
epsilon = 0.000000001;
randomDirection = {100.0, 3 Pi, Exp[1.0]}//N;

faces = {};
rotations = { IdentityMatrix[3] };
generators = { RotationMatrix[2 Pi/5, {phi, 1, 0}], RotationMatrix[Pi, {0, 0, 1}], RotationMatrix[Pi, {1, 0, 0}] };
getAngle[v1_, v2_] := ArcCos[v1.v2/Norm[v1]/Norm[v2]]/Pi*180.;

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
        (* Print[Length[rotations]]; *)
    ];
];

(* the piece part is rounded, the motion part is accurate *)
piecesAndMotion = {{{-2.618033989, -1., -1.6180339890000002}, {-1., 5.538511314930554*^-15, -0.2594638151136062}}, {{-3.2360679770000003, 0., 0.}, {-1., 5.538511314930994*^-15, -0.2594638151136062}}, {{-2.618033989, 1., -1.6180339890000002}, {-0.9999999999999961, 1.022076697164456*^-14, -0.25946381511360955}}, {{-2.618033989, -1., 1.6180339890000002}, {-0.7783914453408208, -0.358570173636287, 0.32071491318186607}}, {{-2.618033989, 1., 1.6180339890000002}, {-0.7783914453408202, 0.3585701736362913, 0.32071491318186235}}, {{-1.6180339890000002, 2.618033989, 1.}, {-0.4198212717045427, 0.9387489019317549, 0.09910635852267992}}, {{-1.6180339890000002, 2.618033989, -1.}, {-0.4198212717045411, 0.9387489019317545, 0.09910635852266715}}, {{0., 3.2360679770000003, 0.}, {-0.419821271704532, 0.9387489019317575, 0.09910635852267703}}, {{-1.6180339890000002, -2.618033989, 1.}, {-0.4198212717045304, -0.9387489019317512, 0.09910635852269456}}, {{-1.6180339890000002, -2.618033989, -1.}, {-0.41982127170452715, -0.9387489019317566, 0.09910635852268292}}, {{0., -3.2360679770000003, 0.}, {-0.4198212717045197, -0.9387489019317528, 0.09910635852268576}}, {{-1., -1.6180339890000002, 2.618033989}, {-0.1982127170453735, -0.5801787282954594, 0.6792850868181401}}, {{-1., 1.6180339890000002, 2.618033989}, {-0.19821271704535318, 0.5801787282954639, 0.6792850868181333}}, {{-1., 1.6180339890000002, -2.618033989}, {-0.06125109806825871, 0.3585701736362884, -0.8396425434090657}}, {{-1., -1.6180339890000002, -2.618033989}, {-0.06125109806825111, -0.3585701736362874, -0.8396425434090724}}, {{0., 0., 3.2360679770000003}, {0.1603574565909181, 3.727568533865629*^-15, 0.9008936414773132}}, {{0., 0., -3.2360679770000003}, {0.16035745659092526, -1.0958413484570056*^-14, -0.9766041623861783}}, {{1., -1.6180339890000002, -2.618033989}, {0.38196601125009777, -0.3585701736362819, -0.8396425434090732}}, {{1., 1.6180339890000002, -2.618033989}, {0.3819660112501064, 0.3585701736362872, -0.8396425434090621}}, {{1., -1.6180339890000002, 2.618033989}, {0.5189276302272103, -0.5801787282954702, 0.6792850868181313}}, {{1.6180339890000002, -2.618033989, -1.}, {0.5189276302272134, -0.5801787282954735, -0.4810723697727668}}, {{1., 1.6180339890000002, 2.618033989}, {0.5189276302272161, 0.5801787282954627, 0.6792850868181347}}, {{1.6180339890000002, -2.618033989, 1.}, {0.5189276302272178, -0.5801787282954685, 0.6792850868181297}}, {{1.6180339890000002, 2.618033989, 1.}, {0.5189276302272186, 0.5801787282954614, 0.6792850868181386}}, {{1.6180339890000002, 2.618033989, -1.}, {0.5189276302272189, 0.580178728295461, -0.48107236977277873}}, {{2.618033989, 1., 1.6180339890000002}, {0.518927630227219, 0.5801787282954596, 0.679285086818133}}, {{2.618033989, -1., 1.6180339890000002}, {0.5189276302272218, -0.5801787282954667, 0.6792850868181364}}, {{2.618033989, -1., -1.6180339890000002}, {0.7405361848863934, -0.22160855465917179, -0.618033988749885}}, {{2.618033989, 1., -1.6180339890000002}, {0.7405361848863973, 0.2216085546591722, -0.6180339887498878}}, {{3.2360679770000003, 0., 0.}, {0.8774978038635038, -2.7705166761654886*^-16, -0.2594638151136097}}};
piecesAndMotion = SortBy[piecesAndMotion, N[#[[1]].randomDirection]&];
flattenMotion = Flatten[Map[Last, piecesAndMotion]];

axis = {-1, 0, phi};

axes = {};
flatMotions = {};
For[rotationId = 1, rotationId <= Length[rotations], rotationId++,
    rotation = rotations[[rotationId]];
    newAxis = axis.rotation;
    AppendTo[axes, newAxis];
    newPiecesAndMotion = Map[#.rotation&, piecesAndMotion, {2}];
    newPiecesAndMotion = SortBy[newPiecesAndMotion, N[#[[1]].randomDirection]&];
    flatPieces = Flatten[Map[First, newPiecesAndMotion]];
    flatMotion = Flatten[Map[Last, newPiecesAndMotion]];
    AppendTo[flatMotions, flatMotion];
];


axesAndMotion = Transpose[ {axes, flatMotions} ];
axesAndMotion = SortBy[axesAndMotion,  N[#[[1]].randomDirection]& ];

(* Map[Print[N[#[[1]]] ]&, axesAndMotion]; *)
(* Map[Print[N[#[[2]]] ]&, axesAndMotion]; *)
dedupedMotions = {};
dedupedAxes = {};
For[axisId = 1, axisId <= Length[axesAndMotion], axisId++,
    If[ Mod[axisId, 5] == 1,
        AppendTo[dedupedMotions, axesAndMotion[[axisId]][[2]]];
        AppendTo[dedupedAxes,    axesAndMotion[[axisId]][[1]]];
    ];
];

sv = SingularValueList[N[dedupedMotions], Tolerance -> 0];
tol = 10^-10 * Max[sv];
numericalRank = Count[sv, _?(# > tol &)];

Print[sv];
Print["numerical Rank: ", numericalRank];
Print[];

dedupedAxes = dedupedAxes//N;
For[axisId = 1, axisId <= Length[dedupedAxes], axisId++,
    Print[IntegerString[axisId, 10, 2], " : ", dedupedAxes[[axisId]]];
];

selectedIndexes = {9, 10, 11, 12};
selectedMotions = dedupedMotions[[selectedIndexes]];
sv = SingularValueList[N[selectedMotions], Tolerance -> 0];
tol = 10^-10 * Max[sv];
numericalRank = Count[sv, _?(# > tol &)];

Print["selected motions: ", selectedIndexes];
Print[sv];
Print["numerical Rank: ", numericalRank];


(* antipodal relation: 1 + 12 == 2 + 11 *)
Print[ Norm[(dedupedMotions[[1]] + dedupedMotions[[12]]) - (dedupedMotions[[2]] + dedupedMotions[[11]]) ]];

Print[ Norm[(dedupedMotions[[1]] + dedupedMotions[[12]]) - (dedupedMotions[[5]] + dedupedMotions[[8]]) ]];


(* phi - scale relation: (8-7)/phi == (9-10); (8-10)/phi = (11-7) *)
Print[ Norm[(dedupedMotions[[8]] - dedupedMotions[[7]]) / phi - (dedupedMotions[[9]] - dedupedMotions[[10]])]];

Print[ Norm[(dedupedMotions[[6]] - dedupedMotions[[5]]) / phi - (dedupedMotions[[9]] - dedupedMotions[[10]])]];

Print[ Norm[(dedupedMotions[[8]] - dedupedMotions[[10]]) / phi - (dedupedMotions[[11]] - dedupedMotions[[7]])]];

Print[ Norm[(dedupedMotions[[10]] - dedupedMotions[[12]]) / phi - (dedupedMotions[[9]] - dedupedMotions[[6]])]];

(* Print["flatPieces"];
Print[flatPieces];
 *)
scalingFactor = 1/Abs[(dedupedMotions[[1]] + dedupedMotions[[12]])[[1]]]*phi/3;
(* Print[scalingFactor];
Print[(Total[dedupedMotions])*scalingFactor]; *)


avgMotion = Total[dedupedMotions]/Length[dedupedMotions];
avgAxis = Total[dedupedAxes]/Length[dedupedAxes];

Print[Norm[(Total[dedupedMotions])*scalingFactor - flatPieces ]];

For[pair1 = 1, pair1 <= Length[dedupedAxes] && pair1 <= 1, pair1++,
    For[pair2 = pair1 + 1, pair2 <= Length[dedupedAxes], pair2++,
        testPair = {pair1, pair2};

        ratioDistance = EuclideanDistance[dedupedMotions[[testPair[[1]]]], dedupedMotions[[testPair[[2]]]]] / EuclideanDistance[dedupedAxes[[ testPair[[1]]]], dedupedAxes[[testPair[[2]]]]]/0.8864342186366941;
        motionAngle = getAngle[dedupedMotions[[testPair[[1]]]], dedupedMotions[[testPair[[2]]]]];
        centralMotionAngle = getAngle[dedupedMotions[[testPair[[1]]]] - avgMotion, dedupedMotions[[testPair[[2]]]] - avgMotion];
        axisAngle = getAngle[ dedupedAxes[[ testPair[[1]]]], dedupedAxes[[testPair[[2]]]] ];
        Print[testPair, " ratioDistance: ", ratioDistance, " motionAngle: ", motionAngle, " centralMotionAngle: ", centralMotionAngle, " axisAngle: ", axisAngle];
    ];
];

(* 
01 : {-1.618033988749895, -1., 0.}
02 : {-1.618033988749895, 1., 0.}
03 : {-1., 0., -1.618033988749895}
04 : {-1., 0., 1.618033988749895}
05 : {-1.1102230246251565*^-16, -1.618033988749895, -0.9999999999999999}
06 : {-1.1102230246251565*^-16, -1.618033988749895, 0.9999999999999999}
07 : {-1.1102230246251565*^-16, 1.618033988749895, -0.9999999999999999}
08 : {-1.1102230246251565*^-16, 1.618033988749895, 0.9999999999999999}
09 : {0.9999999999999998, -1.6653345369377348*^-16, -1.6180339887498947}
10 : {0.9999999999999998, -1.6653345369377348*^-16, 1.6180339887498947}
11 : {1.6180339887498947, -1., 2.220446049250313*^-16}
12 : {1.6180339887498947, 1., 2.220446049250313*^-16}

Two relations contribute to reduction of dimensions:
- antipodal relation: For every antipodal pair of motions, the sum = uniform expansion. So, 1 + 12 == 2 + 11
- phi-scale relation: For two axes where the diff vectors are parallel but scaled by phi, the motion is also linearly dependent:
  - (8-7)/phi == (9-10); 
  - (8-10)/phi = (11-7)

Phi-scale relation implies antipodal relation.
(8-7)/phi = 9 - 10
(6-5)/phi = 9 - 10
=> 8-7 = 6-5
=> 8+5 = 6+7

Phi-scale relation suffices to prove the dimension <= 4
Start from 9, 10, 11, 12, full rank, dim = 4 (separate from the relation)
8 is linearly dependent because (10-11)/phi = 9 - 8
7 is linearly dependent because (8-10)/phi = 11 - 7
6 is linearly dependent because (10-12)/phi = 9 - 6
now we have antipodal points 6,7. All other points are linearly dependent.

Study how deduped motion vectors (90D) are embeded into a 3D space. Feasible motions: 4D region. If we cut the 4D region by a hyperplane Total[v]==1, we get a 3D region. The 3D region has 12 extreme points, forming an icosahedron. If we associate the axis (dedupedAxis) to the icosahedron region of the dedupedMotion, the relation is the isomorphism between an icosahedron and a great icosahedron

{1, 2} ratioDistance: 1.6180339887498907 motionAngle: 30.74113539383537 centralMotionAngle: 116.56505117707806 axisAngle: 63.43494882292201
{1, 3} ratioDistance: 1.6180339887498907 motionAngle: 30.74113539383537 centralMotionAngle: 116.56505117707765 axisAngle: 63.43494882292201
{1, 4} ratioDistance: 1.618033988749892 motionAngle: 30.74113539383539 centralMotionAngle: 116.56505117707752 axisAngle: 63.43494882292201
{1, 5} ratioDistance: 1.6180339887498962 motionAngle: 30.74113539383548 centralMotionAngle: 116.5650511770783 axisAngle: 63.43494882292201
{1, 6} ratioDistance: 1.6180339887498811 motionAngle: 30.741135393835194 centralMotionAngle: 116.56505117707697 axisAngle: 63.43494882292201
{1, 7} ratioDistance: 0.6180339887498997 motionAngle: 18.85699116878547 centralMotionAngle: 63.43494882292272 axisAngle: 116.56505117707799
{1, 8} ratioDistance: 0.6180339887498931 motionAngle: 18.856991168785257 centralMotionAngle: 63.43494882292211 axisAngle: 116.56505117707799
{1, 9} ratioDistance: 0.618033988749884 motionAngle: 18.85699116878496 centralMotionAngle: 63.43494882292114 axisAngle: 116.56505117707799
{1, 10} ratioDistance: 0.618033988749898 motionAngle: 18.85699116878543 centralMotionAngle: 63.43494882292269 axisAngle: 116.56505117707799
{1, 11} ratioDistance: 0.6180339887498899 motionAngle: 18.856991168785175 centralMotionAngle: 63.434948822921605 axisAngle: 116.56505117707799
{1, 12} ratioDistance: 1. motionAngle: 36.31108348286407 centralMotionAngle: 180. axisAngle: 179.99999879258172

*)
