#!/usr/bin/env wolframscript

numVertices[degTally_] := Total[degTally[[All, 2]]];

countDeg[degTally_, k_] :=
  Total @ Cases[degTally, {k, n_} :> n];

maxDeg[degTally_] :=
  Max[degTally[[All, 1]]];

vertexDegrees[p_] :=
  Sort @ Tally @ Counts @ Flatten @ PolyhedronData[p, "EdgeIndices"];

interestingPolyhedra = {#, vertexDegrees[#], numVertices[vertexDegrees[#]]} & /@
  Select[
    PolyhedronData["Convex"],
    Function[p,
      Module[{dt, n3, nOther},
        Quiet @ Check[
          dt = vertexDegrees[p];
          n3 = countDeg[dt, 3];
          nOther = numVertices[dt] - n3;
          n3 >= nOther &&
          numVertices[dt] <= 32 &&
          maxDeg[dt] > 3,
          False
        ]
      ]
    ]
  ];

Map[Print, interestingPolyhedra];

(* 

{AugmentedDodecahedron, {{3, 15}, {4, 5}, {5, 1}}, 21}
{AugmentedHexagonalPrism, {{3, 8}, {4, 5}}, 13}
{AugmentedPentagonalPrism, {{3, 6}, {4, 5}}, 11}
{AugmentedTruncatedCube, {{3, 16}, {4, 12}}, 28}
{BilinskiDodecahedron, {{3, 8}, {4, 6}}, 14}
{ElongatedDodecahedron, {{3, 16}, {4, 2}}, 18}
{ElongatedTriangularPyramid, {{3, 4}, {4, 3}}, 7}
{Gyrobifastigium, {{3, 4}, {4, 4}}, 8}
{PentagonalCupola, {{3, 10}, {4, 5}}, 15}
{PentagonalRotunda, {{3, 10}, {4, 10}}, 20}
{{Pyramid, 4}, {{3, 4}, {4, 1}}, 5}
{{Pyramid, 5}, {{3, 5}, {5, 1}}, 6}
{RhombicDodecahedron, {{3, 8}, {4, 6}}, 14}
{RhombicTriacontahedron, {{3, 20}, {5, 12}}, 32}
{SmallTriakisOctahedron, {{3, 8}, {8, 6}}, 14}
{SquareCupola, {{3, 8}, {4, 4}}, 12}
{TetragonalAntiwedge, {{3, 4}, {4, 2}}, 6}
{TrapezoRhombicDodecahedron, {{3, 8}, {4, 6}}, 14}
{TriakisIcosahedron, {{3, 20}, {10, 12}}, 32}
{TriakisTetrahedron, {{3, 4}, {6, 4}}, 8}
{TriangularCupola, {{3, 6}, {4, 3}}, 9}
{TridiminishedIcosahedron, {{3, 6}, {4, 3}}, 9}

 *)