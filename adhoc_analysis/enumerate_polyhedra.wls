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
          (* n3 >= nOther && *)
          (* numVertices[dt] > 32 && *)
          n3 < nOther &&
          n3 >= 0.5 nOther &&
          maxDeg[dt] > 3,
          False
        ]
      ]
    ]
  ];

Map[Print, interestingPolyhedra];

(* 
numVertices <= 32
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

numVertices > 32
{AugmentedTruncatedDodecahedron, {{3, 50}, {4, 15}}, 65}
{Conext21Polyhedron, {{3, 56}, {4, 24}}, 80}
{MetabiaugmentedTruncatedDodecahedron, {{3, 40}, {4, 30}}, 70}
{ParabiaugmentedTruncatedDodecahedron, {{3, 40}, {4, 30}}, 70}
{PentagonalHexecontahedron, {{3, 80}, {5, 12}}, 92}
{PentagonalIcositetrahedron, {{3, 32}, {4, 6}}, 38}
{RhombicEnneacontahedron, {{3, 60}, {5, 12}, {6, 20}}, 92}
{TridiminishedRhombicosidodecahedron, {{3, 30}, {4, 15}}, 45}


n3 < nOther && n3 >= 0.5 nOther:

{AugmentedTridiminishedIcosahedron, {{3, 4}, {4, 6}}, 10}
{AugmentedTruncatedTetrahedron, {{3, 6}, {4, 9}}, 15}
{{Dipyramid, 3}, {{3, 2}, {4, 3}}, 5}
{ElongatedPentagonalCupola, {{3, 10}, {4, 15}}, 25}
{ElongatedPentagonalPyramid, {{3, 5}, {4, 5}, {5, 1}}, 11}
{ElongatedPentagonalRotunda, {{3, 10}, {4, 20}}, 30}
{ElongatedSquareCupola, {{3, 8}, {4, 12}}, 20}
{ElongatedSquarePyramid, {{3, 4}, {4, 5}}, 9}
{ElongatedTriangularCupola, {{3, 6}, {4, 9}}, 15}
{GyrateBidiminishedRhombicosidodecahedron, {{3, 20}, {4, 30}}, 50}
{MetabiaugmentedDodecahedron, {{3, 10}, {4, 10}, {5, 2}}, 22}
{MetabidiminishedRhombicosidodecahedron, {{3, 20}, {4, 30}}, 50}
{ParabiaugmentedDodecahedron, {{3, 10}, {4, 10}, {5, 2}}, 22}
{ParabidiminishedRhombicosidodecahedron, {{3, 20}, {4, 30}}, 50}
{RhombicIcosahedron, {{3, 10}, {4, 10}, {5, 2}}, 22}
{TriaugmentedTruncatedDodecahedron, {{3, 30}, {4, 45}}, 75}


 *)