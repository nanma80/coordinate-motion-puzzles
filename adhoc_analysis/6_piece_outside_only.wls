#!/usr/bin/env wolframscript

phi = (Sqrt[5] + 1)/2;

(* group 1: front/back over up/down *)
faceFrontUpOutside = {0, -1, 1};
faceFrontUpInside = {0, -1, -1};

faceBackUpOutside = {0, 1, 1};
faceBackUpInside = {0, 1, -1};

faceFrontDownOutside = {0, -1, -1};
faceFrontDownInside = {0, -1, 1};

faceBackDownOutside = {0, 1, -1};
faceBackDownInside = {0, 1, 1};

(* group 2: up/down over left/right *)
faceUpLeftOutside = {-1, 0, 1};
faceUpLeftInside = {1, 0, 1};

faceDownLeftOutside = {-1, 0, -1};
faceDownLeftInside = {1, 0, -1};

faceUpRightOutside = {1, 0, 1};
faceUpRightInside = {-1, 0, 1};

faceDownRightOutside = {1, 0, -1};
faceDownRightInside = {-1, 0, -1};

(* group 3: left/right over front/back*)
faceLeftFrontOutside = {-1, -1, 0};
faceLeftFrontInside = {-1, 1, 0};

faceRightFrontOutside = {1, -1, 0};
faceRightFrontInside = {1, 1, 0};

faceLeftBackOutside = {-1, 1, 0};
faceLeftBackInside = {-1, -1, 0};

faceRightBackOutside = {1, 1, 0};
faceRightBackInside = {1, -1, 0};

solutions = FindInstance[
    xFront == 1.0
    (* group 1 *)
    && ({xFront, yFront, zFront} - {xUp, yUp, zUp}).faceFrontUpOutside >= 0
    (* && ({xFront, yFront, zFront} - {xUp, yUp, zUp}).faceFrontUpInside >= 0 *)
    && ({xFront, yFront, zFront} - {xDown, yDown, zDown}).faceFrontDownOutside >= 0
    (* && ({xFront, yFront, zFront} - {xDown, yDown, zDown}).faceFrontDownInside >= 0 *)
    && ({xBack, yBack, zBack} - {xUp, yUp, zUp}).faceBackUpOutside >= 0
    (* && ({xBack, yBack, zBack} - {xUp, yUp, zUp}).faceBackUpInside >= 0 *)
    && ({xBack, yBack, zBack} - {xDown, yDown, zDown}).faceBackDownOutside >= 0
    (* && ({xBack, yBack, zBack} - {xDown, yDown, zDown}).faceBackDownInside >= 0 *)
    (* group 2 *)
    && ({xUp, yUp, zUp} - {xLeft, yLeft, zLeft}).faceUpLeftOutside >= 0
    (* && ({xUp, yUp, zUp} - {xLeft, yLeft, zLeft}).faceUpLeftInside >= 0 *)
    && ({xUp, yUp, zUp} - {xRight, yRight, zRight}).faceUpRightOutside >= 0
    (* && ({xUp, yUp, zUp} - {xRight, yRight, zRight}).faceUpRightInside >= 0 *)
    && ({xDown, yDown, zDown} - {xLeft, yLeft, zLeft}).faceDownLeftOutside >= 0
    (* && ({xDown, yDown, zDown} - {xLeft, yLeft, zLeft}).faceDownLeftInside >= 0 *)
    && ({xDown, yDown, zDown} - {xRight, yRight, zRight}).faceDownRightOutside >= 0
    (* && ({xDown, yDown, zDown} - {xRight, yRight, zRight}).faceDownRightInside >= 0 *)
    (* group 3 *)
    && ({xLeft, yLeft, zLeft} - {xFront, yFront, zFront}).faceLeftFrontOutside >= 0
    (* && ({xLeft, yLeft, zLeft} - {xFront, yFront, zFront}).faceLeftFrontInside >= 0 *)
    && ({xLeft, yLeft, zLeft} - {xBack, yBack, zBack}).faceLeftBackOutside >= 0
    (* && ({xLeft, yLeft, zLeft} - {xBack, yBack, zBack}).faceLeftBackInside >= 0 *)
    && ({xRight, yRight, zRight} - {xFront, yFront, zFront}).faceRightFrontOutside >= 0
    (* && ({xRight, yRight, zRight} - {xFront, yFront, zFront}).faceRightFrontInside >= 0 *)
    && ({xRight, yRight, zRight} - {xBack, yBack, zBack}).faceRightBackOutside >= 0
    (* && ({xRight, yRight, zRight} - {xBack, yBack, zBack}).faceRightBackInside >= 0 *)
    (* center of gravity *)
    && xFront + xBack + xUp + xDown + xLeft + xRight == 0
    && yFront + yBack + yUp + yDown + yLeft + yRight == 0
    && zFront + zBack + zUp + zDown + zLeft + zRight == 0
    ,
    {xFront, yFront, zFront, xBack, yBack, zBack, xUp, yUp, zUp, xDown, yDown, zDown, xLeft, yLeft, zLeft, xRight, yRight, zRight}, Reals
];

Print[solutions];

(* 

this is probably too much freedom

xFront -> 1., yFront -> 0., zFront -> 0., 
xBack -> -1., yBack -> 0., zBack -> 0., 
xUp -> 0., yUp -> 0., zUp -> 0., 
xDown -> 0., yDown -> 0., zDown -> 0., 
xLeft -> 0., yLeft -> 1., zLeft -> 0., 
xRight -> 0., yRight -> -1., zRight -> 0.


 *)

