#!/usr/bin/env wolframscript

(* term: U/D -> F/B -> L/R *)
(* template *)
faceUFoverUR = {1, 1, 2};
(* 2-fold rotation by UF, U<->F, L<->R, B<->D *)
faceUFoverFL = {2, 1, 1};
(* 4-fold rotation by U, U,D remain, F->R->B->L->F *)
faceURoverUB = {-1, 1, 2};
faceURoverFR = {-1, 2, 1};
faceUBoverUL = {-1, -1, 2};
faceUBoverBR = {-2, -1, 1};
faceULoverUF = {1, -1, 2};
faceULoverBL = {1, -2, 1};

(* 4-fold rotation by F, F,B remain, U->L->D->R->U *)
faceFLoverUL = {1, 1, 2};
faceFLoverDF = {2, 1, 1};
faceULoverBL = {-1, -2, 1};
faceULoverUF = {-1, -1, 2};
faceBLoverDL = {-1, -2, -1};
faceBLoverUB = {-2, -1, -1};
faceDLoverFL = {1, -2, -1};
faceDLoverDB = {1, -1, -2};

faceDFoverDL = {1, -2, 1};
faceDFoverFR = {2, -1, 1};
faceDLoverDB = {-1, -1, -2};
faceDLoverLF = {-1, -2, -1};
faceDBoverDR = {-1, 1, -2};
faceDBoverBL = {-2, 1, -1};
faceDRoverDF = {1, 1, -2};
faceDRoverRB = {1, 2, -1};

faceRFoverRD = {1, -1, -2};
faceRFoverFU = {2, -1, -1};
faceRDoverRB = {-1, 2, -1};
faceRDoverDF = {-1, 1, -2};
faceRBoverRU = {-1, 2, 1};
faceRBoverBD = {-2, 1, 1};
faceRUoverRF = {1, 2, 1};
faceRUoverUB = {1, 1, 2};



solutions = FindInstance[
    xFront == 1.0
    (* group 1 *)
    && ({xFront, yFront, zFront} - {xUp, yUp, zUp}).faceFrontUpOutside >= 0
    && ({xFront, yFront, zFront} - {xUp, yUp, zUp}).faceFrontUpInside >= 0
    && ({xFront, yFront, zFront} - {xDown, yDown, zDown}).faceFrontDownOutside >= 0
    && ({xFront, yFront, zFront} - {xDown, yDown, zDown}).faceFrontDownInside >= 0
    && ({xBack, yBack, zBack} - {xUp, yUp, zUp}).faceBackUpOutside >= 0
    && ({xBack, yBack, zBack} - {xUp, yUp, zUp}).faceBackUpInside >= 0
    && ({xBack, yBack, zBack} - {xDown, yDown, zDown}).faceBackDownOutside >= 0
    && ({xBack, yBack, zBack} - {xDown, yDown, zDown}).faceBackDownInside >= 0
    (* group 2 *)
    && ({xUp, yUp, zUp} - {xLeft, yLeft, zLeft}).faceUpLeftOutside >= 0
    && ({xUp, yUp, zUp} - {xLeft, yLeft, zLeft}).faceUpLeftInside >= 0
    && ({xUp, yUp, zUp} - {xRight, yRight, zRight}).faceUpRightOutside >= 0
    && ({xUp, yUp, zUp} - {xRight, yRight, zRight}).faceUpRightInside >= 0
    && ({xDown, yDown, zDown} - {xLeft, yLeft, zLeft}).faceDownLeftOutside >= 0
    && ({xDown, yDown, zDown} - {xLeft, yLeft, zLeft}).faceDownLeftInside >= 0
    && ({xDown, yDown, zDown} - {xRight, yRight, zRight}).faceDownRightOutside >= 0
    && ({xDown, yDown, zDown} - {xRight, yRight, zRight}).faceDownRightInside >= 0
    (* group 3 *)
    && ({xLeft, yLeft, zLeft} - {xFront, yFront, zFront}).faceLeftFrontOutside >= 0
    && ({xLeft, yLeft, zLeft} - {xFront, yFront, zFront}).faceLeftFrontInside >= 0
    && ({xLeft, yLeft, zLeft} - {xBack, yBack, zBack}).faceLeftBackOutside >= 0
    && ({xLeft, yLeft, zLeft} - {xBack, yBack, zBack}).faceLeftBackInside >= 0
    && ({xRight, yRight, zRight} - {xFront, yFront, zFront}).faceRightFrontOutside >= 0
    && ({xRight, yRight, zRight} - {xFront, yFront, zFront}).faceRightFrontInside >= 0
    && ({xRight, yRight, zRight} - {xBack, yBack, zBack}).faceRightBackOutside >= 0
    && ({xRight, yRight, zRight} - {xBack, yBack, zBack}).faceRightBackInside >= 0
    (* center of gravity *)
    && xFront + xBack + xUp + xDown + xLeft + xRight == 0
    && yFront + yBack + yUp + yDown + yLeft + yRight == 0
    && zFront + zBack + zUp + zDown + zLeft + zRight == 0
    ,
    {xFront, yFront, zFront, xBack, yBack, zBack, xUp, yUp, zUp, xDown, yDown, zDown, xLeft, yLeft, zLeft, xRight, yRight, zRight}, Reals
];

Print[solutions];

(* 
Found the solution to split in two halves
xFront -> 1., yFront -> -1.0000000000000038, zFront -> 1.0000000000000004, 
xUp -> 1.0000000000000009, yUp -> -1.000000000000003, zUp -> 1.0000000000000002, 
xRight -> 1.0000000000000016, yRight -> -0.9999999999999999, zRight -> 1.0000000000000016

xBack -> -1.0000000000000149, yBack -> 1.000000000000009, zBack -> -0.9999999999999978, 
xDown -> -1.0000000000000047, yDown -> 0.9999999999999988, zDown -> -1.0000000000000038, 
xLeft -> -0.9999999999999998, yLeft -> 1.0000000000000004, zLeft -> -1.0000000000000038, 


 *)

