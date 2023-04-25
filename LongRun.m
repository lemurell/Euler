(* ::Package:: *)

(* Load files. *)
SetDirectory["Z:\\math\\EulerProj"];
<< Lfunctions3.15.m


runList = {};
For[weight=2,weight<=2,weight+=2,
For[Rleft=4,Rleft<=4,Rleft++,
AppendTo[runList,{weight, Rleft}];
];
];

For[Rleft=1,Rleft<=Length[runList],Rleft++,
Rtuple = {runList[[Rleft,2]]}; 
weight = runList[[Rleft,1]];
fileNameBase="Runs/HoloxMaass/NewHoloxMaass_Odd_" <> ToString[weight] <> "_" <> ToString[Rtuple[[1]]];

Ltype = "HoloxMaass";
level = 1;
charvalue = "No";
parity = {1, weight};
OMEGA = I^weight *(-1)^parity[[1]];
realOrImaginary = 1;
coefLimit = getDegree[Ltype];
knownCoef = {};

sideLengths = getSideLengths[Ltype, 1];
RstepList = getRstep[Ltype, 1/20];
sideCounts = getSideCounts[sideLengths, RstepList];

minWidth = 48 / (realOrImaginary + 1);
minMult = 1/200;
maxPower = 4;
initNN = 15;

(* For step 2 *)
 minWidthList = {48} / (realOrImaginary + 1);
 initNNList = {25,30,35,50};
 minMultList = {1/200};



<< eu3.15Total.m
];


NN


sSeq


paraSeq


plist


highplist


slist


paralist//Length


nrOfEquations
