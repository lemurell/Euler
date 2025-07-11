(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[Rleft=3,Rleft<=6,Rleft++,
	AppendTo[runList,{Rleft}];
];

Ltype = "MaassChar";
level = 5;
charvalue = Table[DirichletCharacter[level, 3, x],{x,level-1}];
OMEGA ={Cos[ALPHA], Sin[ALPHA]};
CHARSTEP = 0.1;
extraEq = {bb1[level]^2+bb2[level]==1};
extraEq = {};
realOrImaginary = 0;
coefLimit = getDegree[Ltype];
knownCoef = {};


sideLengths = getSideLengths[Ltype, 1];
RstepList = getRstep[Ltype, 1/10];
sideCounts = getSideCounts[sideLengths, RstepList];

testfunctiontypeSweep = "Classic";
testfunctiontypeZoom = "Classic";
initNN = 15;
minWidth = 38 / (realOrImaginary + 1);
minMult = 1/200;
bwidthStart = 3/2;

maxPower = 4;

For[Rleft=1,Rleft<=Length[runList],Rleft++,
	Rtuple = runList[[Rleft]]; 
	parity = {0, 0};
	fileNameBase="Runs/TestMaassLevel5AlphaChar_" <> ToString[Rtuple[[1]]];

	<< TotalChar.m
];



