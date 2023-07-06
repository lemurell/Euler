(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[a1 = 6, a1<=6, a1++,
	For[a2 = 2,a2<=2, a2++,
		AppendTo[runList,{a1, a2}];
	];
];

Ltype = "SP4";
level = 1;
charvalue = "No";
realOrImaginary = 1;
parity = {0, 0};
coefLimit = getDegree[Ltype];
knownCoef = {};

sideLengths = getSideLengths[Ltype, 1];
RstepList = getRstep[Ltype, 1/10];
sideCounts = getSideCounts[sideLengths, RstepList];

testfunctiontypeSweep = "DS";
testfunctiontypeZoom = "DS";
initNN = 20;
minWidth = 40 / (realOrImaginary + 1);
minMult = 1/200;
bwidthStart = 3/2;

maxPower = 4;

For[Rleft=1,Rleft<=Length[runList],Rleft++,	
	Rtuple = runList[[Rleft]]; 
	OMEGA = (-1)^(parity[[1]] + parity[[2]]);
	
	fileNameBase="Runs/SP4_" <> ToString[parity[[1]]] <> "_" <> ToString[parity[[2]]]<> "_" <> ToString[Rtuple[[1]]] <> "_"<> ToString[Rtuple[[2]]] ;

	<< Total.m
];
