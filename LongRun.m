(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[weight=32,weight<=40,weight+=2,
	For[Rleft=0,Rleft<=1,Rleft++,
		AppendTo[runList,{weight, Rleft}];
	];
];

Ltype = "HoloxMaass";
level = 1;
charvalue = "No";
realOrImaginary = 1;
OddOrEven = 0;
coefLimit = getDegree[Ltype];
knownCoef = {};

sideLengths = getSideLengths[Ltype, 1];
RstepList = getRstep[Ltype, 1/100];
sideCounts = getSideCounts[sideLengths, RstepList];

testfunctiontypeSweep = "DS";
testfunctiontypeZoom = "DS";
initNN = 20;
minWidth = 40 / (realOrImaginary + 1);
minMult = 1/200;
bwidthStart = 3/2;

maxPower = 4;

For[Rleft=1,Rleft<=Length[runList],Rleft++,	
	Rtuple = {runList[[Rleft,2]]}; 
	weight = runList[[Rleft,1]];

	parity = {OddOrEven, weight};
	OMEGA = I^weight *(-1)^parity[[1]];
	
	fileNameBase="Runs/HoloxMaass/DS32HoloxMaass_Even_" <> ToString[weight] <> "_" <> ToString[Rtuple[[1]]];

	<< Total.m
];
