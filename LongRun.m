(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[weight=8,weight<=8,weight+=2,
	For[Rleft=4,Rleft<=4,Rleft++,
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
RstepList = getRstep[Ltype, 1/20];
sideCounts = getSideCounts[sideLengths, RstepList];

testfunctiontypeSweep = "Classic";
testfunctiontypeZoom = "DS";
initNN = 20;
minWidth = 40 / (realOrImaginary + 1);
minMult = 1/200;
maxPower = 4;

For[Rleft=1,Rleft<=Length[runList],Rleft++,	
	Rtuple = {runList[[Rleft,2]]}; 
	weight = runList[[Rleft,1]];

	parity = {OddOrEven, weight};
	OMEGA = I^weight *(-1)^parity[[1]];
	
	fileNameBase="Runs/TestNewVersion";

	<< Total.m
];
