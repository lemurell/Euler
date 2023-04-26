(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[weight=12,weight<=12,weight+=2,
	For[Rleft=95/10,Rleft<=95/10,Rleft++,
		AppendTo[runList,{weight, Rleft}];
	];
];

Ltype = "HoloxMaass";
level = 1;
charvalue = "No";
realOrImaginary = 1;
OddOrEven = 1;
coefLimit = getDegree[Ltype];
knownCoef = {};

sideLengths = getSideLengths[Ltype, 1];
RstepList = getRstep[Ltype, 1/20];
sideCounts = getSideCounts[sideLengths, RstepList];
testfunctiontype = "Classic";
If[testfunctiontype == "Classic",
	initNN = 15;
	minWidth = 40 / (realOrImaginary + 1);
	minMult = 1/200;
,
	truncationList = {5, 7, 10, 12, 14, 16, 18};
];	
maxPower = 4;

For[Rleft=1,Rleft<=Length[runList],Rleft++,	
	Rtuple = {runList[[Rleft,2]]}; 
	weight = runList[[Rleft,1]];

	parity = {OddOrEven, weight};
	OMEGA = I^weight *(-1)^parity[[1]];
	
	fileNameBase="Runs/TestNewVersion";

	<< Total.m
];
