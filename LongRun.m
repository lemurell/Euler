(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[weight=2,weight<=2,weight+=2,
	For[Rleft=4,Rleft<=4,Rleft++,
		AppendTo[runList,{weight, Rleft}];
	];
];

For[Rleft=1,Rleft<=Length[runList],Rleft++,	
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
	Rtuple = {runList[[Rleft,2]]}; 
	weight = runList[[Rleft,1]];
	
	maxPower = 4;
		
	fileNameBase="Runs/HoloxMaass/NewHoloxMaass_Odd_" <> ToString[weight] <> "_" <> ToString[Rtuple[[1]]];

	<< Total.m
];
