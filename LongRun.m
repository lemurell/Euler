(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m


runList = {};
For[weight=2,weight<=2,weight+=2,
	For[Rleft=4,Rleft<=4,Rleft++,
		AppendTo[runList,{weight, Rleft}];
	];
];

Ltype = "HoloxMaass";
level = 1;
charvalue = "No";
realOrImaginary = 1;
coefLimit = getDegree[Ltype];
knownCoef = {};

sideLengths = getSideLengths[Ltype, 1];
RstepList = getRstep[Ltype, 1/20];
sideCounts = getSideCounts[sideLengths, RstepList];
truncationList = {7, 10, 12, 14, 16, 18};	
maxPower = 4;

For[Rleft=1,Rleft<=Length[runList],Rleft++,	
	Rtuple = {runList[[Rleft,2]]}; 
	weight = runList[[Rleft,1]];

	parity = {1, weight};
	OMEGA = I^weight *(-1)^parity[[1]];
	
	fileNameBase="Runs/HoloxMaass/NewHoloxMaass_Odd_" <> ToString[weight] <> "_" <> ToString[Rtuple[[1]]];

	<< Total.m
];
