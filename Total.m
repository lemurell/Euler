(* ::Package:: *)

Clear[result];
st=TimeUsed[];
CurrentStep = 1;

(* Load files. *)
<< Lfunctions.m;

(* Variables to initiate before calling
Ltype = "HoloxMaass";
level = 1;
charvalue = "No";
parity = {1, weight};
OMEGA = I^weight *(-1)^parity[[1]];
realOrImaginary = 1;
coefLimit = getDegree[Ltype];
knownCoef = {};

RstepList = getRstep[Ltype, 1/20];
sideCounts = getSideCounts[sideLengths, RstepList];
Rtuple = {12.90884204, 0.5739285998}

maxPower = 4;

fileNameBase = "Runs/NL/Test";
*)

(* Global precision parameters. *)
PRECISIONMULTIPLE = 4;
TRUNCDIGITS = 7;
DIGITS=TRUNCDIGITS + 6;
PRECISION= Max[PRECISIONMULTIPLE * TRUNCDIGITS, MachinePrecision];
MYPRECISION=PRECISION;

(* Data for L-functions *)
 Rlist = {};
 For[i=-1,i<sideCounts[[1]],i++,
    For[j=-1,j<sideCounts[[2]],j++,
	    For[k=-1,k<sideCounts[[3]],k++,
			AppendTo[Rlist,getRtuple[Ltype, Rtuple, RstepList, {i,j,k}]];
		];
    ];
 ];

(* Global L-function parameters. *)
Ldata = {Ltype, level, charvalue, OMEGA, parity};
Q = QValue[Ltype, level];
klist = getKlist[Ltype];
reslist = {};
polelist = {};
llist = getLlist[Ltype, Rlist[[-1]], parity];

(* Precomputations. *)
logQ = N[Log[Q],PRECISION + 2];
logQlogN = Table[logQ-Log[n],{n,1,1000}];

(* Local testfunction parameters. *)
v=2;
nrOfExtraEquations = 8;
Nlist = getNumberOfCoefficients[OMEGA,klist,llist,parity,v];
NN = Nlist[[TRUNCDIGITS]];
{plist, highplist, nonplist, knownlist}=getUnknowns[Ldata, NN, maxPower, knownCoef];
nrOfUnknowns= (2 - Abs[realOrImaginary]) Length[Union[plist, highplist]];
nrOfEquations = nrOfUnknowns + nrOfExtraEquations;
{slist, paralist,Param} = getSlistDavid[nrOfEquations];

NLmethod = "Secant";
nrOfRuns = 50;
MaxSolutions = 20;
startValuePrec = 0;

incr=2*Pi*v/Log[10]/DIGITS;  (* Not good if v is large *)
M = 10*Sqrt[Log[10]*DIGITS];
Print[N[incr]," M=",N[M]];

expz=Table[Exp[(v+I*k)*logQlogN[[n]]]/(v+I*k),{n,1,NN},{k,-M,M,incr}];

Print["Number of unknowns: ",  nrOfUnknowns];
Print["Unknowns: ", Union[plist, highplist]];
Print["Number of extra equations: ",  nrOfExtraEquations];
Print["Number of R-values: ", Length[Rlist]];
Print["First R-value: ",N[Rlist[[1]],8]];
Print["Filnamn: ",fileNameBase];

Save[fileNameBase ,{DIGITS, PRECISION, TRUNCDIGITS, Ldata, Rlist, g, Param, slist, paralist, NN, minWidth, minWidthList, initNNList, minMultList, maxPower, realOrImaginary}];

RowLength = (sideCounts[[2]] + 1)(sideCounts[[3]] + 1);
If[FileType[ fileNameBase <> "StateData.m"]==File,
	{CurrentStep,lastRowStartValues,RloopStart} = Get[fileNameBase <> "StateData.m"];
,
	Print["Shouldn't happen!"];
	Break;
,
	lastRowStartValues = Table[{},{k,RowLength}];
	RloopStart = 1;
];

For[Rloop = RloopStart,Rloop<=Length[Rlist],Rloop++,
    Rtuple = Rlist[[Rloop]];
    llist = getLlist[Ltype, Rtuple, parity];

	firstIndex = Mod[Rloop,RowLength];
	If[firstIndex == 0, 
		firstIndex = RowLength;
	];

	startValues = lastRowStartValues[[firstIndex]];

	If[sideCounts[[2]]>0 && sideCounts[[3]]==0,
		If[firstIndex > 1,
			For[k=1,k<=Length[lastRowStartValues[[firstIndex-1]]],k++,
				AppendTo[startValues, lastRowStartValues[[firstIndex-1,k]]];
			];
		];
	,
		If[sideCounts[[3]]>0,
			If[Mod[Rloop,sideCounts[[3]]+1] != 1,
				For[k=1,k<=Length[lastRowStartValues[[firstIndex-1]]],k++,
					AppendTo[startValues, lastRowStartValues[[firstIndex-1,k]]];
				];
			];
			If[firstIndex > sideCounts[[3]]+1, 
				For[k=1,k<=Length[lastRowStartValues[[firstIndex-sideCounts[[3]]-1]]],k++,
					AppendTo[startValues, lastRowStartValues[[firstIndex-sideCounts[[3]]-1,k]]];
				];
			];
		];
	];

    result[Rtuple]=solveForOneNL[Ldata, klist, llist, reslist, polelist, slist, paralist, Param, nrOfRuns, NN, M, incr, v, expz, startValues, knownCoef, maxPower, NLmethod, MaxSolutions, startValuePrec, realOrImaginary];
	
    If[Rloop==1,
        Save[fileNameBase <> "Matrix", fullmatrix];
    ];
    
	lastRowStartValues[[firstIndex]]  = result[Rtuple][[1]];

    Save[fileNameBase,result];
    Clear[result];
	Put[{CurrentStep,lastRowStartValues,Rloop+1},fileNameBase <> "StateData.m"];
];
DeleteFile[fileNameBase <> "StateData.m"];

TimeStep1 = (TimeUsed[]-st)/60;
Save[fileNameBase <> "Time.txt",TimeStep1];



<< Step2.m