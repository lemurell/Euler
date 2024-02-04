(* ::Package:: *)

Clear[result];
st=TimeUsed[];

(* Load files. *)
<<Lfunctions.m;

(* Variables to initiate before calling
fileNameBase = "Runs/08/TestVer2.92/SL3_NL" <> ThisFileName ;
maxPower = 4;
Ltype = "GL3";
*)

{CurrentStep,CurrentIndex2,CurrentIndex1,CurrentIndex3} = Get[fileNameBase <> "StateData.m"];

(* Global precision parameters. *)
PRECISIONMULTIPLE = 5;
TRUNCDIGITS=7;
 
nrOfExtraEquations = 8;



Get[fileNameBase];

Switch[getNrOfParameters[Ltype],
	1,   (* ----------- One parameter ----------------------------- *)
Rstep=-Rlist[[1,1]]+Rlist[[2,1]];
margin=1/5;
noConvergence={};
allTogList =  {};
For[i=1,i<Length[Rlist],i++,
		If[Length[result[Rlist[[i]]][[1]]] *Length[result[Rlist[[i+1]]][[1]]] ==0 ,
			AppendTo[noConvergence,N[{Rlist[[i]]}]];
		,
			For[k=1,k<=Length[result[Rlist[[i]]][[1]]],k++,
				{Rapprox, meanR, coefApprox, coef, RMargin, coefMargin}=approximationWithErrorMargin[Rlist[[{i,i+1}]], result, k];

				If[closeProportion[Rapprox, Rlist[[i]], {Rstep *(1+margin)} ] >1/2 && closeProportion[{meanR},Rlist[[i]], {Rstep *(1+margin)} ]==1,
					coef = removePrimeCoefRubbish[coef, coefLimit, realOrImaginary];
					If[Length[coef]>0,
					AppendTo[allTogList,{meanR,coef}];
					];
				];
			];
		];
];
	,2,   (* ----------- Two parameters ----------------------------- *)
row=1;
While[Rlist[[row+1,2]]!=Rlist[[1,2]],
	row++;
];
col=Length[Rlist]/row;
Rstep1=-Rlist[[1,1]]+Rlist[[1+row,1]];
Rstep2=-Rlist[[1,2]]+Rlist[[2,2]];
Rstep = Min[Rstep1, Rstep2];
margin=1/5;
noConvergence={};
allTogList =  {};
For[i=1,i<col,i++,
	For[j=1,j<row,j++,
		ind=j+(i-1)*row;
		If[Length[result[Rlist[[ind]]][[1]]] *Length[result[Rlist[[ind+1]]][[1]]]*Length[result[Rlist[[ind+row]]][[1]]] ==0 ,
			AppendTo[noConvergence,N[{Rlist[[ind]]}]];
		,
			For[k=1,k<=Length[result[Rlist[[ind]]][[1]]],k++,
				{Rapprox, meanR, coefApprox, coef, RMargin, coefMargin}=approximationWithErrorMargin[Rlist[[{ind,ind+1,ind+row}]], result, k];


				If[closeProportion[Rapprox, Rlist[[ind]], {Rstep1 *(1+margin),Rstep2 *(1+margin)} ] >1/2 && closeProportion[{meanR},Rlist[[ind]], {Rstep1 *(1+margin),Rstep2 *(1+margin)} ]==1,
					coef = removePrimeCoefRubbish[coef, coefLimit, realOrImaginary];
					If[Length[coef]>0,
					AppendTo[allTogList,{meanR,coef}];
					Print[{meanR,Rlist[[ind]],Norm[meanR - Rlist[[ind]]]}];
					];
				];
			];
		];
	];
];
	,3,   (* ----------- Three parameters ----------------------------- *)
row=1;
While[Rlist[[row+1,3]]!=Rlist[[1,3]],
	row++;
];
col=2;
While[Rlist[[1+col*row,2]]!=Rlist[[1,2]],
	col++;
];
third=Length[Rlist]/row/col;
RstepListLocal={-Rlist[[1,1]]+Rlist[[1+row*col,1]],-Rlist[[1,2]]+Rlist[[1+row,2]],-Rlist[[1,3]]+Rlist[[2,3]]};
Rstep = Min[RstepListLocal];
margin=1/5;
noConvergence={};
allTogList =  {};
For[i=1,i<third,i++,
	For[j=1,j<col,j++,
		For[m=1,m<row,m++,
			ind=m + (j-1)* row +(i-1)*row*col;
			If[Length[result[Rlist[[ind]]][[1]]] *Length[result[Rlist[[ind+1]]][[1]]]*Length[result[Rlist[[ind+row]]][[1]]]*Length[result[Rlist[[ind+row*col]]][[1]]] ==0 ,
				AppendTo[noConvergence,N[{Rlist[[ind]]}]];
			,
				For[k=1,k<=Length[result[Rlist[[ind]]][[1]]],k++,
					{Rapprox, meanR, coefApprox, coef, RMargin, coefMargin}=approximationWithErrorMargin[Rlist[[{ind,ind+1,ind+row, ind+row*col}]], result, k];
					If[closeProportion[Rapprox, Rlist[[ind]], RstepListLocal * (1+margin)] >1/2 && closeProportion[{meanR},Rlist[[ind]], RstepListLocal * (1+margin) ]==1,
						coef = removePrimeCoefRubbish[coef, coefLimit, realOrImaginary];
						If[Length[coef]>0,
						AppendTo[allTogList,{meanR,coef}];
						];
					];
				];
			];
		];
	];
];
];
Print["Rstep: " , Rstep];


tv=Sort[allTogList];
sameLimit= Rstep;
Print["CandList1: ", {Length[tv],tv}];
For[j=1,j < Length[tv],j++,
	k = j + 1;
	While[k <= Length[tv],
		If[Abs[tv[[k,2,1]]-tv[[j,2,1]]] < 10 * sameLimit && Norm[tv[[k,1]]-tv[[j,1]]] < sameLimit,
			tv = Delete[tv, k];
		,
			k++;
		];
	];
];
Print["CandList1: ", {Length[tv],tv}];

candidates = tv[[All,1]];
candStartV = tv[[All,2]];
candStep2 = candidates;
Save[fileNameBase <> "Candidates.txt", candStep2];


(* Parameters  *)
 RstepStart = getRstep[Ltype, 10^(-TRUNCDIGITS + 2)];
 Rlimit = 1; (*10 RstepList[[1]];*)
 ZoomSteps = 3;
 TRUNCDIGITSstart = TRUNCDIGITS;
 NLmethod = "Secant";
 nrOfRuns = 5;
 MaxSolutions = 1;
 extraAtCenter = 3;
 sameNN = True;
 nrOfExtraEquations = 8;
 successLimit = 10^(-2);
 useSameTruncation = True;
 bwidth = 2;

 st = TimeUsed[];
 Nlist = getNumberOfCoefficients[OMEGA,klist,llist,parity,v,bwidth];
 Print[TimeUsed[]-st, Nlist];
 RstartList=SetPrecision[candidates, PRECISIONMULTIPLE TRUNCDIGITSstart];
 startValueList = SetPrecision[candStartV,PRECISIONMULTIPLE TRUNCDIGITSstart];
 For[Rloop=CurrentIndex2,Rloop<=Length[RstartList],Rloop++,
    fileName=fileNameBase <> "Zoom" <> ToString[Rloop];
    Rinit=RstartList[[Rloop]];
	If[Length[startValueList]>=Rloop && Length[startValueList[[Rloop]]]>0,
		startValues = {startValueList[[Rloop]]};
	,
		startValues = {};
	];
    << ZoomOne.m;
    PrependTo[answer,Ldata];
    If[success,
        saveCandidate[fileNameBase <> "Success",answer];
    ,
        saveCandidate[fileNameBase <> "Fail",answer];
    ];
	Put[{CurrentStep,Rloop+1,CurrentIndex1,CurrentIndex3},fileNameBase <> "StateData.m"];
 ];

If[CurrentStep == 2,
	CurrentStep = 3;
	CurrentIndex2 = Rloop+1;
	CurrentIndex3 = 1;
	Put[{CurrentStep,CurrentIndex2,CurrentIndex1, CurrentIndex3},fileNameBase <> "StateData.m"];
];

TimeStep2 = (TimeUsed[]-st)/60;
Save[fileNameBase <> "Time.txt", TimeStep2];
st=TimeUsed[];
 
{CurrentStep,CurrentIndex2,CurrentIndex1, CurrentIndex3} = Get[fileNameBase <> "StateData.m"];





succList=Get[ fileNameBase <> "Success"];
convLimit=10;
runAgain={};
startV={};
RMarginList = {};
For[k=1,k<=Length[succList],k++,
	If[Length[succList[[k,3,3]]]>0,
		If[Min[Table[Norm[succList[[k,3,3,j]]],{j,Length[succList[[k,3,3]]]}]]<convLimit,
			minIdx=minIndex[succList[[k,3]]];
			AppendTo[runAgain,succList[[k,2]]];
			AppendTo[startV,succList[[k,3,1,minIdx]]];
			AppendTo[RMarginList,Norm[succList[[k,3,3,minIdx]]]];
		];
	];
];
Print["Nr in RunAgain: ", Length[runAgain]];

candidates={};
candStartV={};
candRMargin = {};
sameLimit= Rstep / 2;
For[k=1,k<=Length[runAgain],k++,
	j = 1;
	While[j <= Length[candidates] + 1,
		If[j > Length[candidates],
			AppendTo[candidates,runAgain[[k]]];
			AppendTo[candStartV,startV[[k]]];		
			AppendTo[candRMargin,RMarginList[[k]]];
			j = Infinity;		
		,
			If[Abs[startV[[k,1]]-candStartV[[j,1]]] < sameLimit && Norm[runAgain[[k]]-candidates[[j]]] < sameLimit,
				If[RMarginList[[k]] < candRMargin[[j]],
					candidates[[j]] = runAgain[[k]];
					candStartV[[j]] = startV[[k]];
					candRMargin[[j]] = RMarginList[[k]];
				];
				j = Infinity;
			];
		];
		j++;
	];
];
Print["Candidates: ", candidates];
candStep3 = candidates;
Save[fileNameBase <> "Candidates.txt", candStep3];

(* Parameters  *)
 RstepStart = RstepStart/10000;
 Rlimit = Rlimit;
 ZoomSteps = 2;
 TRUNCDIGITSstart = TRUNCDIGITSstart + 6;
 NLmethod = "Secant";
 nrOfRuns = 5;
 MaxSolutions = 1;
 EndPrecision=5;
 extraAtCenter = 4;
 nrOfExtraEquations = 8;
 fileNameBase =  fileNameBase <> "XX"
 successLimit = 10^(-6);
 useSameTruncation = False;
 
 sameNN = False;
 RstartList=SetPrecision[candidates, PRECISIONMULTIPLE TRUNCDIGITSstart];
 startValueList = SetPrecision[candStartV,PRECISIONMULTIPLE TRUNCDIGITSstart];
 For[Rloop=CurrentIndex3,Rloop<=Length[RstartList],Rloop++,
    fileName=fileNameBase <> "Zoom" <> ToString[Rloop];
    Rinit=RstartList[[Rloop]];
	If[Length[startValueList]>=Rloop,
		startValues = {startValueList[[Rloop]]};
	,
		startValues = {};
	];
    << ZoomOne.m;
    PrependTo[answer,Rinit];
    If[success,
        saveCandidate[fileNameBase <> "Success",answer];
    ,
        saveCandidate[fileNameBase <> "Fail",answer];
    ];
	Put[{CurrentStep,CurrentIndex2,CurrentIndex1, Rloop1+1},StringDrop[fileNameBase,-2] <> "StateData.m"];
 ];

TimeStep3 = (TimeUsed[]-st)/60;
Save[StringDrop[fileNameBase, -2] <> "Time.txt", TimeStep3];
