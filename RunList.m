(* ::Package:: *)

<< Lfunctions.m;
Clear[result];

(* Parameters to initiate before calling *)
 Ltype = "SP4";
 level = 1;
 charvalue = "No";
 (* parity = {OddOrEven, weight};
 OMEGA = (-1)^OddOrEven * I^weight;*)
 OMEGA = (-1)^(parity[[1]]+parity[[2]]);
 Ldata = {Ltype, level, charvalue, OMEGA, parity};
 realOrImaginary = 1;
 coefLimit = 4;
 RstepStart = getRstep[Ltype, 1/10^5];
 (*RstepStart = {10^(-7),10^(-6),10^(-4)};*)
 Rlimit = 1/100;
 Rinit = candidates[[-1]];
 ZoomSteps = 5;
 EndPrecision = 6;
 PRECISIONMULTIPLE = 5;
 TRUNCDIGITSstart = 11;
 PRECISION= Max[PRECISIONMULTIPLE * TRUNCDIGITSstart, MachinePrecision];
 SDIFF = 1/3;
 NLmethod = "Secant";
 nrOfRuns = 20;
 MaxSolutions = 1;
 sameNN = False;
 useSameTruncation = False;
 testfunctiontypeZoom = "DS";
 v=2;
 bwidth = 2;
 Mfactor = 7;
 klist = getKlist[Ltype];
 Q = QValue[Ltype, level];
 logQ = N[Log[Q],PRECISION + 2];
 logQlogN = Table[logQ-Log[n],{n,1,1000}];
 llist = getLlist[Ltype, Rinit, parity];
 Nlist = getNumberOfCoefficients[OMEGA,klist,llist,parity,v,bwidth];
 Print[Nlist];
 Print[{OMEGA,klist,llist,parity,v,bwidth}];
 maxPower = 4;
 nrOfExtraEquations = 8;
 startValues = {}; 
 knownCoef = {{2,1.0484624522346050608`20.0205528820314- 0.37523963887138327`17.57430870932808 I},{5,0.766177577020120230416453-0.1471014533310I}};
 knownCoef={};
 fileNameBase = fileNameBase <>"XX";

 RstartList=SetPrecision[candidates, PRECISIONMULTIPLE TRUNCDIGITSstart];
 startValueList = SetPrecision[candStartV, PRECISIONMULTIPLE TRUNCDIGITSstart];
 Print["Rlist: ", N[RstartList,5]];
 For[Rloop=firstR,Rloop<=Length[RstartList],Rloop++,
    fileName=fileNameBase <> "Zoom" <> ToString[Rloop];
    Rinit=RstartList[[Rloop]];
	If[Length[startValueList]>=Rloop,
		startValues = {startValueList[[Rloop]]};
	];
	If[Length[candParity]>=Rloop,
		parity = candParity[[Rloop]];
        Ldata = {Ltype, level, charvalue, OMEGA, parity};
Print["Parity: ", parity];
	];
	If[Length[candOMEGA]>=Rloop,
		OMEGA = candOMEGA[[Rloop]];
        Ldata = {Ltype, level, charvalue, OMEGA, parity};
Print["OMEGA: ", OMEGA];
	];

    << ZoomOne.m;

    PrependTo[answer,Rinit];
    PrependTo[answer,Ldata];
    If[success,
        saveCandidate[fileNameBase <> "Success",answer];
    ,
        saveCandidate[fileNameBase <> "Fail",answer];
    ];
 ]






