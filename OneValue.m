(* ::Package:: *)

<< Lfunctions.m;
Clear[result];

(* Parameters to initiate before calling *)
 Ltype = "GL3";
 level = 4;
 charvalue = "No";
 parity = {1, 1, 0};
 OMEGA = {1/2, -Sqrt[3]/2};
 Ldata = {Ltype, level, charvalue, OMEGA, parity};
 realOrImaginary = 0;
 coefLimit = 3;
 PRECISIONMULTIPLE = 5;
 TRUNCDIGITS = 6;
 startValuePrec = 4;
 SDIFF = 1/3;
 NLmethod = "Secant";
 nrOfRuns = 100;
 MaxSolutions = 1;
 startNN= 30;
 minWidth = 40;
 minMult = 1/200;
 maxPower = 4;
 unknownsAtStart = 0;
 nrOfExtraEquations = 1;
 extraEq = {bb1[2] + 1/4};
 startValues = {}; 
 knownCoef = {{2, 1.04846245223460506080-0.375239638871383270 I }};
 knownCoef = {};
 fileName = "Runs/NewOMEGA_GL3Level4_c1_110.txt"
 
 RstartList=SetPrecision[candidates, PRECISIONMULTIPLE TRUNCDIGITS];
 startValueList = SetPrecision[candStartV, PRECISIONMULTIPLE TRUNCDIGITS];
 Print["Rlist: ", N[RstartList,5]];
 For[Rloop=1,Rloop<=Length[RstartList],Rloop++,
    Rtuple=RstartList[[Rloop]];
	If[Length[startValueList]>=Rloop,
		startValues = {startValueList[[Rloop]]};
	];
(* Global L-function parameters. *)
Print[Ldata];
Q = QValue[Ltype, level];
klist = getKlist[Ltype];
reslist = {};
polelist = {};
llist = getLlist[Ltype, Rtuple, parity];

(* Local testfunction parameters. *)
 Param = Table[{1/25, 2/5 + (i-1)/2, 0}, {i, 2}];
 v = 2;

    Print[fileName];

        (* Global precision parameters. *)
        DIGITS = TRUNCDIGITS + 8;
		PRECISION = Max[PRECISIONMULTIPLE * TRUNCDIGITS, MachinePrecision];
		MYPRECISION =PRECISION;

	   (* Precomputations. *)
		logQ = N[Log[Q],PRECISION + 2];
		logQlogN = Table[logQ-Log[n],{n,1,1000}];

        (* Local testfunction parameters. *)
		{NN, sSeq, paraSeq} = getFuncEqParameters[{1,0},klist,llist,reslist,polelist,Param,v,minWidth,startNN,minMult,realOrImaginary];
		{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata, NN ,maxPower, knownCoef];
		stillUnknowns = Union[plist,highplist];
		nrOfUnknowns = (2-Abs[realOrImaginary]) Length[stillUnknowns];
        nrOfEquations = nrOfUnknowns + nrOfExtraEquations;
		{slist, paralist} = getSlist[sSeq, paraSeq, nrOfEquations,SDIFF,realOrImaginary];
		Print["Unknowns:  ", stillUnknowns , " " ,Length[stillUnknowns], "st"];
        Print["Size of system: ", nrOfEquations, " x ", nrOfUnknowns];

        (* Precomputations. *)
        incr=2*Pi*v/Log[10]/DIGITS;  (* Not good if v is large *)
        M=Sqrt[Log[10]*DIGITS/Param[[1,1]]];
        expz=Table[Exp[(v+I*k)*logQlogN[[n]]]/(v+I*k),{n,1,NN},{k,-M,M,incr}];

        llist = getLlist[Ltype, Rtuple, parity];
		result[Rtuple]=solveForOneByStep[Ldata, klist, llist, reslist, polelist, slist, paralist, Param, nrOfRuns, NN, M, incr, v, expz, startValues, knownCoef, maxPower, NLmethod, MaxSolutions, startValuePrec, unknownsAtStart,realOrImaginary,extraEq];

        Save[fileName,result];
		Clear[result];
 ]






