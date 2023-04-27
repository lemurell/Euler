(* ::Package:: *)

(* Load files. *)
<< Lfunctions.m;

(* Parameters to initiate before calling 
 fileName
 Rinit = Rstartlist[[Rloop]];
 Ltype = "GL3";
 level = 1;
 charvalue = 1;
 OMEGA = 1;
 parity = 0;
 Ldata = {Ltype, level, charvalues, OMEGA, parity};
 realOrImaginary = 0;
 RstepStart = {1/10^8, 1/10^8};
 Rlimit = 3/10;
 ZoomSteps =2;
 EndPrecision=8;
 PRECISIONMULTIPLE = 3;
 TRUNCDIGITSstart=0;
 SDIFF = 1/3;
 NLmethod = "Secant";
 nrOfRuns = 100;
 minWidthList = {40, 60};
 minMultList = {1/25,1/200};
 MaxSolutions = 5;
 maxPower = 4;
 nrOfExtraEquations = 8;
 startValues = {};
 knownCoef = {{2, 1.04846245223460506080-0.375239638871383270 I }};
 sameNN = False;
*)


(* Global L-function parameters. *)
Print[Ldata];
Q = QValue[Ltype, level];
klist = getKlist[Ltype];
reslist = {};
polelist = {};
llist = getLlist[Ltype, Rinit, parity];

(* Local testfunction parameters. *)
 (*Param = Table[{1/25, 2/5 + (i-1)/2, 0}, {i, 2}];*)
 v = 2;

    (* Precision parameter. *)

    Rtuple = Rinit;
    RstepList = RstepStart;
    Print[fileName];
	TRUNCDIGITS = TRUNCDIGITSstart;

	Rlist = {Rtuple};
	meanR = Rtuple;
	RMargin = RstepList;
    intermed={};
	startValuePrec = -Log[10, Max[RstepList]];
    
    success = True;
	AlreadyExtra = False;

    For[zoomloop=1,zoomloop<=ZoomSteps,zoomloop++,

        (* Global precision parameters. *)
        DIGITS = TRUNCDIGITS + 12;
		PRECISION = Max[PRECISIONMULTIPLE * TRUNCDIGITS, MachinePrecision];
		MYPRECISION =PRECISION;

	   (* Precomputations. *)
		logQ = N[Log[Q],PRECISION + 2];
		logQlogN = Table[logQ-Log[n],{n,1,1000}];

        (* Local testfunction parameters. *)
		(*If[zoomloop>1 || !sameNN,
			minWidth = minWidthList[[Min[zoomloop,Length[minWidthList]]]];
			startNN = initNNList[[Min[zoomloop,Length[initNNList]]]];
			minMult = minMultList[[Min[zoomloop,Length[minMultList]]]];
			{NN, sSeq, paraSeq} = getFuncEqParameters[1,klist,llist,reslist,polelist,Param,v,minWidth,startNN,minMult,realOrImaginary];
		];*)
		
		If[testfunctiontypeZoom == "Classic",
			Param=Table[{1/25,2/5+(i-1)2/5,0},{i,2}];
			{NN, sSeq, paraSeq}  = getFuncEqParameters[1,klist,llist,reslist,polelist,Param,v,minWidth,NN,minMult,realOrImaginary];
		,
			NN = Nlist[[TRUNCDIGITS]];
		];
		If[AlreadyExtra,
			NN += 5;
		];
		startNN = Ceiling[NN * 11/10] ;
		{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata, NN ,maxPower, knownCoef];
		stillUnknowns = Union[plist,highplist];
		nrOfUnknowns = (2-Abs[realOrImaginary]) Length[stillUnknowns];
        nrOfEquations = nrOfUnknowns + nrOfExtraEquations;
		If[testfunctiontypeZoom == "Classic",
			{slist, paralist} = getSlist[sSeq, paraSeq, nrOfEquations,SDIFF,realOrImaginary];
		,
			{slist, paralist, Param} = getSlistDavid[nrOfEquations];
		];
Print["Unknowns:  ", stillUnknowns , " " ,Length[stillUnknowns], "st"];
Print["Size of system: ", nrOfEquations, " x ", nrOfUnknowns];

        (* Precomputations. *)
        incr=2*Pi*v/Log[10]/DIGITS;  (* Not good if v is large *)
        M=5 Sqrt[Log[10]*DIGITS];
        expz=Table[Exp[(v+I*k)*logQlogN[[n]]]/(v+I*k),{n,1,NN},{k,-M,M,incr}];

        llist = getLlist[Ltype, Rtuple, parity];
		result[Rtuple]=solveForOneNL[Ldata, klist, llist, reslist, polelist, slist, paralist, Param, nrOfRuns, NN, M, incr, v, expz, startValues, knownCoef, maxPower, NLmethod, MaxSolutions, startValuePrec, realOrImaginary];
		minIdx = minIndex[result[Rtuple]];
		If[minIdx == 0,
			error = Infinity;
        	counter = Infinity;
		,
        	error=Norm[result[Rtuple][[3,minIdx]]];
            AppendTo[intermed,{Rtuple,error,result[Rtuple][[2, minIdx]]}];
        	counter = 1;
			MaxSolutions = 1;
			startValues={SetPrecision[removePrimeCoefRubbish[result[Rtuple][[1,minIdx]], coefLimit ,realOrImaginary],PRECISION]}; 
		];
        While[(error > 10^(-TRUNCDIGITS+2) || (counter==1 && zoomloop==ZoomSteps)  ) && counter<6,  
			Catch[ 		  
				minError = Infinity;
				For[k=1, k<=getNrOfParameters[Ltype],k++,
					AppendTo[Rlist, Rtuple + RstepList * Table[DiscreteIndicator[i,k,Table[j,{j,getNrOfParameters[Ltype]}]],{i,getNrOfParameters[Ltype]}]  ];
				];
				For[k=-getNrOfParameters[Ltype],k<=-1,k++,
                	llist = getLlist[Ltype, Rlist[[k]], parity];
					result[Rlist[[k]]] = solveForOneNL[Ldata, klist, llist, reslist, polelist, slist, paralist, Param, nrOfRuns, NN, M, incr, v, expz, startValues, knownCoef, maxPower, NLmethod, MaxSolutions, startValuePrec, realOrImaginary];
					If[Length[result[Rlist[[k]]][[1]]]==0,
						error=Infinity;
						Throw["NoConvergence"];
					];
					minError = Min[minError, Norm[result[Rlist[[k]]][[3,1]]]];
            	];
				{Rapprox, meanR, coefApprox, meanCoef, RMargin, coefMargin} = SetPrecision[approximationWithErrorMargin[Rlist[[Table[i,{i,-getNrOfParameters[Ltype]-1,-1}]]], result, minIdx],PRECISION];
Print["meanR: ", meanR];
Print["meanR margin: ",  RMargin];
            	Rtuple = SetPrecision[meanR,PRECISION];          

            	If[Max[Abs[Rtuple-Rinit]]>Rlimit,
                	zoomloop=ZoomSteps;
                	success = False;
                	Break[];
            	];

			    coefAcc=Max[coefMargin[[1]], coefMargin[[(1-Abs[realOrImaginary])Length[coefMargin]/2+1]]];
				If[coefAcc==0,
					coefAcc=Max[1,10^(-Accuracy[coefAcc])];
Print["coefAcc: ", coefAcc];
				];
				
 			   startValuePrec = SetPrecision[Max[2, -Log[10, coefAcc] +2 ],PRECISION];
Print["startValuePrec: ", startValuePrec];

            	AppendTo[Rlist,Rtuple];
            	llist = getLlist[Ltype, Rtuple, parity];
                result[Rtuple]=solveForOneNL[Ldata, klist, llist, reslist, polelist, slist, paralist, Param, nrOfRuns, NN, M, incr, v, expz, startValues, knownCoef, maxPower, NLmethod, MaxSolutions, startValuePrec, realOrImaginary];
				minIdx = minIndex[result[Rtuple]];
				If[minIdx==0,
					error=Infinity;
					Throw["NoConvergence"];
				];
				error=Norm[result[Rtuple][[3,minIdx]]];
Print["Error: ", error, " Coef: ", result[Rtuple][[2]]];
            	AppendTo[intermed,{Rtuple,error,result[Rtuple][[2, minIdx]]}];
			];
			If[error >  10^(4-counter),
				counter = Infinity;
Print["Not fast enough convergence."];
			,
				startValues={SetPrecision[removePrimeCoefRubbish[result[Rtuple][[1,minIdx]], coefLimit, realOrImaginary],PRECISION]}; 
            	RstepList = Table[Min[RMargin[[i]]/10, RstepList[[i]]],{i,Length[RstepList]}];
Print["Rstep: ", RstepList];
            	counter++;
			];
        ];
    
        If[error > 10^(-TRUNCDIGITS+3),  (* OBS \[CapitalADoubleDot]NDRAT TILL SN\[CapitalADoubleDot]LLARE. *)
            zoomloop = ZoomSteps;
            success = False;
        ,
            TRUNCDIGITS += 2;
            RstepList = Table[Max[10^(-TRUNCDIGITS-1),Min[RstepList[[i]],RMargin[[i]]/10]],{i,Length[RstepList]}];
        ];

		If[zoomloop==ZoomSteps && Max[RMargin]>10^(-EndPrecision)/2 && !AlreadyExtra  && !sameNN,
Print["Doing one extra lap to get enough precision."];
			zoomloop--;
			AlreadyExtra=True;
		];
        answer={Rtuple,result[Rtuple],RMargin};
    ];

	newR = { meanR, meanR+RMargin };
    Save[fileName ,{Ldata, Rlist, result, intermed, newR, RMargin, coefMargin, maxPower, realOrImaginary}];
    Clear[result];



