(* ::Package:: *)

Clear[result];

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
	Print["TRUNCDIGITS: ", TRUNCDIGITS];
	
    intermed={};
	startValuePrec = -Log[10, Max[RstepList]];
    
    success = True;
    successLimit = 10^(-8);
	MaxSolutions = 1;
	noImprovement = 0;
	minRMargin = Infinity;

    For[zoomloop=1,zoomloop<=ZoomSteps,zoomloop++,

        (* Global precision parameters. *)
        DIGITS = TRUNCDIGITS + 8;
		PRECISION = Max[PRECISIONMULTIPLE * TRUNCDIGITS, MachinePrecision];
		MYPRECISION =PRECISION;

	   (* Precomputations. *)
		logQ = N[Log[Q],PRECISION + 2];
		logQlogN = Table[logQ-Log[n],{n,1,1000}];

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
        
		Catch[ 	
			Rlist = {Rtuple};	  
			For[k=1, k<=getNrOfParameters[Ltype],k++,
				AppendTo[Rlist, Rtuple + RstepList * Table[DiscreteIndicator[i,k,Table[j,{j,getNrOfParameters[Ltype]}]],{i,getNrOfParameters[Ltype]}]  ];
			];
			
			For[k=1, k<=getNrOfParameters[Ltype] + 1, k++,
                llist = getLlist[Ltype, Rlist[[k]], parity];
				result[Rlist[[k]]] = solveForOneNL[Ldata, klist, llist, reslist, polelist, slist, paralist, Param, nrOfRuns, NN, M, incr, v, expz, startValues, knownCoef, maxPower, NLmethod, MaxSolutions, startValuePrec, realOrImaginary];
				If[Length[result[Rlist[[k]]][[1]]]==0,
					zoomloop=ZoomSteps;
					success = False;
					Throw["NoConvergence"];
				];
				If[k==1,
					startValues={SetPrecision[removePrimeCoefRubbish[result[Rtuple][[1,1]], coefLimit ,realOrImaginary],PRECISION + 2 PRECISIONMULTIPLE]}; 
				];
            ];
			{Rapprox, meanR, coefApprox, meanCoef, RMargin, coefMargin} = SetPrecision[approximationWithErrorMargin[Rlist, result, 1],PRECISION + 2 PRECISIONMULTIPLE];
Print["meanR: ", meanR];
Print["meanR margin: ",  RMargin];
			answer = {Rtuple, result[Rtuple], RMargin, coefMargin};
            AppendTo[intermed,{Rtuple, RMargin, coefMargin, result[Rtuple][[2, 1]]}];      
            Rtuple = meanR;
			
            If[Max[Abs[Rtuple-Rinit]]>Rlimit,
	            zoomloop=ZoomSteps;
                success = False;
            ];

		    coefAcc=Max[coefMargin[[1]], coefMargin[[(1-Abs[realOrImaginary])Length[coefMargin]/2+1]]];
			If[coefAcc==0,
				coefAcc=Max[1,10^(-Accuracy[coefAcc])];
Print["coefAcc: ", coefAcc];
			];
			
			startValuePrec = SetPrecision[Max[2, -Log[10, coefAcc] +2 ], PRECISION + 2 PRECISIONMULTIPLE];
Print["startValuePrec: ", startValuePrec];
			
		];
		
		If[Min[RMargin] < minRMargin,
			noImprovement = 0;
			minRMargin = Min[RMargin];
		,
			noImprovement ++;
			If[noImprovement == 3,
				zoomloop = ZoomSteps;
			]
		];
		
		If[zoomloop < ZoomSteps,
			TRUNCDIGITS += 2;
Print["TRUNCDIGITS: ", TRUNCDIGITS];
			RstepList = Table[Max[RstepList[[i]]/1000, RMargin[[i]]/10],{i,Length[RstepList]}];
Print["RstepList: ", RstepList];
		];
    ];
    
	If[success && Max[RMargin]>successLimit,
		success = False;
	];
	
	newR = { meanR, meanR+RMargin };
    Save[fileName ,{Ldata, Rlist, result, intermed, newR, RMargin, coefMargin, maxPower, realOrImaginary}];
    Clear[result];



