(* ::Package:: *)

addEquations[intpart_, polepart_, realOrImaginary_:0, OMEGA_:{1, 0}]:=Module[{realRow, imagRow, i, j, answer, maxAbs, realRow2 , imagRow2 },
    answer={};
	maxAbs = Max[Abs[intpart]];
    For[i=1,i<Length[polepart],i++,
		If[realOrImaginary == 0,
		    realRow = Re[intpart[[1,1]] - intpart[[i+1,1]]] + OMEGA[[1]](Re[intpart[[1,2]] - intpart[[i+1,2]]]) - OMEGA[[2]](Im[intpart[[1,2]] - intpart[[i+1,2]]]);
		    realRow[[1]] = realRow[[1]] + Re[polepart[[1]] - polepart[[i+1]]];
		    imagRow = Im[intpart[[1,1]] - intpart[[i+1,1]]] + OMEGA[[1]](Im[intpart[[1,2]] - intpart[[i+1,2]]]) + OMEGA[[2]](Re[intpart[[1,2]] - intpart[[i+1,2]]]);
		    imagRow[[1]] = imagRow[[1]] + Im[polepart[[1]] - polepart[[i+1]]];
		    realRow2 = Re[intpart[[1,1]] - intpart[[i+1,1]]] - OMEGA[[1]](Re[intpart[[1,2]] - intpart[[i+1,2]]]) + OMEGA[[2]](Im[intpart[[1,2]] - intpart[[i+1,2]]]);
		    realRow2[[1]] = realRow2[[1]] + Re[polepart[[1]] - polepart[[i+1]]];
		    imagRow2 = Im[-intpart[[1,1]] + intpart[[i+1,1]]] + OMEGA[[1]](Im[intpart[[1,2]] - intpart[[i+1,2]]]) + OMEGA[[2]](Re[intpart[[1,2]] - intpart[[i+1,2]]]);
		    imagRow2[[1]] = imagRow2[[1]] + Im[polepart[[1]] - polepart[[i+1]]];
		,
			If[OMEGA == {-1, 0},
		        realRow = Im[intpart[[1,1]] - intpart[[i+1,1]]] + OMEGA[[1]](Im[intpart[[1,2]] - intpart[[i+1,2]]]) + OMEGA[[2]](Re[intpart[[1,2]] - intpart[[i+1,2]]]);
		        realRow[[1]] = realRow[[1]] + Im[polepart[[1]] - polepart[[i+1]]];
		    ,
			    realRow = Re[intpart[[1,1]] - intpart[[i+1,1]]] + OMEGA[[1]](Re[intpart[[1,2]] - intpart[[i+1,2]]]) - OMEGA[[2]](Im[intpart[[1,2]] - intpart[[i+1,2]]]);
			    realRow[[1]] = realRow[[1]] + Re[polepart[[1]] - polepart[[i+1]]];
		    ,
			    realRow = Re[intpart[[1,1]] - intpart[[i+1,1]]] + OMEGA[[1]](Re[intpart[[1,2]] - intpart[[i+1,2]]]) - OMEGA[[2]](Im[intpart[[1,2]] - intpart[[i+1,2]]]);
			    realRow[[1]] = realRow[[1]] + Re[polepart[[1]] - polepart[[i+1]]];
			];
		];
        If[realOrImaginary == 1,
	       AppendTo[answer, realRow];
	    ,
	       AppendTo[answer, alternateLists[realRow, imagRow2]];
	       AppendTo[answer, alternateLists[imagRow, realRow2]];
		];
	];

    answer/maxAbs
 ]

all2Unknowns[Ldata_, coef_,maxPower_:4, removeRubbish_:False,  maxNrOfCoef_:200, knownCoef_:{}]:=Module[{pcoef, k, idx, plist, highplist, nonplist, knownlist, unknowns},
	pcoef={};
	idx = 1;
	{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata,maxNrOfCoef,maxPower, knownCoef];
	unknowns = Union[plist, highplist];
	For[k=2,k<=Length[coef]/2,k++,
		If[MemberQ[unknowns,k],
			pcoef = Insert[pcoef, coef[[2k-1]],idx];
			AppendTo[pcoef, coef[[2k]]];
			idx++;
		];
	];
	If[removeRubbish,
		pcoef=removePrimeCoefRubbish[pcoef, getDegree[Ldata[[1]]]];
	];
	pcoef
]

alternateLists[list1_,list2_]:=Module[{answer, i},
  answer = Table[{}, {i, 2 Length[list1]}];
  For[i = 1, i <= Length[list1], i++, 
      answer[[2i - 1]] = list1[[i]];  
      answer[[2i]] = list2[[i]];
  ];
  answer
]

approximationWithErrorMargin[Rlist_, result_,resultIdx_:1]:=Module[{nA,signCoef, RApprox, meanR, coefApprox, meanCoef,RMargin,coefMargin, nrOfExtraEquations,minInd, k,r,s,t, dim },
	signCoef = 3;
	minInd[1] = resultIdx;
	dim = Length[Rlist]-1;
	For[k=2,k<= dim+1 ,k++,
		minInd[k] = minDistIndex[result[Rlist[[1]]][[1,resultIdx]], result[Rlist[[k]]], signCoef];
	];
	If[ Product[minInd[k],{k,dim+1}] ==0,
		{}
	,
		nrOfExtraEquations=Length[result[Rlist[[1]]][[3,1]]];
		Switch[dim,
			1,
				nA=Table[newApprox[Rlist,Table[result[Rlist[[k]]][[3,minInd[k],{r}]],{k,dim+1}],Table[result[Rlist[[k]]][[1,minInd[k]]],{k,dim+1}]],{r,1,nrOfExtraEquations}];
			,2,
				nA=Table[newApprox[Rlist,Table[result[Rlist[[k]]][[3,minInd[k],{r,s}]],{k,dim+1}],Table[result[Rlist[[k]]][[1,minInd[k]]],{k,dim+1}]],{r,1,nrOfExtraEquations-1},{s,r+1,nrOfExtraEquations}];
			,3,
				nA=Table[newApprox[Rlist,Table[result[Rlist[[k]]][[3,minInd[k],{r,s,t}]],{k,dim+1}],Table[result[Rlist[[k]]][[1,minInd[k]]],{k,dim+1}]],{r,1,nrOfExtraEquations-2},{s,r+1,nrOfExtraEquations-1},{t,s+1,nrOfExtraEquations}];
		];
(*Print[nA];*)
		RApprox=Flatten[nA,dim-1][[All,1]];
(*Print["RApprox: ",RApprox];*)
(*Save["RApprox",RApprox];*)
		meanR=Mean[RApprox];
(*Print["meanR: ",meanR];*)
		meanR=Table[Mean[Sort[RApprox[[All,k]]][[Range[Ceiling[Length[RApprox]/4], Ceiling[3 Length[RApprox]/4]]]]] ,{k,Length[Rlist[[1]]]}];
(*Print["meanR: ",meanR];*)
		coefApprox=Flatten[nA,dim-1][[All,2]];
		meanCoef=Mean[coefApprox];
		RMargin = Table[Max[Table[Abs[RApprox[[k,j]]-meanR[[j]]],{k,Length[RApprox]}]],{j,Length[RApprox[[k]]]}];
(*Print["Rmargin: ",RMargin];*)
		For[k=1,k<=Length[Rlist[[1]]],k++,
			If[RMargin[[k]] == 0,
Print["Rmargin == 0"];
				RMargin[[k]] = Max[Max[RMargin], Norm[Rlist[[1]]-Rlist[[2]]]];
			];
		];
		coefMargin=Table[ Max[Table[Abs[coefApprox[[k,j]]-meanCoef[[j]]],{k,Length[RApprox]}]],{j,Length[meanCoef]}];
		{RApprox, meanR, coefApprox, meanCoef,RMargin,coefMargin}
	]
]

closeInListIndex[val_,list_,sameLimit_]:=Module[{i, ans},
	i=1;
	ans=0;
	While[i<=Length[list],
		If[Norm[val-list[[i]]] < sameLimit,
			ans=i;
			Break[];
		];
		i++;
	];
	ans
]

(* Gives the proportion of vectors in approx which are within +/- limit of
exact. *)
closeProportion[approx_, exact_, limit_]:=Module[{ctr, k, diff, j, dimCtr},
	ctr=0;
	For[k=1,k<=Length[approx],k++,
		diff=Abs[approx[[k]]-exact];
		dimCtr = 0;
		For[j=1,j<=Length[exact],j++,
			If[diff[[j]]< limit[[j]],
				dimCtr++;
			,
				Break[];
			];
		];
		If[dimCtr == Length[exact],
			ctr++;
		];
	];
	If[Length[approx] == 0,
		1
	,
		ctr/Length[approx]
	]
]

(* L1= {R1,R2,coef,newR} *)
compareLfuntions[L1_,L2_]:=Module[{k, pcoef, R1diff, R2diff, prec1, prec2, R1, R2, R1prec1, R2prec1, R1prec2, R2prec2, bestL, bestR1prec, bestR2prec, bestCoef},
	R1diff = Abs[L1[[1]]-L2[[1]]];
	R2diff = Abs[L1[[2]]-L2[[2]]];
	R1prec1 = Max[Table[Abs[L1[[4, i, 1]]-L1[[4, j, 1]]],{i,4},{j,i+1,4}]];
	R2prec1 = Max[Table[Abs[L1[[4, i, 2]]-L1[[4, j, 2]]],{i,4},{j,i+1,4}]];
	R1prec2 = Max[Table[Abs[L2[[4, i, 1]]-L2[[4, j, 1]]],{i,4},{j,i+1,4}]];
	R2prec2 = Max[Table[Abs[L2[[4, i, 2]]-L2[[4, j, 2]]],{i,4},{j,i+1,4}]];
	If[Max[R1prec1, R2prec1] < Max[R1prec2, R2prec2],
		bestL = L1;
		bestR1prec = -Ceiling[Log[10, 2 R1prec1]];
		bestR2prec = -Ceiling[Log[10, 2 R2prec1]];
	,
		bestL = L2;
		bestR1prec = -Ceiling[Log[10, 2 R1prec2]];
		bestR2prec = -Ceiling[Log[10, 2 R2prec2]];
	];
	
	R1 = SetPrecision[bestL[[1]], bestR1prec + Ceiling[Log[10,bestL[[1]]]]];
	R2 = SetPrecision[bestL[[2]], bestR2prec + Ceiling[Log[10,bestL[[2]]]]];
	Print["Parameters: ", {R1, R2}, " Precision: ", {bestR1prec, bestR2prec}];
	Print["Difference from other: ",  {R1diff, R2diff}];

	bestCoef = Table[SetPrecision[bestL[[3,k]], -Ceiling[Log[10, 2 Abs[L1[[3,k]]- L2[[3,k]]]]] + If[L1[[3,k]]==0,0,Ceiling[Log[10, Abs[L1[[3,k]]]]]]] , {k, Min[Length[L1[[3]]], Length[L2[[3]]]]}];
	pcoef = {};
	For[k=1,k<= Length[bestCoef]/2,k++,
		If[PrimeQ[k],
			AppendTo[pcoef, {k, bestCoef[[2k-1]] + I bestCoef[[2k]]}];
		];
	];
	Print["Prime coefficients: ", pcoef];

	{R1, R2, bestCoef} 
]

(* Computes the two rows for the the equation (real and imaginary part)
that corresponds to the supplied data. In particular s is the complex
point at which the computation is made.
 --------------------------------------------------------------------*)
computeEquation[OMEGA_,klist_,llist_,reslist_,polelist_,Param_,v_,NN_ ,incr_, M_,expz_,s_,realOrImaginary_:0]:=Module[{    exps1,exps2,n,gamma1,gamma2,i,j,k,Aloop,glist1,glist2,answer,polepart,matrixrows},
 answer=Table[{},{j,Length[Param]}];
    polepart=Table[0,{j,Length[Param]}];

        exps1=Table[Exp[s*logQlogN[[n]]],{n,1,NN}];
        exps2=Table[Exp[(1-s)*logQlogN[[n]]],{n,1,NN}];

        gamma1=Table[Product[ N[Gamma[ klist[[j]]((v+I*k)+s)+llist[[j]] ],PRECISION],{j,1,Length[klist]}],{k,-M,M,incr}];

        gamma2=Table[Product[ N[Gamma[ klist[[j]]((v+I*k)+1-s)+Conjugate[llist[[j]]] ],PRECISION],{j,1,Length[klist]}],{k,-M,M,incr}];
        For[Aloop=1,Aloop<=Length[Param],Aloop++,
            glist1=Table[ N[g[s+v+I*k , Param[[Aloop]], s],PRECISION],{k,-M,M,incr}];
            glist2=Table[ N[g[s-v-I*k , Param[[Aloop]], s],PRECISION],{k,-M,M,incr}];

            answer[[Aloop]]= {incr * Table[ 
                      exps1[[n]] * Sum[
                        gamma1[[k]] * glist1[[k]] * expz[[n,k]], {k,1,Length[gamma1]}] ,{n,NN}]/2/Pi, 
                        incr * Table[ exps2[[n]] *  Sum[gamma2[[k]] * glist2[[k]] * expz[[n,k]],
                    {k,1,Length[gamma1]}] ,{n,NN}] / 2/Pi};

            polepart[[Aloop]]=Sum[reslist[[i]] * g[polelist[[i]],Param[[Aloop]], s]/(s-polelist[[i]]),{i,Length[polelist]}];
        ];

        matrixrows = addEquations[answer, polepart, realOrImaginary, OMEGA];
		matrixrows = matrixrows / Max[Abs[matrixrows]]

]


(* Computes the value of the L-function with coefficients in coeflist and
 functional equation given by the data in Q,w
 --------------------------------------------------------------------*)
computeL[s_,coeflist_,Q_,w_,klist_,llist_,reslist_,polelist_,A_,truncdigits_,v_ ]:=Module[{digits, incr, M, NN, g, exps1, exps2, expz, gamma1, gamma2, glist1, glist2, answer, polepart, precisionmultiple},
	digits = truncdigits +8;
	precisionmultiple=4;
    PRECISION = Max[precisionmultiple * truncdigits, MachinePrecision];
    MYPRECISION = PRECISION;
    incr=2*Pi*v/Log[10]/digits;  (* Not good if v is large *)
    M=Sqrt[Log[10]*digits/A];
    NN = Length[coeflist];
    g[z_,b_]:=E^(b*(z-s)^2 );
    logQ=N[Log[Q],PRECISION];
    logQlogN=Table[logQ-Log[n],{n,1,NN}];
    exps1=Table[Exp[s*logQlogN[[n]]],{n,1,NN}];
    exps2=Table[w*Exp[(1-s)*logQlogN[[n]]],{n,1,NN}];
    expz=Table[Exp[(v+I*k)*logQlogN[[n]]]/(v+I*k),{n,1,NN},{k,-M,M,incr}];
    gamma1=Table[Product[ N[Gamma[ klist[[j]]((v+I*k)+s)+llist[[j]] ],PRECISION],{j,1,Length[klist]}],{k,-M,M,incr}];
    gamma2=Table[Product[ N[Gamma[ klist[[j]]((v+I*k)+1-s)+Conjugate[llist[[j]]] ],PRECISION],{j,1,Length[klist]}],{k,-M,M,incr}];
    glist1=Table[ N[g[s+v+I*k , A],PRECISION],{k,-M,M,incr}];
    glist2=Table[ N[g[s-v-I*k , A],PRECISION],{k,-M,M,incr}];
    answer= incr * Table[ 
			 exps1[[n]] * Sum[
					 gamma1[[k]] * glist1[[k]] * expz[[n,k]], {k,1,Length[gamma1]}]
		      + exps2[[n]] *  Sum[gamma2[[k]] * glist2[[k]] * expz[[n,k]],
			    {k,1,Length[gamma1]}] ,{n,NN}] / 2/Pi;
	

        polepart=Sum[reslist[[i]] * g[polelist[[i]],A]/(s-polelist[[i]]),{i,Length[polelist]}];

	(Sum[coeflist[[i]] * answer[[i]],{i,Length[coeflist]}]
		      + polepart)/Q^s/Product[ 
					      Gamma[ klist[[i]] s+llist[[i]] ],
		                      {i,1,Length[klist]}]
]

computeEquationMatrix[Ldata_,klist_,llist_,reslist_,polelist_,phaseFactor_,slist_, paralist_,Param_,NN_,M_,incr_,v_,expz_,realOrImaginary_:0]:=Module[{OMEGA,errorSize,sloop,s,exps1,exps2,n,gamma1,gamma2,i,j,k,Aloop,glist1,glist2,answer,polepart,matrixrows,locParam},
(*st=TimeUsed[];
Print["This is version 2.99 and the level is: ", level];*)
	OMEGA = Ldata[[4]];
    fullmatrix= {};
    answer=Table[{},{j,Length[Param]}];
    polepart=Table[0,{j,Length[Param]}];
	locParam = Param;

    For[sloop=1,sloop<=Length[slist],sloop++,
        s=slist[[sloop]];
		locParam[[All,1]] = paralist[[sloop]];

        exps1=Table[Exp[s*logQlogN[[n]]],{n,1,NN}];
        exps2=Table[Exp[(1-s)*logQlogN[[n]]],{n,1,NN}];
        gamma1=Table[Product[ N[Gamma[ klist[[j]]((v+I*k)+s)+llist[[j]] ],PRECISION],{j,1,Length[klist]}],{k,-M,M,incr}];
        gamma2=Table[Product[ N[Gamma[ klist[[j]]((v+I*k)+1-s)+Conjugate[llist[[j]]] ],PRECISION],{j,1,Length[klist]}],{k,-M,M,incr}];
        For[Aloop=1,Aloop<=Length[Param],Aloop++,
            glist1=Table[ N[g[s+v+I*k , locParam[[Aloop]], s],PRECISION],{k,-M,M,incr}];
            glist2=Table[ N[g[s-v-I*k , locParam[[Aloop]], s],PRECISION],{k,-M,M,incr}];

            answer[[Aloop]]= {N[phaseFactor] * incr * Table[ 
                      exps1[[n]] * Sum[
                        gamma1[[k]] * glist1[[k]] * expz[[n,k]], {k,1,Length[gamma1]}] ,{n,NN}]/2/Pi, 
                        Conjugate[N[phaseFactor]] * incr * Table[ exps2[[n]] *  Sum[gamma2[[k]] * glist2[[k]] * expz[[n,k]],
                    {k,1,Length[gamma1]}] ,{n,NN}] / 2/Pi};

            polepart[[Aloop]]=Sum[reslist[[i]] * g[polelist[[i]],locParam[[Aloop]], s]/(s-polelist[[i]]),{i,Length[polelist]}];
        ];

        matrixrows = addEquations[answer, polepart, realOrImaginary, OMEGA];
        For[j=1,j<=Length[matrixrows],j++,
(*Print["matrixrows: ", Precision[matrixrows[[1]]], "  nr: ", sloop];*)
            AppendTo[fullmatrix, matrixrows[[j]]];
        ];
    ];
(* Print["Smallest element: ", Max[Abs[fullmatrix[[All,-1]]]]];*)
	fullmatrix
]

(* Assume u>=v, converts Bian/Booker to our parameters. *)
convBB2FKL[u_,v_]:={(2v+u)/6,(u-v)/6}


(* Assume R1>=R2, converts our parameters to Bian/Booker . *)
convFKL2BB[Rpair_]:={2(Rpair[[1]]+2 Rpair[[2]]),2(Rpair[[1]]-Rpair[[2]])}

data2Lcalc[Ldata_,eig_,coefList_,fileName_,newR_,prec_:0]:=Module[{ Ltype, level, OMEGA, parity, pcoef, i,j, k, newRerror,len, locPrec, shift},
	Ltype = Ldata[[1]];
	level = Ldata[[2]];
	OMEGA = Ldata[[4]];
	If[Length[Ldata]>4,
		parity = Ldata[[5]];
	,
		parity = {};
	];
Print["Parity: ",parity];
	If[prec==0,
		newRerror = Max[Flatten[Table[ Max[ Table[ Abs[newR[[i,k]] - newR[[j,k]]],{k,Min[2,Length[newR[[i]]]]}]] ,{i,Length[newR]-1}, {j,i+1,Length[newR]}]]];
		If[newRerror == 0,		
			locPrec=Accuracy[newRerror];
		,
			locPrec =Floor[-Log[10,newRerror]-Log[10,2]];
		];
	,
		locPrec = prec;
	];
    str=OpenWrite[fileName];
	If[Ltype == "SP4" || Ltype == "SP6" || Ltype == "SP4Std" || Ltype == "HoloxMaass" || Ltype == "Maass",
	    Write[str,2];
	,
	    Write[str,3];
	];
    Write[str,0];
    Write[str,Length[coefList]/2];
    Write[str,0];
	Switch[Ltype,
		"GL3",
			If[parity=={1},
				shift = "0.5 ";
			,
				shift = "0 ";
			];
		    Write[str,3];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[Chop[eig[[2]],10^(-8)], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-(eig[[1]]+eig[[2]]), locPrec]],"\n"];
		,"Maass",
			If[parity[[1]]==1,
				shift = "0.5 ";
			,
				shift = "0 ";
			];
		    Write[str,2];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[eig[[1]]/2, locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[-eig[[1]]/2, locPrec]],"\n"];
		,"HoloxMaass",
			If[parity[[1]]==1,
				shift = "0.5 ";
			,
				shift = "0 ";
			];
		    Write[str,3];
		    Write[str,1];
		    WriteString[str,ToString[N[(parity[[2]]-1)/2, locPrec]]," 0","\n"];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[-eig[[1]], locPrec]],"\n"];
		,"GL3Holo",
			If[parity[[1]]==1,
				shift = "0.5 ";
			,
				shift = "0 ";
			];
		    Write[str,2];
		    Write[str,1];
		    WriteString[str,ToString[N[(parity[[2]]-1)/2, locPrec]],ToString[N[-eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,shift,ToString[N[eig[[1]], locPrec]],"\n"];
		,"SP4",
			shift = {"0 ", "0 "};
			If[Length[parity]==2,
				If[parity[[1]]==1,
					shift[[1]] = "0.5 ";
				];
				If[parity[[2]]==1,
					shift[[2]] = "0.5 ";
				];
			];
		    Write[str,4];
		    Write[str,0.5];
		    WriteString[str, shift[[1]], ToString[N[eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str, shift[[2]],ToString[N[eig[[2]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str, shift[[1]],ToString[N[-eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str, shift[[2]],ToString[N[-eig[[2]], locPrec]],"\n"];
		,"SP4Std",
		    Write[str,5];
		    Write[str,0.5];
		    WriteString[str,"0 0","\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[2]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-eig[[2]], locPrec]],"\n"];
		,"GL4",
		    Write[str,4];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[2]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[3]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-(eig[[1]]+eig[[2]]+eig[[3]]), locPrec]],"\n"];
		,"SP6",
		    Write[str,6];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[2]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[eig[[3]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-eig[[1]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-eig[[2]], locPrec]],"\n"];
		    Write[str,0.5];
		    WriteString[str,"0 ",ToString[N[-eig[[3]], locPrec]],"\n"];
	];
	Write[str,N[QValue[Ltype,level]]];
    WriteString[str, Re[N[OMEGA,14]], " ", Im[N[OMEGA,14]], "\n"];
    Write[str,0];
	If[Ltype == "SP4" || Ltype == "SP4Std" || Ltype == "SP6" || Ltype == "HoloxMaass" || Ltype == "Maass",
	    For[i=1,i<=Length[coefList]/2,i++,
	        WriteString[str,ToString[N[Chop[coefList[[2i-1]],10^(-5)], locPrec]],"\n"];
	    ];
    ,
	    For[i=1,i<=Length[coefList]/2,i++,
	        WriteString[str,ToString[N[Chop[coefList[[2i-1]],10^(-5)], locPrec]]," ",ToString[N[Chop[coefList[[2i]],10^(-5)], locPrec]],"\n"];
	    ];
	];
	Close[str];
	locPrec
]

(* Finds which of the coefficientlists in l1 and l2 match best.
Returns a list of pairs {[i,j}, ... }, where this means that coefficients
nr i in l1 match best with coefficients j in l2.  ----------------------*)
divideSolutions[l1_,l2_,nr_:3]:=Module[{diff, diffOrder, ans, i, j, k},
	diff=Table[Norm[l1[[i,Range[1,nr]]]-l2[[j,Range[1,nr]]]]+Norm[l1[[i, Range[1+Length[l1[[i]]]/2,nr+Length[l1[[i]]]/2]]]-l2[[j,Range[1+Length[l2[[j]]]/2,nr+Length[l2[[j]]]/2]]]],{i,Length[l1]},{j,Length[l2]}];
	diffOrder=Sort[Flatten[Table[{diff[[i,j]],i,j},{i,Length[l1]},{j,Length[l2]}],1]];
	ans={};
	While[Length[diffOrder]>0,
		AppendTo[ans,diffOrder[[1,{2,3}]]];
		k=2;
		While[k<=Length[diffOrder],
			If[diffOrder[[k,2]]==diffOrder[[1,2]] ||diffOrder[[k,3]]==diffOrder[[1,3]],
				diffOrder=Delete[diffOrder,k];
			,
				k++; 
			];
		];
		diffOrder=Delete[diffOrder,1];
	];
	ans
]

(* Returns a floating point value that is value with the number of 
significant digits equal to digits. (Use mainly to get smaller datafiles
when using high precision when high precision not needed in the final answer.
 --------------------------------------------------------------------- *)
fixedNrOfDigits[value_,digits_]:=Module[{intexp},
	If[value==0,
		0
	,
	    intexp = -Round[Log[10, value]] + digits - 1;
        N[Round[value * 10^intexp]/10^intexp]
	]
]



(* Global test function.  *)
g[z_, b_, s_] := (z - s)^b[[3]]*E^( b [[1]]*(z - s)^2 -Sign[Im[s]] b[[2]] I (z - s));

getArchiveFileIndex[param_,targetFileBase_,prec_:7]:=Module[{content, Rtuple, coef, newR, error, k, ans},
	ans = 0;
	k=1;
	While[FileType[ targetFileBase <> ToString[k]<> ".txt"]==File,
		content = Get[targetFileBase <> ToString[k] <> ".txt"];
		Rtuple = content[[1]];
		If[Max[Abs[param[[1]] - Rtuple[[1]]],Abs[param[[Min[2,Length[Rtuple]]]] - Rtuple[[Min[2,Length[Rtuple]]]]]] < 10^(-prec),
			ans = k;
			Break[];
		];
		k++;
	]; 
	ans
]

getArchiveFileName[Ldata_,archiveFileType_]:=Module[{Ltype, level, OMEGA, targetFileBase, listFile, listFilePart, listFileAgain,parity},
	Ltype = Ldata[[1]];
	level = Ldata[[2]];
	OMEGA = Ldata[[4]];
	parity = Ldata[[5]];
	Switch[Ltype,
	"GL3",
		If[parity=={1},
			targetFileBase = "SL3Final/Odd/SL3Odd_Maass";
			listFile ="SL3Final/Odd/SL3Odd_All.txt";
			listFilePart ="SL3Final/Odd/SL3Odd_AllPart.txt";
			listFileAgain ="SL3Final/Odd/SL3Odd_AllAgain.txt";
		,
			If[level==1,
				targetFileBase = "SL3Final/SL3_MaassTEST";
				listFile ="SL3Final/SL3_AllTEST.txt";
				listFilePart ="SL3Final/SL3_AllPartTEST.txt";
				listFileAgain ="SL3Final/SL3_AllAgainTEST.txt";
			,
				targetFileBase = "SL3Final/Level/SL3_Level" <> ToString[level] <> "OMEGA" <> OMEGAcode[OMEGA];
				listFile ="SL3Final/Level/SL3_Level" <> ToString[level] <> "OMEGA" <> OMEGAcode[OMEGA]<>".txt";
				listFilePart ="SL3Final/Level/SL3_PartLevel" <> ToString[level] <> "OMEGA" <> OMEGAcode[OMEGA]<>".txt";
				listFileAgain ="SL3Final/Level/SL3_AgainLevel" <> ToString[level] <> "OMEGA" <> OMEGAcode[OMEGA]<>".txt";
			];
		];
	,"SP4",
		If[level==1,
			targetFileBase = "SP4Final/SP4_Maass";
			listFile ="SP4Final/SP4_All.txt";
			listFilePart ="SP4Final/SP4_AllPart.txt";
			listFileAgain ="SP4Final/SP4_AllAgain.txt";
		,
			targetFileBase = "SP4Final/Level/SP4_Level" <> ToString[level]<> "OMEGA" <> OMEGAcode[OMEGA];
			listFile ="SP4Final/Level/SP4_Level" <> ToString[level]<> "OMEGA" <> OMEGAcode[OMEGA]<>".txt";
			listFilePart ="SP4Final/Level/SP4_PartLevel" <> ToString[level]<> "OMEGA" <> OMEGAcode[OMEGA]<>".txt";
			listFileAgain ="SP4Final/Level/SP4_AgainLevel" <> ToString[level]<> "OMEGA" <> OMEGAcode[OMEGA]<>".txt";
		];
	,"SP4Std",
			targetFileBase = "SP4StdFinal/SP4Std_Maass";
			listFile ="SP4StdFinal/SP4Std_All.txt";
			listFilePart ="SP4StdFinal/SP4Std_AllPart.txt";
			listFileAgain ="SP4StdFinal/SP4Std_AllAgain.txt";
	,"HoloxMaass",
			targetFileBase = "HoloxMaassFinal/HoloxMaass_Maass";
			listFile ="HoloxMaassFinal/HoloxMaass_All.txt";
			listFilePart ="HoloxMaassFinal/HoloxMaass_AllPart.txt";
			listFileAgain ="HoloxMaassFinal/HoloxMaass_AllAgain.txt";
	,"GL3Holo",
			targetFileBase = "GL3HoloFinal/GL3Holo_Maass";
			listFile ="GL3HoloFinal/GL3Holo_All.txt";
			listFilePart ="GL3HoloFinal/GL3Holo_AllPart.txt";
			listFileAgain ="GL3HoloFinal/GL3Holo_AllAgain.txt";
	,"GL4Holo",
			targetFileBase = "GL4HoloFinal/GL4Holo_Maass";
			listFile ="GL4HoloFinal/GL4Holo_All.txt";
			listFilePart ="GL4HoloFinal/GL4Holo_AllPart.txt";
			listFileAgain ="GL4HoloFinal/GL4Holo_AllAgain.txt";
	,"GL4",
			targetFileBase = "GL4Final/GL4_Maass";
			listFile ="GL4Final/GL4_All.txt";
			listFilePart ="GL4Final/GL4_AllPart.txt";
			listFileAgain ="GL4Final/GL4_AllAgain.txt";
	,"SP6",
			targetFileBase = "SP6Final/SP6_Maass";
			listFile ="SP6Final/SP6_All.txt";
			listFilePart ="SP6Final/SP6_AllPart.txt";
			listFileAgain ="SP6Final/SP6_AllAgain.txt";
	];
	Switch[archiveFileType,
	"Base",
		targetFileBase
	,"High",
		listFile
	,"Low",
		listFilePart
	,"Again",
		listFileAgain
	, _,
		Print["Wrong argument for archiveFileType in getArchiveFileName"];
		""
	]
]

getCandidatesFromZeroLines[zeroLines_,indIdx_, sameLimit_]:=Module[{nrEnoughClose, idx, nr, k,err, i, cd, j, nrOfInd, R1zero, R2zero, coefZero, cl, tbl, closeList},
	nrOfInd = Length[zeroLines[[1,1,4]]];
	cl = Table[{},{k,nrOfInd}];
	nrEnoughClose = 2;
	For[nr=1,nr<=Length[zeroLines],nr++,
		For[idx=1,idx<=nrOfInd,idx++,
			err=Table[{zeroLines[[nr,k,1]]//N, zeroLines[[nr,k,4,idx]], zeroLines[[nr,k,2]]},{k,Length[zeroLines[[nr]]]}];
			If[Length[err]>1,
				If[err[[1,2]](2err[[1,2]]-err[[2,2]]) <= 0,
					R1zero=secantzero[err[[1,1,1]],err[[2,1,1]],err[[1,2]],err[[2,2]]];
					R2zero=secantzero[err[[1,1,2]],err[[2,1,2]],err[[1,2]],err[[2,2]]];
					coefZero = interpolCoef[R1zero, err[[All,1,1]], err[[All,3]]];
					AppendTo[cl[[idx]],{{R1zero,R2zero}, coefZero}];
				];
				For[j=1,j<Length[err],j++,
					If[err[[j,2]]err[[j+1,2]] <= 0,
						R1zero=secantzero[err[[j,1,1]],err[[j+1,1,1]],err[[j,2]],err[[j+1,2]]];
						R2zero=secantzero[err[[j,1,2]],err[[j+1,1,2]],err[[j,2]],err[[j+1,2]]];
						coefZero = interpolCoef[R1zero, err[[All,1,1]], err[[All,3]]];
						AppendTo[cl[[idx]],{{R1zero,R2zero}, coefZero}];
					];
				];
				If[err[[j,2]](2err[[j,2]]-err[[j-1,2]]) <= 0,
					R1zero=secantzero[err[[j,1,1]],err[[j-1,1,1]],err[[j,2]],err[[j-1,2]]];
					R2zero=secantzero[err[[j,1,2]],err[[j-1,1,2]],err[[j,2]],err[[j-1,2]]];
					coefZero = interpolCoef[R1zero, err[[All,1,1]], err[[All,3]]];
					AppendTo[cl[[idx]],{{R1zero,R2zero}, coefZero}];
				];
			,
				AppendTo[cl[[idx]],{err[[1,1]], err[[1,3]]}];
			];
		];
	];
	cl=Delete[cl,indIdx];
	tbl=Table[0,{i,Length[cl]},{j,Length[cl[[i]]]}];
	closeList={};
	For[i=1,i<=Length[cl],i++,
		For[j=1,j<=Length[cl[[i]]],j++,
			tbl[[i,j]]=Sort[Table[Min[Table[Norm[cl[[i,j,1]]-cl[[k,m,1]]],{m,Length[cl[[k]]]}]],{k,Length[cl]}]];
			k=Length[cl];
			While[tbl[[i,j,k]]>10 sameLimit,
				k--;
			];
Print[cl[[i,j,1]],"   ", tbl[[i,j,4]]/sameLimit, "  ", k];
			If[k>nrEnoughClose,
				AppendTo[closeList,{cl[[i,j,1]],tbl[[i,j,4]]/sameLimit,k, cl[[i,j,2]]}];
			]
		];
	];
	cd=closeList;
	i=1;
	While[i<Length[cd],
		j=i+1;
		While[j<=Length[cd],
			If[Norm[cd[[j,1]]-cd[[i,1]]] < 10 sameLimit,
				cd=Delete[cd,j];
			,
				j++;
			];
		];
		i++;
	];
	{cd, closeList, cl}
]


(* Returns the degree of the L-function of the given type. *)
getDegree[Ltype_]:=
	Switch[Ltype
	,"Maass",
		2
	,"Holomorph",
		2
	,"GL3",
		3
	,"GL3Holo",
		3
	,"GL4Holo",
		4
	,"SiegelPara",
		4
	,"SiegelSpin",
		4
	,"SP4",
		4
	,"GL4",
		4
	,"HoloxMaass",
		4
	,"SP4Std",
		5
	,"SP6",
		6
	]



(* Returns the list of deltas from the Gamma-factors to eliminate cancellation.
 See (30) in Rubinsteinn. 
 --------------------------------------------------------------------- *)
getDeltaList[s_,digits2_, klist_,llist_]:=Module[{t, j, arg, answer, c},
    answer={};
    c=digits2 * Log[10]/Sum[klist[[i]],{i,Length[klist]}];
    For[j=1,j<=Length[klist],j++,
	t=Im[klist[[j]] * s + llist[[j]]];
	If[Abs[t]<=2*c/Pi,
	   arg = Pi/2;
	,
	   arg=Abs[c/t];
	];
	AppendTo[answer, Exp[I * Sign[t] * (Pi/2-arg)]];
    ];
    answer
]




(* Returns the number of coefficients to use (NN), the list of s values (sSeq) and parameters
for the quadratic term in the test function (paraSeq) *) 
getFuncEqParameters[OMEGA_,klist_,llist_,reslist_,polelist_,Param_,v_,minWidth_:40,initNN_:30,minMult_:1/200,realOrImaginary_:0]:=Module[{expz,  incr, M, continue,  sStep,  s, n, k, matrixrows, sStart, sEnd, NN, paraList, paraIndexStart, paraIndex, paraSeq, sSeq, error, prec, thisIndex, locParam, raiseDigits, hasRaisedPara, extraPrec},
	NN=initNN;
	locParam = Param;
	Print[NN];
	paraList = {1/3,1/5,1/8,1/12,1/16,1/25,1/40,1/80,1/200};
	paraIndexStart = 6;
(*	paraList = {1/25};
	paraIndexStart = 1;*)
    sStep = 4;
	extraPrec = 2;
	sSeq = {};
	paraSeq = {};
    incr=2*Pi*v/Log[10]/DIGITS;  (* Not good if v is large *)
    M=Sqrt[Log[10]*DIGITS/Param[[1,1]]];
	expz=Table[Exp[(v+I*k)*logQlogN[[n]]]/(v+I*k),{n,1,NN},{k,-M,M,incr}];

	s=1/2 - I;
	continue = True;
	paraIndex = paraIndexStart;
	thisIndex = 0;
	hasRaisedPara = False;
	raiseDigits = 0;
	While[continue,
		locParam[[All,1]] = paraList[[paraIndex]];
		matrixrows = computeEquation[OMEGA,klist,llist,reslist,polelist,locParam,v,NN,incr, M,expz,s];
		error = Max[Abs[matrixrows[[1,{-2,-1}]]]];
		prec = Max[Abs[matrixrows[[2]]]];
(*Print["s: ", s, "  Parameter: ", paraList[[paraIndex]]];
Print["Precision: ", prec];
Print["Error: ", error];*)
		If[ prec > 10 ^(-TRUNCDIGITS - extraPrec )/2 ,
			continue = False;
			raiseDigits = 1;
		,
			If[error > 10 ^(-TRUNCDIGITS )/2,
				If[thisIndex==0,
					If[paraIndex < Length[paraList]  && paraList[[paraIndex+1]]>=minMult,
						paraIndex ++;
						hasRaisedPara = True;
					,
						continue=False;
					];
				,
					PrependTo[sSeq, Im[s]];
					PrependTo[paraSeq, paraList[[thisIndex]]];
					s -=I sStep;
					paraIndex = thisIndex;
					thisIndex = 0;
				];
			,
				thisIndex = paraIndex;
				If[error < 10 ^(-TRUNCDIGITS-1)/2 && paraIndex>1  && !hasRaisedPara,
					paraIndex --;
				,
					PrependTo[sSeq, Im[s]];
					PrependTo[paraSeq, paraList[[thisIndex]]];
					s -=I sStep;
					hasRaisedPara = False;
					thisIndex = 0;
				];
			];
		];
	];

	If[realOrImaginary==1,  (* Selfdual so only positive half *)
		continue=False;
	,
		s=1/2 + I;
		continue = True;
		paraIndex = paraIndexStart;
		hasRaisedPara = False;
		thisIndex = 0;
	];
	While[continue,
		locParam[[All,1]] = paraList[[paraIndex]];
		matrixrows = computeEquation[OMEGA,klist,llist,reslist,polelist,locParam,v,NN,incr, M,expz,s];
		error = Max[Abs[matrixrows[[1,{-2,-1}]]]];
		prec = Max[Abs[matrixrows[[2]]]];
(*Print["s: ", s, "  Parameter: ", paraList[[paraIndex]]];
Print["Precision: ", prec];
Print["Error: ", error];*)
		If[ prec > 10 ^(-TRUNCDIGITS  - extraPrec)/2 ,
			continue = False;
			raiseDigits++;
		,
			If[error > 10 ^(-TRUNCDIGITS)/2,
				If[thisIndex==0,
					If[paraIndex < Length[paraList]  && paraList[[paraIndex+1]]>=minMult,
						paraIndex ++;
						hasRaisedPara = True;
					,
						continue=False;
					];
				,
					AppendTo[sSeq, Im[s]];
					AppendTo[paraSeq, paraList[[thisIndex]]];
					s +=I sStep;
					paraIndex = thisIndex;
					thisIndex = 0;
				];
			,
				thisIndex = paraIndex;
				If[error < 10 ^(-TRUNCDIGITS-1)/2 && paraIndex>1  && !hasRaisedPara,
					paraIndex --;
				,
					AppendTo[sSeq, Im[s]];
					AppendTo[paraSeq, paraList[[thisIndex]]];
					s +=I sStep;
					hasRaisedPara = False;
					thisIndex = 0;
				];
			];
		];
	];

	If[Length[sSeq]>0,
Print["Length: ", sSeq[[-1]]-sSeq[[1]], "  Interval: ", sSeq[[1]], ",",  sSeq[[-1]]];

		If[sSeq[[-1]]-sSeq[[1]]<minWidth,
			If[raiseDigits == 2,
				DIGITS += 3;
				Print["Raising DIGITS to ",DIGITS];
			];
			{NN, sSeq, paraSeq}=getFuncEqParameters[OMEGA,klist,llist,reslist,polelist,Param,v,minWidth,NN+5,minMult,realOrImaginary];
		];
	,
		{NN, sSeq, paraSeq}=getFuncEqParameters[OMEGA,klist,llist,reslist,polelist,Param,v,minWidth,NN+5,minMult,realOrImaginary];
	];

	{NN, sSeq, paraSeq}
]


getHeckeFactor[Ldata_,coef_,n_,RealPart_,isComplexPair_:True,maxPower_:4]:=Module[{Ltype,level,ch,fi, pp,k,p,r,a,b,c,d},
	Ltype = Ldata[[1]];
	level = Ldata[[2]];
	ch = Ldata[[3]];
	If[n==1,
		If[RealPart,
			1
		,
			0
		]
	,
		fi=FactorInteger[n];
		If[Length[fi]>1,  (* At least two primes, so multiplicative *)
			pp = fi[[1,1]]^fi[[1,2]];
			If[RealPart,
				getHeckeFactor[Ldata,coef,pp,True,isComplexPair,maxPower]*getHeckeFactor[Ldata,coef,n/pp,True,isComplexPair,maxPower]-getHeckeFactor[Ldata,coef,pp,False,isComplexPair,maxPower]*getHeckeFactor[Ldata,coef,n/pp,False,isComplexPair,maxPower]
			,
				getHeckeFactor[Ldata,coef,pp,True,isComplexPair,maxPower]*getHeckeFactor[Ldata,coef,n/pp,False,isComplexPair,maxPower]+getHeckeFactor[Ldata,coef,pp,False,isComplexPair,maxPower]*getHeckeFactor[Ldata,coef,n/pp,True,isComplexPair,maxPower]
			]

		,    (* A prime power *)
			{p,r}=fi[[1]];
			If[r>maxPower,
				If[RealPart,
					coef[[2p^r-1]]
				,
					coef[[2p^r]]
				]
			,
			{a,b}=coef[[{2p-1,2p}]];
			If[Length[coef]>=2p^2,
				{c,d}=coef[[{2p^2-1,2p^2}]];
			];
			Switch[getDegree[Ltype],
			2,
				If[RealPart,
					If[GCD[p,level]>1,
							Switch[r,
								1,
									a
								,2,
									a^2-b^2
								,3,
(*Print["Is doing level: ",level];*)
									a^3 - 3a b^2
								,4,
									a^4  - 6 a^2 b^2 + b^4 
								, _,
									coef[[2p^r-1]]
								]
					,
						Switch[r,
							1,
								a
							,2,
								a^2-b^2-1
							,3,
								a^3-2a-3a b^2
							,4,
								1 -3 a^2+a^4 + 3 b^2 -6 a^2 b^2 + b^4
							, _,
								coef[[2p^r-1]]
						]
					]
				,
					If[GCD[p,level]>1,
							Switch[r,
								1,
									b
								,2,
									2 a b
								,3,
									3 a^2 b  - b^3
								,4,
									4 a^3 b - 4 a b^3 
								, _,
									coef[[2p^r]]
							]
					,
						Switch[r,
							1,
								b
							,2,
								2a b
							,3,
								-b^3+3a^2b - 2b
							,4,
								4 a^3 b - 4 a b^3 - 6 a b
							, _,
								coef[[2p^r]]
						]
					]
				]
			,3,
				If[RealPart,
					If[GCD[p,level]>1,
						If[ch=="No",
							Switch[r,
								1,
									a
								,2,
									c
								,3,
(*Print["Is doing level: ",level];*)
									2 a c -a^3 + 3a b^2 -2 b d 
								,4,
									a^2 - 2 a^4  - b^2 + 12 a^2 b^2 - 2 b^4 - c + 3 a^2 c  - 3 b^2 c  -  6 a b d 
								, _,
									coef[[2p^r-1]]
								]
						,
							Switch[r,
								1,
									a
								,2,
									-ch+a^2-b^2
								,3,
(*Print["Is doing level: ",level];*)
									-2 a ch +a^3 - 3a b^2 
								,4,
									ch^2 - 3a^2 ch+ a^4+3b^2 ch-6a^2 b^2+b^4						
								, _,
									coef[[2p^r-1]]
								]
						]
					,
						Switch[r,
							1,
								a
							,2,
								a^2-b^2-a
							,3,
								a^3-2a^2-2b^2-3a b^2+1
							,4,
								2 a+a^2-3 a^3+a^4-b^2-3 a b^2-6 a^2 b^2+b^4
							, _,
								coef[[2p^r-1]]
						]
					]
				,
					If[GCD[p,level]>1,
						If[ch=="No",
							Switch[r,
								1,
									b
								,2,
									d
								,3,
									2 b c - 3 a^2 b  + b^3 + 2 a d
								,4,
									2 a b - 8 a^3 b + 8 a b^3  + 6 a b c  - d + 3 a^2 d - 3 b^2 d
								, _,
									coef[[2p^r]]
							]
						,
							Switch[r,
								1,
									b
								,2,
									2a b
								,3,
									-2 b ch + 3 a^2 b  - b^3
								,4,
									6 a b ch+ 4 a^3 b -4 a b^3 
								, _,
									coef[[2p^r]]
							]
						]
					,
						Switch[r,
							1,
								b
							,2,
								b+2b a
							,3,
								-b^3+3a^2b
							,4,
								2 b-2a b-3a^2b+4a^3b-3b^3-4a b^3
							, _,
								coef[[2p^r]]
						]
					]
				]
			,4,
				If[RealPart,
					If[Divisible[level,p],
						coef[[2p^r-1]]
(*Print["Non-trivial level hasn't been implemented for degree 4"];*)
					,
						Switch[r,
							1,
								a
							,2,
								c
							,3,
								a-a^3+3a b^2+2a c-2b d
							,4,
								-1 + 2*a^2 - a^4  + 2*b^2 + 6*a^2*b^2 - b^4 + a^2*c  - b^2*c + c^2  -  2*a*b*d  - d^2
							, _,
								coef[[2p^r-1]]
						]
					]
				,
					If[Divisible[level,p],
						coef[[2p^r]]
(*Print["Non-trivial level hasn't been implemented for degree 4"];*)
					,
						Switch[r,
							1,
								b
							,2,
								d
							,3,
								-b-3b a^2+b^3+2b c+2a d
							,4,
								- 4*a^3*b +  4*a*b^3 + 2*a*b*c + a^2*d - b^2*d + 2*c*d 
							, _,
								coef[[2p^r]]
						]
					]
				]
			,5,
				If[RealPart,
					If[Divisible[level,p],
Print["Non-trivial level hasn't been implemented for degree 5"];
					,
						Switch[r,
							1,
								a
							,2,
								c
							,3,
								a^2-a^3-b^2+3 a b^2-c+2 a c-2 b d
							,4,
								-a+2 a^3-a^4+2 a b^2+6 a^2 b^2-b^4-2 a c+a^2 c-b^2 c+c^2-2 b d-2 a b d-d^2
							, _,
								coef[[2p^r-1]]
						]
					]
				,
					If[Divisible[level,p],
Print["Non-trivial level hasn't been implemented for degree 5"];
					,
						Switch[r,
							1,
								b
							,2,
								d
							,3,
								-2 a b-3 a^2 b+ b^3+2 b c+ d+2 a d
							,4,
								 b-2 a^2 b-4 a^3 b-2 b^3+4 a b^3-2 b c+2 a b c+2 a d+ a^2 d- b^2 d+2 c d
							, _,
								coef[[2p^r]]
						]
					]
				]
			,6,
				If[RealPart,
					coef[[2p^r-1]]
				,
					coef[[2p^r]]
				]
			]
			]
		]
	]
]


getKlist[Ltype_]:=Module[{},
	Switch[Ltype
	,"Maass",
		{1/2,1/2}
	,"Holomorph",
		{1}
	,"GL3",
		{1/2,1/2,1/2}
	,"GL3Holo",
		{1/2, 1/2, 1/2}
	,"GL4Holo",
		{1/2, 1/2, 1/2, 1/2}
	,"SiegelSpin",
		{1,1}
	,"SiegelPara",
		{1,1}
	,"SP4",
		{1/2,1/2,1/2,1/2}
	,"GL4",
		{1/2,1/2,1/2,1/2}
	,"HoloxMaass",
		{1,1/2,1/2}
	,"SP4Std",
		{1/2,1/2,1/2,1/2,1/2}
	,"SP6",
		{1/2,1/2,1/2,1/2,1/2,1/2}
	]
]

getLlist[Ltype_,param_,parity_:{}]:=Module[{},
	Switch[Ltype
	,"Maass",
		{(I*param[[1]] + parity[[1]])/2,  (-I*param[[1]]+ parity[[1]])/2}
	,"Holomorph",
		{(param[[1]]-1)/2}
	,"GL3",
		{I*param[[1]]+ parity[[1]]/2, I*param[[2]]+ parity[[1]]/2, -I*(param[[1]]+param[[2]])}
	,"GL3Holo",
		{(parity[[2]]-1)/4-I*param[[1]]/2,(parity[[2]]+1)/4-I*param[[1]]/2, I*param[[1]]+parity[[1]]/2}
	,"GL4Holo",
		{(parity[[3]]-1)/4-I*(param[[1]]+param[[2]])/2,(parity[[3]]+1)/4-I*(param[[1]]+param[[2]])/2, I*param[[1]]+parity[[1]]/2, I*param[[2]]+parity[[2]]/2}
	,"SiegelPara",
		{1/2,1/2}
	,"SiegelSpin",
		{1/2, param[[1]]-3/2}
	,"SP4",
		{I*param[[1]]+parity[[1]]/2, I*param[[2]]+parity[[2]]/2, -I*param[[1]]+parity[[1]]/2, -I*param[[2]]+parity[[2]]/2}
	,"GL4",
		{I*param[[1]], I*param[[2]], I*param[[3]], -I*(param[[1]] + param[[2]]+ param[[3]])}
	,"SP4Std",
		{parity[[3]]/2,I*param[[1]]+parity[[1]]/2, I*param[[2]]+parity[[2]]/2, -I*param[[1]]+parity[[1]]/2, -I*param[[2]]+parity[[2]]/2}
	,"HoloxMaass",
		{(parity[[2]]-1)/2, I*param[[1]]+parity[[1]]/2, -I*param[[1]]+parity[[1]]/2}
	,"SP6",
		{I*param[[1]], I*param[[2]], I*param[[3]], -I*param[[1]], -I*param[[2]], -I*param[[3]]}
	]
]

getMultipleStartValues[startValues_, nrOfRuns_, nrOfUnknowns_, realOrImaginary_:0]:=Module[{i, k, ans, locStart},
	ans = {};
	For[i=1,i<=nrOfRuns+1,i++,
			locStart = startValues;
			While[Length[locStart]/(2-Abs[realOrImaginary]) < nrOfUnknowns,
				If[realOrImaginary==0,
					k = Length[locStart]/2 +1;
					locStart=Insert[locStart, SetPrecision[2 Random[]-1, PRECISION], k];
				];
				AppendTo[locStart,SetPrecision[2 Random[]-1, PRECISION]];
			];
			AppendTo[ans, locStart];
	];
	ans
]

(* Returns a list of the number of coefficients to use to have a truncation error at most 10^k where k takes all values from minTruncation to maxTruncation *) 
getNumberOfCoefficients[OMEGA_,klist_,llist_,parity_,v_,bwidth_:2,quadfact_:0, maxTruncation_:20,minTruncation_:1, maxNN_:70]:=Module[{digits, expz, incr, M, s, n, k, matrixrows, NN, param, truncation, Nlist, significantrow, precisonrow},
	param={{quadfact,0,0},{quadfact,Min[bwidth, 3/2],0}};
	digits = 20;
    incr=2*Pi*v/Log[10]/digits;  
    If[quadfact>0,
        M=Sqrt[Log[10]*digits/quadfact];
    ,
	    M=7 Sqrt[Log[10]*digits];
    ];         
	expz=Table[Exp[(v+I*k)*logQlogN[[n]]]/(v+I*k),{n,1,maxNN},{k,-M,M,incr}];
	s=1/2 + I;
	matrixrows = computeEquation[OMEGA,klist,llist,{},{},param,v,maxNN,incr, M,expz,s];
	If[OMEGA == {-1, 0},
		significantrow = matrixrows[[2]];
		precisonrow = matrixrows[[1]];
	,
		significantrow = matrixrows[[1]];
		precisonrow = matrixrows[[2]];
	,
		significantrow = matrixrows[[1]];
		precisonrow = matrixrows[[2]];
	];
	Print["Precision: ", Max[Abs[precisonrow[[{-2,-1}]]]]];
	Print["Max truncation: ", Max[Abs[significantrow[[{-2,-1}]]]]];
	
	truncation = minTruncation;
	Nlist = {};
	NN=1;
	While[truncation<=maxTruncation && 2NN<=Length[significantrow],
		If[Max[Abs[significantrow[[{2 NN-1,2NN}]]]]<10^-truncation,
			AppendTo[Nlist,NN+1];
			truncation++;
		,
			NN ++;
		];
	];

	If[truncation <=maxTruncation,
		Print["Warning, the value of maxNN was too small. Using estimate for all truncation greater than ", truncation-1];
		While[truncation <=maxTruncation,
			AppendTo[Nlist, 2 Nlist[[-1]] - Nlist[[-2]]+1];
			truncation++;
		];
	];
	Nlist
]

getNrOfParameters[Ltype_]:=
	Switch[Ltype
	,"Maass",
		1
	,"Holomorph",
		1
	,"GL3",
		2
	,"GL3Holo",
		1
	,"GL4Holo",
		2
	,"SiegelPara",
		1
	,"SiegelSpin",
		1
	,"SP4",
		2
	,"GL4",
		3
	,"SP4Std",
		2
	,"HoloxMaass",
		1
	,"SP6",
		3
	]





getParameterListIndex[param_,list_,prec_:7]:=Module[{k, ans},
	ans = 0;
	For[k=1, k<=Length[list],k++,	
		If[Max[Abs[param[[1]] - list[[k,1]]],Abs[param[[Min[2,Length[param]]]] - list[[k,Min[2,Length[param]]]]]] < 10^(-prec),
			ans = k;
			Break[];
		];
	]; 
	ans
]

getPhaseFactor[Ltype_, param_]:=
	If[Ltype == "GL3Holo",
		1
	,
		1
	]
	
getRstep[Ltype_,step_]:= Table[step,{i,getNrOfParameters[Ltype]}]

getRtuple[Ltype_,Rtuple_, RstepList_, iter_]:=
	Switch[getNrOfParameters[Ltype]
	,1,
		{Rtuple[[1]] + iter[[1]]*RstepList[[1]]}
	,2,
		{Rtuple[[1]] + iter[[1]]*RstepList[[1]], Rtuple[[2]] + iter[[2]]*RstepList[[2]]}
	,3,
		{Rtuple[[1]] + iter[[1]]*RstepList[[1]], Rtuple[[2]] + iter[[2]]*RstepList[[2]], Rtuple[[3]] + iter[[3]]*RstepList[[3]]}
	]

getSideCounts[sideLengths_, RstepList_]:=Module[{i,answer},
	answer = Table[0,{i,Length[sideLengths]}];
	For[i=1,i<=Length[sideLengths],i++,
		If[sideLengths[[i]]>0,
			answer[[i]] = Ceiling[sideLengths[[i]]/RstepList[[i]]]+2;
		];
	];
	answer
]


getSideLengths[Ltype_,length_]:=
	Switch[getNrOfParameters[Ltype]
	,1,
		{length, 0 , 0}
	,2,
		{length, length, 0}
	,3,
		{length, length, length}
	]



(* Returns a list of nrOfEquations equally spaced s-values in the intervals specified by
sSeq and the corresponding quadratic multiplier specified in paraSeq. 
------------------------------------------------------------------------*)
getSlist[sSeq_, paraSeq_, nrOfEquations_, SDIFF_:1/3, realOrImaginary_:0]:=Module[{sStep, slist, paralist, idx, i, j, parts, sStart, sEnd},
	(* Set the proportion of points to different parts of the line.
		Each of the parts get 1/Length[parts] of the s values. *)
	If[realOrImaginary==0,
		parts = {15/100,20/100,30/100,20/100,15/100};
	,
		parts = {41/100,33/100,26/100};
	];

	(* Get the list of s values to use. *)
	slist = {};
	For[i=1,i<=Length[parts],i++,
		sStep= Length[parts]*parts[[i]]*(sSeq[[-1]]-sSeq[[1]])/nrOfEquations;
(*Print["sStep: ", N[sStep]];*)
		sStart = sSeq[[1]]-SDIFF + Sum[parts[[j]],{j,i-1}] * (sSeq[[-1]]-sSeq[[1]]);
		sEnd = sSeq[[1]]-SDIFF + Sum[parts[[j]],{j,i}] * (sSeq[[-1]]-sSeq[[1]]);
		If[i<= Length[parts]-Mod[nrOfEquations, Length[parts]],
			sEnd -= sStep;
		];
		slist = Flatten[AppendTo[slist, Table[1/2+I j,{j,sStart, sEnd, sStep}]]];
	];
Print["nrOfEquations ", nrOfEquations, " nrOfSvalues ", Length[slist]];
	
	(* Make sure the number of s values is the correct one. 
		Shouldn't ever occur.*)
	While[Length[slist] < nrOfEquations,
		AppendTo[slist, (slist[[-(nrOfEquations - Length[slist])-1]] + slist[[-(nrOfEquations - Length[slist])-1]])/2];
	];
	slist = slist[[Range[1,nrOfEquations]]];

	(* Move any point that is too close to the central point. *)
	For[i=1,i<=nrOfEquations,i++,
		If[Abs[Im[slist[[i]]]] < 1/10,
			slist[[i]] = 1/2 + I Sign[Im[slist[[i]]]]/10;
		];
	];

	(*Get the quadratic multiplier in the test function to use for each s value. *)
	paralist = {};
	idx = 1;
	For[i=1,i<=nrOfEquations,i++,
		While[idx<Length[sSeq] && Abs[Im[slist[[i]]]-sSeq[[idx+1]]]<Abs[Im[slist[[i]]]-sSeq[[idx]]],
			idx++;
		];
		AppendTo[paralist, paraSeq[[idx]]];
	];
	
	{slist, paralist}
]


getSlistDavid[nrOfEquations_,width_:2]:=Module[{up, slist, down, i, blist,continue},
	slist = {1/2+I};
	up=1;
	down=2;
	blist={{0,0,0}};
	For[i=1,i<=nrOfEquations,i++,
		AppendTo[blist,{0,up/down,0}];
		continue=True;
		While[continue,
			up ++;
			If[up>=width down,
				up=1;
				down++;
				continue=False;
			,
				If[GCD[up,down]==1,
					continue=False;
				];
			];
		];
	];
	{slist,{0}, blist}
]


getUnknowns[Ldata_,NN_,maxPower_,knownCoef_:{}]:=Module[{Ltype,plist, highplist, nonplist, knownlist, p, fi, charvalue, level},
	Ltype = Ldata[[1]];
	level = Ldata[[2]];
	charvalue = Ldata[[3]];
	plist={};
	nonplist={};
	highplist={};
	knownlist={};
	p=1;
	While[p<=NN,
		If[MemberQ[knownCoef[[All,1]],p],
			AppendTo[knownlist,p];
		,
			fi=FactorInteger[p];
			If[(Length[fi]==1 && GCD[level,p]>1 &&  getDegree[Ltype]>3) || PrimeQ[p] || (PrimeQ[Sqrt[p]] && (getDegree[Ltype]>3 || (GCD[level,p]>1 && charvalue =="No"))) || (PrimeQ[p^(1/3)] && Ltype=="SP6") ,
				AppendTo[plist,p];
			,
				If[Length[fi]==1 && Max[fi[[All,2]]]>maxPower,
					AppendTo[highplist,p];
				,
					AppendTo[nonplist,p];
				];
			];
		];
		p++;
	];
	{plist, highplist, nonplist, knownlist}
]


interpolCoef[x_,xlist_,coeflist_]:=Module[{flist, intOrder},
	intOrder = Min[2, Length[xlist]-1];
	flist=Table[Interpolation[Table[{xlist[[i]],coeflist[[i,k]]},{i,Length[xlist]}],InterpolationOrder-> intOrder],{k,Length[coeflist[[1]]]}];
	Table [flist[[k]][x],{k,Length[flist]}]
];

isInParameterList[param_,list_,prec_:7]:= getParameterListIndex[param,list,prec]>0

minIndex[result_]:=Module[{k, errorSize},
	If[Length[result[[3]]]>0,
		errorSize = Table[Norm[result[[3,k]]] ,{k,Length[result[[3]]]}];
		Ordering[errorSize][[1]]
	,
		0
	]
]

minDistIndex[aplist_,res_,nr_]:=Module[{distSize, k},
	If[Length[aplist]>0 && Length[res[[1]]]>0,
		distSize = Table[Norm[aplist[[Range[1,nr]]]-res[[1,j]][[Range[1,nr]]]]+Norm[aplist[[Range[1+Floor[Length[aplist]/2],nr+Floor[Length[aplist]/2]]]]-res[[1,j]][[Range[1+Floor[Length[aplist]/2],nr+Floor[Length[aplist]/2]]]]],{j,Length[res[[1]]]}];
		Ordering[distSize][[1]]
	,
		0
	]
]

(* J_ij = f_i (x+delta_j *)
newApprox[Rlist_,flist_,coeflist_:{}]:=Module[{delta,i,J,idx,f0, dim,R1, coef1},
	dim=Length[Rlist]-1;
(*Print["flist: ",flist];*)
	delta=Table[0,{i, dim}];
	idx=Table[0,{i, dim+1}];
	For[i=2,i<=dim+1,i++,
		idx[[i]]=Ordering[Table[Abs[Rlist[[1,j]]-Rlist[[i,j]]],{j,Length[Rlist[[1]]]}]][[-1]];
		delta[[idx[[i]]]]=Rlist[[i,idx[[i]]]]-Rlist[[1,idx[[i]]]];
	];
	f0 = Table[flist[[1,idx[[j]]]],{j,2,dim+1}];
(*Print["f0: ",f0];*)
	J=Check[Table[(flist[[i+1,idx[[j]]]]-flist[[1,idx[[j]]]])/delta[[idx[[i+1]]]],{j,2,dim+1},{i,dim}],IdentityMatrix[dim]];
	R1=Rlist[[1]]-PseudoInverse[J] . f0;
	If[Length[coeflist]>0,
		coef1=Check[Table[coeflist[[1,k]] +Sum[( coeflist[[i,k]]-coeflist[[1,k]])((R1-Rlist[[1]])[[idx[[i]]]])/delta[[idx[[i]]]],{i,2,dim+1}],{k,Length[coeflist[[1]]]}], coeflist[[1]]];
	,
		coef1 = {};
	];
	{ R1, coef1}
]

OMEGAcode[OMEGA_]:=Module[{},
	Switch[OMEGA,
		Cos[(2 \[Pi])/3]+I Sin[(2 \[Pi])/3],
			"root31"
		,Cos[(4 \[Pi])/3]+I Sin[(4 \[Pi])/3],
			"root32"
		, _,
			ToString[OMEGA]
	]
]


prime2All[Ldata_,pcoef_,maxPower_:4,complexCoef_:True,maxNrOfCoef_:200]:=Module[{coef, k, continue, plist, highplist, nonplist, knowns, knownlist},
	coef=pcoef;
	continue = True;
	k=1;
	{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata,maxNrOfCoef,maxPower];
	knowns = Union[plist, highplist];
	While[continue,
		If[MemberQ[knowns,k],
			If[k>Length[coef]/2,
				continue = False;
			];
		,
			If[complexCoef,
				coef = Insert[coef, getHeckeFactor[Ldata,coef,k,True,True,maxPower], 2k-1];
				coef = Insert[coef, getHeckeFactor[Ldata,coef,k,False,True,maxPower], 2k];
			,
				(* This hasn't been implemented yet
					coef = Insert[coef, getHeckeFactor[Ldata,coef,k,True,True,maxPower], k];*)
				Print["WARNING: This hasn't been implemented yet in prime2All"];
			];
		];
		k++;
	];
	coef
]

QValue[Ltype_,level_]:=Module[{},
	Switch[Ltype
	,"Maass",
		Sqrt[level]/Pi
	,"Holomorph",
		Sqrt[level]/2/Pi
	,"GL3",
		Sqrt[level]/Pi^(3/2)
	,"GL3Holo",
		Sqrt[level]/Pi^(3/2)
	,"GL4Holo",
		Sqrt[level]/Pi^2
	,"SiegelPara",
		Sqrt[level]/4/Pi^2
	,"SiegelSpin",
		Sqrt[level]/4/Pi^2
	,"SP4",
		Sqrt[level]/Pi^2
	,"HoloxMaass",
		Sqrt[level]/2/Pi^2
	,"GL4",
		Sqrt[level]/Pi^2
	,"SP4Std",
		Sqrt[level]/Pi^(5/2)
	,"SP6",
		Sqrt[level]/Pi^3
	]
]

removePrimeCoefRubbish[pcoef_, limit_:3, realOrImaginary_:0]:=Module[{loc},
	loc=pcoef;
	If[realOrImaginary==0,
		While[Length[loc]>0 && Max[Abs[loc[[-1]]] , Abs[loc[[Length[loc]/2]]]]>limit,
			loc  = Delete[loc,Length[loc]/2];
			loc  = Delete[loc,-1];
		];
	,
		While[Length[loc]>0 && Abs[loc[[-1]]]>limit,
			loc  = Delete[loc,-1];
		];
	];
	loc
]

removeCoefRubbish[Ldata_,coef_, maxPower_:4, maxNrOfCoef_:200, limit_:4,realOrImaginary_:0]:=Module[{unknowns, loc, plist, highplist, nonplist, knownlist, k},
	loc=coef;
	{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata,maxNrOfCoef,maxPower];
	unknowns = Union[plist, highplist];
	For[k=3-Abs[realOrImaginary],k<=Length[coef]/(2-Abs[realOrImaginary]),k++,
		If[realOrImaginary==0,
			If[PrimeQ[k] && Max[Abs[coef[[{2k-1,2k}]]]]>limit,
				loc=loc[[Range[1,2(k-1)]]];
				Break[];
			];
		,
			If[PrimeQ[k] && Abs[coef[[k]]]>limit,
				loc=loc[[Range[1, k-1]]];
				Break[];
			];
		];
	];
	loc
]

saveCandidate[fileName_,saveData_]:=Module[{},
    dataList=Get[fileName];
    If[dataList==$Failed,
        dataList={};
    ];
    AppendTo[dataList,saveData];
    Put[dataList,fileName];
]

saveFinal[directory_:"Upload"]:=Module[{sourceFiles,Rtuple, coef, error, pLI, finalLimit, oKLimit, prec, idx, dataListAgain, listFileAgain, dataListPart, listFilePart, saveData, dataList, listFile, k, targetFileBase,  pcoef, len,  res, parity},
	sourceFiles = FileNames["*XXZoom*",directory];
	Get[sourceFiles[[1]]];  (* Get Ldata *)
	finalLimit = 10^(-10)/2;
	oKLimit = 10^(-7)/2;
	
	If[Length[Ldata]>4,
		parity = Ldata[[5]];
	,
		parity = {};
	];
	targetFileBase = getArchiveFileName[Ldata, "Base"];
	listFile = getArchiveFileName[Ldata, "High"];
	listFilePart = getArchiveFileName[Ldata, "Low"];
	listFileAgain = getArchiveFileName[Ldata, "Again"];

		dataList=Get[listFile];
        If[dataList==$Failed,
             dataList={};
		,
			Print["SHOULDN'T HAPPEN"];
		,
			Print["Backup"];
			Save[listFile <> "backup", dataList];
	    ];
	    dataListPart=Get[listFilePart];
	    If[dataListPart==$Failed,
	        dataListPart={};
		,
			Print["SHOULDN'T HAPPEN"];
		,
			Save[listFilePart <> "backup", dataListPart];
	    ];
	    dataListAgain=Get[listFileAgain];
	    If[dataListAgain==$Failed,
	        dataListAgain={};
		,
			Print["SHOULDN'T HAPPEN"];
		,
			Save[listFileAgain <> "backup", dataListAgain];
	    ];
	For[idx=1, idx<=Length[sourceFiles],idx++,
		If[FileType[sourceFiles[[idx]]]==File,
Print[sourceFiles[[idx]]];
			Get[sourceFiles[[idx]]];
			If[Length[Ldata]>4,
				parity = Ldata[[5]];
Print["parity: ", parity];
			];
			Rtuple = Rlist[[-1]];
			res = result[Rtuple];
			If[Length[res[[1]]]<1,
				len = 0;
			,
				len= Length[res[[1,1]]];
			];
			If[Abs[realOrImaginary]==1,
Print["Real coefficients!"];
				For[i=1,i<=len,i++,
					AppendTo[res[[1,1]],0];
				];
				len = 2*len;
			];
			If[len==0,   (* No solution. *)
				coef = {};
			,
				pcoef=alternateLists[res[[1,1, Range[1, len/2]]], res[[1,1, Range[len/2+1, len]]]];
				coef =prime2All[Ldata,pcoef, maxPower];
				error = res[[3,1]];
			];
			prec = Max[ Table[ Max[ Table[ Abs[newR[[i,k]] - newR[[j,k]]],{k,Min[2,Length[newR[[i]]]]}]] ,{i,Length[newR]-1}, {j,i+1,Length[newR]}]];
			If[prec > oKLimit || prec == -Infinity || len==0,
				Print["Not enough precision (", N[prec,2], "): ", Rtuple];
				If[! (isInParameterList[Rtuple, dataList[[All,1]]] || isInParameterList[Rtuple, dataListPart[[All,1]]] ),
					pLI = getParameterListIndex[Rtuple, dataListAgain[[All,1]]];
					If[pLI==0 || dataListAgain[[pLI, 3]]>prec,
						saveData={Rtuple, coef, prec, parity};
						AppendTo[dataListAgain, saveData];
						Print["Append to again"];
						If[pLI>0,
							dataListAgain = Delete[dataListAgain, pLI];
							Print["Delete from again"];
						];
					];
				];
			,
				If[isInParameterList[Rtuple, dataList[[All,1]]] && prec > finalLimit,
					Print["Already in final and now not enough precision (", N[prec,2], "): ", Rtuple];
					pLI = getParameterListIndex[Rtuple, dataListAgain[[All,1]]];
					If[pLI>0,
						dataListAgain = Delete[dataListAgain, pLI];
						Print["Delete from again: ", Rtuple];
					];
				,
					If[isInParameterList[Rtuple, dataList[[All,1]]] || isInParameterList[Rtuple, dataListPart[[All,1]]] ,
						k = getArchiveFileIndex[Rtuple, targetFileBase];
					,
						k=1;
						While[FileType[ targetFileBase <> ToString[k]<> ".txt"]==File,
							k++;
						];						
					];
					saveData = {Rtuple, coef, newR, error, parity};
					PutAppend[saveData, targetFileBase <> ToString[k] <> ".txt"];
					pLI = getParameterListIndex[Rtuple, dataListAgain[[All,1]]];
					If[pLI>0,
						dataListAgain = Delete[dataListAgain, pLI];
						Print["From again to better! (", N[prec,2], "): ", Rtuple];
					];
					If[prec < finalLimit,
						pLI = getParameterListIndex[Rtuple, dataListPart[[All,1]]];
						If[pLI>0,
							dataListPart = Delete[dataListPart, pLI];
							AppendTo[dataList,saveData];
							Print["From part to final (", N[prec,2], "): ", Rtuple];
						,
							pLI = getParameterListIndex[Rtuple, dataList[[All,1]]];
							If[pLI>0,
								If[Max[Abs[dataList[[pLI, 3,1]]-dataList[[pLI, 3,2]]]]>Max[Abs[saveData[[ 3,1]]-saveData[[ 3,2]]]],
									dataList[[pLI]] = saveData;
									Print["Updated final (", N[prec,2], "): ", Rtuple];
								,
									Print["Already better precision in final (", N[prec,2], "): ", Rtuple];
								];
							,
								AppendTo[dataList,saveData];
								Print["Appended to final (", N[prec,2], "): ", Rtuple];
							];
						];

					,
						pLI = getParameterListIndex[Rtuple, dataListPart[[All,1]]];
(*Print["PLI: ", pLI,saveData[[3]]];*)
						If[pLI>0,
							If[Max[Abs[dataListPart[[pLI, 3,1]]-dataListPart[[pLI, 3,2]]]]>Max[Abs[saveData[[ 3,1]]-saveData[[ 3,2]]]],
								dataListPart[[pLI]] = saveData;
								Print["Updated part (", N[prec,2], "): ", Rtuple];
							,
								Print["Already better precision in part (", N[prec,2], "): ", Rtuple];
							];
						,
							AppendTo[dataListPart,saveData];
							Print["Appended to part (", N[prec,2], "): ", Rtuple];
						];
					];
				];
			];
		];
	];
    Put[dataList,listFile];
    Put[dataListPart,listFilePart];
    Put[dataListAgain,listFileAgain];
]

(* Gives the x-value where the line through {x1,y1} and {x2,y2} crosses
the x-axis.
-------------------------------------------------------------------*)
secantzero[x1_,x2_,y1_,y2_]:=(x1 y2-x2 y1)/(y2-y1)

solveForOneNL[Ldata_,klist_,llist_,reslist_,polelist_,phaseFactor_,slist_, paralist_,Param_,nrOfRuns_,NN_,M_,incr_,v_,expz_,startValues_:{},knownCoef_:{},maxPower_:4,NLmethod_:"Secant",MaxSolutions_:Infinity, startValuePrec_:0,realOrImaginary_:0, extraEq_:{}]:=Module[{errorSize,allSol, aplist,errorList,answer},
    fullmatrix = computeEquationMatrix[Ldata,klist,llist,reslist,polelist,phaseFactor,slist, paralist,Param,NN,M,incr,v,expz,realOrImaginary];
Print["Time: ", TimeUsed[]-st];
	{allSol, aplist,errorList}=solveNL[Ldata,fullmatrix,nrOfRuns,startValues,knownCoef,maxPower,NLmethod,MaxSolutions, startValuePrec,realOrImaginary,extraEq];
Print["Time: ", TimeUsed[]-st];
	errorSize = Table[Norm[errorList[[k]]],{k,Length[errorList]}];
    Print[ " Rtuple: ", N[Im[llist]], "  ap:  ", N[aplist,4],  "  Errorsize:  ", N[errorSize,3]];

   {allSol, aplist,errorList}
]

solveForOneByStep[Ldata_,klist_,llist_,reslist_,polelist_,phaseFactor_,slist_, paralist_,Param_,nrOfRuns_,NN_,M_,incr_,v_,expz_,startValues_:{},knownCoef_:{},maxPower_:4,NLmethod_:"Secant",MaxSolutions_:Infinity, startValuePrec_:0,unknownsAtStart_:0,realOrImaginary_:0, extraEq_:{}]:=Module[{errorSize,allSol, aplist,errorList,answer, unknownsPresent, nrOfCoef, nrOfEquations, plist, highplist, nonplist, knownlist, stillUnknowns, locstartValues},
    fullmatrix = computeEquationMatrix[Ldata,klist,llist,reslist,polelist,slist, paralist,Param,NN,M,incr,v,expz,realOrImaginary];
	locstartValues = startValues;
	nrOfCoef=Length[fullmatrix[[1]]]/(2-Abs[realOrImaginary]);
	nrOfEquations=Length[fullmatrix]/(2-Abs[realOrImaginary]);
	{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata,nrOfCoef,maxPower,knownCoef];
	stillUnknowns=Union[plist,highplist];
	If[unknownsAtStart==0,
		unknownsPresent = Length[stillUnknowns];
	,
		unknownsPresent = unknownsAtStart;
	];
(*Print[unknownsPresent, " ",stillUnknowns];*)
	While[unknownsPresent<=Length[stillUnknowns],
		If[unknownsPresent==Length[stillUnknowns],
			{allSol, aplist,errorList}=solveNL[Ldata,fullmatrix ,nrOfRuns,locstartValues,knownCoef,maxPower,NLmethod,MaxSolutions, startValuePrec,realOrImaginary,extraEq];
		,
			{allSol, aplist,errorList}=solveNL[Ldata,fullmatrix[[Range[1, (4- 3 Abs[realOrImaginary]) unknownsPresent],Range[1,(2-Abs[realOrImaginary])*stillUnknowns[[unknownsPresent]]]]] ,nrOfRuns,locstartValues,knownCoef,maxPower,NLmethod,MaxSolutions, startValuePrec, realOrImaginary, extraEq];
		];
(*Print["Time: ", TimeUsed[]-st];*)
		errorSize = Table[Norm[errorList[[k]]],{k,Length[errorList]}];
        Print[ " Unknowns: ", stillUnknowns[[Range[1,unknownsPresent]]], "  ap:  ", N[aplist,4]];
		locstartValues = {removePrimeCoefRubbish[allSol[[1]],4,realOrImaginary]};
(*Print[N[locstartValues,5]];*)
		unknownsPresent += 1;
	];

   {allSol, aplist,errorList}
]

solveNL[Ldata_,fullmatrix_,nrOfRuns_,startValues_:{},knownCoef_:{},maxPower_:4,NLmethod_:"Secant",MaxSolutions_:Infinity, startValuePrec_:0,realOrImaginary_:0, extraEq_:{}]:=Module[{nrOfSignParam, stepVector, interpolLine, maxStep, locStartValues, nrOfSolutions, locStart, nrOfEquations,eqnsInCheck,allEqns,errorList,checkEqns,eqnsInSystem,stepFactor,normError,nr,start,sol,aplist,allSol,j,k,eqns,stillUnknowns,heckeCond,highplist,nonplist,plist,knownlist,unknowns,nrOfCoef,realUnknowns,imagUnknowns,OMEGA,extraGoal,idx,place,allUnknowns,placesToRemove,pos,myStartValues},
	OMEGA = Ldata[[4]];
	If[NumberQ[N[OMEGA[[1]]]] && NumberQ[N[OMEGA[[2]]]],
		nrOfSignParam = 0;
	,
		nrOfSignParam=2;
	];
	myStartValues = startValues;
	nrOfCoef=Length[fullmatrix[[1]]]/(2-Abs[realOrImaginary]);
	nrOfEquations=Length[fullmatrix]/(2-Abs[realOrImaginary]);
	{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata,nrOfCoef,maxPower, {}];
	allUnknowns=Union[plist,highplist];
	{plist, highplist, nonplist, knownlist} = getUnknowns[Ldata,nrOfCoef,maxPower, knownCoef];
	stillUnknowns=Union[plist,highplist];
	stepFactor = nrOfEquations / ((2-Abs[realOrImaginary]) Length[stillUnknowns] - Length[extraEq] + nrOfSignParam);
	extraGoal = 2;
	Switch[realOrImaginary,
	0,
		realUnknowns=Table[bb1[i],{i,nrOfCoef}];
		imagUnknowns=Table[bb2[i],{i,nrOfCoef}];
		unknowns=alternateLists[realUnknowns,imagUnknowns];
		heckeCond=Union[Table[realUnknowns[[nonplist[[k]]]]->getHeckeFactor[Ldata,unknowns,nonplist[[k]],True,True,maxPower],{k,Length[nonplist]}],Table[imagUnknowns[[nonplist[[k]]]]->getHeckeFactor[Ldata,unknowns,nonplist[[k]],False,True,maxPower],{k,Length[nonplist]}]];

		(* Add the known coefficients to the list of conditions to apply. *)
		For[k=1,k<=Length[knownlist],k++,
			idx = Position[knownCoef[[All,1]],knownlist[[k]]][[1,1]];
			AppendTo[heckeCond, realUnknowns[[knownlist[[k]]]]->Re[knownCoef[[idx,2]]]];
			AppendTo[heckeCond, imagUnknowns[[knownlist[[k]]]]->Im[knownCoef[[idx,2]]]];
		];
		(* Remove the known coefficients from the startvalues. *)
		For[j=1,j<=Length[myStartValues],j++,
(*Print[myStartValues[[j]]];*)
			placesToRemove={};
			For[k=1,k<=Length[knownlist],k++,
				pos = Position[allUnknowns, knownlist[[k]]][[1,1]];
				If[pos < Length[myStartValues[[j]]]/2,
					AppendTo[placesToRemove, {pos}];
					AppendTo[placesToRemove, {pos + Length[myStartValues[[j]]]/2}];
				];
			];
(*Print[placesToRemove];*)
			If[Length[placesToRemove]>0,
				myStartValues[[j]] = Delete[myStartValues[[j]],placesToRemove];
(*Print[myStartValues[[j]]];*)
			];
		];
	,1,
		unknowns=Table[bb1[i],{i,nrOfCoef}];
		imagUnknowns=Table[0,{i,nrOfCoef}];
		heckeCond=Table[unknowns[[nonplist[[k]]]]->getHeckeFactor[Ldata,alternateLists[unknowns,imagUnknowns],nonplist[[k]],True,True,maxPower],{k,Length[nonplist]}];

		(* Add the known coefficients to the list of conditions to apply. *)
		For[k=1,k<=Length[knownlist],k++,
			idx = Position[knownCoef[[All,1]],knownlist[[k]]][[1,1]];
			AppendTo[heckeCond, unknowns[[knownlist[[k]]]]->Re[knownCoef[[idx,2]]]];
		];
		(* Remove the known coefficients from the startvalues. *)
		For[j=1,j<=Length[myStartValues],j++,
(*Print[myStartValues[[j]]];*)
			placesToRemove={};
			For[k=1,k<=Length[knownlist],k++,
				pos = Position[allUnknowns, knownlist[[k]]][[1,1]];
				If[pos < Length[myStartValues[[j]]],
					AppendTo[placesToRemove, {pos}];
				];
			];
(*Print[placesToRemove];*)
			If[Length[placesToRemove]>0,
				myStartValues[[j]] = Delete[myStartValues[[j]],placesToRemove];
(*Print[myStartValues[[j]]];*)
			];
		];
	];
(*Print[heckeCond];*)
	If[OMEGA == {-1, 0}, 
		allEqns=ReplaceRepeated[fullmatrix[[Table[(2-Abs[realOrImaginary]) * k,{k,nrOfEquations}],Range[1,(2-Abs[realOrImaginary]) nrOfCoef]]] . unknowns , heckeCond];
	,
		allEqns=ReplaceRepeated[fullmatrix[[Table[(2-Abs[realOrImaginary]) * k-1 + Abs[realOrImaginary],{k,nrOfEquations}],Range[1,(2-Abs[realOrImaginary]) nrOfCoef]]] . unknowns , heckeCond];
	,
		allEqns=ReplaceRepeated[fullmatrix[[Table[(2-Abs[realOrImaginary]) * k-1 + Abs[realOrImaginary],{k,nrOfEquations}],Range[1,(2-Abs[realOrImaginary]) nrOfCoef]]] . unknowns , heckeCond];
	];
Print[Dimensions[allEqns]];
	eqnsInSystem = Table[Floor[k],{k,1,nrOfEquations,stepFactor}];
Print[eqnsInSystem, " ",stepFactor, " ", nrOfEquations];
	eqnsInCheck = Complement[Table[k,{k,nrOfEquations}],eqnsInSystem];
	eqns=allEqns[[eqnsInSystem]];
	eqns = Join[eqns, extraEq];
	checkEqns = allEqns[[eqnsInCheck]];
	nrOfSolutions = 0; 
	allSol={};
	aplist={};
	errorList={};
	maxStep = 1;
	If[MaxSolutions == 1 && Length[myStartValues] == 1,
		locStartValues = getMultipleStartValues[myStartValues[[1]], nrOfRuns, Length[stillUnknowns],realOrImaginary];
	,
		locStartValues = myStartValues;
	];
	If[startValuePrec>0 && Length[myStartValues]>0,
		interpolLine = Interpolation[{{1,startValuePrec},{Max[2,Length[myStartValues[[1]]]/(2-Abs[realOrImaginary])], -Log[10,maxStep]}},InterpolationOrder->1];
		stepVector =Table[10^(Ceiling[-interpolLine[k]]),{k,Length[myStartValues[[1]]]/(2-Abs[realOrImaginary])}];
		While[Length[stepVector]<Length[stillUnknowns],
			AppendTo[stepVector,maxStep];
		];
	];
	For[nr=1,nr<=nrOfRuns + Length[myStartValues],nr++,
(*Print["L\[ODoubleDot]sning: ", nr];*)
		If[nr<=Length[locStartValues],
			locStart = locStartValues[[nr]];
			While[Length[locStart]/(2-Abs[realOrImaginary])< Length[stillUnknowns],
				If[realOrImaginary==0,
					k = Length[locStart]/2 +1;
					locStart=Insert[locStart, 0, k];
				];
				AppendTo[locStart,0];
			];
			If[startValuePrec>0,
				If[realOrImaginary==0,
					start=Union[Table[{realUnknowns[[stillUnknowns[[k]]]],locStart[[k]] - stepVector[[k]]/2 ,locStart[[k]] + stepVector[[k]]/2 },{k,Length[stillUnknowns]}],Table[{imagUnknowns[[stillUnknowns[[k]]]],locStart[[k+Length[stillUnknowns]]]- stepVector[[k]]/2 ,locStart[[k+Length[stillUnknowns]]]+ stepVector[[k]]/2 },{k,Length[stillUnknowns]}]];
				,
					start=Table[{unknowns[[stillUnknowns[[k]]]],locStart[[k]] - stepVector[[k]]/2 ,locStart[[k]] + stepVector[[k]]/2 },{k,Length[stillUnknowns]}];
				];		
			,
				If[realOrImaginary==0,
					start=Union[Table[{realUnknowns[[stillUnknowns[[k]]]],locStart[[k]]  },{k,Length[stillUnknowns]}],Table[{imagUnknowns[[stillUnknowns[[k]]]],locStart[[k+Length[stillUnknowns]]] },{k,Length[stillUnknowns]}]];
				,
					start=Table[{unknowns[[stillUnknowns[[k]]]],locStart[[k]]  },{k,Length[stillUnknowns]}];
				];
			];
		,
			If[realOrImaginary==0,
				start=Union[Table[{realUnknowns[[stillUnknowns[[k]]]],SetPrecision[2 Random[]-1, PRECISION]},{k,Length[stillUnknowns]}],Table[{imagUnknowns[[stillUnknowns[[k]]]],SetPrecision[2 Random[]-1, PRECISION]},{k,Length[stillUnknowns]}]];
			,
				start=Table[{unknowns[[stillUnknowns[[k]]]],SetPrecision[2 Random[]-1, PRECISION]},{k,Length[stillUnknowns]}];
			];
		];
		iter=0;
		eval=0;
		If[nrOfSignParam == 2,
			If[Length[start[[1]]]==2,
				AppendTo[start, {OMEGA[[1]], 1}];
				AppendTo[start, {OMEGA[[2]], 0}];
			,
				AppendTo[start, {OMEGA[[1]], -1, 1}];
				AppendTo[start, {OMEGA[[2]], -1, 1}];
			];
		];
(*Print[start];
Print["Still unknowns: ",stillUnknowns];
Print["Antal eqn: ",Length[eqns],"  ",eqns[[-1]]];*)
		sol=FindRoot[eqns, start,  WorkingPrecision->PRECISION, AccuracyGoal->TRUNCDIGITS + extraGoal, PrecisionGoal->TRUNCDIGITS + extraGoal, Method -> NLmethod, StepMonitor:>iter++,EvaluationMonitor:>eval++];
(*Print["Iterationer: ", iter, " Evaluations: ", eval];*)
		normError=Norm[eqns /. sol];
		If[normError<10^(-TRUNCDIGITS+2),	
			nrOfSolutions++;
			If[nrOfSolutions>=MaxSolutions,	
				nr=Infinity;
			];
			For[k=1,k<=Length[allSol],k++,
				If[Max[Abs[allSol[[k,1]]-sol[[1,2]]],Abs[allSol[[k,2]]-sol[[2,2]]]]<10^(-TRUNCDIGITS+3),
					aplist[[k,2]]++;
					k=Infinity;
				];
			];
			If[k!=Infinity,  (* New solution *)
				(* Add the known coefficients to solution *)
				For[k=1,k<=Length[knownlist],k++,
					idx = Position[knownCoef[[All,1]],knownlist[[k]]][[1,1]];
(*Print["index: ", idx];*)
					place = 1;
					While[place<=Length[stillUnknowns],
						If[knownlist[[k]]<stillUnknowns[[place]],
(*Print["sol: ", sol];*)
							If[realOrImaginary==0,
								sol = Insert[sol, realUnknowns[[knownlist[[k]]]]->Re[knownCoef[[idx,2]]],place];
								stillUnknowns = Insert[stillUnknowns, knownlist[[k]],place];
								sol = Insert[sol, imagUnknowns[[knownlist[[k]]]]->Im[knownCoef[[idx,2]]],place + (Length[sol]+1)/2];
							,
								sol = Insert[sol, unknowns[[knownlist[[k]]]]->Re[knownCoef[[idx,2]]],place];
								stillUnknowns = Insert[stillUnknowns, knownlist[[k]],place];
							];
(*Print["sol: ", sol];*)
							place = Infinity;
						,
							place++;
						];
					]
				];

				AppendTo[allSol,sol[[All,2]]];
				AppendTo[aplist,{{sol[[1,2]],sol[[2,2]]},1}];
				AppendTo[errorList,checkEqns /. sol];
			];
		];
	];
(*Print["Prec errorList: ", Precision[errorList]];
Print["Prec checkEqns: ", Precision[checkEqns]];*)
	SetPrecision[errorList,PRECISION];
	{allSol, aplist,errorList}
]

(* Returns x as a string truncated to n decimal places. *)
truncate[x_,n_]:= Module[{real,imag},
	real=Re[x];
	imag = Im[x];
	If[imag==0,
		ToString[N[IntegerPart[(10^n)*real]/(10^n),Ceiling[Log[10,Abs[real]]]+n]]
	,
		ToString[N[IntegerPart[(10^n)*real]/(10^n),Ceiling[Log[10,Abs[real]]]+n]+N[IntegerPart[(10^n)*imag]/(10^n),Ceiling[Log[10,Abs[imag]]]+n]*I]
	]
]

writeAllToLcalc[Ldata_]:=Module[{locLdata, targetFolder, eigenvalueFile, lcalcBase, Ltype, level, OMEGA, parity, nrOfPrimes, coeflist, streamList, high, low, i, j, eig, coef, newR, ctr},
	Ltype = Ldata[[1]];
	level = Ldata[[2]];
	OMEGA = Ldata[[4]];
	If[Length[Ldata]>4,
		parity = Ldata[[5]];
	,
		parity = {};
	];
	locLdata = Ldata;
	If[Length[locLdata]==4,
		AppendTo[locLdata,{}];
	];
	nrOfPrimes = 9;
	coeflist = Table[{},{k,nrOfPrimes + 1}];
	high=Sort[Get[getArchiveFileName[Ldata, "High"]]];
	low=Sort[Get[getArchiveFileName[Ldata, "Low"]]];
	Switch[Ltype,
	"GL3",
		If[parity == {1},
			targetFolder = "Publish/SL3/Odd/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/sl3MaassOdd";
		,
			If[level==1,
				targetFolder = "Publish/SL3/Level1/";
				eigenvalueFile = targetFolder <> "EigenvalueList.txt";
				lcalcBase = targetFolder <> "Lcalc/sl3Maass";
			,
				targetFolder = "Publish/SL3/Level" <> ToString[level] <> "/";
				eigenvalueFile = targetFolder <> "EigenvalueListOMEGA" <> OMEGAcode[OMEGA] <> ".txt";
				lcalcBase = targetFolder <> "Lcalc/sl3MaassLevel" <> ToString[level] <> "OMEGA" <> OMEGAcode[OMEGA] <> "_";
			];
		];
	,"GL4",
			targetFolder = "Publish/GL4/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/sl4Maass";
	,"SP4",
			targetFolder = "Publish/SP4/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/sp4Maass";
	,"SP4Std",
			targetFolder = "Publish/SP4Std/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/sp4stdMaass";
	,"HoloxMaass",
			targetFolder = "Publish/HoloxMaass/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/HoloxMaass";
	,"GL3Holo",
			targetFolder = "Publish/GL3Holo/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/GL3Holo";
	,"GL4Holo",
			targetFolder = "Publish/GL4Holo/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/GL4Holo";
	,"SP6",
			targetFolder = "Publish/SP6/";
			eigenvalueFile = targetFolder <> "EigenvalueList.txt";
			lcalcBase = targetFolder <> "Lcalc/SP6Maass";
	];
	streamList = OpenWrite[eigenvalueFile];
Print[eigenvalueFile];
	ctr=1;
	For[i=1,i<=Length[high],i++,
		eig = high[[i,1]];
		coef = removeCoefRubbish[Ldata, high[[i,2]]];
		newR = high[[i,3]];
		If[Length[high[[i]]]>4,
			locLdata[[5]] = high[[i,5]];  (* Change parity *)
		];
		data2Lcalc[locLdata, eig, coef, lcalcBase <> ToString[ctr] <> ".txt",newR,0];
		WriteString[streamList, ToString[ctr], ". ", N[eig, 15],"\n"];
		For[j=1,j<=nrOfPrimes,j++,
			AppendTo[coeflist[[j]], fixedNrOfDigits[coef[[2 Prime[j]-1]],15] + I fixedNrOfDigits[coef[[2 Prime[j]]], 15]];
		];
		j = nrOfPrimes+1;
		While[Prime[j]< Length[coef]/2 - 5,
			AppendTo[coeflist[[nrOfPrimes+1]], fixedNrOfDigits[coef[[2 Prime[j]-1]],15] + I fixedNrOfDigits[coef[[2 Prime[j]]], 15]];
			j++;
		];
		ctr++;
	];
	For[i=1,i<=Length[low],i++,
		eig = low[[i,1]];
		coef = removeCoefRubbish[Ldata,low[[i,2]]];
		newR = low[[i,3]];
		If[Length[low[[i]]]>4,
			locLdata[[5]] = low[[i,5]];  (* Change parity *)
		];
		data2Lcalc[locLdata, eig,coef,lcalcBase <> ToString[ctr] <> ".txt", newR, 0];
		WriteString[streamList, ToString[ctr], ". ",  N[eig, 15],"\n"];
Print["Done with data2lcalc 2"];
		For[j=1,j<=nrOfPrimes,j++,
			AppendTo[coeflist[[j]], fixedNrOfDigits[coef[[2 Prime[j]-1]],15] + I fixedNrOfDigits[coef[[2 Prime[j]]], 15]];
		];
Print["Done with data2lcalc 3"];
		j = nrOfPrimes+1;
		While[Prime[j]< Length[coef]/2 - 5,
			AppendTo[coeflist[[nrOfPrimes+1]], fixedNrOfDigits[coef[[2 Prime[j]-1]],15] + I fixedNrOfDigits[coef[[2 Prime[j]]], 15]];
			j++;
		];
		ctr++;
	];
	For[j=1,j<=nrOfPrimes,j++,
		Put[coeflist[[j]], targetFolder <> "Coef" <> ToString[Prime[j]] <>".txt"];
	];
	Put[coeflist[[nrOfPrimes+1]], targetFolder <> "CoefHigher.txt"];

	Close[streamList];
]


(* Format of Lcalc file
3  (complex coefficients)[
0  (unknown type)
16 (# of coef)
0  (not periodic)
3  (# of Gamma-factors)
.5 (Gamma-factors)
0 R1
.5
0 R2
.5
0 -(R1+R2)
1/Pi^(3/2)  (Q)
1 0  (omega)
0  (# of poles)
ReA1 ImA1  (coefficients)
ReA2 ImA2
....
ReA16 ImA16
*)
