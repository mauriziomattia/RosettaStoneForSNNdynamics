(* ::Package:: *)

(* ::Title:: *)
(*LIFLibrary*)


(* ::Subtitle:: *)
(*Written by Maurizio Mattia, 2018*)


(* ::Subsubtitle:: *)
(*Library of functions computing eigenvalues and eigenfunctions of the Fokker-Plank operator L for the leaky-IF (LIF) neuron. Formulas are from (Brunel & Hakim, Neural Comput 1999) and (Deniz & Rotter, Phys Rev E 2017). *)
(*Special thanks to Gianni Valerio Vinci for helpful discussions, pointing out an error in the normalization factor of eigenfunctions Subscript[\[CapitalPhi], n] expressed as a combination of parabolic cylinder functions and contributing to understand bifurcation of DD \[Lambda]s*)


(* ::Section:: *)
(*Numerical expression for the gain function \[CapitalPhi] of LIF neurons*)


(* ::Text:: *)
(*To conclude we can put together all these results to compose a fast and precise numerical implementation of the gain function \[CapitalPhi] for the LIF neuron.*)


(* ::Input::Initialization:: *)
Psi[w_]=If[w>=-3,E^w^2 (Erf[w]+1),2/Sqrt[\[Pi]] NIntegrate[E^-y^2 E^(2y w),{y,0,5}]];


(* ::Subsubsection:: *)
(*For large |a| and |b|:*)


(* ::Input:: *)
(*MFPT1[a_,b_]:=Sqrt[\[Pi]]NIntegrate[Psi[w],{w,a,b}];*)


(* ::Subsubsection:: *)
(*For relatively small |a| and |b|:*)


(* ::Input::Initialization:: *)
MFPTfrom0[b_]=Sqrt[\[Pi]]DawsonF[b]Psi[b]-2\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(b\)]\(DawsonF[x] \[DifferentialD]x\)\);
MFPTto0[a_]=-MFPTfrom0[a];
MFPT2[a_,b_]=MFPTfrom0[b]+MFPTto0[a];


(* ::Subsubsection:: *)
(*For large b\[GreaterGreater]10:*)


(* ::Text:: *)
(*At very strong subthreshold regime we expect to have huge mean FPTs. Setting as upper limit for the mean FPT of 10^100 times the decay constant \[Tau]V, the upper bound is b<15. After this safely can be assumed MFTP[a,b]->\[Infinity].  *)


(* ::Input:: *)
(*FindRoot[Log[MFPTfrom0[b]]==100Log[10],{b,1}][[1]]*)


(* ::Subsubsection:: *)
(*For a<=b<=0:*)


(* ::Input::Initialization:: *)
MFPT3[a_,b_]:=NIntegrate[((E^(2 b w)-E^(2 a w))/w)Exp[-w^2],{w,0,\[Infinity]}]


(* ::Subsection:: *)
(*Here is the optimal estimate of mean FPT in unit of \[Tau]V*)


(* ::Input::Initialization:: *)
MFPT[a_,b_]:=If[b<=0,MFPT3[a,b],If[b>15,\[Infinity],If [a<=0,MFPT3[a,0]+MFPT2[0,b],MFPT2[a,b]]]]


(* ::Input::Initialization:: *)
Phi[\[Mu]_,\[Sigma]_,\[Tau]_,H_,\[Theta]_,Tarp_]:=1/(Tarp+\[Tau] MFPT[(H-\[Mu])/\[Sigma],(\[Theta]-\[Mu])/\[Sigma]]);


(* ::Text:: *)
(*Here is the resulting plot to be compared to the one produced with the first approximated expression for \[CapitalPhi] (see at the beginning of the notebook). *)


(* ::Input:: *)
(*Plot[{Phi[\[Mu],0.25,0.010,10,20,0.002],Phi[\[Mu],1,0.010,10,20,0.002],Phi[\[Mu],4,0.010,10,20,0.002]},{\[Mu],-10,60},PlotStyle->{Red,Green,Blue},PlotRange->All]*)


(* ::Subsection:: *)
(*Other useful functions to fix \[Mu] and \[Sigma]*)


(* ::Text:: *)
(*Given the firing rate \[Nu] and the mean input current \[Mu], which is the value of \[Sigma]?*)


(* ::Input::Initialization:: *)
GetSigma[\[Nu]_,\[Mu]_,\[Tau]_,Vr_,Vt_,Tarp_]:=\[Sigma]/.FindRoot[Phi[\[Mu],\[Sigma],\[Tau],Vr,Vt,Tarp]==\[Nu],{\[Sigma],10.0}];


(* ::Text:: *)
(*In the limit \[Sigma]->0 what is the maximum value of \[Mu] at fixed \[Nu]?*)


(* ::Input:: *)
(*V[t]/.Last[DSolve[{\[Tau] V'[t]==-V[t]+\[Mu] \[Tau],V[0]==Vr},V[t],t]]*)


(* ::Input:: *)
(*Solve[%==Vt,t]*)


(* ::Input:: *)
(*Solve[\[Tau] Log[(-Vr+\[Mu] \[Tau])/(-Vt+\[Mu] \[Tau])]+Tarp==1/\[Nu],\[Mu]]*)


(* ::Input::Initialization:: *)
MaxMu[\[Nu]_,\[Tau]_,Vr_,Vt_,Tarp_]:=(-Vr+E^((1-Tarp \[Nu])/(\[Nu] \[Tau])) Vt)/((-1+E^((1-Tarp \[Nu])/(\[Nu] \[Tau]))) \[Tau]);


(* ::Section:: *)
(*Implementing (Denis & Rotter, 2017) formulas (and few alternatives)*)


(* ::Text:: *)
(*Independent solutions of the following Sturm-Liouville problem [Eqs. (A2) and (A3)]:*)


(* ::Input::Initialization:: *)
\[Phi]1DR[x_,\[Lambda]_]:=Hypergeometric1F1[1/2-\[Lambda]/2,1/2,-x^2];


(* ::Input::Initialization:: *)
\[Phi]2DR[x_,\[Lambda]_]:=Gamma[\[Lambda]/2]/Gamma[(\[Lambda]+1)/2] Hypergeometric1F1[(1-\[Lambda])/2,1/2,-x^2]+2x Hypergeometric1F1[1-\[Lambda]/2,3/2,-x^2];


(* ::Text:: *)
(*Characteristic equation (CE) whose solution gives the eigenvalues of the Fokker-Planck (FP) operator [Eq. (A8)]:*)


(* ::Input::Initialization:: *)
CE[\[Lambda]_,Xr_,Xt_]:= \[Phi]2DR[Xt,\[Lambda]]E^Xt^2-  \[Phi]2DR[Xr,\[Lambda]]E^Xr^2;


(* ::Text:: *)
(*This version of the CE is rather robust to find the real eigenvalues of the FP at strong drift-dominated regime:*)


(* ::Input::Initialization:: *)
CEDD[\[Lambda]_,Xr_,Xt_]:= E^\[Lambda] (E^-Xr^2/Gamma[\[Lambda]/2] \[Phi]2DR[Xt,\[Lambda]]-E^-Xt^2/Gamma[\[Lambda]/2] \[Phi]2DR[Xr,\[Lambda]]);


(* ::Text:: *)
(*Another independent solution resulting from Mathematica and the its related CE (not used in the following):*)


(* ::Input::Initialization:: *)
\[Phi]2DRH[x_,\[Lambda]_]:=-E^-x^2 HermiteH[-\[Lambda],x];


(* ::Input::Initialization:: *)
CEH[\[Lambda]_,Xr_,Xt_]:= \[Phi]2DRH[Xt,\[Lambda]]E^Xt^2-  \[Phi]2DRH[Xr,\[Lambda]]E^Xr^2;


(* ::Text:: *)
(*Coefficients a[\[Lambda]], b[\[Lambda]] and d[\[Lambda]] from Eqs. (A9) and (A10) to compose the eigenfunctions of the FP operator:*)


(* ::Input::Initialization:: *)
aDR[\[Lambda]_,Xr_]:=E^Xr^2 \[Phi]2DR[Xr,\[Lambda]];


(* ::Text:: *)
(*note that due to the characteristic equation aDR[\[Lambda],Xr] == aDR[\[Lambda],Xt].*)


(* ::Input::Initialization:: *)
bDR[\[Lambda]_,Xt_]:=-E^Xt^2\[Phi]1DR[Xt,\[Lambda]];


(* ::Input::Initialization:: *)
dDR[\[Lambda]_,Xr_,Xt_]:=E^Xr^2 \[Phi]1DR[Xr,\[Lambda]]-E^Xt^2 \[Phi]1DR[Xt,\[Lambda]];


(* ::Text:: *)
(*Here are the eigenfunctions \[Phi][\[Lambda]] for any \[Lambda] != 0:*)


(* ::Input::Initialization:: *)
\[Phi]DR[x_,\[Lambda]_,Xr_,Xt_]:=If[x>=Xr,aDR[\[Lambda],Xr]\[Phi]1DR[x,\[Lambda]]+bDR[\[Lambda],Xt]\[Phi]2DR[x,\[Lambda]],dDR[\[Lambda],Xr,Xt]\[Phi]2DR[x,\[Lambda]]];


(* ::Input::Initialization:: *)
Der\[Phi]DR[x_,\[Lambda]_,Xr_,Xt_]=-(1/2)(aDR[\[Lambda],Xr]\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\ \(\[Phi]1DR[x, \[Lambda]]\)\)+bDR[\[Lambda],Xt]\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\ \(\[Phi]2DR[x, \[Lambda]]\)\));
flux[\[Mu]\[Tau]_,\[Sigma]\[Tau]_,\[Lambda]\[Tau]_]:=Module[{Xr,Xt},
Xr=(Vr-\[Mu]\[Tau])/\[Sigma]\[Tau]; Xt=(Vt-\[Mu]\[Tau])/\[Sigma]\[Tau];Der\[Phi]DR[Xt,\[Lambda]\[Tau],Xr,Xt]];


(* ::Text:: *)
(*while these are the eigenfunctions \[Psi][\[Lambda]] of the adjoint FP operator (B10):*)


(* ::Input::Initialization:: *)
\[Psi]DRNN[x_,\[Lambda]_]:=E^x^2 \[Phi]2DR[x,\[Lambda]]; (* This is the not normalized version. *)
RelXmin=5; (* This set the lower bound of the integral to compute the inner products. *)
\[Psi]DR[x_,\[Lambda]_,Xr_,Xt_]:=\[Psi]DRNN[x,\[Lambda]]/NIntegrate[\[Psi]DRNN[y,\[Lambda]]\[Phi]DR[y,\[Lambda],Xr,Xt],{y,Xr-RelXmin(Xt- Xr),Xt}];


(* ::Text:: *)
(*Here is the stationary probability density, the eigenfunction associated to \[Lambda]=0. It is an elaboration of Eq. (3.10) in (Brunel & Hakim, 1999):*)


(* ::Input::Initialization:: *)
\[Phi]BH0[x_,Xr_,Xt_]:=Sqrt[\[Pi]] /MFPT[Xr,Xt] E^-x^2 (Erfi[Xt]-If[x<=Xr,Erfi[Xr],Erfi[x]])


(* ::Input::Initialization:: *)
flux0[\[Mu]\[Tau]_,\[Sigma]\[Tau]_]:=Module[{Xr,Xt},
Xr=(Vr-\[Mu]\[Tau])/\[Sigma]\[Tau]; Xt=(Vt-\[Mu]\[Tau])/\[Sigma]\[Tau];1/MFPT[Xr,Xt]];


(* ::Subsection:: *)
(*Formulas based on parabolic cylinder function Subscript[D, \[Nu]](z)*)


(* ::Input::Initialization:: *)
\[Phi]1DRD[x_,\[Lambda]_]:=2^(\[Lambda]/2-1)/Sqrt[\[Pi]] Gamma[(\[Lambda]+1)/2]E^(-(x^2/2)) (ParabolicCylinderD[-\[Lambda],Sqrt[2]x]+ParabolicCylinderD[-\[Lambda],-Sqrt[2]x]);
\[Phi]2DRD[x_,\[Lambda]_]:=2^(\[Lambda]/2)/Sqrt[\[Pi]] Gamma[\[Lambda]/2]E^(-(x^2/2)) ParabolicCylinderD[-\[Lambda],-Sqrt[2]x];


(* ::Input::Initialization:: *)
CEDorig[\[Lambda]_,Xr_,Xt_]:= \[Phi]2DRD[Xt,\[Lambda]]E^Xt^2-  \[Phi]2DRD[Xr,\[Lambda]]E^Xr^2;


(* ::Input::Initialization:: *)
CED[\[Lambda]_,Xr_,Xt_]:= 1/Sin[(\[Pi] \[Lambda])/2] (-(\[Lambda]/E))^(\[Lambda]/2) (E^(Xt^2/2) ParabolicCylinderD[-\[Lambda],-Sqrt[2]Xt]-E^(Xr^2/2) ParabolicCylinderD[-\[Lambda],-Sqrt[2]Xr]);


(* ::Input::Initialization:: *)
aDRD[\[Lambda]_,Xr_]:=E^Xr^2 \[Phi]2DRD[Xr,\[Lambda]];
bDRD[\[Lambda]_,Xt_]:=-E^Xt^2\[Phi]1DRD[Xt,\[Lambda]];
dDRD[\[Lambda]_,Xr_,Xt_]:=E^Xr^2 \[Phi]1DRD[Xr,\[Lambda]]-E^Xt^2 \[Phi]1DRD[Xt,\[Lambda]];


(* ::Input::Initialization:: *)
\[Phi]DRD[x_,\[Lambda]_,Xr_,Xt_]:=If[x>=Xr,aDRD[\[Lambda],Xr]\[Phi]1DRD[x,\[Lambda]]+bDRD[\[Lambda],Xt]\[Phi]2DRD[x,\[Lambda]],dDRD[\[Lambda],Xr,Xt]\[Phi]2DRD[x,\[Lambda]]];


(* ::Input::Initialization:: *)
Der\[Phi]DRD[x_,\[Lambda]_,Xr_,Xt_]=-(1/2)(aDRD[\[Lambda],Xr]\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\ \(\[Phi]1DRD[x, \[Lambda]]\)\)+bDRD[\[Lambda],Xt]\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\ \(\[Phi]2DRD[x, \[Lambda]]\)\));


(* ::Input::Initialization:: *)
\[Psi]DRDNN[x_,\[Lambda]_]:=E^x^2 \[Phi]2DRD[x,\[Lambda]]; (* This is the not normalized version. *)
RelXmin=5; (* This set the lower bound of the integral to compute the inner products. *)
\[Psi]DRD[x_,\[Lambda]_,Xr_,Xt_]:=\[Psi]DRDNN[x,\[Lambda]]/NIntegrate[\[Psi]DRDNN[y,\[Lambda]]\[Phi]DRD[y,\[Lambda],Xr,Xt],{y,Xr-RelXmin(Xt- Xr),Xt}];


(* ::Section:: *)
(*Compute eigenvalues \[Lambda]*)


(* ::Subsection:: *)
(*Given a guess for the \[Lambda]s returns their exact values*)


(* ::Input::Initialization:: *)
GetLambdasFromGuess[Xr_,Xt_,\[Lambda]\[Tau]0_]:=
Module[{\[Lambda]\[Tau],\[Lambda]\[Tau]1,ndx,NumericalZero=10^-10},
If[Depth[\[Lambda]\[Tau]0]==1,\[Lambda]\[Tau]1={\[Lambda]\[Tau]0},\[Lambda]\[Tau]1=\[Lambda]\[Tau]0];
\[Lambda]\[Tau]=Table[
If[Im[\[Lambda]\[Tau]1[[n]]]==0,
If[Xt==-Xr,
\[Lambda]\[Tau]/.FindRoot[CE[\[Lambda]\[Tau],Xr,Xt]==0,{\[Lambda]\[Tau],Re[\[Lambda]\[Tau]1[[n]]]}],
\[Lambda]\[Tau]/.FindRoot[CEDD[\[Lambda]\[Tau],Xr,Xt]==0,{\[Lambda]\[Tau],Re[\[Lambda]\[Tau]1[[n]]]}]
],
\[Lambda]\[Tau]/.FindRoot[CED[\[Lambda]\[Tau],Xr,Xt]==0,{\[Lambda]\[Tau],\[Lambda]\[Tau]1[[n]]}]
],
{n,Length[\[Lambda]\[Tau]1]}];

ndx=Flatten[Position[\[Lambda]\[Tau],x_/;Re[x]>0]];
If[Length[ndx]>0,Print["Warning: ", Length[ndx]," eigenvalues result to have Re[\[Lambda]]>0."]];

ndx=Flatten[Position[Table[
Min[Abs[CE[\[Lambda]\[Tau][[n]],Xr,Xt]],Abs[CED[\[Lambda]\[Tau][[n]],Xr,Xt]],Abs[CEDD[\[Lambda]\[Tau][[n]],Xr,Xt]]],
{n,Length[\[Lambda]\[Tau]]}],x_/;x>NumericalZero]];
If[Length[ndx]>0,Print["Warning: ", Length[ndx]," eigenvalues do not satisfy CE[\[Lambda]]==0 ",ndx,"."]];

(* \[Lambda]\[Tau]=SortBy[\[Lambda]\[Tau],-Re[#]&] *)
If[Length[\[Lambda]\[Tau]]==1,\[Lambda]\[Tau]=\[Lambda]\[Tau][[1]]];
\[Lambda]\[Tau]
]


(* ::Subsection:: *)
(*Return the DD \[Lambda]s at bifurcation, i.e. Xr == -Xt*)


(* ::Input::Initialization:: *)
GetDDLambdasAtBifurcation[Xt_,NumOfLambdas_]:=
Module[{\[Lambda]\[Tau]0},
(* The theoretical guesses. *)
\[Lambda]\[Tau]0=Table[(-3 n^2 \[Pi]^2+3 Xt^2-Xt^4)/(6 Xt^2),{n,NumOfLambdas}];
(* Refines numerically the above guesses. *)
GetLambdasFromGuess[-Xt,Xt,\[Lambda]\[Tau]0] 
];


(* ::Subsection:: *)
(*Guess for real DD \[Lambda]s (OBSOLETE)*)


(* ::Text:: *)
(*Previous guess is correct only under the assumption of \[Mu] \[Tau] much greater than Vt. If this is not the case we keep only the real part and for the imaginary part we take the following expression:*)


(* ::Input::Initialization:: *)
\[Alpha][k_,Xr_,Xt_]:=((2\[Pi] k)/(Xt-Xr))^2;
\[Beta][Xr_,Xt_]:=((Xt+Xr)/2)^2;
\[Gamma][k_,Xr_,Xt_]:=(-4 k^2 \[Pi]^2+2 Xr (Xr-Xt) Log[1+Sqrt[1-E^(-2 Xr (Xr-Xt))]]+Log[1+Sqrt[1-E^(-2 Xr (Xr-Xt))]]^2)/(2((Xr-Xt)^2) );
LambdaDDGuess[k_,Xr_,Xt_]:=If[Xt==-Xr,1/2-(k^2 \[Pi]^2)/(2 Xt^2)-Xt^2/6,If[k>1,\[Gamma][k,Xr,Xt]+I Sqrt[\[Alpha][k,Xr,Xt] \[Beta][Xr,Xt]],2\[Gamma][k,Xr,Xt]+I Sqrt[\[Alpha][k,Xr,Xt] \[Beta][Xr,Xt]]]];
LambdaDD[k_,Xr_,Xt_]:=\[Lambda]/.FindRoot[CED[\[Lambda],Xr,Xt]==0,{\[Lambda],LambdaDDGuess[k,Xr,Xt]}];
