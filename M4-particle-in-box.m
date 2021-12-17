(* ::Package:: *)

(*Particle-in-box-M4*)
(*Jared Frazier*)
(*Description:

Plot the probability density functions for the time-independent wave function
in one and two dimensions. Derive the expectation of position and momentum. 
Determine the normalization constants.

*)
(*Due: 10-12-2020*)

Clear["Global`*"];

(*-----------------------------*)
   (*1-1 => Plotting 1D PDFs*)
(*-----------------------------*)

Print["---------------------------------------------
-----------------Question 1-----------------
---------------------------------------------"
]

(*P(x) = Subscript[\[Psi], n](x) for n[1,3]*)
waveN1 := 2*(Sin[Pi*x])^2;
waveN2 := 2*(Sin[2*Pi*x])^2;
waveN3 := 2*(Sin[3*Pi*x])^2;

(*Plot these wave functions from a = 0, a = L where L = 1*)
Print["Probability Density Functions For |\!\(\*SubscriptBox[\(\[Psi]\), \(n\)]\)(x)\!\(\*SuperscriptBox[\(|\), \(2\)]\) In 1D"]
Show[
	Plot[ 
		waveN1, 
		{x, 0, 1},
		PlotLabel->"|\!\(\*SubscriptBox[\(\[Psi]\), \(n\)]\)(x)\!\(\*SuperscriptBox[\(|\), \(2\)]\) for n = 1,2,3" ,
		AxesLabel->{"Length x", "|\!\(\*SubscriptBox[\(\[Psi]\), \(n\)]\)(x)\!\(\*SuperscriptBox[\(|\), \(2\)]\)"},
		PlotStyle->{Red, Dashed},
		Filling->Axis,
		PlotLegends->{"\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(x) = 2\!\(\*SuperscriptBox[\(sin\), \(2\)]\)\[Pi]x"}
	],
	Plot[
		waveN2, 
		{x, 0, 1},
		PlotStyle->{Green,Thick},
		Filling->Axis,
		PlotLegends->{"\!\(\*SubscriptBox[\(\[Psi]\), \(2\)]\)(x) = 2\!\(\*SuperscriptBox[\(sin\), \(2\)]\)2\[Pi]x"}
	],
	Plot[
		waveN3,
		{x,0,1},
		PlotRange->Full,
		Filling->Axis,
		PlotLegends->{"\!\(\*SubscriptBox[\(\[Psi]\), \(3\)]\)(x) = 2\!\(\*SuperscriptBox[\(sin\), \(2\)]\)3\[Pi]x"}
	]
]

(*-----------------------------*)
   (*1-2 => Plotting 2D PDFs*)
(*-----------------------------*)
Print["Density Functions for \!\(\*SubscriptBox[\(\[Psi]\), \(n\)]\)(x)\!\(\*SubscriptBox[\(\[Psi]\), \(n\)]\)(y)"]
(*2D wave functions for n=1,2,3,*)
waveN1X := Sqrt[2]*Sin[Pi*x];
waveN2X := Sqrt[2]*Sin[2*Pi*x];
waveN3X := Sqrt[2]*Sin[3*Pi*x];
waveN1Y := Sqrt[2]*Sin[Pi*y];
waveN2Y := Sqrt[2]*Sin[2*Pi*y];
waveN3Y := Sqrt[2]*Sin[2*Pi*y];
(*2D system*)
DensityPlot[
	waveN1X*waveN1Y,
	{x,0,1},
	{y,0,1},
	PlotLabel->"\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(x)\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(y) Density Plot"
]
DensityPlot[
	waveN2X*waveN2Y,
	{x,0,1},
	{y,0,1},
	PlotLabel->"\!\(\*SubscriptBox[\(\[Psi]\), \(2\)]\)(x)\!\(\*SubscriptBox[\(\[Psi]\), \(2\)]\)(y) Density Plot"
]
DensityPlot[
	waveN3X*waveN3Y,
	{x,0,1},
	{y,0,1},
	PlotLabel->"\!\(\*SubscriptBox[\(\[Psi]\), \(3\)]\)(x)\!\(\*SubscriptBox[\(\[Psi]\), \(3\)]\)(y)Density Plot"
]
(*-----------------------------*)
 (*2 => Expectation Derivation*)
(*-----------------------------*)
Print["---------------------------------------------
-----------------Question 2-----------------
---------------------------------------------"
]
(*Derive <x>*)
Print["Derive <x>"]
Print["(1) Steps 2-3 are the human intuition steps needed before Mathematica:"]
Print["(2) \!\(\*SuperscriptBox[\(sin\), \(2\)]\)(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\)) =\!\(\*FractionBox[\(1\\\  - \\\ cos \((\*FractionBox[\(n\[Pi]x\), \(L\)])\)\), \(2\)]\) in the integrand of <x> = \!\(\*FractionBox[\(2\), \(L\)]\)\[Integral]\!\(\*SuperscriptBox[\(xsin\), \(2\)]\)(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx 
from the double angle identity"]
Print["(3) Then the sum of two integrals can be used to solve for <x> 
such that <x> =\!\(\*FractionBox[\(1\), \(L\)]\)[\[Integral]xdx + \[Integral]cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx]"]
Print["(4) Given steps 2-3, I will use Mathematica to evaluate the integrals
\[Integral]xdx and \[Integral]cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx] separately if 'n' is any positive integer on the interval
x = 0 to x = L"]
firstIntegral = Integrate[x, {x,0,L}]; (*Compute both integrals from 0 to L*)
secondIntegral = Assuming[
					Element[n, PositiveIntegers],
					Integrate[Cos[(n*Pi*x)/L], {x, 0, L}]
				  ];
Print["(5) From mathematica, '\!\(\*FractionBox[\(1\), \(L\)]\) * ", firstIntegral, " + ", secondIntegral, 
"' is <x>."]
Print["(6) <x> = \!\(\*FractionBox[\(L\), \(2\)]\)"]
Print[]

(*Derive <px>*)
Print["Derive <\!\(\*SubscriptBox[\(p\), \(x\)]\)>"]
Print["(1) Since \!\(\*SubscriptBox[OverscriptBox[\(p\), \(^\)], \(x\)]\) = -i\[HBar]\!\(\*FractionBox[\(d\), \(dx\)]\), then <\!\(\*SubscriptBox[\(p\), \(x\)]\)> = -i\[HBar]\[Integral]\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))\!\(\*FractionBox[\(d\), \(dx\)]\)\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx
simplifies to <\!\(\*SubscriptBox[\(p\), \(x\)]\)> = -i\[HBar]\!\(\*FractionBox[\(2  n\[Pi]\), SuperscriptBox[\(L\), \(2\)]]\)\[Integral]sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx"]
expectationOfMomentum = Assuming [ (*Compute expectation of momentum*)
							Element[n, PositiveIntegers],
							Integrate[Sin[(n*Pi*x)/L]*Cos[(n*Pi*x)/L], {x, 0, L}]
						];
Print["(2) Using mathematica to evaluate the second integral from
step 1 on 
x = 0 to x = L where n is an element of the set of positive integers, then 
<\!\(\*SubscriptBox[\(p\), \(x\)]\)> = -i\[HBar]\!\(\*FractionBox[\(2  n\[Pi]\), SuperscriptBox[\(L\), \(2\)]]\)\[Integral]sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\)) = -i\[HBar]\!\(\*FractionBox[\(2  n\[Pi]\), SuperscriptBox[\(L\), \(2\)]]\)*", expectationOfMomentum,
"\n(3) Therefore, <\!\(\*SubscriptBox[\(p\), \(x\)]\)> = 0"]
Print[]

(*Derive <x^2> for Subscript[\[Sigma]^2, x]*)
Print["Derive <\!\(\*SuperscriptBox[\(x\), \(2\)]\)> For \!\(\*SubscriptBox[SuperscriptBox[\(\[Sigma]\), \(2\)], \(x\)]\) "]
Print["(1) \!\(\*SubscriptBox[\(\[Sigma]\), \(x\)]\) = \!\(\*SqrtBox[SubscriptBox[SuperscriptBox[\(\[Sigma]\), \(2\)], \(x\)]]\) = \!\(\*SqrtBox[\(\(<\)\*SuperscriptBox[\(x\), \(2\)]\(>\)\(\\\ \)\(-\\\ \(\(<\)\(x\)\*SuperscriptBox[\(>\), \(2\)]\)\)\)]\). Therefore, I must derive <\!\(\*SuperscriptBox[\(x\), \(2\)]\)> first
before calculating \!\(\*SubscriptBox[\(\[Sigma]\), \(x\)]\)."]
Print["(2) <\!\(\*SuperscriptBox[\(x\), \(2\)]\)> = \!\(\*FractionBox[\(2\), \(L\)]\)\[Integral]\!\(\*SuperscriptBox[\(x\), \(2\)]\)\!\(\*SuperscriptBox[\(sin\), \(2\)]\)(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx, just like with <x> except now \!\(\*SuperscriptBox[\(x\), \(2\)]\)
is in the integrand."]
Print["(3) \!\(\*SuperscriptBox[\(sin\), \(2\)]\)(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\)) = \!\(\*FractionBox[\(1\\\  - \\\ cos \((\*FractionBox[\(n\[Pi]x\), \(L\)])\)\), \(2\)]\) so the integral in step 2 expands
to \!\(\*FractionBox[\(1\), \(L\)]\)[\[Integral]\!\(\*SuperscriptBox[\(x\), \(2\)]\)dx + \[Integral]cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx]"]
Print["(4) From step 5 of the derivation of <x>, \[Integral]cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx is zero
if n is an element of the set of positive integers."]
Print["(5) \[Integral]\!\(\*SuperscriptBox[\(x\), \(2\)]\)dx = ", Integrate[x^2, {x,0,L}], " so \!\(\*FractionBox[\(1\), \(L\)]\) * ", Integrate[x^2, {x,0,L}]]
Print["(6) <\!\(\*SuperscriptBox[\(x\), \(2\)]\)> = \!\(\*FractionBox[SuperscriptBox[\(L\), \(2\)], \(3\)]\)"]
Print["(7) \!\(\*SubscriptBox[\(\[Sigma]\), \(x\)]\) = \!\(\*SqrtBox[SuperscriptBox[SubscriptBox[\(\[Sigma]\), \(x\)], \(2\)]]\) = \!\(\*SqrtBox[\(\(<\)\*SuperscriptBox[\(x\), \(2\)]\(>\)\(\\\ \)\(-\\\ \(\(<\)\(x\)\*SuperscriptBox[\(>\), \(2\)]\)\)\)]\) = ", Sqrt[L^2/3 - L^2/4]] (*Subscript[\[Sigma], x]*)

(*Derive <Subscript[p, x]^2> for Subscript[\[Sigma]^2, Subscript[p, x]]*)
Print[]
Print["Derive <\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(x\)], \(2\)]\)> For \!\(\*SubscriptBox[SuperscriptBox[\(\[Sigma]\), \(2\)], SubscriptBox[\(p\), \(x\)]]\)"]
Print["(1) \!\(\*SubscriptBox[\(\[Sigma]\), SubscriptBox[\(p\), \(x\)]]\) = \!\(\*SqrtBox[SuperscriptBox[SubscriptBox[\(\[Sigma]\), SubscriptBox[\(p\), \(x\)]], \(2\)]]\) = \!\(\*SqrtBox[\(\(<\)\*SuperscriptBox[SubscriptBox[\(p\), \(x\)], \(2\)]\(>\)\(\\\ \)\(-\\\ \(\(<\)\*SubscriptBox[\(p\), \(x\)]\*SuperscriptBox[\(>\), \(2\)]\)\)\)]\) and some human intuition is required
to set-up the proper integral, as in previous steps."]
Print["(2) This is the general expression
for  <\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(x\)], \(2\)]\)> = \[Integral]\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))(-i\[HBar]\!\(\*FractionBox[\(d\), \(dx\)]\)(-i\[HBar]\!\(\*FractionBox[\(d\), \(dx\)]\)\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))))dx."]
Print["(3) Several constants in the expression can be pulled 
out of the integrand: \!\(\*FractionBox[\(\(2\) \*SuperscriptBox[\(i\), \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)]\(\\\ \)\), \(L\)]\)\[Integral]sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))\!\(\*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]]\)sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx"]
Print["(4) I will evaluate \!\(\*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]]\)sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\)) in the integrand using Mathematica: ", 
Assuming[Element[n, PositiveIntegers], D[Sin[(n*Pi*x)/L] ,{x, 2}]]]
Print["(5) Then step 3 simplifies to -\!\(\*FractionBox[\(\(2\) \*SuperscriptBox[\(i\), \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)] \*SuperscriptBox[\(n\), \(2\)] \*SuperscriptBox[\(\[Pi]\), \(2\)]\(\\\ \)\), SuperscriptBox[\(L\), \(3\)]]\)\[Integral]sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))sin(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))dx, the
integrand of which was already computed in step 5 of 'Derive <x>'
by using the double angle identity"]
Print["(6) Step 5, therefore, simplifies to -\!\(\*FractionBox[\(\(2\) \*SuperscriptBox[\(i\), \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)] \*SuperscriptBox[\(n\), \(2\)] \*SuperscriptBox[\(\[Pi]\), \(2\)]\(\\\ \)\), SuperscriptBox[\(L\), \(3\)]]\)[\[Integral]xdx - \[Integral]cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\))], so
evaluating the second integral from x=0 to x=L using Mathematica, 
\[Integral]cos(\!\(\*FractionBox[\(n\[Pi]x\), \(L\)]\)) = ", Assuming[Element[n, PositiveIntegers],
Integrate[Cos[(n*Pi*x)/L], {x, 0, L}]]]
Print["(7)Then all that's left is -\!\(\*FractionBox[\(\(2\) \*SuperscriptBox[\(i\), \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)] \*SuperscriptBox[\(n\), \(2\)] \*SuperscriptBox[\(\[Pi]\), \(2\)]\(\\\ \)\), SuperscriptBox[\(L\), \(3\)]]\)\[Integral]xdx for x=0 to x=L,
and \[Integral]xdx is simply ", Integrate[x, {x,0,L}]]
Print["(8) <\!\(\*SuperscriptBox[SubscriptBox[\(p\), \(x\)], \(2\)]\)> = ",  (2(\[HBar]^2) (n^2) (\[Pi]^2) )/L^3 * Integrate[x, {x,0,L}]]
Print["(9) Then \!\(\*SubscriptBox[\(\[Sigma]\), SubscriptBox[\(p\), \(x\)]]\) = \!\(\*SqrtBox[SuperscriptBox[SubscriptBox[\(\[Sigma]\), SubscriptBox[\(p\), \(x\)]], \(2\)]]\) = \!\(\*SqrtBox[\(\(<\)\*SuperscriptBox[SubscriptBox[\(p\), \(x\)], \(2\)]\(>\)\(\\\ \)\(-\\\ \(\(<\)\*SubscriptBox[\(p\), \(x\)]\*SuperscriptBox[\(>\), \(2\)]\)\)\)]\) = ", Sqrt[(2(\[HBar]^2) (n^2) (Pi^2) )/L^3 * Integrate[x, {x,0,L}]]]

(*-----------------------------*)
(*3 => Acceptable wave function*)
(*-----------------------------*)
Print["---------------------------------------------
-----------------Question 3-----------------
---------------------------------------------"
]
(*Find normalization constant*)
Print["Find Normalization Constant 'A'"]
Print["(1) f(x) is normal if \[Integral]f(x)f(x)dx = 1. Therefore,
I can solve for A from \[Integral]\!\(\*SuperscriptBox[\(A\), \(2\)]\)[x(x-L)(x-L)]dx = 1 since 'A' is a constant."]
Print["(2) \!\(\*SuperscriptBox[\(A\), \(2\)]\)" , Integrate[(x*(x-L))^2, {x, 0, L}], " = 1"]
Print["(3) f(x) is normalized if A = ", Sqrt[1/(Integrate[(x*(x-L))^2, {x, 0, L}])]]
Print["(4) A = \!\(\*SqrtBox[FractionBox[\(30\), SuperscriptBox[\(L\), \(5\)]]]\)"]

Print[]
(*Plot f(x) with exact ground state wave function assuming L = 1*)
Print["Plot the PDFs assuming L = 1"]
L = 1;
groundStateWaveFunc = Sqrt[2/L]*Sin[(Pi*x)/L];
fx = Sqrt[30/L^5]*x*(L-x);
Show[
	Plot[ (*Plot fx squared*)
		(fx)^2,
		{x, 0, L},
		PlotRange->All,
		PlotLabel->"PDFs of |f(x)\!\(\*SuperscriptBox[\(|\), \(2\)]\), |\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(x)\!\(\*SuperscriptBox[\(|\), \(2\)]\), and f(x)\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(x)",
		AxesLabel->{"Length x","|\[Phi](x)\!\(\*SuperscriptBox[\(|\), \(2\)]\)"},
		PlotLegends->{"(f(x)\!\(\*SuperscriptBox[\()\), \(2\)]\) = (\!\(\*SqrtBox[FractionBox[\(30\), SuperscriptBox[\(L\), \(5\)]]]\)x(L-x)\!\(\*SuperscriptBox[\()\), \(2\)]\)"},
		PlotStyle->{Red, Dashed},
		ImageSize->Large
	],
	Plot[ (*Plot the ground state wave function squared*)
		(groundStateWaveFunc)^2,
		{x, 0, L},
		PlotLegends->{"(\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(x)\!\(\*SuperscriptBox[\()\), \(2\)]\) = (\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\!\(\*FractionBox[\(\[Pi]x\), \(L\)]\))\!\(\*SuperscriptBox[\()\), \(2\)]\)"},
		PlotStyle->{Orange, Bold}
	],
	Plot[ (*Plot difference between fx and ground state then square it*)
		(groundStateWaveFunc-fx)^2,
		{x, 0, L},
		PlotLegends->{"\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)(x)-f(x)\!\(\*SuperscriptBox[\()\), \(2\)]\) = (\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\!\(\*FractionBox[\(\[Pi]x\), \(L\)]\))-(\!\(\*SqrtBox[FractionBox[\(30\), SuperscriptBox[\(L\), \(5\)]]]\)x(L-x)\!\(\*SuperscriptBox[\()\), \(2\)]\)"},
		PlotStyle->{Green, Thick}
	]
]

(*Probability calculations for each PDF*)
Print["Calculate Probability Of x = 0 To x = \!\(\*FractionBox[\(L\), \(3\)]\) For Each PDF |\[Phi](x)\!\(\*SuperscriptBox[\(|\), \(2\)]\)"]
probGroundState = Integrate[groundStateWaveFunc^2, {x, 0, L/3}];
probFx = Integrate[fx^2, {x, 0, L/3}];
probProduct = Integrate[(groundStateWaveFunc-fx)^2, {x, 0, L/3}];
Print["(i)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), FractionBox[\(L\), \(3\)]]\)",groundStateWaveFunc^2, " = ", probGroundState]
Print["(ii)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), FractionBox[\(L\), \(3\)]]\)",fx^2, " = ", probFx]
Print["(iii)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), FractionBox[\(L\), \(3\)]]\)", (groundStateWaveFunc-fx)^2, " = ", probProduct]
Print[]
(*-----------------------------*)
 (*4 => Superposition*)
(*-----------------------------*)
Print["---------------------------------------------
-----------------Question 4-----------------
---------------------------------------------"
]
Print["Solve for normalization constant A of initial PDF"]
Print["(1) \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)\!\(\*SuperscriptBox[\(\[Psi]\[Psi]\), \(*\)]\) = 1, then \!\(\*SuperscriptBox[\(A\), \(2\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)(\!\(\*SubscriptBox[\(c\), \(1\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)+\!\(\*SubscriptBox[\(c\), \(2\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(2\)]\))(\!\(\*SubscriptBox[\(c\), \(1\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(\[Psi]\), \(*\)], \(1\)]\)+(\!\(\*SubscriptBox[\(c\), \(2\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(\[Psi]\), \(*\)], \(2\)]\))dx = 1."]
Print["(2) Expanding the integral in step 2, and using the fact that the
dot product of two orthonormal eigenstates is 0,
\!\(\*SuperscriptBox[\(A\), \(2\)]\)[(\!\(\*FractionBox[\(1\), \(4\)]\)\!\(\*SuperscriptBox[\()\), \(2\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(\[Psi]\), \(*\)], \(1\)]\)dx + 0 + 0 + (\!\(\*FractionBox[\(3\), \(4\)]\)\!\(\*SuperscriptBox[\()\), \(2\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(2\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(\[Psi]\), \(*\)], \(2\)]\)dx] = 1"]
Print["(3) Since \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(x\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(\[Psi]\), \(*\)], \(x\)]\)dx = 1, then \!\(\*SuperscriptBox[\(A\), \(2\)]\)(\!\(\*FractionBox[\(1\), \(16\)]\) +\!\(\*FractionBox[\(9\), \(16\)]\)) = 1. 
Therefore, A = \!\(\*SqrtBox[FractionBox[\(16\), \(10\)]]\)= \!\(\*FractionBox[\(4\), SqrtBox[\(10\)]]\)"]
Print[]
(*Plot Initial PDFs*)
normConstA = 4/Sqrt[10];
Print["Plot initial PDF Of |\[Psi]\!\(\*SuperscriptBox[\(|\), \(2\)]\) Where n=1, m=2, \!\(\*SubscriptBox[\(c\), \(1\)]\)=0.25, \!\(\*SubscriptBox[\(c\), \(2\)]\)=0.75, and A = \!\(\*FractionBox[\(4\), SqrtBox[\(10\)]]\)"]
Plot[
	(normConstA*(0.25*waveN1X+0.75*waveN2X))^2, 
	{x, 0, L}, 
	PlotRange->All, 
	PlotLabel->"|\[Psi](x,t=0)\!\(\*SuperscriptBox[\(|\), \(2\)]\)= (A(0.25\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(\[Pi]x)+0.75\!\(\*SqrtBox[FractionBox[\(2\), \(L\)]]\)sin(2\[Pi]x)\!\(\*SuperscriptBox[\()\), \(2\)]\) vs Length x",
	AxesLabel->{"Length x","|\[Psi](x,t=0)\!\(\*SuperscriptBox[\(|\), \(2\)]\)"},
	ImageSize->Large
]



