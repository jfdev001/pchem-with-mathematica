(* ::Package:: *)

(*H-Atom-M5*)
(*Jared Frazier*)
(*Description:

Reproducing radial distribution functions, spherical harmonics, and probability density functions
for the hydrogen atom from Quantum Chemistry 2nd Edition (McQuarrie).

*)
(*Date: 10/25/2020*)

Clear["Global`*"]

(*----------------------------*)
(*Question 1: Radial Functions*)
(*----------------------------*)

Print["(*----------------------------*)
(*Question 1: Radial Functions*)
(*----------------------------*)"]

(*Part a -- table form*)
Print["(a) Radial Wave Functions \!\(\*SubscriptBox[\(R\), \(nl\)]\)(r) through n = 2 where  0 <= l <= n-1 for n \[Element] \!\(\*SuperscriptBox[\(\[DoubleStruckCapitalZ]\), \(+\)]\):"]
tableEle = {{"\!\(\*SubscriptBox[\(R\), \(10\)]\)(r) = 2(\!\(\*FractionBox[\(Z\), SubscriptBox[\(a\), \(0\)]]\)\!\(\*SuperscriptBox[\()\), \(3/2\)]\)\!\(\*SuperscriptBox[\(e\), \(-\[Rho]\)]\)"}, {"\!\(\*SubscriptBox[\(R\), \(20\)]\)(r) = (\!\(\*FractionBox[\(Z\), \(2 \*SubscriptBox[\(a\), \(0\)]\)]\)\!\(\*SuperscriptBox[\()\), \(3/2\)]\)(2-\[Rho])\!\(\*SuperscriptBox[\(e\), \(\(-\[Rho]\)/2\)]\)"}, {"\!\(\*SubscriptBox[\(R\), \(21\)]\)(r) = \!\(\*FractionBox[\(1\), \(\[Sqrt]3\)]\)(\!\(\*FractionBox[\(Z\), \(2 \*SubscriptBox[\(a\), \(0\)]\)]\)\!\(\*SuperscriptBox[\()\), \(3/2\)]\)\!\(\*SuperscriptBox[\(\[Rho]e\), \(\(-\[Rho]\)/2\)]\)"}};
Grid[tableEle, Frame->All]
Print["where Z is the nuclear charge and \[Rho] = Zr/\!\(\*SubscriptBox[\(a\), \(0\)]\) and SubscriptBox[a,0\]is the Bohr radius."]
a0 = 1;
Print["Z = 1, \[Rho] = -r/\!\(\*SubscriptBox[\(a\), \(0\)]\), and \!\(\*SubscriptBox[\(a\), \(0\)]\) = 5.29E-11m = 1 bohrs \n"]

(*Part b -- plot radial distribution functions pg328 MQ2e*)
Print["(b) Plot the corresponding radial distributions (\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(R\), \(nl\)], \(2\)]\)(r)), also in table form:"]

(*Radial functions of r to plot for part b*)
R1s = r^2*(2/a0^(3/2)*Exp[-r/a0])^2;
R2s = r^2*(1/(2*a0^(3/2))*(2-(r/a0))*Exp[-r/(2*a0)])^2;
R2p = r^2*(1/(Sqrt[3]*(2*a0)^(3/2))*(r/a0)*Exp[-r/(2*a0)])^2;

(*Plots for each radial dist functions for table -- range for r from pg 329 in MQ2e*)
r1 = Plot[R1s, {r,0,5}, PlotRange->All, PlotStyle->{Orange}, PlotLabel->"Probability density for 1s \!\(\*SubscriptBox[\(H\), \(orbital\)]\)",
          AxesLabel->{"r/\!\(\*SubscriptBox[\(a\), \(0\)]\)","\!\(\*SuperscriptBox[\(r\), \(2\)]\)[\!\(\*SubscriptBox[\(R\), \(10\)]\)(r)\!\(\*SuperscriptBox[\(]\), \(2\)]\)/\!\(\*SubscriptBox[\(a\), \(0\)]\)"}, ImageSize->250];
r2 = Plot[R2s, {r,0,15}, PlotRange->All, PlotStyle->{Red}, PlotLabel->"Probability density for 2s \!\(\*SubscriptBox[\(H\), \(orbital\)]\)",
		  AxesLabel->{"r/\!\(\*SubscriptBox[\(a\), \(0\)]\)","\!\(\*SuperscriptBox[\(r\), \(2\)]\)[\!\(\*SubscriptBox[\(R\), \(20\)]\)(r)\!\(\*SuperscriptBox[\(]\), \(2\)]\)/\!\(\*SubscriptBox[\(a\), \(0\)]\)"}, ImageSize->250];
r3 = Plot[R2p, {r,0,15}, PlotRange->All, PlotStyle->{Blue}, PlotLabel->"Probability density for 2p \!\(\*SubscriptBox[\(H\), \(orbital\)]\)",
		  AxesLabel->{"r/\!\(\*SubscriptBox[\(a\), \(0\)]\)","\!\(\*SuperscriptBox[\(r\), \(2\)]\)[\!\(\*SubscriptBox[\(R\), \(21\)]\)(r)\!\(\*SuperscriptBox[\(]\), \(2\)]\)/\!\(\*SubscriptBox[\(a\), \(0\)]\)"}, ImageSize->250];
radialDistTable = {r1,r2,r3};
TableForm[radialDistTable, TableDirections->Row]

(*Part c -- max probability radius for 1s state*)
R1sFunc[r_] := r^2*(2/a0^(3/2)*Exp[-r/a0])^2;
R1sDerivative = D[R1s, {r, 1}];         (*First derivative*)
critPt = NSolve[R1sDerivative == 0, r]; (*Critical points*)
Print["\n(c) Maximum probability radius for 1s state is when r = ", r/.critPt[[2]]]
Print["Therefore \!\(\*SuperscriptBox[\(r\), \(2\)]\)[\!\(\*SubscriptBox[\(R\), \(10\)]\)(1)\!\(\*SuperscriptBox[\(]\), \(2\)]\)/\!\(\*SubscriptBox[\(a\), \(0\)]\) = ", R1sFunc[r/.critPt[[2]]]]

(*Part d -- 1s wave function is orthogonal to that of of the 2s*)
Print["\n(d) Show that 1s radial function is orthogonal to that of 2s"]
R10 = 2(1/a0)^(3/2) Exp[-r/a0];
R20 = (1/(2*a0))^(3/2)*(2-(r/a0))*Exp[-r/(2*a0)];
orthogonalRadFuncs = Integrate[(r^2)*R10*R20, {r, 0, Infinity}];
Print["From Mathematica calculations, <\!\(\*SubscriptBox[\(R\), \(10\)]\) | \!\(\*SubscriptBox[\(R\), \(20\)]\)> = ", orthogonalRadFuncs]

(*-------------------------------*)
(*Question 2: Spherical Harmonics*)
(*-------------------------------*)

(*Pg MQ2e Pg 323*)

Print["\n(*-------------------------------*)
(*Question 2: Spherical Harmonics*)
(*-------------------------------*)"]

(*Spherical Harmonic Elements of Table*)
shR1 = SphericalHarmonicY[0, 0, \[Theta], \[CurlyPhi]];
shR2 = SphericalHarmonicY[1, 0, \[Theta], \[CurlyPhi]];
shR3 = SphericalHarmonicY[1, 1, \[Theta], \[CurlyPhi]];
shR4 = SphericalHarmonicY[1, -1, \[Theta], \[CurlyPhi]];
shR5 = SphericalHarmonicY[2, 0, \[Theta], \[CurlyPhi]];
shR6 = SphericalHarmonicY[2, 1, \[Theta], \[CurlyPhi]];
shR7 = SphericalHarmonicY[2, -1, \[Theta], \[CurlyPhi]];
shR8 = SphericalHarmonicY[2, 2, \[Theta], \[CurlyPhi]];
shR9 = SphericalHarmonicY[2, -2, \[Theta], \[CurlyPhi]];

(*Print grid -- part 1*)
Print["Part 1 -- Table of \!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(\[ScriptL]\)], SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]]\)(\[Theta], \[CurlyPhi]) for 0 \[LessEqual] \[ScriptL] \[LessEqual] 2 and -\[ScriptL] \[LessEqual] \!\(\*SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]\) \[LessEqual] +\[ScriptL]"]
Grid[
	{
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(\[ScriptL]\)], SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]]\)(\[Theta], \[CurlyPhi]) = ", "Definition"},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(0\)], \(0\)]\)(\[Theta], \[CurlyPhi]) = ",  shR1},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(1\)], \(0\)]\)(\[Theta], \[CurlyPhi]) = ", shR2},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(1\)], \(1\)]\)(\[Theta], \[CurlyPhi]) = ", shR3},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(1\)], \(-1\)]\)(\[Theta], \[CurlyPhi]) = ", shR4},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(2\)], \(0\)]\)(\[Theta], \[CurlyPhi]) = ", shR5},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(2\)], \(1\)]\)(\[Theta], \[CurlyPhi]) = ", shR6},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(2\)], \(-1\)]\)(\[Theta], \[CurlyPhi]) = ", shR7},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(2\)], \(2\)]\)(\[Theta], \[CurlyPhi]) = ", shR8},
		{"\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(2\)], \(-2\)]\)(\[Theta], \[CurlyPhi]) = ", shR9}
	},
	Frame->All
]

(*Print spherical plot*)
Print["Part 2 -- Visualizing the Spherical Harmonics"]
Print["---------------"]
Print["\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(0\)], \(0\)]\)(\[Theta], \[CurlyPhi]) plot"]
Print["---------------"]
SphericalPlot3D[shR1, {\[Theta], 0, Pi}, {\[CurlyPhi], 0, 2*Pi}]
Print["----------------"]
Print["\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(1\)], SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]]\)(\[Theta], \[CurlyPhi]) plots"]
Print["----------------"]
Table[
	SphericalPlot3D[
		{Abs[SphericalHarmonicY[1, m, \[Theta], \[CurlyPhi]]]},
		{\[Theta], 0, Pi},
		{\[CurlyPhi], 0, 2*Pi}
	],
	{m, 0, 1}
]
Print["----------------"]
Print["\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(2\)], SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]]\)(\[Theta], \[CurlyPhi]) plots"]
Print["----------------"]
Table[
	SphericalPlot3D[
		{Abs[SphericalHarmonicY[2, m, \[Theta], \[CurlyPhi]]]},
		{\[Theta], 0, Pi},
		{\[CurlyPhi], 0, 2*Pi}
	],
	{m, 0, 2}
]

(*------------------*)
(*Question 3: Energy*)
(*------------------*)
Print["\n(*------------------*)
(*Question 3: Energy*)
(*------------------*)"]

Print["(a) For \!\(\*SubscriptBox[\(\[Psi]\), SubscriptBox[\(\[ScriptN]\[ScriptL]\[ScriptM]\), \(\[ScriptL]\)]]\)(r, \[Theta], \[Phi]) = \!\(\*SubscriptBox[\(R\), \(\[ScriptN]\[ScriptL]\)]\)(r)\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(\[ScriptL]\)], SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]]\)(\[Theta], \[Phi]) then for 1s function \[ScriptN] = 1, \[ScriptL] = 0, \!\(\*SubscriptBox[\(\[ScriptM]\), \(\[ScriptL]\)]\) = 0"]
Print["\!\(\*SubscriptBox[\(\[Psi]\), \(1  s\)]\)(r, \[Theta], \[Phi]) = \!\(\*SubscriptBox[\(\[Psi]\), \(100\)]\)(r, \[Theta], \[Phi]) = \!\(\*SubscriptBox[\(R\), \(10\)]\)(r)\!\(\*SuperscriptBox[SubscriptBox[\(Y\), \(0\)], \(0\)]\)(\[Theta], \[Phi]) = \!\(\*FractionBox[\(1\), SqrtBox[\(\[Pi]\)]]\)(\!\(\*FractionBox[\(1\), SubscriptBox[\(a\), \(0\)]]\)\!\(\*SuperscriptBox[\()\), \(3/2\)]\)(-r/\!\(\*SubscriptBox[\(a\), \(0\)]\))"]
Print["\tExplanation: For the 1s function to satisfy the radial Schrodinger equation,\!\(\*SubscriptBox[\(\\\ \), \(\(r\)\*OverscriptBox[\(\[RightArrow]\), \(lim\)]\(\[Infinity]\)\(\\\ \)\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(1  s\)]\)(r, \[Theta], \[Phi]) = 0.
This is because the number of radial nodes -- radial sections in the spherical coordinate plane
where \!\(\*SubscriptBox[\(R\), \(\[ScriptN]\[ScriptL]\)]\)(r) = 0 -- that is predicted from \[ScriptN] - \[ScriptL] - 1 should also be described by \!\(\*SubscriptBox[\(\[Psi]\), \(1  s\)]\)(r, \[Theta], \[Phi]). Since \[ScriptN], the principle
quantum number, is 1 and \[ScriptL], the angular quantum number, is 0, \[ScriptN] - \[ScriptL] - 1 predicts that number of radial nodes is 0.
\!\(\*SubscriptBox[\(\\\ \), \(\(r\)\*OverscriptBox[\(\[RightArrow]\), \(lim\)]\(\[Infinity]\)\(\\\ \)\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(1  s\)]\)(r, \[Theta], \[Phi]) = 0 predicts the same thing. Therefore, \!\(\*SubscriptBox[\(\[Psi]\), \(1  s\)]\)(r, \[Theta], \[Phi] satisfies the radial Schrodinger equation."]
Print["The energy of the electron in the 1s state is \!\(\*SubscriptBox[\(E\), \(n\)]\) = -\!\(\*FractionBox[\(13.6\\\ eV\), SuperscriptBox[\(n\), \(2\)]]\), or simply \!\(\*SubscriptBox[\(E\), \(1\)]\) = -13.6 eV]"]
 
 Print["\n(b) Deriving the energy for two trial functions \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) and \!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\)"]
 phi1 = Exp[-\[Alpha]*r];
 phi2 = Exp[-\[Alpha]*r^2];
 Print["Part 1 --> \[Epsilon][\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\)] for \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) = \!\(\*SuperscriptBox[\(e\), \(-\[Alpha]r\)]\)"]
 Print["To get \[Epsilon][\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\)] = \!\(\*FractionBox[\(\(\(<\)\(\\\ \)\*SubscriptBox[\(\[Phi]\), \(1\)]\)\\\  | \\\ \(\(\*OverscriptBox[\(H\), \(^\)] \*SubscriptBox[\(\[Phi]\), \(1\)]\)\(\\\ \)\(>\)\)\), \(\(\(<\)\(\\\ \)\*SubscriptBox[\(\[Phi]\), \(1\)]\)\\\  | \\\ \(\*SubscriptBox[\(\[Phi]\), \(1\)]\(\\\ \)\(>\)\)\)]\), I will first get < \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) | \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) >"]
 Print["(1) < \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) | \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) > = \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(\[Phi]\), \(1\)], \(2\)]\)dr, or ", 
 Normal[Assuming[Element[\[Alpha], PositiveReals], Integrate[r^2 * phi1^2, {r, 0, Infinity}]]]
 ]
 Print["(2) Then < \[Phi] | \!\(\*OverscriptBox[\(H\), \(^\)]\)\[Phi] > for \!\(\*OverscriptBox[\(H\), \(^\)]\) = \!\(\*FractionBox[\(-\*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(2 \*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(r\), \(2\)]\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)(\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)) - \!\(\*FractionBox[SuperscriptBox[\(e\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)] r\)]\)."]
 Print["(3) Substituting \!\(\*OverscriptBox[\(H\), \(^\)]\) into < \[Phi] | \!\(\*OverscriptBox[\(H\), \(^\)]\)\[Phi] > => 
< \[Phi] | \!\(\*FractionBox[\(-\*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(2 \*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(r\), \(2\)]\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)(\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)) - \!\(\*FractionBox[SuperscriptBox[\(e\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)] r\)]\) | \[Phi] >"]
Print["(4) This is equivalent to: 
- < \[Phi] | \!\(\*FractionBox[\(-\*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(2 \*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(r\), \(2\)]\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)(\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)) | \[Phi] > - < \[Phi] | \!\(\*FractionBox[SuperscriptBox[\(e\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)] r\)]\) | \[Phi] >"]
lhsIntegration = Assuming[
	Element[\[Alpha], PositiveReals], 
	Integrate[
		Exp[-\[Alpha]*r]*(-2*\[Alpha]*r*Exp[-\[Alpha]*r]+\[Alpha]^2*r^2*Exp[-\[Alpha]*r]),
		{r, 0, Infinity}]
];
Print["(5) The left side of the equation is -\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(-\[Alpha]r\)]\)(\!\(\*FractionBox[\(-\*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(2 \*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(r\), \(2\)]\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)(\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\))\!\(\*SuperscriptBox[\(e\), \(-\[Alpha]r\)]\))dr.
Integrating using mathematica, this is: \!\(\*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2 \*SubscriptBox[\(m\), \(e\)]\)]\)*", lhsIntegration]
rhsIntegration = Assuming[
	Element[\[Alpha], PositiveReals], 
	Integrate[
		r*Exp[-2*\[Alpha]*r],
		{r, 0, Infinity}]
];
Print["(6) The right side of the equation is -\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(-\[Alpha]r\)]\)((\!\(\*FractionBox[SuperscriptBox[\(e\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)] r\)]\))\!\(\*SuperscriptBox[\(e\), \(-\[Alpha]r\)]\))dr, or
\!\(\*FractionBox[\(\(-1\) \*SuperscriptBox[\(\[ScriptE]\), \(2\)]\), \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)]\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)(\!\(\*SuperscriptBox[\(re\), \(-\[Alpha]r\)]\)*\!\(\*SuperscriptBox[\(e\), \(-\[Alpha]r\)]\))dr, or \!\(\*FractionBox[\(-1\), \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)]\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(re\), \(\(-2\) \[Alpha]r\)]\)dr. 
Integrating using mathematica: \!\(\*FractionBox[\(-1\), \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)]\)]\) * ", rhsIntegration]
Print["(7) Combining the results from step 1 (the denominator), step 5 and 6 (the numerator),
\[Epsilon][\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\)] =\!\(\*FractionBox[\(\([\*FractionBox[\(\[HBar]\), \(2 \*SubscriptBox[\(m\), \(\[ScriptE]\)]\)]*\*FractionBox[\(1\), \(4  \[Alpha]\)]]\)\(-\)\([\*FractionBox[\(1 \*SuperscriptBox[\(\[ScriptE]\), \(2\)]\), \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)]\)]*\*FractionBox[\(1\), \(4 \*SuperscriptBox[\(\[Alpha]\), \(2\)]\)]]\)\(\\\ \)\), \([\*FractionBox[\(1\), \(4 \*SuperscriptBox[\(\[Alpha]\), \(3\)]\)]]\)]\)"]

Print["\nPart 2 --> \[Epsilon][\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\)] for \!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\) = \!\(\*SuperscriptBox[\(e\), \(-\*SuperscriptBox[\(\[Alpha]r\), \(2\)]\)]\)"]
Print["(1) Separating the numerator out and denominator out as before, 
the denominator is < \!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\) | \!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\) > = \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[SubscriptBox[\(\[Phi]\), \(2\)], \(2\)]\)dr, or ", 
Normal[Assuming[Element[\[Alpha], PositiveReals], Integrate[r^2 * phi2^2, {r, 0, Infinity}]]]]
lhsIntegrationPhi2 = Assuming[
	Element[\[Alpha], PositiveReals], 
	Integrate[
		Exp[-\[Alpha]*r^2]*(3*\[Alpha]*r^2*Exp[-\[Alpha]*r^2]-2*\[Alpha]^2r^4*Exp[-\[Alpha]*r^2]),
		{r, 0, Infinity}
	]
];
Print["(2) The left side of the numerator is \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(-\*SuperscriptBox[\(\[Alpha]r\), \(2\)]\)]\)(\!\(\*FractionBox[\(-\*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(2 \*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(r\), \(2\)]\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)(\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\))\!\(\*SuperscriptBox[\(e\), \(-\*SuperscriptBox[\(\[Alpha]r\), \(2\)]\)]\))dr,
so the result is -\!\(\*FractionBox[\(\[HBar]\), SubscriptBox[\(m\), \(e\)]]\)*", lhsIntegrationPhi2]
rhsIntegrationPhi2 = Assuming[
	Element[\[Alpha], PositiveReals], 
	Integrate[
		r*Exp[-2*\[Alpha]*r^2],
		{r, 0, Infinity}]
];
Print["(3) The right side of the numerator is -\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(-\*SuperscriptBox[\(\[Alpha]r\), \(2\)]\)]\)((\!\(\*FractionBox[SuperscriptBox[\(\[ScriptE]\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)] r\)]\))\!\(\*SuperscriptBox[\(e\), \(-\*SuperscriptBox[\(\[Alpha]r\), \(2\)]\)]\))dr, 
so the result is -\!\(\*FractionBox[SuperscriptBox[\(\[ScriptE]\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)]\)]\)* ", rhsIntegrationPhi2]
Print["(4) Putting the results of step 1, 2, and 3 together:
\[Epsilon][\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\)] =\!\(\*FractionBox[\(\(-\([\*FractionBox[\(\[HBar]\), SubscriptBox[\(m\), \(e\)]]*\*FractionBox[\(3\\\ \*SqrtBox[FractionBox[\(\[Pi]\), \(2\)]]\), \(16\\\ \*SqrtBox[\(\[Alpha]\)]\)]]\)\)\(\\\ \)\(-\)\(\\\ \)\([\*FractionBox[SuperscriptBox[\(\[ScriptE]\), \(2\)], \(4 \*SubscriptBox[\(\[Pi]\[Epsilon]\), \(0\)]\)]*\*FractionBox[\(1\), \(4\\\ \[Alpha]\)]]\)\(\\\ \)\), \([\*FractionBox[SqrtBox[FractionBox[\(\[Pi]\), \(2\)]], \(8\\\ \*SuperscriptBox[\(\[Alpha]\), \(3/2\)]\)]]\)]\)"]



