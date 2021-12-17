(* ::Package:: *)

(*Jared Frazier*)
(*Variational-method-m6*)
(*Description:

Reproduction of variational principle findings from Physical Chemistry: A Molecular Approach (1997)
with demonstrated trial functions.

*)
(*Date: 11/18/2020*)
Clear["Global`*"]
(*-----------------------------------------------------*)
(*Question 1: Approximate Solution to Particle in 1D-box*)
(*-----------------------------------------------------*)

Print["
(*-----------------------------------------------------*)
(*Question 1: Approximate Solution to Particle in 1D-box*)
(*----------------------------------------------------*)"]
Print["(1) Find the normalization constant A and assume L = 1"]
Print["(2) Normalization constant can be be found by 
\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)|\[Phi]\!\(\*SuperscriptBox[\(|\), \(2\)]\)dx = 1 where \[Phi] = Ax(x-L)"]
Print["(3) Substituting \[Phi] into step (2), 
\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)(Ax(x-L)\!\(\*SuperscriptBox[\()\), \(2\)]\)dx = \!\(\*SuperscriptBox[\(A\), \(2\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)\!\(\*SuperscriptBox[\(x\), \(2\)]\)(x-L\!\(\*SuperscriptBox[\()\), \(2\)]\)dx = 1"]

(*Definite integral calculation*)
L = 1; (*Length of box*)
phi = A*x*(x-L); (*Phi trial function without normalization const*)
closedIntegralPhiSquare = Integrate[phi^2, {x, 0, L}]; (*Result of integral*)
(*/Definite integral calculation*)

Print["(4) Using Mathematica to calculate the definite integral
where L = 1, ", closedIntegralPhiSquare, " = 1"]
Print["(5) Normalization constant A therefore is A = +-\!\(\*SqrtBox[\(30\)]\)"]
Print["(6)The Hamiltonian for a harmonic oscillator (particle
in a 1D-box) is given by \!\(\*OverscriptBox[\(H\), \(^\)]\) = -\!\(\*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)]\)\!\(\*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]]\) + V(x) where V(x) = 0
in the bounds of the box 0 to L (L = 1)."]
Print["(7) \!\(\*SubscriptBox[\(E\), \(\[Phi]\)]\) = \!\(\*FractionBox[\(\(\(<\)\(\\\ \)\(\[Phi]\)\)\\\  | \\\ \*OverscriptBox[\(H\), \(^\)]\\\  | \\\ \(\(\[Phi]\)\(\\\ \)\(>\)\)\), \(\(\(<\)\(\\\ \)\(\[Phi]\)\)\\\  | \\\ \(\(\[Phi]\)\(\\\ \)\(>\)\)\)]\) = \!\(\*FractionBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\(\*SuperscriptBox[\(\[Phi]\), \(*\)] \*OverscriptBox[\(H\), \(^\)] \(\[Phi]dx\)\(\\\ \)\)\), \(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\(\(|\)\(\[Phi]\)\*SuperscriptBox[\(|\), \(2\)]\(dx\)\(\\\ \)\)\)]\) =
= \!\(\*FractionBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]Ax \((x - L)\) - \*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)] \*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]] Ax \((x - L)\) dx\), \(1\)]\) =  \!\(\*FractionBox[\(\(-\*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)]\) \(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]Ax \((x - L)\) \*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]] Ax \((x - L)\) dx\)\), \(1\)]\)"]

(*Second derivative of \[Phi] with respect to x and integral of phi*)
secondDerivPhi = D[phi, {x, 2}];
closedIntegralPhi = Integrate[phi, {x, 0, L}];
(*/Second derivative of \[Phi] with respect to x*)

Print["(8) The second derivative with respect to x of the integral
in the denominator of step (6) is ", secondDerivPhi, " 
then \!\(\*SubscriptBox[\(E\), \(\[Phi]\)]\) = \!\(\*FractionBox[\(\(-2\) A \*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)] \(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]Ax \((x - L)\) dx\)\), \(1\)]\) where the integral is ", 
closedIntegralPhi, " then \!\(\*SubscriptBox[\(E\), \(\[Phi]\)]\) =\!\(\*FractionBox[\(\*SuperscriptBox[\(A\), \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(6  m\)]\) = \!\(\*FractionBox[\(30 \*SuperscriptBox[\(h\), \(2\)]\), \(24 \*SuperscriptBox[\(m\[Pi]\), \(2\)]\)]\)"]

(*Calculation of Subscript[E, \[Phi]] and Subscript[E, o]*)
h = Quantity["PlanckConstant"]; (*J/s*)
electronMass = Quantity["ElectronMass"]; (*kg*) 
ePhi = (30*(h)^2)/(24*(electronMass)*Pi^2);
eNaught = (h)^2/(8*(electronMass));
(*/Calculation of Subscript[E, \[Phi]] and Subscript[E, 0]*)

Print["(9) Using mathematica, \!\(\*SubscriptBox[\(E\), \(\[Phi]\)]\) = ", N[ePhi, 4]]

Print["(10) This is a reasonable answer because \!\(\*SubscriptBox[\(E\), \(\[Phi]\)]\) > \!\(\*SubscriptBox[\(E\), \(0\)]\), where \!\(\*SubscriptBox[\(E\), \(0\)]\) = \!\(\*FractionBox[\(\*SuperscriptBox[\(\[Pi]\), \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)]\), \(2 \*SuperscriptBox[\(mL\), \(2\)]\)]\) = ", N[eNaught, 4] ]

(*Calculation of percent error*)
percentError =N[( ePhi - eNaught)/eNaught * 100, 4];
(*/Calculation of percent error*)

Print["(11) Therefore, \!\(\*SubscriptBox[\(E\), \(\[Phi]\)]\) = ", ePhi, " with ", percentError, "% error."]

(*---------------------------------------*)
(*Question 2: One parameter trial function*)
(*---------------------------------------*)
Print[
"\n(*---------------------------------------*)
(*Question 2: One parameter trial function*)
(*---------------------------------------*)"
]
Print["(1) First normalizing the function \[Phi](\[Alpha]) = A(\!\(\*SuperscriptBox[\(L\), \(2\)]\)- \!\(\*SuperscriptBox[\(x\), \(2\)]\))(\!\(\*SuperscriptBox[\(L\), \(2\)]\)-\!\(\*SuperscriptBox[\(\[Alpha]x\), \(2\)]\)) and using
\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)|\[Phi]\!\(\*SuperscriptBox[\(|\), \(2\)]\)dx = 1, therefore \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)(A(\!\(\*SuperscriptBox[\(L\), \(2\)]\)- \!\(\*SuperscriptBox[\(x\), \(2\)]\))(\!\(\*SuperscriptBox[\(L\), \(2\)]\)-\!\(\*SuperscriptBox[\(\[Alpha]x\), \(2\)]\))\!\(\*SuperscriptBox[\()\), \(2\)]\)dx = 1"]

(*Normalization integral*)
phiAlpha = A*(L^2-x^2)(L^2-\[Alpha]*x^2); (*Without A*)
normIntegralPhiAlpha = Integrate[phiAlpha^2, {x, 0, L}];
(*/Normalization integral*)

Print["(2) Using Mathematica, the result of integration is
", normIntegralPhiAlpha, " = 1"]

(*Solve for A*)
normConstant = Solve[normIntegralPhiAlpha == 1, A];
positiveNormConstant = A/.normConstant[[2]];
phiAlphaOnly = positiveNormConstant*(L^2-x^2)*(L^2-\[Alpha]*x^2); (*In terms of alpha, x, and L only*)
(*/Solve for A*)

Print["(3) Using Mathematica, the normalization constant
is A = ", positiveNormConstant]
Print["(4) Knowing the normalization constant,
E(\[Alpha]) = \!\(\*FractionBox[\(\(\(<\)\(\\\ \)\(\[Phi] \((\[Alpha])\)\)\)\\\  | \\\ \*OverscriptBox[\(H\), \(^\)]\\\  | \\\ \(\(\[Phi] \((\[Alpha])\)\)\(\\\ \)\(>\)\)\), \(\(\(<\)\(\\\ \)\(\[Phi] \((\[Alpha])\)\)\)\\\  | \\\ \(\(\[Phi] \((\[Alpha])\)\)\(\\\ \)\(>\)\)\)]\) = \!\(\*FractionBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\(\*SuperscriptBox[\((\[Phi] \((\[Alpha])\))\), \(*\)] \*OverscriptBox[\(H\), \(^\)] \(\[Phi]\) \((\[Alpha])\) \(dx\)\(\\\ \)\)\), \(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\(\(|\)\(\[Phi] \((\[Alpha])\)\)\*SuperscriptBox[\(|\), \(2\)]\(dx\)\(\\\ \)\)\)]\) = \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)(\[Phi](\[Alpha])\!\(\*SuperscriptBox[\()\), \(*\)]\)\!\(\*OverscriptBox[\(H\), \(^\)]\)\[Phi](\[Alpha])dx"]
Print["(5) \!\(\*OverscriptBox[\(H\), \(^\)]\) = -\!\(\*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)]\)\!\(\*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]]\) + V(x) but V(x) = 0 in the 1D box, so substituting
\!\(\*OverscriptBox[\(H\), \(^\)]\) and \[Phi](\[Alpha]) into step (4), E(\[Alpha]) = -\!\(\*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)", phiAlphaOnly," \!\(\*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]]\)", phiAlphaOnly, " dx"]

(*Second derivative of \[Phi] with respect to x and integral of phiAlpha*)
secondDerivPhiAlpha = D[phiAlphaOnly, {x, 2}];
(*/Second derivative of \[Phi] with respect to x and integral of phiAlpha*)

Print["(5) Second derivative inside the integrand is 
\!\(\*FractionBox[SuperscriptBox[\(d\), \(2\)], SuperscriptBox[\(dx\), \(2\)]]\)", phiAlphaOnly, " = ",  Simplify[secondDerivPhiAlpha], " The integral now becomes
E(\[Alpha]) = -\!\(\*FractionBox[SuperscriptBox[\(\[HBar]\), \(2\)], \(2  m\)]\)\!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(L\)]\)", Simplify[phiAlphaOnly*secondDerivPhiAlpha], "dx"]

(*Solution to integral*)
phiAlphaIntegrand = phiAlphaOnly*secondDerivPhiAlpha;
integralPhiAlpha = Integrate[phiAlphaIntegrand, {x, 0, L}];
(*/Solution to integral*)

Print["(6) The result of that integral is E(\[Alpha]) = ", Simplify[-h^2/(8*electronMass*Pi^2)* integralPhiAlpha]]

(*Derivative ePhi with respect to \[Alpha] set to 0*)
eAlpha = -h^2/(8*electronMass*Pi^2) * integralPhiAlpha;
eAlphaFunc[\[Alpha]_] = -h^2/(8*electronMass*Pi^2) * integralPhiAlpha;
derivEAlpha = D[eAlpha, {\[Alpha], 1}];
optimizedEAlpha = NSolve[derivEAlpha == 0, \[Alpha]];
alphaOne = \[Alpha]/.optimizedEAlpha[[1]];
alphaTwo = \[Alpha]/.optimizedEAlpha[[2]];
(*/Derivative ePhi with respect to \[Alpha] set to 0*)

Print["(7) Then optimizing the parameter \[Alpha], \!\(\*FractionBox[\(\[PartialD]E \((\[Alpha])\)\), \(\[PartialD]\[Alpha]\)]\) = ", Simplify[derivEAlpha], " = 0" ]
Print["(8) \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\) = ", N[alphaOne, 3], " or \!\(\*SubscriptBox[\(\[Alpha]\), \(2\)]\) = ", N[alphaTwo, 3]];
Print["(9) For \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\) = ", N[alphaOne, 3], " E(\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)) = ", eAlphaFunc[alphaOne]]
Print["    For \!\(\*SubscriptBox[\(\[Alpha]\), \(2\)]\) = ", N[alphaTwo, 3], " E(\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)) = ", eAlphaFunc[alphaTwo]]

(*%error both instances*)
percentErrorAlpha1 = (eAlphaFunc[alphaOne] - eNaught)/eNaught*100;
(*/%error both instances*)

Print["(10) The optimal parameter value for this function must therefore be
\!\(\*SubscriptBox[\(\[Alpha]\), \(\(1\)\(\\\ \)\)]\)= ", alphaOne, " with corresponding energy E(\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)) = ", eAlphaFunc[alphaOne]]

(*Residual Plot*)
exactPsi = Sqrt[2]*Sin[Pi*x];
phiOptimalAlphaOnly = phiAlphaOnly/.\[Alpha]->alphaOne;
Plot[
	{exactPsi-phiOptimalAlphaOnly, phiOptimalAlphaOnly, 
	exactPsi},
	{x, 0, L}, PlotLabel->"Residual plot",
	AxesLabel->{"Length"}, PlotRange->Full, ImageSize->Large,
	PlotLegends->{"\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\) - \[Phi](\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\),x)", "\[Phi](\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\),x)", "\!\(\*SubscriptBox[\(\[Psi]\), \(1\)]\)"}
]
(*----------------------------------*)
(*Question 3: Gaussian Trial Function*)
(*Pg 382 MCQ 2e*)
(*----------------------------------*)
Print[
"\n(*----------------------------------*)
(*Question 3: Gaussian Trial Function*)
(*----------------------------------*)"
]
Print["(1) Find normaliztion constant A of \[Phi](r,\[Alpha]) = \!\(\*SuperscriptBox[\(Ae\), \(-\*SuperscriptBox[\(\[Alpha]r\), \(2\)]\)]\) using procedures
in previous steps (i.e. Denominator is < \[Phi] | \[Phi] > =\!\(\*SubsuperscriptBox[\(\[Integral]\), \(-\[Infinity]\), \(\[Infinity]\)]\)|\[Phi]\!\(\*SuperscriptBox[\(|\), \(2\)]\)\!\(\*SuperscriptBox[\(r\), \(2\)]\)dr = 1"]

(*Normalize the Gaussian*)
gaussian = A*Exp[-1*\[Alpha]*r^2]; (*Gaussian trial function*)
gaussianNormalizationIntegral = Normal[Integrate[r^2*gaussian^2, {r, 0, Infinity}]];
gaussianNormalizationConst = Simplify[A/.(Solve[gaussianNormalizationIntegral == 1, A])[[2]]];
normalizedGaussian = gaussian/.A->gaussianNormalizationConst;
(*/Normalize the Gaussian*)

Print["(2) The normalization constant A therefore is A = ", gaussianNormalizationConst]
Print["(3) Numerator of E(\[Alpha]) = \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)\!\(\*SuperscriptBox[\(\[Phi]\), \(*\)]\)\!\(\*OverscriptBox[\(H\), \(^\)]\)\!\(\*SuperscriptBox[\(\[Phi]r\), \(2\)]\)dr where \!\(\*OverscriptBox[\(H\), \(^\)]\) = -\!\(\*FractionBox[\(1\), \(2\)]\)(\!\(\*FractionBox[\(d\), SuperscriptBox[\(dr\), \(2\)]]\)+\!\(\*FractionBox[\(2\), \(r\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\))-\!\(\*FractionBox[\(1\), \(r\)]\)"]

(*Hamiltonian operator*)
secondDerivGaussian = D[normalizedGaussian, {r, 2}];
firstDerivGaussian = D[normalizedGaussian, {r, 1}];
divRGaussian = normalizedGaussian/r;
hamiltonianOnGaussian = -1/2*(secondDerivGaussian+2/r*firstDerivGaussian) - divRGaussian;
(*/Hamilitonian operator*)

Print["(4) \!\(\*OverscriptBox[\(H\), \(^\)]\)\[Phi] = -\!\(\*FractionBox[\(A\), \(2\)]\)(\!\(\*FractionBox[\(d\), SuperscriptBox[\(dr\), \(2\)]]\)\[Phi] + \!\(\*FractionBox[\(2\), \(r\)]\)\!\(\*FractionBox[\(d\), \(dr\)]\)\[Phi]) - \!\(\*FractionBox[\(\[Phi]\), \(r\)]\) = ", Simplify[hamiltonianOnGaussian]]
Print["(5) Then E(\[Alpha]) = \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)", Simplify[normalizedGaussian*hamiltonianOnGaussian], " \!\(\*SuperscriptBox[\(r\), \(2\)]\)dr"]

(*Evaluate Hamiltonian with L = Infinity for hydrogen*)
eHAtom = Simplify[Normal[Integrate[normalizedGaussian*hamiltonianOnGaussian*r^2, {r, 0, Infinity}]]];
(*/Evaluate Hamiltonian with L = Infinity for hydrogen*)

Print["(6) E(\[Alpha]) = ", eHAtom]

(*Optimize E(\[Alpha])*)
firstDerivEHAtom = D[eHAtom, {\[Alpha], 1}];
optimizedEHAtom = NSolve[firstDerivEHAtom == 0, \[Alpha]];
alphaEHAtom = \[Alpha]/.optimizedEHAtom[[1]];
functionEHAtom[a_] = -2*Sqrt[2/Pi]*Sqrt[a]+(3*a/2);
(*/Optimize E(\[Alpha])*)
Print["(7) Optimization of E(\[Alpha]) => \!\(\*FractionBox[\(\[PartialD]E \((\[Alpha])\)\), \(\[PartialD]\[Alpha]\)]\) = ", firstDerivEHAtom, " = 0"]
Print["(8) Optimal parameter for \[Alpha] = ", alphaEHAtom]
Print["(9) Minimum energy is E(\[Alpha]) = ", functionEHAtom[alphaEHAtom], ". This
is a reasonable value since \!\(\*SubscriptBox[\(E\), \(min\)]\) in the textbook is the same constant (-0.424) times
\!\(\*FractionBox[\(\*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(e\), \(4\)]\), \(16 \*SuperscriptBox[\(\[Pi]\), \(2\)] \*SuperscriptBox[SubscriptBox[\(\[Epsilon]\), \(0\)], \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)]\)]\) and \!\(\*SubscriptBox[\(E\), \(0\)]\)'s constant is -0.500 times the \!\(\*FractionBox[\(\*SubscriptBox[\(m\), \(e\)] \*SuperscriptBox[\(e\), \(4\)]\), \(16 \*SuperscriptBox[\(\[Pi]\), \(2\)] \*SuperscriptBox[SubscriptBox[\(\[Epsilon]\), \(0\)], \(2\)] \*SuperscriptBox[\(\[HBar]\), \(2\)]\)]\) for Gaussian trial function
for H-atom"]



