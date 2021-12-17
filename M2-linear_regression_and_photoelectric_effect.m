(* ::Package:: *)

(*Fitting-M2*)
(*Jared Frazier*)
(*09-20-20*)
(*Description: 
First derive the line of best fit for a linear equation Overscript[y, ^] = Subscript[\[Beta], 0] + Subscript[\[Beta], 1]X, that is 
derive the coffecients Subscript[\[Beta], 0] and Subscript[\[Beta], 1]. Once these two coefficients are derived, use the following observations
from the Cathode Ray experiment to derive the variables associated with the equation modeling stopping 
voltage and frequency of light:

[1] The photoelectric effect is the emission of electrons when electromagnetic radiation (e.g., light),
hits a material.
[2] The expected emission of electrons was essentially continuous, meaning as the intensity (average power transfer over 
one period of the wave) of light increased, it was expected that the kinetic energy of emitted electrons would be
proportional. This expectation was refuted with experimental results.
[3] Incident light (radiation) had to have a frequency exceeding a certain threshold before electrons were ejected.
[4] The stopping potential V is the voltage at which no electrons that are emitted from the metal plate subject
to incident radiation are able to reach another plate (due to electron-electron repulsion). The stopping voltage
was measured for different frequencies \[Nu] of the incident beam.
[5] Lastly, it is observed from Einstein and the Cathode Ray Experiment that the kinetic energy of the ejected  electrons
must be Subscript[E, k] = h\[Nu] - \[Phi] = h\[Nu] - Subscript[h\[Nu], 0] where the units of Subscript[E, k] are the electron charge e and the voltage V. Solving for 
voltage, the stopping voltage becomes a simple linear equation whose parameters can be determined by minimzizing the 
sum of squared errors w.r.t Subscript[\[Beta], 0] and Subscript[\[Beta], 1] such that V(\[Nu]) = h/e\[Nu] - h/eSubscript[\[Nu], 0] = Subscript[\[Beta], 1]\[Nu] + Subscript[\[Beta], 0].

Refs:
https://aklectures.com/lecture/photoelectric-effect-compton-effect-wave-particle-duality/stopping-voltage-and-work-function
https://en.wikipedia.org/wiki/Photoelectric_effect*)

Clear["Global`"]

(*Import Notation package and define symbols*)
<<Notation`
Clear["Global`*"];
Notation[ParsedBoxWrapper[SubscriptBox[OverscriptBox["\[Beta]", "^"], "0"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["betaNaught"]]
Notation[ParsedBoxWrapper[SubscriptBox[OverscriptBox["\[Beta]", "^"], "1"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["betaOne"]]
Notation[ParsedBoxWrapper[SubscriptBox[OverscriptBox["Y", "^"], "k"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["yHat"]]
Notation[ParsedBoxWrapper[SubscriptBox["Y", "k"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["yData"]]
Notation[ParsedBoxWrapper[SubscriptBox["X", "k"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["xData"]]
Notation[ParsedBoxWrapper[SubscriptBox["Y", "avg"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["yBar"]]
Notation[ParsedBoxWrapper[SubscriptBox["X", "avg"]] \[DoubleLongLeftRightArrow] ParsedBoxWrapper["xBar"]]

(*Deriving the slope (Subscript[Overscript[\[Beta], ^], 0]) and intercept(Subscript[Overscript[\[Beta], ^], 1])*)
(*General form of linear regression equation, betaNaught and betaOne are unknown constants currently*)
yHat = betaNaught + betaOne * xData;

(*Define the general form for the sum of square errors*)
sse = Sum[(yData - yHat)^2, {k, 1, N}]; 
sseSummand = (yData - yHat)^2;

(*Display logic*)
Print["Part 1 - Derivation of slope 'a' (\!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\)) and intercept 'b' (\!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\)):"]
Print["(1) Linear relationship between stopping voltage and frequency of incident light
\tV(\[Nu]) = a\[Nu] + b can be modeled with the same equation \!\(\*SubscriptBox[OverscriptBox[\(Y\), \(^\)], \(k\)]\) = ", yHat, " where \!\(\*SubscriptBox[OverscriptBox[\(Y\), \(^\)], \(k\)]\) is the
\tpredicted response variable (V(\[Nu])) and \!\(\*SubscriptBox[\(X\), \(k\)]\) is the kth explanatory variable (\[Nu])."]
Print["(2) The sum of squared errors, also known as the mean squared error, is
\tthe sum of the difference between all response values in the data set (\!\(\*SubscriptBox[\(Y\), \(k\)]\)) and 
\tpredicted response values (\!\(\*SubscriptBox[OverscriptBox[\(Y\), \(^\)], \(k\)]\)) and is modeled by the equation SSE = \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)(\!\(\*SubscriptBox[\(Y\), \(k\)]\)-\!\(\*SubscriptBox[OverscriptBox[\(Y\), \(^\)], \(k\)]\)\!\(\*SuperscriptBox[\()\), \(2\)]\) 
\tfor the \!\(\*SuperscriptBox[\(k\), \(th\)]\) term in a data set with N terms."]
Print["(3) Substituting the definition of \!\(\*SubscriptBox[OverscriptBox[\(Y\), \(^\)], \(k\)]\) into SSE yields an equation which,
\tonce differentiated with respect to \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\) and \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\) and set equal to 0, may be solved
\tsimultaneously to yield the minimum values of \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\) and \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\) to minimize the 'distance' 
\tbetween the prediction 'line' and values of \!\(\*SubscriptBox[\(Y\), \(k\)]\): SSE = \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)(\!\(\*SubscriptBox[\(Y\), \(k\)]\)-(\!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\) + \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\)\!\(\*SubscriptBox[\(X\), \(k\)]\))\!\(\*SuperscriptBox[\()\), \(2\)]\)"]
(*Differentiate with respect to Subscript[Overscript[\[Beta], ^], 0] and Subscript[Overscript[\[Beta], ^], 1], derivative of summation is sum of derivatives*)
summandDerivBetaNaught = D[sseSummand, {betaNaught, 1}];
summandDerivBetaOne = D[sseSummand, {betaOne, 1}];
Print["(4)  \!\(\*FractionBox[\(\[PartialD]SSE\), \(\[PartialD]\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\)]\)\!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)(\!\(\*SubscriptBox[\(Y\), \(k\)]\)-(\!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\) + \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\)\!\(\*SubscriptBox[\(X\), \(k\)]\))\!\(\*SuperscriptBox[\()\), \(2\)]\) = \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)", summandDerivBetaNaught, 
"\tand \!\(\*FractionBox[\(\[PartialD]SSE\), \(\[PartialD]\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\)]\)\!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)(\!\(\*SubscriptBox[\(Y\), \(k\)]\)-(\!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(0\)]\) + \!\(\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\)\!\(\*SubscriptBox[\(X\), \(k\)]\))\!\(\*SuperscriptBox[\()\), \(2\)]\) = \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)", summandDerivBetaOne]

(*Solve partial derivative with respect to betaNaught \[Equal]  0*)
soln1 = Solve[{summandDerivBetaNaught==0}, {betaNaught}];
storeSoln1 = betaNaught/.soln1[[1]];

(*Summation of a constant*)
summationRule = Sum[betaNaught, {k, 1, N}];

(*Display solution for Subscript[Overscript[B, ^], 0] in terms of Subscript[Overscript[B, ^], 1]*)
Print["(5) The solution for ", summationRule, " = ", storeSoln1, " for N values in the  
\tdata is equivalent to \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(0\)]\) = \!\(\*OverscriptBox[\(Y\), \(_\)]\) - \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(1\)]\)\!\(\*OverscriptBox[\(X\), \(_\)]\) where \!\(\*OverscriptBox[\(Y\), \(_\)]\) and \!\(\*OverscriptBox[\(X\), \(_\)]\) are the average of the X and 
\tY values in the data."]

(*Equivalent solution for Subscript[Overscript[B, ^], 0]*)
betaNaughtAvg = yBar - betaOne*xBar;

(*Solve for Subscript[Overscript[B, ^], 1] by substituting Subscript[Overscript[B, ^], 0] in terms of Subscript[Overscript[B, ^], 1] into \[PartialD]SSE/\[PartialD]Subscript[Overscript[\[Beta], ^], 1]*)
betaOneDerivOneVarSummand = xData * (-betaNaughtAvg - betaOne*xData + yData);

(*Display new summation*)
Print["(6) Substituting \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(0\)]\) in terms of \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(1\)]\) into \!\(\*FractionBox[\(\[PartialD]SSE\), \(\[PartialD]\*SubscriptBox[OverscriptBox[\(\[Beta]\), \(^\)], \(1\)]\)]\), \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)" betaOneDerivOneVarSummand, 
" which simplifies to \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)", FullSimplify[betaOneDerivOneVarSummand] ]
(*Self expand*)
Print["(7) Expanding the summation using summation rules,
\t\!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)\!\(\*SubscriptBox[\(X\), \(k\)]\)(\!\(\*SubscriptBox[\(Y\), \(k\)]\)-\!\(\*SubscriptBox[\(Y\), \(avg\)]\)) - \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(1\)]\)\!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)\!\(\*SubscriptBox[\(X\), \(k\)]\)(\!\(\*SubscriptBox[\(X\), \(k\)]\)-\!\(\*SubscriptBox[\(X\), \(avg\)]\))=0. Mathematica cannot recognize that there are 
\ttwo variables \!\(\*SubscriptBox[\(Y\), \(k\)]\) and \!\(\*SubscriptBox[\(X\), \(k\)]\) in the simplified equation. It treats all of them as 
\tconstants and the syntax x[[k]] or y[[k]] cannot be used either to indicate 
\tthe \!\(\*SuperscriptBox[\(k\), \(th\)]\) variable of a predefined list (requires actual data)."]
Print["(8) Therefore, \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(1\)]\) = \!\(\*FractionBox[\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\*SubscriptBox[\(X\), \(k\)] \((\*SubscriptBox[\(Y\), \(k\)] - \*SubscriptBox[\(Y\), \(avg\)])\)\), \(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\*SubscriptBox[\(X\), \(k\)] \((\*SubscriptBox[\(X\), \(k\)] - \*SubscriptBox[\(X\), \(avg\)])\)\)]\) and \!\(\*SubscriptBox[OverscriptBox[\(B\), \(^\)], \(0\)]\) = ", betaNaughtAvg]
Print[]

(*Function to convert wavelength to frequency*)
toFrequency[wave_] := (3*10^8)/(wave*10^(-9));

(*x and y values in photoelectric data set*)
xList = {toFrequency[365], toFrequency[405], toFrequency[436], 
         toFrequency[546], toFrequency[577]}; 
yList = {1.855, 1.366, 1.143, 0.618, 0.498};
numEle = 5;

(*Use derived equations to calc linear equation *)
Print["Part 2 - Calculate slope and intercept of the equation for stopping voltage as a function of frequency V(\[Nu]):"]

(*Calculate averages of x and y lists*)
xAvg = Sum[xList[[i]], {i,1,numEle}]/numEle;
yAvg = Sum[yList[[i]], {i,1,numEle}]/numEle;
calcBetaOne = Sum[xList[[i]]*(yList[[i]]-yAvg), {i, 1, numEle}]/Sum[xList[[i]]*(xList[[i]]-xAvg), {i, 1, numEle}];
calcBetaNaught = yAvg-calcBetaOne*xAvg;

(*Line of best fit*)
linEq[\[Nu]_] := calcBetaNaught + calcBetaOne*\[Nu];
Print["(1) My line of best fit is \!\(\*OverscriptBox[\(V\), \(^\)]\)(\[Nu]) = ", linEq[\[Nu]], " and Planck's constant
\th is ", calcBetaOne, " eV*s which is close to the actual value
\t4.1357 \[Times] \!\(\*SuperscriptBox[\(10\), \(-15\)]\) eV*s. Therefore, %error is ", 100*(Abs[calcBetaOne-4.13157*10^(-15)]/(4.13157*10^(-15))), "%"]

(*Standard Deviation and Correlation coefficient r*)
sx = (Sum[(xList[[i]] - xAvg)^2, {i,1,numEle}]/(numEle-1))^(1/2);
sy = (Sum[(yList[[i]] - yAvg)^2, {i,1,numEle}]/(numEle-1))^(1/2);
r = 1/(numEle-1)*Sum[(xList[[i]] - xAvg)/sx*(yList[[i]]-yAvg)/sy, {i,1,numEle}];

(*Coefficient of determination*)
r2 = r^2;
Print["(2) The coefficient of determination \!\(\*SuperscriptBox[\(R\), \(2\)]\) = ", r2, " from R = \!\(\*FractionBox[\(1\), \(N - 1\)]\)*\!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\)(\!\(\*FractionBox[\(\*SubscriptBox[\(X\), \(k\)] - \*OverscriptBox[\(X\), \(_\)]\), SubscriptBox[\(S\), \(x\)]]\))(\!\(\*FractionBox[\(\*SubscriptBox[\(Y\), \(k\)] - \*OverscriptBox[\(Y\), \(_\)]\), SubscriptBox[\(S\), \(y\)]]\))"]

(*How well the line predicts values*)
eqError = Sqrt[Sum[(yList[[i]]-linEq[xList[[i]]])^2, {i, 1, numEle}]/(numEle-2)];
Print["(3) The error associated with predicted value \!\(\*OverscriptBox[\(V\), \(^\)]\)(\[Nu]) is modeled by the equation
\tS = \!\(\*SqrtBox[FractionBox[\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(N\)]\*SuperscriptBox[\((V \*SubscriptBox[\((\[Nu])\), \(k\)]\\\  - \\\ \*OverscriptBox[\(V\), \(^\)] \*SubscriptBox[\((\[Nu])\), \(k\)])\), \(2\)]\), \(N\\\  - \\\ 2\)]]\) = \[PlusMinus]", eqError, " Volts" ]

(*Plot data with my line of best fit*)
Print[]
xy2dList = {{toFrequency[365], 1.855},{ toFrequency[405], 1.366},{toFrequency[436],1.143 },
            {toFrequency[546], 0.618}, {toFrequency[577], 0.498}};
Print["Part 3 - Plot the data with my line of best fit and mathematicas line of best fit:"]
Print["(1) My line of best fit is \!\(\*OverscriptBox[\(V\), \(^\)]\)(\[Nu]) = ", linEq[\[Nu]], " and the plot is shown below:"]
Show[{Plot[linEq[x], {x,toFrequency[577],toFrequency[365]}, PlotLabel->"My Line of Best Fit for Photelectric Effect",
AxesLabel->{"Frequency \[Nu] (\!\(\*SuperscriptBox[\(s\), \(-1\)]\))", "Stopping Voltage V"}, ImageSize->700],  ListPlot[xy2dList]}]

(*Mathematica's line of best fit formula*)
mathematicaLinEq = LinearModelFit[xy2dList, \[Nu], \[Nu]];
functionMathematicaLinEq[\[Nu]_] = Normal[mathematicaLinEq];
Print["(2) The equation mathematica finds for the line of best fit is
\tV(\[Nu]) = ", Normal[mathematicaLinEq], " which directly matches my line of best of fit, and the plot is shown below:"]
Show[{Plot[functionMathematicaLinEq[x], {x,toFrequency[577],toFrequency[365]}, PlotLabel->"Mathematica's Line of Best Fit for Photelectric Effect",
AxesLabel->{"Frequency \[Nu] (\!\(\*SuperscriptBox[\(s\), \(-1\)]\))", "Stopping Voltage V"}, ImageSize->700],  ListPlot[xy2dList]}] 
