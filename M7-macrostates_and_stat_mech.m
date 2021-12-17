(* ::Package:: *)

(*Jared Frazier*)
(*Microstates-M7*)
(*Description:

Relate the multiplicity of macrostates to principles of thermodynamics.

*)
(*Date: 02/15/2021*)
Clear["Global`*"]

(**************************************************************************)
                (*Part 1 -- The Multiplicity of Macrostates*)
(**************************************************************************)



(*Subroutines in for multiplicity function*)
totalStatesM[N_] := 2 ^ N 
multiplicity[N_, n_] := N! / (n! * (N-n)!)
probability[\[CapitalOmega]_, M_] := \[CapitalOmega] / M
entropy[\[CapitalOmega]_] := 1 * Log[\[CapitalOmega]]

(*
:Purpose: Compute multiplicity table with variable macrostates.

:param macrostatesN_: Number of macrostates in the system.
:param state1_: Define a state a particle can have (e.g. "H", 1, "L", etc.).
:param state2_: Define another state a particle can have.
:param q3_: Should be True for question 3 and changes output for microstates column.

:return: <Table object> representing multiplicity
*)
computeMultiplicityTable[macrostatesN_, state1_, state2_, q3_] := (
	multiplicityTable = Table[0, macrostatesN+1, 5];
	pos = 1;
	For[n = macrostatesN, n >= 0, n--,  (*'n' the number of heads for a given macrostate varies*)
		(*Initial macrostate*)
		macroList = List[];
	
		(*Write "H" to list based on 'n'*)
		For[i = 1, i <= n, i++,
			AppendTo[macroList, state1];
		];
	
		(*Complete list with "T"s*)
		For[i = 1, i <= macrostatesN-n, i++,
			AppendTo[macroList, state2];
		];
	
		(*List of microstates corresponding to that macrostate*)
		microList = Permutations[macroList];
	
		(*Concatenated microList*)
		concatenatedMicroList = List[];
		microListDims = Dimensions[microList];
		For[row = 1, row <= microListDims[[1]], row++,
			AppendTo[concatenatedMicroList, StringRiffle[microList[[row]], ""]];
		];
		
		(*Set microstate str for table*)
		omega = multiplicity[macrostatesN, n];
		If [q3,
			(*TRUE*) 
			on = 0;
			off = 0;
			(*Iterate through macrostate string and count number of 'on' and 'off' states*)
			For [i = 1, i <= macrostatesN, i++,
				If [ToString[microList[[1, i]]] == state1, on++, off++]
			];
			(*Set the on/off string*)
			resultString = ToString[on] <> "/" <> ToString[off];,
			
			(*FALSE*)
			(*Create column microstate string*)
			resultString = "";
			For[i = 1, i <= omega, i++,
				resultString = StringJoin[resultString, ToString[Part[concatenatedMicroList, i]]];
				If[i != omega, resultString = StringJoin[resultString, "\n"]];
			];
		]
		
		(*Create the table*)
		capitalM = totalStatesM[macrostatesN];
		multiplicityTable[[pos, 1]] = n;                               (*Macrostate*)
		multiplicityTable[[pos,2]] = resultString;             (*Microstates*)
		multiplicityTable[[pos, 3]] = omega;                           
		multiplicityTable[[pos, 4]] = probability[omega, capitalM];    (*Probabiltity*)
		multiplicityTable[[pos, 5]] = SetPrecision[entropy[omega], 4]; (*Entropy*)
		pos++;
	];

	(*Return the table*)
	headerTable = Prepend[multiplicityTable, {"Macrostate", "Microstates", "Multiplicity", "Probability", "Entropy"}];
	Return[finalTable = Transpose[headerTable]];
)


Print["(**************************************************************************)
                (*Part 1 -- The Multiplicity of Macrostates*)
(**************************************************************************)"]
tenCoinTable = computeMultiplicityTable[10, "H", "T", False];
Grid[tenCoinTable, Frame->All]
Print["Is the 10 coin distribution more peaked than the 4 coin distribution?\n
Yes, the 10 coin distribution is significantly more peaked as 10 coins lends
to a greater number of microstates per macrostate."]


(**************************************************************************)
      (*Part 2 -- Evolution of Macrostates for a non-isolated system*)
(**************************************************************************)
(*Coin flip function*)
coinFlip := RandomChoice[{"H","T"}]

(*Lists for results*)
yNumHeads = List[];       (*Counts number of heads after a given toss*)
yNumTails = List[];       (*Counts number of tails after a given toss*)
xNumTosses = Range[100];  (*1 to 100th toss list*)


Print["(**************************************************************************)
      (*Part 2 -- Evolution of Macrostates for a non-isolated system*)
(**************************************************************************)"]

Print["What macrostate will be reached in the process of flipping N ordered
coins randomly starting from the NH state?"]

(*Do 100 tosses*)
For [toss = 1, toss <= 100, toss++,
	cntHeads = 0;
	cntTails = 0;
	(*Do 100 coin flips*)
	For [flip = 1, flip <= 100, flip++,
		If [coinFlip == "H", cntHeads++, cntTails++];
	];
	(*Append to parallel lists*)
	AppendTo[yNumHeads, cntHeads];
	AppendTo[yNumTails, cntTails];
]

(*Graph results*)
headData = Transpose[{xNumTosses, yNumHeads}]; 
tailData = Transpose[{xNumTosses, yNumTails}];
ListPlot [
	{headData, tailData},
	PlotRange->Full,
	ImageSize->{750, 750},
	PlotLegends->{"Heads",
				 "Tails"},
	PlotLabel->"Number of Coin Faces Appearing from 100 Flips vs. Number of Simulations (Coin Tosses)",
	AxesLabel->{"Number of Coin Tosses", "Number of Face"}
]

(*Calculate entropy*)
yEntropyList = List[];
For [i = 1, i <= Length[yNumHeads], i++,
	AppendTo[yEntropyList, entropy[multiplicity[100, Part[yNumHeads, i]]]];
]


(*Graph entropy*)
entropyData = Transpose[{xNumTosses, yEntropyList}];
ListPlot[
	entropyData,
	ImageSize->Large,
	PlotLabel->"Entropy vs. Number of Simulations (Coin Tosses)",
	AxesLabel->{"Number of Coin Tosses", "S = k ln(\[CapitalOmega])"}
]

Print["What pattern can you find during the process (how does the macrostate evolve)? 
Is there a tendency toward equilibrium?\nThe macrostate generally oscillates around the 
expected value for the macrostate which is 50.\n\nWhat happens in terms of entropy?\nThe entropy
also seems to remain in a constrained region on the graph, which is to be expected if multiplicity
is somewhat consistent.
"]



(**************************************************************************)
     (*Part 3 -- Dependence on Statistical Multiplicity on the volume*)
(**************************************************************************)
Print["(**************************************************************************)
     (*Part 3 -- Dependence of Statistical Multiplicity on the volume*)
(**************************************************************************)"]
Print["Compute L/R Multiplicity for N=4"]
volumeTable = computeMultiplicityTable[4, "L", "R", True];
Grid[volumeTable, Frame->All]
Print["Which macrostate is the largest multiplicity?
n=2 macrostate has the largest multiplicity of 6."]
Print["\nDo your findings support that particles occupy the entire volume after removing a restraint?
Yes, my findings support the statistical nature of the observation that
particles spread out in the available volume."]


(**************************************************************************)
                        (*Part 4 -- Gas diffusion*)
(**************************************************************************)
(*Initial Lists*)
initMacrostate = {"A", "A", "A", "A", "B", "B", "B", "B"};
allPermutations = Permutations[initMacrostate];

(*This will be a list of strings holding all {?A, ?B/?A, ?B} format*)
macrostateList = List[];

(*Each row in allPermutations is a list with 8 columns containing permutations of AAAABBBB*)
For [row = 1, row <= Length[allPermutations], row++,
	(*Counting macrostates
	OnLeft is defined as columns 1-4
	OnRight is defined as columns 5-8
	*)
	cntAOnLeft = 0;    
	cntAOnRight = 0;
	cntBOnLeft = 0;
	cntBOnRight = 0;
	
	(*Macrostate string of "?A, ?B/?A, ?B" format*)
	macrostateString = "";
	
	(*Count particles on left side*)
	For [letter = 1, letter <= 4, letter++,
		If [ToString[allPermutations[[row, letter]]] == "A", 
		cntAOnLeft++, 
		cntBOnLeft++
		];
	];
	
	(*Count particles on right side*)
	For [letter = 5, letter <= 8, letter++,
		If [ToString[allPermutations[[row, letter]]] == "A", 
			cntAOnRight++, 
			cntBOnRight++
		];
	];
	
	(*Build the macrostateList*)
	toAdd = ToString[cntAOnLeft] <> "A, " <> ToString[cntBOnLeft] <> "B/" <> ToString[cntAOnRight] <> "A, " <> ToString[cntBOnRight] <> "B";
	AppendTo[macrostateList, toAdd];
]

(*Count multiplicity*)
noDuplicateMacrostatesList = DeleteDuplicates[macrostateList];
multiplicityCounterList = List[];                               (*Parallel to noDuplicateMacrostatesLists*)
(*Iterate through noDuplicateMacrostatesList and the total macrostateList to count multiplicity*)
For [i = 1, i <= Length[noDuplicateMacrostatesList], i++,
	multiplicityCounter = 0;
	For [mstate = 1, mstate <= Length[macrostateList], mstate++,
		(*If the current macrostateList element matches the current noDuplicateMacrostatesList element, multiplicityCounter increments*)
		If [Part[macrostateList, mstate] == Part[noDuplicateMacrostatesList, i],
			multiplicityCounter++;
		];
	];
	AppendTo[multiplicityCounterList, multiplicityCounter];
]

Print["(**************************************************************************)
                        (*Part 4 -- Gas diffusion*)
(**************************************************************************)"]

(*Format and display table*)
question4Table = Prepend[Transpose[{noDuplicateMacrostatesList, multiplicityCounterList}], {"Macrostate", "Multiplicity"}];
Grid[question4Table, Frame->All]

Print["The equilibrium composition is 2A, 2B/2A, 2B based on the table above."]


(**************************************************************************)
                          (*Part 5 -- Heat Flow*)
(**************************************************************************)
Print["(**************************************************************************)
                          (*Part 5 -- Heat Flow*)
(**************************************************************************)"]
Print["What are the total energies and multiplicities of the isolated systems B {0111} and A {0011}?"]
Print["The blue and gray columns represent system B and A, respectively. The total energy
of system B (n=3) is 3 and the multiplicity is 4. The total energy of system
A (n=2) is 2 and the multiplicity is 6.\n"]
isolatedSystems = computeMultiplicityTable[4, "1", "0", False];
Grid[isolatedSystems, Frame->All, Background -> {{None, None, LightBlue, LightGray, None, None}}]
Print["Calculate the multiplicity of the overall system with total energy of 5 units?
What is the expected direction of heat flow?"]
Print["The brown column represents the combined system with total energy of 5. 
The multiplicity of this system is 56 from the table. The direction of heat flow is
such that the total energy of the system is 4. This means heat from either A or B  must be 
lost to the surroundings such that the overall system only has 4 units of energy."]
combinedSystems = computeMultiplicityTable[8, "1", "0", False];
Grid[combinedSystems, Frame->All, Background->{{None, None, None,None, LightBrown}}]



