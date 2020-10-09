(* $Id: QuCalc.m.pre 2.2 2000-06-22 10:14:12-04 dumais Exp $ *)

(*
Auteur:
    Paul Dumais
    dumais@iro.umontreal.ca

    Département d'informatique et recherche opérationnelle
    Université de Montréal
    CP 6128, succ. Centre-Ville
    Montréal, QC
    H3C 3J7
    Canada


"\"Mathematica\" is a registered trademark of \"Wolfram Research\".  QuCalc\nmay and must be distributed freely.  It must be distributed without\nmodifications, including the name of the authors and this message."
*)


BeginPackage["QuCalc`", {"LinearAlgebra`Orthogonalization`"}]
  rcsRevision = StringReplace["$Revision: 2.2 $", "$" -> ""];
  rcsDate     = StringReplace["$Date: 2000-06-22 10:14:12-04 $", "$" -> ""];
  Print[
    "Welcome to QuCalc, the quantum computation package for Mathematica 4.0.\n\nAuthor: Paul Dumais, dumais@iro.umontreal.ca\nHelp messages: Paul Dumais and Hugo Touchette",
    "\n\n",
    "\"Mathematica\" is a registered trademark of \"Wolfram Research\".  QuCalc\nmay and must be distributed freely.  It must be distributed without\nmodifications, including the name of the authors and this message.",
    "\n\n",
    rcsRevision,
    "\n",
    rcsDate,
    "\n\n",
    "Type \"?intro\" or \"?list\" for help."
  ];

  $PrePrint = MatrixForm;

  intro::usage = "General remarks relative to QuCalc\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nQuCalc is a Mathematica package whose pupose is to simulate and solve\nproblems related to quantum computing.  QuCalc has been developed by\nPaul Dumais at Laboratoire d'informatique th\[EAcute]orique et quantique\n(Montr\[EAcute]al University) and at the Crypto and Quantum Info Lab (McGill\nUniversity).\n\nA minimal knowledge of Mathematica is required in order to use the\nQuCalc package effectively. It is useful to recall that Mathematica\ndoes not incorporate in its kernel a data type for matrices: instead,\none must use a \"list of lists\" to represent matrices. It is strongly\nrecommended to carry out and analyze the examples suggested in the\nhelp sections of the this package. These examples form an essential\ncomplement to the main help text.\n\nThe paradigm adopted for the QuCalc Mathematica package is that of\nconverting states and operators into vectors and matrices. Kets are\nrepresented as column vectors written in the canonical basis (i.e. the\nbasis into which the Pauli matrix \\sigma_z is diagonal); unitary\ntransformations are represented as square matrices; mixed states\ncorrespond to density matrices; etc.\n\nThe identifiers defined when QuCalc is loaded in Mathematica all begin\nwith lowercase letters. Following standard usage, identifiers\nbeginning with a capital letter are those of Mathematica.  Some\nidentifiers of Mathematica are also overloaded in QuCalc such as\n\"CircleTimes\".\n\nWhenever you don't know how to interpret a result x returned by\nQuCalc, try FullForm[x]. This will tell you how QuCalc \"interprets\"\nthe result.  If x appears more complicated than it should be, test on\nit several simplification tools included in Mathematica, such as\nSimplify[x], FullSimplify[x], TrigReduce[x], ComplexExpand[x], ...\n\nIn general, functions do not test their inputs. When invalid inputs\nare supplied to a function, the result is not defined (error message\nfrom Mathematica, bizarre answer, etc.).\n\nTo obtain help pertaining to a specific function f, type \"?f\". To\nobtain the list of the subjects for which help is supplied, type\n\"?list\".\n\nErrors or bugs found in the QuCalc package, as well as questions or\ncomments should be addressed to Paul Dumais, dumais@iro.umontreal.ca.";
  list::usage  = "General information:\n\tintro           list\n\nMathematica identifiers overloaded in QuCalc:\n\tCircleDot       CirclePlus      CircleTimes     Dot\n\tOverBar         OverTilde       Power           SuperDagger\n\tVee             Wedge\n\nData types:\n\tens             ket             schmidt         sqo\n\tstate           supop           unit\n\nTests:\n\tensQ            ketQ            sqoQ            stateQ\n\tsupopQ          unitQ\n\nConstants:\n\tbb84            cnot            knot            mm\n\tmmm             not             phim            phip\n\tpsim            psip            sigx            sigy\n\tsigz            wh              xm              xp\n\tym              yp              zm              zp\n\nOther functions:\n\tanc             band            bits            block\n\tbnot            bor             bscal           bxor\n\tcircuit         cycle           ctrl            dag\n\tdotexp          eigen           eigenVal        eigenVect\n\tentropy         fgate           fidelity        fourier\n\tgate            id              kron            krondiv\n\tkronexp         ktrl            maxmix          phase\n\trandomUnit      rotx            roty            rotz\n\tswap            trout           unvec           vec";

  CircleDot::usage   = "bscal[s1, s2]\nCircleDot[s1, s2]\n\tReturns the boolean scalar product of the binary lists s1 and\n\ts2.\n\n\tThis function can be called using the Mathematica infix\n\toperator \"CircleDot\", by typing: \"ESC\" \"c\" \".\" \"ESC\".\n\n\tTry:\n\t\tbscal[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tband, bits, bnot, bor, bxor";
  CirclePlus::usage  = "bxor[s1, s2, ...]\nCirclePlus[s1, s2, ...]\n\tReturns the result of a bitwise \"exclusive or\" operation\n\tperformed over bit lists.\n\n\tThis function may be called in infix notation using the\n\tMathematica operator \"CirclePlus\" by typing: \"ESC\" \"c\" \"+\"\n\t\"ESC\".\n\n\tTry:\n\t\tbxor[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tband, bits, bnot, bor, bscal";
  CircleTimes::usage = "kron[x, y, ...]\nCircleTimes[x, y, ...]\n\tReturns the Kronecker product of x,y,... The arguments can be\n\tmatrices, \"ket\", \"unit\", \"state\", etc., as long as the product is\n\tmeaningful.\n\n\tThis function in Mathematica can be called in infix notation\n\tusing the operator \"CircleTimes\" by typing: \"ESC\" \"c\" \"*\"\n\t\"ESC\".  Note that this operator has a weaker precedence than the\n\t\"Dot\" product, which denoted with \".\".\n\n\tTry:\n\t\tkron[zp, xp]\n\t\tkron[wh, cnot]\n\t\tkron[zm, maxmix[2]]\n\t\tkron[mm, mm]\n\t\tkron[mm, wh]\n\n\tSee also:\n\t\tDot, circuit, krondiv, kronexp";
  Dot::usage         = "Dot[x, y, ...]\nx.y. ...\n\tDot gives the product (composition, application) of its arguments.\n\tThe arguments may be in the form of matrices, \"ket\", \"unit\",\n\t\"state\", etc, as long as the product is defined.\n\n\tNote that if u is a \"unit\" and r is a \"state\", then u.r returns the\n\tdensity matrix resulting from the application of u on r. Do not type\n\tu.r.dag[u].\n\n\tTry:\n\t\tdag[xp] . zp\n\t\tsigz . zm\n\t\tsigz . sigz\n\t\twh . state[zm]\n\t\tmmm . xp\n\t\twh . mm . wh\n\n\tSee also:\n\t\tcircuit, dotexp, kron";
  OverBar::usage     = "OverBar[x]\n\tReturns the complex conjugate of x. It is equivalent to the\n\tMathematica function \"Conjugate\".\n\n\tThis function can be entered using 2D notation by typing:\n\t\"CTRL-&\" \"_\".\n\n\tTry:\n\t\tOverBar[I]\n\n\tSee also:\n\t\tdag";
  OverTilde::usage   = "bnot[s]\nOverTilde[s]\n\tReturns the binary complement (negation) of a bit list.\n\n\tThis function can be entered in 2D notation by typing:\n\t\"CTRL-&\" \"~\".\n\n\tTry:\n\t\tbnot[{0,1,0}]\n\n\tSee also:\n\t\tband, bits, bor, bscal, bxor";
  \[Placeholder];
  Power::usage       = "dotexp[x, n]\nPower[x, n]\nx^n\n\tThese operators return the result of the n-fold product x.x. ... .x.\n\tThe argument x can be of type \"unit\", \"supop\", \"sqo\", etc., as long\n\tas the result is meaningful.\n\n\tTry:\n\t\tSimplify[rotx[t]^3]\n\n\tSee also:\n\t\tDot, kronexp";
  SuperDagger::usage = "dag[x]\nSuperDagger[x]\n\tReturns the conjugate transpose of x.\n\n\tThis function can be called in 2D notation by typing:\n\t\"CTRL-^\" \"ESC\" \"d\" \"g\" \"ESC\".\n\n\tTry:\n\t\tdag[wh] == wh";
  Vee::usage         = "bor[s1, s2, ...]\nVee[s1, s2, ...]\n\tPerforms a bitwise \"or\" operation over bit lists.\n\n\tThis functions can be called in \"infix\" notation using the\n\tMathematica operator \"Vee\" by typing: \"ESC\" \"v\" \"ESC\".\n\n\tTry:\n\t\tbor[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tband, bits, bnot, bscal, bxor";
  Wedge::usage       = "band[s1, s2, ...]\nWedge[s1, s2, ...]\n\tPerforms a bitwise \"and\" operation over binary strings, i.e.,\n\tlists of 0's and 1's.\n\n\tThis function may be entered in infix notation using the\n\toperator \"Wedge\" of Mathematica, by typing: \"ESC\" \"^\" \"ESC\".\n\n\tTry:\n\t\tband[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tbits, bnot, bor, bscal, bxor";

  anc::usage        = "anc[nA, nB, nC]\n\tReturns a data structure of type \"sqo\" representing the\n\toperation of adding a constant ancillary state |0> of\n\tdimension nB in between two Hilbert spaces of dimensions nA\n\tand nC.\n\n\tTry:\n\t\tr = anc[2,2,2] . ket[\"11\"]\n\t\teigen[r]\n\n\tSee also:\n\t\tsqo, trout";
  band::usage       = "band[s1, s2, ...]\nWedge[s1, s2, ...]\n\tPerforms a bitwise \"and\" operation over binary strings, i.e.,\n\tlists of 0's and 1's.\n\n\tThis function may be entered in infix notation using the\n\toperator \"Wedge\" of Mathematica, by typing: \"ESC\" \"^\" \"ESC\".\n\n\tTry:\n\t\tband[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tbits, bnot, bor, bscal, bxor";
  bb84::usage       = "bb84\n\tConstant of type \"ens\" representing a statistical ensemble of 4\n\tequiprobable states corresponding to the qubits 0 and 1 written in\n\tthe canonical and diagonal bases.\n\n\tTry:\n\t\tbb84\n\t\twh . bb84\n\t\tstate[bb84]\n\n\tSee also:\n\t\tens, zp";
  bits::usage       = "bits[x, n]\n\tConverts the integer x into its corresponding binary\n\texpression (list of 0's and 1's).\n\n\tTry:\n\t\tket[bits[5,3]]\n\t\tket[bits[5,4]]\n\n\tSee also:\n\t\tband, bnot, bor, bscal, bxor";
  block::usage      = "block[u, v]\n\tBlock product of the unitary transformations u and v.\n\n\tTry:\n\t\tblock[wh, wh]\n\n\tSee also:\n\t\tkron, unit";
  bnot::usage       = "bnot[s]\nOverTilde[s]\n\tReturns the binary complement (negation) of a bit list.\n\n\tThis function can be entered in 2D notation by typing:\n\t\"CTRL-&\" \"~\".\n\n\tTry:\n\t\tbnot[{0,1,0}]\n\n\tSee also:\n\t\tband, bits, bor, bscal, bxor";
  bor::usage        = "bor[s1, s2, ...]\nVee[s1, s2, ...]\n\tPerforms a bitwise \"or\" operation over bit lists.\n\n\tThis functions can be called in \"infix\" notation using the\n\tMathematica operator \"Vee\" by typing: \"ESC\" \"v\" \"ESC\".\n\n\tTry:\n\t\tbor[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tband, bits, bnot, bscal, bxor";
  bscal::usage      = "bscal[s1, s2]\nCircleDot[s1, s2]\n\tReturns the boolean scalar product of the binary lists s1 and\n\ts2.\n\n\tThis function can be called using the Mathematica infix\n\toperator \"CircleDot\", by typing: \"ESC\" \"c\" \".\" \"ESC\".\n\n\tTry:\n\t\tbscal[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tband, bits, bnot, bor, bxor";
  bxor::usage       = "bxor[s1, s2, ...]\nCirclePlus[s1, s2, ...]\n\tReturns the result of a bitwise \"exclusive or\" operation\n\tperformed over bit lists.\n\n\tThis function may be called in infix notation using the\n\tMathematica operator \"CirclePlus\" by typing: \"ESC\" \"c\" \"+\"\n\t\"ESC\".\n\n\tTry:\n\t\tbxor[{0,1,1,1}, {0,0,1,1}]\n\n\tSee also:\n\t\tband, bits, bnot, bor, bscal";
  circuit::usage    = "circuit[m]\n\tThis function performs a series of compositions and Kronecker\n\tproducts associated with the elements of the matrix m. The\n\tlines of m represent the quantum wires which carry qubits from\n\tleft to right.  The columns of m represent the quantum gates\n\tto be applied on qubits.  In other words, the matrix m is a\n\tgate array.\n\n\tFor example, circuit[{{wh,id[]}, {not,not}}] returns a data\n\tstructure of type \"unit\" (a unitary transformation) equivalent\n\tto a circuit of two gates acting on two quantum wires. The\n\tfirst wire is tranformed by a Walsh-Hadamard (wh) operation\n\tfollowed by a negation (not) operation. The second wire, is\n\tleft unchanged by the first identity (id[]) operation, and is\n\tthen negated by a \"not\" operation.\n\n\tTo build a gate acting coherently on a subset of selected\n\tqubits (leaving unchanged the rest of the qubits), place the\n\tdesired transformation on the last qubit of the selected\n\tsubset in question. The other wires are assigned integers\n\tdescribing a \"chained list\" in top-down order.  For example,\n\t\tcircuit[{{id[]},{1},{3},{id[]},{id[]},{ctrl[ctrl[not]]}}]\n\trepresents a Toffoli gate (ctrl[ctrl[not]]), acting on the 2nd, 3rd\n\tand 6th wires of a quantum circuit of 6 qubits.\n\n\tThe function \"circuit\" is very useful when combined with the\n\t2D notation of Mathematica commonly used to create matrices:\n\t\"CRTL-RETURN\" to create the lines, and \"CTRL-,\" to create the\n\tcolumns. In that case, the circuit entered in Mathematica\n\tlooks exactly as in standard gate array notation. Note also\n\tthat the \"Placeholder\" of Mathematica (the little blank square\n\tappearing in 2D notation), stands for the identity\n\ttransformation. It is thus unnecessary to fill out all the\n\tentries of a circuit.\n\n\tTry:\n\t\tcircuit[{{ 1,     1,   1  },\n\t\t         {cnot, knot, cnot}}]\n\n\t\tcircuit[{{ 1,   mm},\n\t\t         {cnot, mm}}]\n\n\tSee also:\n\t\tDot, ctrl, cycle, fgate, gate, kron, ktrl, swap";
  cnot::usage       = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  ctrl::usage       = "ctrl[u]\nctrl[n, u]\nktrl[u]\n\tVarious forms of controlled gates (conditional operations).\n\n\tctrl[u]: controlled u gate;\n\tctrl[n,u]: u gate controlled by n wires;\n\tktrl[u]: inverted controlled u gate.\n\n\tTry:\n\t\tctrl[not] == cnot\n\t\tctrl[2,not]\n\t\tktrl[wh]\n\n\tSee also:\n\t\tcircuit, cnot, knot";
  cycle::usage      = "cycle[n,i,j]\n\tUnitary transformation acting on n qubits which performs a circular\n\tpermutation of the qubits i to j (in top-down order).\n\n\tTry:\n\t\tc = cycle[3,1,3]\n\t\tc . ket[\"010\"] == ket[\"001\"]\n\t\tPower[c,2] == dag[c]\n\n\tSee also:\n\t\tcircuit";
  dag::usage        = "dag[x]\nSuperDagger[x]\n\tReturns the conjugate transpose of x.\n\n\tThis function can be called in 2D notation by typing:\n\t\"CTRL-^\" \"ESC\" \"d\" \"g\" \"ESC\".\n\n\tTry:\n\t\tdag[wh] == wh";
  dotexp::usage     = "dotexp[x, n]\nPower[x, n]\nx^n\n\tThese operators return the result of the n-fold product x.x. ... .x.\n\tThe argument x can be of type \"unit\", \"supop\", \"sqo\", etc., as long\n\tas the result is meaningful.\n\n\tTry:\n\t\tSimplify[rotx[t]^3]\n\n\tSee also:\n\t\tDot, kronexp";
  eigen::usage      = "eigen[r]\neigenVal[r]\neigenVect[r]\n\tReturns the eigenvalues and eigenvectors of a state.\n\n\teigen[r]: returns a statistical ensemble (type \"ens\") formed with\n\tall the pairs eigenvalues-eigenvectors of r.\n\n\teigenVal[r]: returns a diagonal matrix whose elements are the\n\teigenvalues of the state r, in the data type \"state\".\n\n\teigenVect[r]: returns a unitary matrix whose columns are the\n\teigevectors of r.\n\n\tTry:\n\t\teigen[maxmix[2]]\n\t\teigenVect[state[xp]]\n\t\teigenVal[state[xp]]\n\t\teigen[state[yp]] == yp\n\n\tSee also:\n\t\tens, schmidt, state, unit";
  eigenVal::usage   = "eigen[r]\neigenVal[r]\neigenVect[r]\n\tReturns the eigenvalues and eigenvectors of a state.\n\n\teigen[r]: returns a statistical ensemble (type \"ens\") formed with\n\tall the pairs eigenvalues-eigenvectors of r.\n\n\teigenVal[r]: returns a diagonal matrix whose elements are the\n\teigenvalues of the state r, in the data type \"state\".\n\n\teigenVect[r]: returns a unitary matrix whose columns are the\n\teigevectors of r.\n\n\tTry:\n\t\teigen[maxmix[2]]\n\t\teigenVect[state[xp]]\n\t\teigenVal[state[xp]]\n\t\teigen[state[yp]] == yp\n\n\tSee also:\n\t\tens, schmidt, state, unit";
  eigenVect::usage  = "eigen[r]\neigenVal[r]\neigenVect[r]\n\tReturns the eigenvalues and eigenvectors of a state.\n\n\teigen[r]: returns a statistical ensemble (type \"ens\") formed with\n\tall the pairs eigenvalues-eigenvectors of r.\n\n\teigenVal[r]: returns a diagonal matrix whose elements are the\n\teigenvalues of the state r, in the data type \"state\".\n\n\teigenVect[r]: returns a unitary matrix whose columns are the\n\teigevectors of r.\n\n\tTry:\n\t\teigen[maxmix[2]]\n\t\teigenVect[state[xp]]\n\t\teigenVal[state[xp]]\n\t\teigen[state[yp]] == yp\n\n\tSee also:\n\t\tens, schmidt, state, unit";
  ens::usage        = "ens[...]\n\tData type representing a statistical ensemble of states.\n\n\tA valid \"ens\" object contains a list of pairs {p,r} where p is\n\ta real number between 0 and 1 (a probability), and r is an\n\tobject of type \"ket\" or \"state\". The p's must add to unity.\n\n\tAn \"ens\" object containing a trivial list of only one pair {p,r}\n\twith p=1 is automatically converted to a \"state\" or \"ket\" r.\n\n\tTry:\n\t\tbb84\n\t\tFullForm[bb84]\n\t\tmmm . xp\n\t\tFullForm[mmm . xp]\n\t\teigen[maxmix[2]]\n\t\tFullForm[eigen[maxmix[2]]]\n\n\tSee also:\n\t\tbb84, eigen, ensQ, ket, state, sqo";
  ensQ::usage       = "ensQ[x]\nketQ[x]\nsqoQ[x]\nstateQ[x]\nsupopQ[x]\nunitQ[x]\n\tBoolean functions testing the validity of x according to the type\n\tstructure associated with the function name.\n\n\tTry:\n\t\tunitQ[wh]\n\t\tSimplify[unitQ[fourier[3]]]\n\t\tketQ[ket[zp + zm]]\n\n\tSee also:\n\t\tens, ket, sqo, state, supop, unit";
  entropy::usage    = "entropy[r]\n\tNumerical calculation of the von Neumann entropy of the state r.\n\n\tTry:\n\t\tentropy[state[bb84]]\n\n\tSee also:\n\t\tfidelity";
  fgate::usage      = "fgate[m, n, f]\n\tReturns a unitary transformation, acting on m+n qubits, that\n\tis equivalent to the function f acting on classical inputs.\n\n\tThe function f is a Mathematica \"pure function\". It takes\n\ton input a list of m bits and must return a list of n\n\tbits. Thus, f is a (classical) function of m bits to n bits.\n\n\tTry:\n\t\tfgate[1,1,(#)&]\n\t\tfgate[1,2,({#[[1]], #[[1]]})&]\n\n\tSee also:\n\t\tcircuit, gate";
  fidelity::usage   = "fidelity[r1, r2]\n\tNumerical calculation of the fidelity function of two mixed states.\n\n\tTry:\n\t\tfidelity[state[zp], state[xp]]\n\n\tSee also:\n\t\tentropy";
  fourier::usage    = "fourier[m]\n\tm-dimensional Fourier transform.\n\n\tTry:\n\t\tfourier[4]\n\n\tSee also:\n\t\tcircuit";
  gate::usage       = "gate[n, b, f]\ngate[n, f]\n\tReturns a unitary transformation of dimension b^n (2^n if b is\n\tommitted) described by the function f.\n\n\tThe function f is a Mathematica \"pure function\". It takes as\n\tinput n integers between 0 and b-1, and returns a \"ket\". For\n\tdifferent values in its arguments f must return orthogonal kets.\n\n\tTry:\n\t\tgate[2,(ket[{#1, bxor[#1,#2]}])&]\n\n\tSee also:\n\t\tcircuit, fgate";
  id::usage         = "id[n]\nid[]\n\tIdentity transformation on n qubits (1 if n is ommited).\n\n\tTry:\n\t\tid[8]\n\n\tSee also:\n\t\tcircuit";
  ket::usage        = "ket[...]\n\tData type representing a pure state.\n\n\tA well-defined \"ket\" contains a matrix of dimension n x 1, and must\n\tbe of norm 1.\n\n\tNote that in 2D notation, lines of a matrix can be created using\n\t\"CTRL-RETURN\".\n\n\tTo obtain a \"bra\" from a \"ket\", use dag[ket[...]]. There is no\n\tdata type \"bra\" in QuCalc.\n\n\tThe expression ket[0] (or ket[1]) is converted in such a way to\n\trepresent the bit 0 (or 1) in the canonical basis of a 2-dimensional\n\tHilbert space, i.e., a qubit. To obtain a register of many qubits,\n\tone can concatenate several bits into a list or into a character\n\tstring. To obtain a quantum register whose dimension is a power\n\tnot reducible in base 2, one can add as a subscript the value of the\n\tbase in question. See the examples below.\n\n\tIn 2D notation, \"Subscript[x,y]\" may be entered by typing: \"x\"\n\t\"CTRL-_\" \"y\".\n\n\tTry:\n\t\tket[0]\n\t\tFullForm[ket[0]]\n\t\tFullForm[ket[0][[1]]]\n\t\tket[{0,1,0}]\n\t\tket[\"010\"]\n\t\tket[Subscript[2, 3]]\n\t\tket[Subscript[{0,2,0}, 3]]\n\t\tket[Subscript[\"020\", 3]]\n\t\t-zm\n\t\tket[(1/Sqrt[2])(zp+zm)]\n\n\tSee also:\n\t\tketQ, phip, zp";
  ketQ::usage       = "ensQ[x]\nketQ[x]\nsqoQ[x]\nstateQ[x]\nsupopQ[x]\nunitQ[x]\n\tBoolean functions testing the validity of x according to the type\n\tstructure associated with the function name.\n\n\tTry:\n\t\tunitQ[wh]\n\t\tSimplify[unitQ[fourier[3]]]\n\t\tketQ[ket[zp + zm]]\n\n\tSee also:\n\t\tens, ket, sqo, state, supop, unit";
  knot::usage       = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  kron::usage       = "kron[x, y, ...]\nCircleTimes[x, y, ...]\n\tReturns the Kronecker product of x,y,... The arguments can be\n\tmatrices, \"ket\", \"unit\", \"state\", etc., as long as the product is\n\tmeaningful.\n\n\tThis function in Mathematica can be called in infix notation\n\tusing the operator \"CircleTimes\" by typing: \"ESC\" \"c\" \"*\"\n\t\"ESC\".  Note that this operator has a weaker precedence than the\n\t\"Dot\" product, which denoted with \".\".\n\n\tTry:\n\t\tkron[zp, xp]\n\t\tkron[wh, cnot]\n\t\tkron[zm, maxmix[2]]\n\t\tkron[mm, mm]\n\t\tkron[mm, wh]\n\n\tSee also:\n\t\tDot, circuit, krondiv, kronexp";
  krondiv::usage    = "krondiv[v, w]\n\t\"Kronecker division\". krondiv returns a matrix y of dimension\n\tN/n x 1, such that v = kron[w,y], where v is a matrix of\n\tdimension N x 1, and w is matrix of dimension n x 1.\n\n\tThe function krondiv only accepts Mathematica matrices, and no other\n\ttypes of data, such as \"ket\".\n\n\tNote that lines of a column vector can be created in Mathematica\n\tusing the 2D notation: \"CTRL-RETURN\".\n\n\tTry:\n\t\tkrondiv[{{4},{5},{8},{10},{12},{15}}, {{1},{2},{3}}]\n\n\tSee also:\n\t\tkron, vec, unvec";
  kronexp::usage    = "kronexp[x, n]\n\tn-fold Kronecker product.\n\n\tTry:\n\t\tkronexp[wh, 3]\n\t\tkronexp[mm, 3]\n\t\tkronexp[zp, 3] == ket[\"000\"]\n\n\tSee also:\n\t\tkron, dotexp";
  ktrl::usage       = "ctrl[u]\nctrl[n, u]\nktrl[u]\n\tVarious forms of controlled gates (conditional operations).\n\n\tctrl[u]: controlled u gate;\n\tctrl[n,u]: u gate controlled by n wires;\n\tktrl[u]: inverted controlled u gate.\n\n\tTry:\n\t\tctrl[not] == cnot\n\t\tctrl[2,not]\n\t\tktrl[wh]\n\n\tSee also:\n\t\tcircuit, cnot, knot";
  maxmix::usage     = "maxmix[n]\n\tReturns the maximally mixed density matrix of dimension n.\n\n\tTry:\n\t\tmaxmix[3]\n\n\tSee also:\n\t\tstate";
  mm::usage         = "mm\n\tConstant of type \"supop\" representing a 1-qubit measuring process in\n\tthe canonical basis.\n\n\tTry:\n\t\tmm\n\t\tmm . xp\n\t\tmm . zp\n\t\tkron[mm, mm]\n\t\tmm . mm == mm\n\t\twh . mm . wh\n\n\tSee also:\n\t\tmmm, supop";
  mmm::usage        = "mmm\n\tConstant of type \"sqo\" representing a 1-qubit measuring process in\n\tthe canonical basis.\n\n\tTry:\n\t\tmmm\n\t\tmmm . xp\n\t\tmmm . zp\n\t\tkron[mmm, id[]]\n\t\tmmm . mmm == mm\n\n\tSee also:\n\t\tmm, sqo";
  not::usage        = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  phase::usage      = "phase[t]\nrotx[t]\nroty[t]\nrotz[t]\n\tphase[t] is a global phase change with angle t applied to 1 qubit.\n\n\trotx[t], roty[t] and rotz[t] are rotations around the axis\n\tcorresponding to the function name in the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tphase[t]\n\t\trotx[t]\n\t\troty[t]\n\t\trotz[t]\n\n\tSee also:\n\t\tunit, randomUnit";
  phim::usage       = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  phip::usage       = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  psim::usage       = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  psip::usage       = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  randomUnit::usage = "randomUnit[]\n\tReturns a random unitary transformation of dimension 2.  The\n\tdistribution of the outcome is based on a uniform random\n\tchoice of the Euler angles in the Bloch-Poincar\[EAcute] sphere.\n\n\tMore precisely, randomUnit[] is defined as follows:\n\t\tphase[t1] . rotz[t2] . roty[t3] . rotz[t4]\n\twhere t1, t2, t3, and t4 are independent and identically distributed\n\tuniform random variables chosen between 0 and 2 Pi.\n\n\tTry:\n\t\trandomUnit[] . zp\n\n\tSee also:\n\t\tphase, unit";
  rotx::usage       = "phase[t]\nrotx[t]\nroty[t]\nrotz[t]\n\tphase[t] is a global phase change with angle t applied to 1 qubit.\n\n\trotx[t], roty[t] and rotz[t] are rotations around the axis\n\tcorresponding to the function name in the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tphase[t]\n\t\trotx[t]\n\t\troty[t]\n\t\trotz[t]\n\n\tSee also:\n\t\tunit, randomUnit";
  roty::usage       = "phase[t]\nrotx[t]\nroty[t]\nrotz[t]\n\tphase[t] is a global phase change with angle t applied to 1 qubit.\n\n\trotx[t], roty[t] and rotz[t] are rotations around the axis\n\tcorresponding to the function name in the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tphase[t]\n\t\trotx[t]\n\t\troty[t]\n\t\trotz[t]\n\n\tSee also:\n\t\tunit, randomUnit";
  rotz::usage       = "phase[t]\nrotx[t]\nroty[t]\nrotz[t]\n\tphase[t] is a global phase change with angle t applied to 1 qubit.\n\n\trotx[t], roty[t] and rotz[t] are rotations around the axis\n\tcorresponding to the function name in the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tphase[t]\n\t\trotx[t]\n\t\troty[t]\n\t\trotz[t]\n\n\tSee also:\n\t\tunit, randomUnit";
  schmidt::usage    = "schmidt[...]\n\tData type representing a Schmidt decomposition of a pure state\n\tseen as a bipartite state.\n\n\tA well-defined \"schmidt\" contains a list of triplets {p, k1, k2},\n\twhere the p's are the probabilities (Schmidt coefficients) of the\n\tdecomposition. The k1's (the k2's) are objects of type \"ket\", which\n\tform an orthonormal basis (Schmidt basis) of Alice's (Bob's)\n\tsub-space.\n\n\tThe expression schmidt[n1,n2,k], where k is a \"ket\" of dimension\n\tn1*n2 is converted so as to represent a Schmidt decomposition of k\n\twith n1 and n2 as the respective dimensions of the spaces of the two\n\tparties (Alice and Bob).\n\n\tTry:\n\t\ts = schmidt[2,2,phip]\n\t\tFullForm[s]\n\t\tFullForm[s[[1]]]\n\t\tket[s]\n\n\tSee also:\n\t\tket, eigen";
  sigx::usage       = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  sigy::usage       = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  sigz::usage       = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  sqo::usage        = "sqo[...]\n\tData type representing a \"selective quantum operation\" (SQO).\n\n\tMathematically, a SQO is defined as a family of square or\n\trectangular matrices A_{i,j}, 1<=i<=m, 1<=j<=k_i, of dimension\n\ts x r, such that\n\t\t\\sum A_{i,j}^\\dagger A_{i,j}\n\tis equal to the r x r identity matrix. A SQO is a general\n\tquantum operation which includes, as special cases, unitary\n\ttransformations, superoperators, POVMs (Positive Operator\n\tValued Measure), trace-out's, or the effect of adding\n\tancillae. When applied to a mixed state \\rho, of arbitrary\n\tdimension r, a SQO returns with probability p_i the state\n\t\t\\sigma_i = (1/p_i) B_i\n\tof dimension s, where\n\t\tB_i = \\sum_{j} A_{i,j} \\rho A_{i,j}^\\dagger\n\t\tp_i = tr(B_i).\n\n\tIt is very unlikely that you will have to build a \"sqo\"\n\tentirely from scratch, except possibly for the case where a\n\tmeasurement along non-orthogonal axes has to be made. In\n\tgeneral, it is sufficient to combine the sqo's already defined\n\tin QuCalc with some unitary transformations to obtain a general\n\tquantum operation. To get acquainted with the details of a\n\t\"sqo\", refer to the examples given below.\n\n\tTry:\n\t\tmmm\n\t\tFullForm[mmm]\n\t\tq = kron[mmm, mm]\n\t\tFullForm[q]\n\t\tmmm . xp\n\n\tSee also:\n\t\tanc, mmm, sqoQ, supop, trout";
  sqoQ::usage       = "ensQ[x]\nketQ[x]\nsqoQ[x]\nstateQ[x]\nsupopQ[x]\nunitQ[x]\n\tBoolean functions testing the validity of x according to the type\n\tstructure associated with the function name.\n\n\tTry:\n\t\tunitQ[wh]\n\t\tSimplify[unitQ[fourier[3]]]\n\t\tketQ[ket[zp + zm]]\n\n\tSee also:\n\t\tens, ket, sqo, state, supop, unit";
  state::usage      = "state[...]\n\tData type representing (in general) a mixed state.\n\n\tA well-defined \"state\" must contain a positive hermitian square\n\tmatrix with trace equal to one.\n\n\tNote that in 2D notation, \"CTRL-RETURN\" can be used to create\n\tthe lines of a matrix, whereas \"CTRL-,\" is used to create the\n\tcolumns.\n\n\tAn object \"state\" representing a pure state is not automatically\n\tconverted into a \"ket\"; use \"eigen\" to force that conversion.\n\n\tTry:\n\t\tr = state[phip]\n\t\tFullForm[r]\n\t\tFullForm[r[[1]]]\n\t\tstate[(1/2)(state[phip] + state[phim])]\n\n\tSee also:\n\t\teigen, ket, maxmix, stateQ";
  stateQ::usage     = "ensQ[x]\nketQ[x]\nsqoQ[x]\nstateQ[x]\nsupopQ[x]\nunitQ[x]\n\tBoolean functions testing the validity of x according to the type\n\tstructure associated with the function name.\n\n\tTry:\n\t\tunitQ[wh]\n\t\tSimplify[unitQ[fourier[3]]]\n\t\tketQ[ket[zp + zm]]\n\n\tSee also:\n\t\tens, ket, sqo, state, supop, unit";
  supop::usage      = "supop[...]\n\tData type representing a superoperator.\n\n\tA well-defined \"supop\" must contain a list of square matrices A_j\n\tsuch that\n\t\t\\sum_j A_j^\\dagger A_j\n\tis equal to the identity matrix. When applied to a mixed state \\rho,\n\ta \"supop\" returns the mixed state\n\t\t\\sum_{j} A_j \\rho A_j^\\dagger\n\tas long as the dimensions of the matrices are compatible.\n\n\tIn 2D notation, \"CTRL-RETURN\" can be used to create the lines of a\n\tmatrix, whereas \"CTRL-,\" is used to create the columns.\n\n\tTry:\n\t\tmm\n\t\tFullForm[mm]\n\t\tmm . xp\n\n\tSee also:\n\t\tmm, sqo, supopQ";
  supopQ::usage     = "ensQ[x]\nketQ[x]\nsqoQ[x]\nstateQ[x]\nsupopQ[x]\nunitQ[x]\n\tBoolean functions testing the validity of x according to the type\n\tstructure associated with the function name.\n\n\tTry:\n\t\tunitQ[wh]\n\t\tSimplify[unitQ[fourier[3]]]\n\t\tketQ[ket[zp + zm]]\n\n\tSee also:\n\t\tens, ket, sqo, state, supop, unit";
  swap::usage       = "swap[n, i, j]\nswap[]\n\tReturns a unitary transformation acting on n qubits, and whose\n\taction is to swap qubits i and j.\n\n\tswap[] is equivalent to swap[2,1,2].\n\n\tTry:\n\t\tc = swap[4,1,3]\n\t\tc . ket[\"0110\"] == ket[\"1100\"]\n\t\tPower[c,2] == id[4]\n\n\tSee also:\n\t\tcircuit";
  trout::usage      = "trout[nA, nB, nC]\n\tThis functions returns a \"sqo\" data structure representing the\n\toperation of tracing out a Hilbert space of dimension nB in\n\tbetween 2 Hilbert spaces of dimension nA and nC.\n\n\tTry:\n\t\ttrout[2, 2, 1] . phip\n\t\ttrout[1, 2, 2] . phip\n\t\ttrout[2, 4, 2] . ket[\"0110\"] == state[ket[\"00\"]]\n\n\tSee also:\n\t\tsqo, anc";
  unit::usage       = "unit[...]\n\tData type corresponding to unitary transformations.\n\n\tA valid \"unit\" data is a square unitary matrix.\n\n\tNote that a matrix in 2D notation can be built using \"CRTL-RETURN\"\n\tto create lines and \"CTRL-,\" to create columns.\n\n\tTry:\n\t\tFullForm[sigz]\n\t\tFullForm[sigz[[1]]]\n\t\tI wh\n\t\tunit[Exp[(I Pi)/4] wh]\n\t\tunit[(1/Sqrt[3])(sigx + sigy + sigz)]\n\n\tSee also:\n\t\tcnot, fourier, id, knot, not, phase, randomUnit, rotx,\n\t\tsigx, unitQ, wh";
  unitQ::usage      = "ensQ[x]\nketQ[x]\nsqoQ[x]\nstateQ[x]\nsupopQ[x]\nunitQ[x]\n\tBoolean functions testing the validity of x according to the type\n\tstructure associated with the function name.\n\n\tTry:\n\t\tunitQ[wh]\n\t\tSimplify[unitQ[fourier[3]]]\n\t\tketQ[ket[zp + zm]]\n\n\tSee also:\n\t\tens, ket, sqo, state, supop, unit";
  unvec::usage      = "vec[x]\nunvec[v, n]\n\tThis function transforms a matrix x into a vector. The columns of\n\tthe matrix (from left to right) are concatenated in the vector in\n\ttop-down order.\n\n\tunvec[v, n] undo the transformation by re-transforming a vector v\n\tinto a matrix containing n lines.\n\n\tThese functions are defined for Mathematica matrices; they cannot be\n\tapplied to a \"ket\" or any other data types specific to QuCalc.\n\n\tIn 2D notation \"CTRL-RETURN\" can be used to create lines of a\n\tmatrix, whereas \"CTRL-,\" is used to create the columns.\n\n\tTry:\n\t\tv = vec[{{1,2,3},{4,5,6},{7,8,9},{10,11,12}}]\n\t\tFullForm[v]\n\t\tunvec[v, 2]\n\n\tSee also:\n\t\tkrondiv";
  vec::usage        = "vec[x]\nunvec[v, n]\n\tThis function transforms a matrix x into a vector. The columns of\n\tthe matrix (from left to right) are concatenated in the vector in\n\ttop-down order.\n\n\tunvec[v, n] undo the transformation by re-transforming a vector v\n\tinto a matrix containing n lines.\n\n\tThese functions are defined for Mathematica matrices; they cannot be\n\tapplied to a \"ket\" or any other data types specific to QuCalc.\n\n\tIn 2D notation \"CTRL-RETURN\" can be used to create lines of a\n\tmatrix, whereas \"CTRL-,\" is used to create the columns.\n\n\tTry:\n\t\tv = vec[{{1,2,3},{4,5,6},{7,8,9},{10,11,12}}]\n\t\tFullForm[v]\n\t\tunvec[v, 2]\n\n\tSee also:\n\t\tkrondiv";
  wh::usage         = "cnot\nknot\nnot\nsigx\nsigy\nsigz\nwh\n\tConstants of type \"unit\".\n\n\tcnot: controled-not;\n\tknot: inversed controled-not;\n\tnot: negation (of 1 qubit);\n\tsigx, sigy, sigz: Pauli matrices;\n\twh: Walsh-Hadamard transformation.\n\n\tTry:\n\t\tnot == sigx\n\n\tSee also:\n\t\tunit";
  xm::usage         = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  xp::usage         = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  ym::usage         = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  yp::usage         = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  zm::usage         = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";
  zp::usage         = "phim\nphip\npsim\npsip\nxm\nxp\nym\nyp\nzm\nzp\n\tConstants of type \"ket\".\n\n\tphim, phip, psim, psip: Bell's states.\n\n\txm, xp, ym, yp, zm, zp: Pure states corresponding to the poles and\n\tthe 4 cardinal points on the equator of the Bloch-Poincar\[EAcute] sphere.\n\n\tTry:\n\t\tcnot . kron[wh, id[]] . ket[\"00\"] == phip\n\n\tSee also:\n\t\tket";

  Begin["`Private`"]

    (*-----------
      Utile :
        kron (CircleTimes)
          Produit de Kronecker.

        vec
          Mise en vecteur d'une matrice, les colonnes de la matrice,
          prises de gauche à droite, sont concaténées de haut en bas.

        unvec[x, m]
          Remet un vecteur x en matrice de m lignes.

        krondiv[v,w]
          Trouve un vecteur x, de dim. N/n x 1, tel que v==kron[w,x],
          où v est un vecteur de dim. N x 1
          et w est un vecteur de dim. n x 1.

        kronexp
          Produit de Kronecker n fois.

        dotexp
          Produit n fois.

        dag (SuperDagger)
          Conjuguée transposée.

        bits
          Convertit un entier en binaire, i.e. une liste de 0 et 1.

        bnot (OverTilde)
          Complémente une liste de 0 et 1.

        bscal (CircleDot)
          Produit scalaire de deux chaînes de bits, i.e. deux listes de 0 et 1.

        bxor (CirclePlus), band (Wedge) et bor (Vee)
          Somme modulo 2 ("ou exclusif") de deux chaînes de bits, i.e. deux
          listes de 0 et 1, "et" et "ou".

        nullListQ
          Teste si tous les éléments d'une liste sont nuls.

        squareMatrixQ[x]
          Teste si x est une matrice carrée.

      --------------------------------------------------------------------*)

    kron[u_ /; MatrixQ[u], v_ /; MatrixQ[v]] := 
      Module[
        {w},
	w = Outer[Times, u, v]; 
	Partition[Flatten[Transpose[w, {1,3,2,4}]], Dimensions[w][[2]] Dimensions[w][[4]]]
      ];
    SetAttributes[kron, OneIdentity];
    kron[u_, v_, w__] := Fold[kron, u, {v, w}];
    CircleTimes = kron;

    vec[x_ /; MatrixQ[x]] := List /@ Flatten[Transpose[x]];

    unvec[x_ /; MatrixQ[x], m_Integer] := Transpose[Partition[Flatten[x], m]];

    krondiv[v_ /; MatrixQ[v], w_ /; MatrixQ[w]] :=
      Transpose[
        LinearSolve[
          w, Transpose[unvec[v, Dimensions[v][[1]]/Dimensions[w][[1]]]]
        ]
      ];

    kronexp[x_, n_Integer /; n > 0] := 
      Module[
        {i, xx},
        For[i=1; xx=x, i<n, i++, xx = kron[xx, x]];
        xx
      ];

    dotexp[x_, n_Integer /; n > 0] := 
      Module[
        {i, xx},
        For[i=1; xx=x, i<n, i++, xx = xx . x];
        xx
      ];

    dag[u_ /; MatrixQ[u]] := Transpose[Conjugate[u]];
    SuperDagger = dag;
    OverBar = Conjugate;

    band[l1_, l2_] := l1 l2;
    SetAttributes[band, OneIdentity];
    band[u_, v_, w__] := Fold[band, u, {v, w}];
    Wedge = band;

    bits[x_, n_] := IntegerDigits[x, 2, n];

    bnot[l_] := 1-l;
    OverTilde = bnot;

    bor[l1_, l2_] := Mod[l1 + l2 + l1 l2, 2];
    SetAttributes[bor, OneIdentity];
    bor[u_, v_, w__] := Fold[bor, u, {v, w}];
    Vee = bor;

    bscal[l1_, l2_] := Mod[Apply[Plus, l1 l2], 2];
    CircleDot = bscal;

    bxor[l1_, l2_] := Mod[l1 + l2, 2];
    SetAttributes[bxor, OneIdentity];
    bxor[u_, v_, w__] := Fold[bxor, u, {v, w}];
    CirclePlus = bxor;

    nullListQ[x_] := Fold[And, True, (# == 0)& /@ Flatten[x]];

    squareMatrixQ[x_] := MatrixQ[x] && (Dimensions[x][[1]] == Dimensions[x][[2]]);


    (*--------------------------------------------------------------------
      ket[_] : contient un vecteur colonne, matrice nx1, normalisé.
      bra[_] : est synonyme (et en aval de) de dag[ket[_]].
      --------------------------------------------------------------------*)

    ket /: dag[ket[k_]] := bra[k];

    ket /: Format[ket[k_ /; MatrixQ[k]]] := MatrixForm[k];
    bra /: Format[bra[k_ /; MatrixQ[k]]] := MatrixForm[dag[k]];

    ket /: kron[ket[k1_], ket[k2_]] := ket[kron[k1, k2]];
    bra /: kron[bra[k1_], bra[k2_]] := bra[kron[k1, k2]];

    ket /: Dot[bra[k1_], ket[k2_]] := (dag[k1] . k2)[[1,1]];
    ket /: Dot[ket[k1_], bra[k2_]] := sqMatr[k1 . dag[k2]];

    ket /: Dimensions[ket[k_]] := Dimensions[k];
    bra /: Dimensions[bra[k_]] := Dimensions[dag[k]];

    ket /: Times[Complex[0,1], ket[k_]]  := ket[(I) k];
    ket /: Times[-1, ket[k_]]            := ket[(-1) k];
    ket /: Times[Complex[0,-1], ket[k_]] := ket[(-I) k];
    ket /: Times[a_, ket[k_]] := colVect[a k];
    ket /: Plus[ket[k1_], ket[k2_]] := colVect[k1 + k2];
    ket /: Plus[ket[k1_], colVect[k2_]] := colVect[k1 + k2];
    ket /: Plus[colVect[k1_], ket[k2_]] := colVect[k1 + k2];
    colVect /: Plus[colVect[k1_], colVect[k2_]] := colVect[k1 + k2];
    colVect /: Times[a_, colVect[k_]] := colVect[a k];
    ket[colVect[k_]] := ket[k];
    ket[ket[k_]] := ket[k];

    bra /: Times[Complex[0,1], bra[k_]]  := bra[(I) k];
    bra /: Times[-1, bra[k_]]            := bra[(-1) k];
    bra /: Times[Complex[0,-1], bra[k_]] := bra[(-I) k];
    bra /: Times[a_, bra[k_]] := linVect[a k];
    bra /: Plus[bra[k1_], bra[k2_]] := linVect[k1 + k2];
    bra /: Plus[bra[k1_], linVect[k2_]] := linVect[k1 + k2];
    bra /: Plus[linVect[k1_], bra[k2_]] := linVect[k1 + k2];
    linVect /: Plus[linVect[k1_], linVect[k2_]] := linVect[k1 + k2];
    linVect /: Times[a_, linVect[k_]] := linVect[a k];
    bra[linVect[k_]] := bra[k];
    bra[bra[k_]] := bra[k];
    colVect /: Dot[linVect[k1_], colVect[k2_]] := (dag[k1] . k2)[[1,1]];

    ket[List[x__Integer]] := ket[Subscript[{x}, 2]];
    ket[Subscript[List[x__Integer], n_Integer]] := 
      Module[
        {i,l},
        l = {x};
        ket[Table[{If[FromDigits[l, n] == i, 1, 0]}, {i,0,(n^Length[l])-1}]]
      ];

    ket[s_String] := ket[Subscript[s, 2]];
    ket[Subscript[s_String, n_Integer]] :=
      Module[
        {i,l}, 
        l = ToExpression[Characters[s]]; 
        ket[Table[{If[FromDigits[l, n] == i, 1, 0]}, {i,0,(n^Length[l])-1}]]
      ];

    ket[0] := ket[{{1},{0}}];
    ket[1] := ket[{{0},{1}}];
    ket[Subscript[i_, n_Integer] /; ((0 <= i) && (i < n))] := ket[Table[{If[i==j,1,0]}, {j,0,n-1}]];

    psim = ket[(1/Sqrt[2]) {{0}, {1}, {-1}, {0}}];
    psip = ket[(1/Sqrt[2]) {{0}, {1}, {1},  {0}}];
    phim = ket[(1/Sqrt[2]) {{1}, {0}, {0}, {-1}}];
    phip = ket[(1/Sqrt[2]) {{1}, {0}, {0},  {1}}];
    zp = ket[{{1},{0}}];
    zm = ket[{{0},{1}}];
    xp = ket[(1/Sqrt[2]) {{1},  {1}}];
    xm = ket[(1/Sqrt[2]) {{1}, {-1}}];
    yp = ket[(1/Sqrt[2]) {{1},  {I}}];
    ym = ket[(1/Sqrt[2]) {{I},  {1}}];

    ketQ[k_ket] := dag[k] . k == 1;


    (*--------------------------------------------------------------------
      unit[_] : transformation unitaire, contient une matrice carree
      unitaire.
      --------------------------------------------------------------------*)

    unit /: Format[unit[u_ /; MatrixQ[u]]] := MatrixForm[u];

    unitQ[unit[u_]] := squareMatrixQ[u] && (dag[u] . u == IdentityMatrix[Dimensions[u][[1]]]);

    unit /: dag[unit[u_]] := unit[dag[u]];
    unit /: Det[unit[u_]] := Det[u];

    unit /: kron[unit[u_], unit[v_]] := unit[kron[u,v]];
    unit /: kron[unit[u_], v_ /; MatrixQ[v]] := kron[u,v];
    unit /: kron[u_ /; MatrixQ[u], unit[v_]] := kron[u,v];

    unit /: Dot[unit[u_], ket[k_]] := ket[u . k];
    unit /: Dot[bra[k_], unit[u_]] := bra[dag[dag[k] . u]];
    unit /: Dot[unit[u_], unit[v_]] := unit[u . v];
    unit /: Power[u_unit, n_Integer] := dotexp[u,n];
    unit /: Times[Complex[0,1], unit[u_]]  := unit[(I) u];
    unit /: Times[-1, unit[u_]]            := unit[(-1) u];
    unit /: Times[Complex[0,-1], unit[u_]] := unit[(-I) u];
    unit /: Times[a_, unit[u_]] := sqMatr[a u];
    unit /: Plus[unit[u1_], unit[u2_]] := sqMatr[u1 + u2];
    unit /: Plus[unit[u_], sqMatr[m_]] := sqMatr[u + m];
    unit /: Plus[sqMatr[m_], unit[u_]] := sqMatr[m + u];
    unit[sqMatr[u_]] := unit[u];
    unit[unit[u_]] := unit[u];

    unit /: Dot[unit[u_], state[r_]] := state[u . r . dag[u]];
    unit /: Dot[u_unit, e_ens] := Map[({(#)[[1]], u.((#)[[2]])})&, e, {2}];

    unit /: Dimensions[unit[u_]] := Dimensions[u];

    id[n_Integer] := unit[IdentityMatrix[2^n]];
    id[] = id[1];
    wh = unit[(1/Sqrt[2]) {{1,1},{1,-1}}];

    (* Produit en blocs de matrices unitaires. *)
    block[unit[u_], unit[v_]] :=
      Module[
        {d, dd, i, j, t},
        d = Dimensions[u][[1]];
        dd = d + Dimensions[v][[1]];
        t = Table[
          If[i <= d, If[j <= d, u[[i,j]], 0], If[j <= d, 0,  v[[i-d,j-d]]]],
          {i, 1, dd},
          {j, 1, dd}
        ];
        unit[t]
      ];

    ctrl[m_Integer, unit[u_]] :=
      Module[
        {d, dd, i, j, t},
        d = Dimensions[u][[1]];
        dd = 2^m d; 
        t = Table[
          If[(i <= dd - d) || (j <= dd - d), If[i == j, 1, 0], u[[i+d-dd,j+d-dd]]],
          {i, 1, dd},
          {j, 1, dd}
        ];
        unit[t]
      ];
    ctrl[u_unit] := ctrl[1, u];

    not = unit[{{0,1}, {1,0}}];
    cnot = unit[{{1,0,0,0}, {0,1,0,0}, {0,0,0,1}, {0,0,1,0}}];
    knot = unit[{{1,0,0,0}, {0,0,0,1}, {0,0,1,0}, {0,1,0,0}}];
    sigx = unit[{{0,1},{1,0}}];
    sigy = unit[{{0,-I},{I,0}}];
    sigz = unit[{{1,0},{0,-1}}];

    (* cf. BBC+95 *)
    phase[t_] := unit[{{Exp[I t],0}, {0,Exp[I t]}}];
    (* pour "voir" les rotations : règle de la main gauche *)
    rotx[t_]  := unit[{{Cos[t/2], I Sin[t/2]}, {I Sin[t/2], Cos[t/2]}}];
    roty[t_]  := unit[{{Cos[t/2], Sin[t/2]}, {-Sin[t/2], Cos[t/2]}}];
    rotz[t_]  := unit[{{Exp[I t/2], 0}, {0, Exp[-I t/2]}}];
    randomUnit[] := phase[2 Pi Random[]] . rotz[2 Pi Random[]] . roty[2 Pi Random[]] . rotz[2 Pi Random[]];

    gate[n_Integer, f_Function] := gate[n, 2, f];
    gate[n_Integer, base_Integer, f_Function] :=
      Module[
        {d, t, k, i, j},
        d = base^n;
        t = Table[0, {d}, {d}];
        For[k = 0, k < d, k++,
          t += (f @@ IntegerDigits[k, base, n])[[1]] . Table[If[k == j, 1, 0], {i,1,1}, {j,0,d-1}];
        ];
        unit[t]
      ];

    fgate[m_Integer, n_Integer, f_Function] :=
      gate[
	m+n,
	(
	  Module[{xx,bb},
	    xx={##}[[Range[1,m]]];
	    bb={##}[[Range[m+1,m+n]]];
	    ket[Join[xx, bxor[bb, f[xx]]]]
	  ]
	)&
      ];

    swap[n_, i_, j_] :=
      Module[
        {args, k},
        gate[
          n, 
          (
            args = {##};
            ket[
              Flatten[
                {
                  Table[args[[k]], {k,1,i-1}], 
                  args[[j]], 
                  Table[args[[k]], {k,i+1,j-1}], 
                  args[[i]], 
                  Table[args[[k]], {k,j+1,n}]
                }
              ]
            ]
          )&		
        ]
      ];
    swap[] = swap[2,1,2];


    cycle[n_, i_, j_] :=
      Module[
        {args, k},
        gate[
          n, 
          (
            args = {##};
            ket[
              Flatten[
                {
                  Table[args[[k]], {k,1,i-1}], 
                  args[[j]], 
                  Table[args[[k]], {k,i,j-1}], 
                  Table[args[[k]], {k,j+1,n}]
                }
              ]
            ]
          )&		
        ]
      ];


    ktrl[u_unit] := Module[{l}, l = Log[2,Dimensions[u][[1]]] + 1; dag[cycle[l,1,l]] . ctrl[u] . cycle[l,1,l]];


    circuit[cc_List] := 

      (*  cc est une matrice décrivant le circuit tel qu'on le dessine:
       *  les lignes de la matrice représentent les fils du circuit et
       *  les portes sont décrites sur une même colonne de la matrice.
       *  Le résultat est une transformation unitaire, ou un
       *  super-opérateur, ou un SQO équivalent au circuit donné.
       *)

      Module[
        {c, n, p, i, j, t, tj, tjj, k, cyc, ii, oldciplusk},
	c = cc;
	n = Dimensions[c][[1]];
	p = Dimensions[c][[2]];
	t = id[n];
	For[j=1, j<=p, j++,
          (*  "pour chaque porte, de gauche à droite" *)

	  tj = id[0];
	  tjj = id[n];
	  For[i=1, i<=n, i++,
            (* "pour chaque fil, de haut en bas" *)

            (*
            Print["trace i=", i];
            Print["trace t=", MatrixForm[t]];
            Print["trace tj=", MatrixForm[tj]];
            Print["trace c=", c];
            Print["trace tjj=", MatrixForm[tjj]];
            *)

	    If[c[[i,j]] === \[Placeholder],
	      tj = kron[tj, id[1]]
	    (* else *) ,
	      If[Head[c[[i,j]]] === Integer,
                k = c[[i,j]];
	        cyc = cycle[n, i+1, i+k];
		t = cyc . t;
		tjj = tjj . dag[cyc];
		c[[i,j]] = 1;
		oldciplusk = c[[i+k,j]];
		For[ii = i+k, ii > i+1, ii--,
		  c[[ii,j]] = c[[ii-1,j]];
                  If[(Head[c[[ii,j]]]===Integer) && (c[[ii,j]] > i+k+1-ii), 
                    c[[ii,j]]--
                  ]
		];
		c[[i+1,j]] = oldciplusk;
		If[Head[c[[i+1,j]]] === Integer, 
                  c[[i+1,j]] += (k-1)
                ]
              (* else *) ,
	        tj = kron[tj, c[[i,j]]]
	      ]
	    ]
	  ];
	  t =  tjj . tj . t
	];
	t
      ];

    fourier[m_] :=
      unit[
        (1/Sqrt[m]) Module[{j,k,w}, 
          w = Exp[(2 Pi I)/m]; Sum[w^(j k) ket[Subscript[j,m]] . dag[ket[Subscript[k,m]]], {j,0,m-1}, {k,0,m-1}]
        ]
      ];


    (*--------------------------------------------------------------------
      supop[_] : super-opérateur, contient une liste de k matrices
      (A_j, 1 <= j <= k) carrées d'ordre r, telles que :
          Somme(j = 1..k, dag(A_j) A_j) = I_rxr.

      state[_] : contient une matrice de densité.
      --------------------------------------------------------------------*)

    supop /: Format[supop[s_]] := MatrixForm /@ s;

    supop /: kron[supop[s1_], supop[s2_]] := supop[Flatten[Table[kron[s1[[k1]], s2[[k2]]], {k1,1,Length[s1]}, {k2,1,Length[s2]}], 1]];
    supop /: kron[unit[u_], supop[s_]] := supop[Table[kron[u, s[[k]]], {k,1,Length[s]}]];
    supop /: kron[supop[s_], unit[u_]] := supop[Table[kron[s[[k]], u], {k,1,Length[s]}]];

    supop /:
      Dot[supop[s1_], supop[s2_]] :=
        supop[DeleteCases[Flatten[Table[Dot[s1[[k1]], s2[[k2]]], {k1,1,Length[s1]}, {k2,1,Length[s2]}], 1], x_ /; nullListQ[x], 1]];
    supop /:
      Dot[supop[s_], unit[u_]] := 
        supop[DeleteCases[Table[Dot[s[[k]], u], {k,1,Length[s]}], x_ /; nullListQ[x], 1]];
    supop /:
      Dot[unit[u_], supop[s_]] := 
        supop[DeleteCases[Table[Dot[u, s[[k]]], {k,1,Length[s]}], x_ /; nullListQ[x], 1]];

    supopQ[supop[a_]] := Module[{k, r},
        k = Dimensions[a][[1]];
        r = Dimensions[a][[2]];
        (Dimensions[a][[2]] == Dimensions[a][[3]]) && (Sum[dag[a[[j]]] . a[[j]], {j,1,k}] == IdentityMatrix[r])
    ];

    supop[u_unit] := u;
    supop[s_supop] := s;
    supop[{u_ /; MatrixQ[u]}] := unit[u];
    supop[{x___List, a_List , y___List, a_List, z___List}] := supop[{x, Sqrt[2] a, y, z}];

    state /: Format[state[r_ /; MatrixQ[r]]] := MatrixForm[r];
    state[ket[k_]] := state[k . dag[k]];
    kron[state[r1_], state[r2_]] := state[kron[r1, r2]];
    kron[r_state, k_ket] := kron[r, state[k]];
    kron[k_ket, r_state] := kron[state[k], r];

    state /: Times[a_, state[r_]] := sqMatr[a r];
    state /: Plus[state[r1_], state[r2_]] := sqMatr[r1 + r2];
    state /: Plus[state[r1_], sqMatr[r2_]] := sqMatr[r1 + r2];
    state /: Plus[sqMatr[r1_], state[r2_]] := sqMatr[r1 + r2];
    sqMatr /: Plus[sqMatr[r1_], sqMatr[r2_]] := sqMatr[r1 + r2];
    sqMatr /: Times[a_, sqMatr[r_]] := sqMatr[a r];
    state /: Dot[state[r1_], state[r2_]] := sqMatr[r1 . r2];
    state /: Dot[state[r1_], sqMatr[r2_]] := sqMatr[r1 . r2];
    state /: Dot[sqMatr[r1_], state[r2_]] := sqMatr[r1 . r2];
    sqMatr /: Dot[sqMatr[r1_], sqMatr[r2_]] := sqMatr[r1 . r2];
    state[sqMatr[r_]]   := state[r];
    state[state[r_]]    := state[r];

    posMatr /: Times[a_, posMatr[r_]] := sqMatr[a r];
    posMatr /: Plus[posMatr[r1_], posMatr[r2_]] := posMatr[r1 + r2];
    posMatr /: Plus[posMatr[r1_], sqMatr[r2_]] := sqMatr[r1 + r2];
    posMatr /: Plus[sqMatr[r1_], posMatr[r2_]] := sqMatr[r1 + r2];
    posMatr /: Dot[posMatr[r1_], posMatr[r2_]] := sqMatr[r1 . r2];
    posMatr /: Dot[posMatr[r1_], sqMatr[r2_]] := sqMatr[r1 . r2];
    posMatr /: Dot[sqMatr[r1_], posMatr[r2_]] := sqMatr[r1 . r2];
    posMatr[sqMatr[q_]] := posMatr[q];

    state /: Plus[state[r1_], posMatr[r2_]] := posMatr[r1 + r2];
    state /: Plus[posMatr[r1_], state[r2_]] := posMatr[r1 + r2];
    state /: Dot[state[r1_], posMatr[r2_]] := sqMatr[r1 . r2];
    state /: Dot[posMatr[r1_], state[r2_]] := sqMatr[r1 . r2];
    state[posMatr[r_]]  := state[r];

    stateQ[state[r_]] := (r == dag[r]) && (Tr[r] == 1) && (And @@ ((# >= 0)& /@ Eigenvalues[r]));

    state   /: Tr[r_state]     := 1;
    posMatr /: Tr[posMatr[q_]] := Tr[q];
    sqMatr  /: Tr[sqMatr[m_]]  := Tr[m];

    eigenVal[posMatr[q_]]  := posMatr[DiagonalMatrix[Reverse[Eigenvalues[q]]]];
    eigenVect[posMatr[q_]] :=
      unit[Transpose[GramSchmidt[Reverse[Eigensystem[q][[2]]], InnerProduct -> ((Conjugate[#1] . #2)&)]]];
    squRoot[q_posMatr] := Module[{u, d},
        u = N[eigenVect[q][[1]]];
        d = N[eigenVal[q][[1]]];
        posMatr[dag[u] . Sqrt[d] . u]
    ];
    eigenVal[state[r_]]      := state[eigenVal[posMatr[r]]];
    eigenVect[state[r_]]     := eigenVect[posMatr[r]];
    state /: squRoot[state[r_]] := squRoot[posMatr[r]];
    eigenList[r_state] := Transpose[{Tr[eigenVal[r][[1]], List], ket /@ Map[List, Transpose[eigenVect[r][[1]]], {2}]}];
    eigenList[e_ens]   := eigenList[state[e]];
    eigenList[k_ket]   := eigenList[state[k]];
    eigen[r_state] := ens[DeleteCases[eigenList[r], {0,_}]];
    eigen[e_ens]   := eigen[state[e]];
    eigen[k_ket]   := eigen[state[k]];


    entropy[state[r_]]   := Re[N[-(Plus @@ ((If[#==0,0,# Log[2,#]])& /@ Eigenvalues[r]))]];
    state /: Dot[bra[k1_], state[r_], ket[k2_]] := Re[N[(dag[k1] . r . k2)[[1,1]]]];
    fidelity[r1_state, r2_state] := Re[N[(Tr[squRoot[posMatr[squRoot[r1] . r2 . squRoot[r1]]]])^2]];


    maxmix[n_Integer] := state[(1/n) IdentityMatrix[n]];

    supop /: Dot[supop[a_], state[r_]] := state[Sum[a[[j]] . r . dag[a[[j]]], {j,Length[a]}]];
    supop /: Dot[supop[a_], ket[k_]] := Module[{r}, r = k . dag[k]; state[Sum[a[[j]] . r . dag[a[[j]]], {j,Length[a]}]]];
    supop /: Dot[s_supop, e_ens] := Map[({(#)[[1]], s.((#)[[2]])})&, e, {2}];


    (* Mesure complète (dont on ignore le résultat) sur 1 fil, dans la base standard. *)
    mm = supop[Table[If[(i == k) && (j == k), 1, 0], {k,1,2}, {i,1,2}, {j,1,2}]];



    (*---------------------------------------------------------------------------------
      sqo[_] :
        "selective quantum operation", voir la fameuse feuille de Richard Cleve.
        Implantation : liste de m (i=1..m) "sqol".  Chaque "sqol" est une liste de
        k_i (j=1..k_i) matrices A_ij rectangulaires d'ordre (s,r).  On a que
        Somme(dag(A_ij) . A_ij) = I_rxr.  Remarque : les "sqol" internes se comportent
        comme des "supop", mais n'ont pas forcément la condition d'unitarité, et leurs
        matrices composantes ne sont pas forcément carrées.

      ens[_] :
        "ensemble"
        Implantation : liste de m paires {p, r} où p est une probabilité et r est
        soit un "state" ou un "ket".  La somme des p doit être égale à 1.
      ---------------------------------------------------------------------------------*)

    sqolVertDim[sqol[a_]] := Dimensions[a][[2]];
    sqolHorDim[sqol[a_]] := Dimensions[a][[3]];
    sqolSum[sqol[a_]] := Module[{k}, k = Dimensions[a][[1]]; Sum[dag[a[[j]]] . a[[j]], {j,1,k}]];
    sqol /: Format[sqol[a_]] := MatrixForm /@ a;
    sqol /:
      Dot[sqol[a1_], sqol[a2_]] :=
        sqol[DeleteCases[Flatten[Table[Dot[a1[[k1]], a2[[k2]]], {k1,1,Length[a1]}, {k2,1,Length[a2]}], 1], x_ /; nullListQ[x], 1]];
    sqol /: kron[sqol[a1_], sqol[a2_]] := sqol[Flatten[Table[kron[a1[[k1]], a2[[k2]]], {k1,1,Length[a1]}, {k2,1,Length[a2]}], 1]];
    sqol[{x___List, a_List , y___List, a_List, z___List}] := sqol[{x, Sqrt[2] a, y, z}];

    sqo /: Format[sqo[o_]] := o;
   
    sqoQ[sqo[o_]] := Module[{m, r}, m = Length[o]; r = sqolHorDim[o[[1]]]; Sum[sqolSum[o[[i]]], {i,1,m}] == IdentityMatrix[r]];

    sqo[u_unit] := u;
    sqo[s_supop] := s;
    sqo[s_sqo] := s;
    sqo[{sqol[s_ /; Dimensions[s][[2]] == Dimensions[s][[3]]]}] := supop[s];

    sqo /:
      Dot[sqo[s1_], sqo[s2_]] :=
        sqo[DeleteCases[Flatten[Table[Dot[s1[[k1]], s2[[k2]]], {k1,1,Length[s1]}, {k2,1,Length[s2]}], 1], sqol[{}], 1]];
    sqo /: Dot[unit[u_], sqo[s_]] :=
        sqo[DeleteCases[Flatten[Table[Dot[sqol[{u}], s[[k]]], {k,1,Length[s]}], 1], sqol[{}], 1]];
    sqo /: Dot[sqo[s_], unit[u_]] :=
        sqo[DeleteCases[Flatten[Table[Dot[s[[k]], sqol[{u}]], {k,1,Length[s]}], 1], sqol[{}], 1]];
    sqo /: Dot[supop[s1_], sqo[s2_]] :=
        sqo[DeleteCases[Flatten[Table[Dot[sqol[s1], s2[[k]]], {k,1,Length[s2]}], 1], sqol[{}], 1]];
    sqo /: Dot[sqo[s1_], supop[s2_]] :=
        sqo[DeleteCases[Flatten[Table[Dot[s1[[k]], sqol[s2]], {k,1,Length[s1]}], 1], sqol[{}], 1]];

    sqo /: kron[sqo[s1_], sqo[s2_]] :=
        sqo[Flatten[Table[kron[s1[[k1]], s2[[k2]]], {k1,1,Length[s1]}, {k2,1,Length[s2]}], 1]];
    sqo /: kron[unit[u_], sqo[s_]] :=
        sqo[Flatten[Table[kron[sqol[{u}], s[[k]]], {k,1,Length[s]}], 1]];
    sqo /: kron[sqo[s_], unit[u_]] :=
        sqo[Flatten[Table[kron[s[[k]], sqol[{u}]], {k,1,Length[s]}], 1]];
    sqo /: kron[supop[s1_], sqo[s2_]] :=
        sqo[Flatten[Table[kron[sqol[s1], s2[[k]]], {k,1,Length[s2]}], 1]];
    sqo /: kron[sqo[s1_], supop[s2_]] :=
        sqo[Flatten[Table[kron[s1[[k]], sqol[s2]], {k,1,Length[s1]}], 1]];


    sqo /: Dot[sqo[s_], state[r_]] :=
      ens[
        Module[
          {m, si, ri, pi},
          m = Length[s];
          DeleteCases[
            Table[
                si = s[[i, 1]];
                ri = Sum[si[[j]] . r . dag[si[[j]]], {j,1,Length[si]}];
                pi = Tr[ri];
                If[pi != 0, {pi, state[(1/pi) ri]}, 0, {pi, state[(1/pi) ri]}],
              {i,1,m}
            ],
            0,
            1
          ]
        ]
      ];
    sqo /: Dot[s_sqo, k_ket] := s . state[k];
    sqo /: Dot[s_sqo, e_ens] := Map[({(#)[[1]], s . ((#)[[2]])})&, e, {2}];

    (* Mesure complète (dont on n'ignore pas le résultat) sur 1 fil, dans la base standard. *)
    mmm = sqo[{sqol[{{{1,0},{0,0}}}], sqol[{{{0,0},{0,1}}}]}];

    (* "Trace Out" de l'espace B, parmi les espaces (A,B,C) *)
    trout[nA_, nB_, nC_] := sqo[{sqol[Table[kron[IdentityMatrix[nA], {Table[If[j == jj, 1, 0], {jj,1,nB}]}, IdentityMatrix[nC]], {j,1,nB}]]}];

    (* Ajout de l'ancille |0> (de l'espace B) à un état (des espaces A et C) *)
    anc[nA_, nB_, nC_] := sqo[{sqol[{kron[IdentityMatrix[nA], Table[{If[ii==1,1,0]}, {ii,1,nB}], IdentityMatrix[nC]]}]}];

    ens /: Format[ens[e_]] := e;

    state[ens[e_]] := state[Sum[e[[i,1]] state[e[[i,2]]], {i,1,Length[e]}]];
    ens[{{1,r_state}}] := r;
    ens[r_state] := r;
    ens[{{1,k_ket}}] := k;
    ens[k_ket] := k;
    ens[ens[e_]] := ens[e];
    ens[{x___List, {p_, ens[e_List]}, y___List}] := ens[{x, Sequence @@ Map[({p * (#)[[1]], (#)[[2]]})&, e, {1}], y}];
    ens[{x___List, {p1_, r_}, y___List, {p2_, r_}, z___List}] := ens[{x, {p1 + p2, r}, y, z}];

    ensQ[ens[e_]] :=
        (And @@ ((# >= 0)& /@ Transpose[e][[1]])) && ((Plus @@ Transpose[e][[1]]) == 1) && (And @@ ((stateQ[state[#]])& /@ Transpose[e][[2]]));

    bb84 = ens[{{1/4,zp}, {1/4,zm}, {1/4,xp}, {1/4,xm}}];

    kron[e_ens, k_ket]   := Map[({(#)[[1]], kron[(#)[[2]], k]})&, e, {2}];
    kron[k_ket, e_ens]   := Map[({(#)[[1]], kron[k, (#)[[2]]]})&, e, {2}];
    kron[e_ens, r_state] := Map[({(#)[[1]], kron[(#)[[2]], r]})&, e, {2}];
    kron[r_state, e_ens] := Map[({(#)[[1]], kron[r, (#)[[2]]]})&, e, {2}];


    (*---------------------------------------------------------------------------------
      schmidt[nA, nB, k] :
        Calcule une décomposition de schmidt du ket k dans l'espace Alice-Bob.
        nA : nombre de dimensions de l'espace d'Alice,
        nB : nombre de dimensions de l'espace de Bob.
        Retourne:
          schmidt[{{v_1,a_1,b_1}, {v_2,a_2,b_2}, ...}]
        où
          v_i : les valeurs propres des matrices marginales,
          a_i : les vecteurs de la base de schmidt côté Alice,
          b_i : les vecteurs de la base de schmidt côté Bob.
      ---------------------------------------------------------------------------------*)

    schmidt[nA_, nB_, k_ket] := 
      schmidt[
        Apply[
          (List[Sqrt[#1], #2, ket[krondiv[(1/Sqrt[#1])kron[state[#2][[1]], IdentityMatrix[nB]] . k[[1]], (#2)[[1]]]]])&,
          DeleteCases[eigenList[trout[nA, nB, 1] . k], {0,_}],
          1
        ]
      ];

    schmidt /: Format[schmidt[s_]] := s;
    ket[schmidt[s_]] := ket[Plus @@ Apply[((#1) kron[(#2), (#3)])& , s, {1}]];
    state[s_schmidt] := state[ket[s]];

  End[];
EndPackage[];
