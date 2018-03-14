%{
	/*---
	Flow*: A Verification Tool for Cyber-Physical Systems.
	Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
	Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

	The code is released as is under the GNU General Public License (GPL).
	---*/


	#include "modelParser.h"
	#include "MyLogger.h"

	extern int yyerror(const char *);
	extern int yyerror(string);
	extern int yylex();
	extern int yyparse();
	extern void parseMpfr(string *str, mpfr_t num);
	bool err;
%}

%union
{
	double dblVal;
	string *identifier;
	vector<Interval> *intVec;
	vector<int> *iVec;
	vector<double> *dVec;
	vector<Monomial> *monoVec;
	vector<Polynomial> *polyVec;
	Monomial *mono;
	Polynomial *poly;
	TaylorModelVec *tmVec;
	Matrix *mat;
	vector<vector<double> > *dVecVec;
	vector<PolynomialConstraint> *vecConstraints;
	ResetMap *resetMap;
	Flowpipe *pFlowpipe;
	TaylorModel *ptm;
	vector<TaylorModel> *ptmVec;
	Interval *pint;
	vector<string> *strVec;
	TreeNode *pNode;
}


%token<dblVal> NUM
%token<identifier> IDENT
%token<identifier> MPFRNUM
%token STATEVAR TMVAR TM EQ GEQ LEQ ASSIGN END
%token MODE INIT BELONGSTO
%token POLYODE1 POLYODE2 POLYODE3
%token VISUALIZE PARAAGGREG INTAGGREG TMAGGREG
%token OUTPUT
%token CONTINUOUS HYBRID
%token SETTING
%token FIXEDST FIXEDORD ADAPTIVEST ADAPTIVEORD
%token MIN MAX
%token REMEST
%token INTERVAL OCTAGON GRID
%token QRPRECOND IDPRECOND COMPIDPRECOND SHRINRWRAPPING REM NOPRECOND
%token QRPRECOND1 QRPRECOND2 QRPRECOND3
%token TIME
%token MODES JUMPS INV GUARD RESET START MAXJMPS
%token PRINTON PRINTOFF UNSAFESET
%token CONTINUOUSFLOW HYBRIDFLOW
%token TAYLOR_PICARD TAYLOR_REMAINDER TAYLOR_POLYNOMIAL NONPOLY_CENTER
%token EXP SIN COS LOG SQRT
%token NPODE_TAYLOR CUTOFF PRECISION
%token GNUPLOT MATLAB COMPUTATIONPATHS
%token LINEARODE PAR
%token METHOD
%token ALGORITHM ALG_FLOW ALG_SIMPLE_IMPL ALG_SIMPLE_COMP ALG_SMALL_COMP FLOW_IMPL
%token DECOMPOSITION 
%token NODECOMPOSITION
%token MYMODEL
%token MYMODELS
%token MYMONO
%token MYHORNERFORMS
%token MYINTGVEC
%token MYINTRVEC
%token FLOWPIPES
%token VARS

%type <poly> polynomial
%type <poly> ODEpolynomial
%type <poly> constraint_polynomial
%type <poly> reset_polynomial
%type <poly> interval_polynomial
%type <poly> linear_polynomial
%type <tmVec> ode
%type <tmVec> linear_ode
%type <iVec> orders
%type <pFlowpipe> init
%type <resetMap> reset
%type <dVec> real_valued_vector
%type <dVecVec> real_valued_vectors
%type <dVec> vector_components
%type <vecConstraints> polynomial_constraints
%type <vecConstraints> linear_constraints
%type <tmVec> taylor_model
%type <tmVec> interval_taylor_model
%type <intVec> taylor_model_domain
%type <intVec> intervals
%type <intVec> remainders
// %type <intVec> local_remainders
%type <ptm> non_polynomial_rhs_picard
%type <pint> non_polynomial_rhs_remainder
%type <poly> non_polynomial_rhs_no_remainder
%type <identifier> non_polynomial_rhs_string
%type <identifier> non_polynomial_rhs_center
%type <strVec> npode
%type <pNode> computation_path
%type <iVec> compVarIds
%type <poly> my_poly
%type <ptm> my_taylor_model
%type <ptmVec> my_taylor_models
%type <mono> my_mono
%type <iVec> my_integer_vector
%type <intVec> my_interval_vector
%type <polyVec> my_polys
%type <pint> mpfr_interval




%left GEQ LEQ EQ 
%left '+' '-'
%left '*' '/'
%nonassoc uminus
%right '^'

%start model

%%

model: CONTINUOUS '{' continuous '}'
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}
  
  OutputWriter writer(continuousProblem.outputFileName);
  continuousProblem.settings->writer = &writer;
  
	clock_t begin, end;
	begin = clock();
	mlog1(sbuilder() << "print: " << continuousProblem.bPrint);
	continuousProblem.run();
	end = clock();
	printf("%ld flowpipes computed4.\n", continuousProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	
	continuousProblem.bSafetyChecking = false;


  
  #ifdef no_output
    cout << "not creating flow output" << endl;
  #else
  	printf("Preparing for plotting and dumping...\n");
  	continuousProblem.composition();
  	printf("Done.\n");
  #endif
	
  /* FLOWSTAR OUTPUT REMOVED
	continuousProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, continuousProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	continuousProblem.dump(fpDumping);
	printf("Done.\n");
	fclose(fpDumping);
	*/
  logger.reset();
  
  //tprint("fl_part");
  tprint("fl_int");
  
  double integrTime = double(end - begin) / CLOCKS_PER_SEC;
  mlog1(sbuilder() << "computation time: " << integrTime);
  writer.info.push_back(sbuilder() << "computation time: " << integrTime);
  
  #ifdef no_output
    ;
  #else
    cout << "creating flow output" << endl;
    logger.log(continuousProblem.flowpipesCompo.size());
    writer.fromFlowstar(continuousProblem.flowpipesCompo,
        continuousProblem.domains);
    writer.writeCSV();
    writer.writeInfo();
  #endif
  
  
  if(pSerializer != NULL)
    pSerializer->serialize();
}
|
CONTINUOUS '{' continuous '}' unsafe_continuous
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	continuousProblem.run();
	end = clock();
	printf("%ld flowpipes computed1.\n", continuousProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	continuousProblem.bSafetyChecking = true;

	printf("Preparing for plotting and dumping...\n");
	continuousProblem.composition();
	printf("Done.\n");

	continuousProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, continuousProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	continuousProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = continuousProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
|
HYBRID '{' hybrid '}'
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	hybridProblem.run();
	end = clock();
	printf("%ld flowpipes computed2.\n", hybridProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	hybridProblem.bSafetyChecking = false;

	printf("Preparing for plotting and dumping...\n");
	printf("Done.\n");

	hybridProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, hybridProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}

	printf("Dumping the Taylor model flowpipes...\n");
	hybridProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);
}
|
HYBRID '{' hybrid '}' unsafe_hybrid
{
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	hybridProblem.run();
	end = clock();
	printf("%ld flowpipes computed3.\n", hybridProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	hybridProblem.bSafetyChecking = true;

	printf("Preparing for plotting and dumping...\n");
	printf("Done.\n");

	hybridProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, hybridProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	hybridProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = hybridProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
|
stateVarDecls plotting CUTOFF NUM OUTPUT IDENT unsafe_continuous CONTINUOUSFLOW '{' tmVarDecls continuous_flowpipes '}'
{
	clock_t begin, end;

	if($4 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$4,$4);
	continuousProblem.cutoff_threshold = I;

	strcpy(continuousProblem.outputFileName, $6->c_str());

	continuousProblem.plot_2D();

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = continuousProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
|
stateVarDecls plotting CUTOFF NUM OUTPUT IDENT CONTINUOUSFLOW '{' tmVarDecls continuous_flowpipes '}'
{
	if($4 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$4,$4);
	continuousProblem.cutoff_threshold = I;
	
	strcpy(continuousProblem.outputFileName, $6->c_str());

	continuousProblem.plot_2D();
}
|
stateVarDecls modeDecls COMPUTATIONPATHS '{' computation_paths '}' plotting OUTPUT IDENT unsafe_hybrid HYBRIDFLOW '{' hybrid_flowpipes '}'
{
	clock_t begin, end;
	strcpy(hybridProblem.outputFileName, $9->c_str());
	generateNodeSeq(hybridProblem.traceNodes, hybridProblem.traceTree);

	hybridProblem.plot_2D();

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = hybridProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}

	delete $9;
}
|
stateVarDecls modeDecls COMPUTATIONPATHS '{' computation_paths '}' plotting OUTPUT IDENT HYBRIDFLOW '{' hybrid_flowpipes '}'
{
	strcpy(hybridProblem.outputFileName, $9->c_str());
	generateNodeSeq(hybridProblem.traceNodes, hybridProblem.traceTree);
	hybridProblem.plot_2D();

	delete $9;
}
|
TAYLOR_PICARD '{' non_polynomial_rhs_picard '}'
{
	$3->getExpansion(parseResult.expansion);
	parseResult.remainder = $3->getRemainder();
	delete $3;
}
|
TAYLOR_REMAINDER '{' non_polynomial_rhs_remainder '}'
{
	parseResult.remainder = (*$3);
	delete $3;
}
|
TAYLOR_POLYNOMIAL '{' non_polynomial_rhs_no_remainder '}'
{
	parseResult.expansion = (*$3);
	delete $3;
}
|
NONPOLY_CENTER '{' non_polynomial_rhs_center '}'
{
	parseResult.strExpansion = (*$3);
	delete $3;
}
|
MYMODEL '{' my_taylor_model '}'
{
  mlog1("tm");
  parseResult.model = TaylorModel(*($3));
	delete $3;
}
|
models_wrapper
|
MYMONO '{' my_mono '}'
{
  Monomial m(*$3);
	parseResult.mono = m;
	delete $3;
}
|
MYHORNERFORMS '{' my_polys '}'
{
  parseResult.polys = $3;
}
|
MYINTGVEC '<' my_integer_vector '>'
{
	parseResult.integerVec = *$3;
}
|
MYINTRVEC '<' my_interval_vector '>'
{
	parseResult.intervalVec = *$3;
}
|
VARS '{' parsing_vars  '}' FLOWPIPES '{' multiple_models '}' {
}
;

parsing_vars:
IDENT {
  parseSetting.clear();
	parseSetting.addVar(*$1);
} | 
parsing_vars ',' IDENT {
	parseSetting.addVar(*$3);
};

multiple_models: models_wrapper {
  parseResult.pipes.push_back(parseResult.tmv);
  if(parseResult.name.empty() == false)
    parseResult.names.push_back(parseResult.name);
} | 
multiple_models ',' models_wrapper  {
  parseResult.pipes.push_back(parseResult.tmv);
  if(parseResult.name.empty() == false)
    parseResult.names.push_back(parseResult.name);
};

models_wrapper: MYMODELS '{' my_taylor_models '}' {
  //mlog1("models wrapper");
  parseResult.tmv = TaylorModelVec(*$3);
  //mlog("tmv", parseResult.tmv);
} | IDENT ':' MYMODELS '{' my_taylor_models '}' {
  //mlog1("models wrapper");
  parseResult.name = *$1;
  parseResult.tmv = TaylorModelVec(*$5);
  //mlog("tmv", parseResult.tmv);
};

continuous_flowpipes: continuous_flowpipes '{' interval_taylor_model taylor_model_domain '}'
{
	continuousProblem.flowpipesCompo.push_back(*$3);
	continuousProblem.domains.push_back(*$4);

	delete $3;
	delete $4;
}
|
'{' interval_taylor_model taylor_model_domain '}'
{
	continuousProblem.flowpipesCompo.push_back(*$2);
	continuousProblem.domains.push_back(*$3);

	delete $2;
	delete $3;
}
;

modeDecls: modeDecls IDENT '{' '{' CUTOFF NUM '}' '{' polynomial_constraints '}' '}'
{
	TaylorModelVec tmvDummy;

	if($6 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$6,$6);
	mode_local_setting.cutoff_threshold = I;
	
	hybridProblem.declareMode(*$2, tmvDummy, *$9, 0, mode_local_setting);

	delete $2;
	delete $9;
}
|
IDENT '{' '{' CUTOFF NUM '}' '{' polynomial_constraints '}' '}'
{
	TaylorModelVec tmvDummy;

	if($5 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$5,$5);
	mode_local_setting.cutoff_threshold = I;

	hybridProblem.declareMode(*$1, tmvDummy, *$8, 0, mode_local_setting);

	delete $1;
	delete $8;
}
;

hybrid_flowpipes: hybrid_flowpipes IDENT '{' tmVarDecls continuous_flowpipes '}'
{
	int id = hybridProblem.getIDForMode(*$2);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.modeIDs.push_back(id);
	hybridProblem.flowpipesCompo.push_back(continuousProblem.flowpipesCompo);
	hybridProblem.domains.push_back(continuousProblem.domains);

	continuousProblem.flowpipesCompo.clear();
	continuousProblem.domains.clear();
	continuousProblem.tmVarTab.clear();
	continuousProblem.tmVarNames.clear();
}
|
IDENT '{' tmVarDecls continuous_flowpipes '}'
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.modeIDs.push_back(id);
	hybridProblem.flowpipesCompo.push_back(continuousProblem.flowpipesCompo);
	hybridProblem.domains.push_back(continuousProblem.domains);

	continuousProblem.flowpipesCompo.clear();
	continuousProblem.domains.clear();
	continuousProblem.tmVarTab.clear();
	continuousProblem.tmVarNames.clear();
}
;

computation_paths: computation_paths computation_path ';'
{
}
|
computation_path ';'
{
}
;

computation_path: computation_path '(' NUM ',' '[' NUM ',' NUM ']' ')' '-' '>' IDENT
{
	int id = hybridProblem.getIDForMode(*$13);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$13).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	list<TreeNode *>::iterator iter = $$->children.begin();
	bool found = false;
	for(; iter!=$$->children.end(); ++iter)
	{
		if((*iter)->jumpID == $3 && (*iter)->modeID == id)
		{
			$$ = *iter;
			found = true;
			break;
		}
	}

	if(!found)
	{
		Interval I($6, $8);
		TreeNode *tmp = new TreeNode((int)$3, id, I);
		tmp->parent = $$;
		$$->children.push_back(tmp);
		$$ = tmp;
	}

	delete $13;
}
|
IDENT
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(hybridProblem.traceTree == NULL)
	{
		Interval intZero;
		hybridProblem.traceTree = new TreeNode(0, id, intZero);
		$$ = hybridProblem.traceTree;
	}
	else
	{
		if(hybridProblem.traceTree->modeID == id)
		{
			$$ = hybridProblem.traceTree;
		}
		else
		{
			parseError("Invalid computation path.", lineNum);
			exit(1);
		}
	}

	delete $1;
}
;

print: PRINTON
{
	continuousProblem.bPrint = true;
	hybridProblem.bPrint = true;
}
|
PRINTOFF
{
	continuousProblem.bPrint = false;
	hybridProblem.bPrint = false;
}
;

unsafe_continuous: UNSAFESET '{' polynomial_constraints '}'
{
	continuousProblem.unsafeSet = *$3;
	delete $3;
}
;

unsafe_hybrid: UNSAFESET '{' hybrid_constraints '}'
{
}
;

hybrid_constraints: hybrid_constraints IDENT '{' polynomial_constraints '}'
{
	int id = hybridProblem.getIDForMode(*$2);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.unsafeSet[id] = *$4;
	hybridProblem.bVecUnderCheck[id] = true;
	delete $4;
}
|
{
	vector<PolynomialConstraint> vecEmpty;
	for(int i=0; i<hybridProblem.modeNames.size(); ++i)
	{
		hybridProblem.unsafeSet.push_back(vecEmpty);
		hybridProblem.bVecUnderCheck.push_back(false);
	}
}
;

polynomial_constraints: polynomial_constraints constraint_polynomial LEQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval B($4);
	PolynomialConstraint pc(*$2, B);
	$$->push_back(pc);

	delete $2;
}
|
polynomial_constraints constraint_polynomial GEQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval I(-1);
	$2->mul_assign(I);

	Interval B(-$4);
	PolynomialConstraint pc(*$2, B);
	$$->push_back(pc);

	delete $2;
}
|
polynomial_constraints constraint_polynomial EQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval B($4);
	PolynomialConstraint pc1(*$2, B);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	Interval mB(-$4);
	PolynomialConstraint pc2(*$2, mB);
	$$->push_back(pc2);

	delete $2;
}
|
polynomial_constraints constraint_polynomial BELONGSTO '[' NUM ',' NUM ']'
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*$2, $7);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	PolynomialConstraint pc2(*$2, -$5);
	$$->push_back(pc2);

	delete $2;
}
|
polynomial_constraints constraint_polynomial LEQ IDENT
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*$4);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$4);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	PolynomialConstraint pc(*$2, range);
	$$->push_back(pc);

	delete $2;
	delete $4;
}
|
polynomial_constraints constraint_polynomial GEQ IDENT
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*$4);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$4);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	Interval I(-1);
	$2->mul_assign(I);
	range *= I;
	PolynomialConstraint pc(*$2, range);
	$$->push_back(pc);

	delete $2;
	delete $4;
}
|
polynomial_constraints constraint_polynomial EQ IDENT
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*$4);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$4);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	PolynomialConstraint pc1(*$2, range);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*$2, range);
	$$->push_back(pc2);

	delete $2;
	delete $4;
}
|
polynomial_constraints constraint_polynomial BELONGSTO '[' IDENT ',' IDENT ']'
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*$7);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$7);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$7).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*$2, range);
	$$->push_back(pc1);

	id = continuousProblem.getIDForPar(*$5);

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$5);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval I(-1);
	$2->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*$2, range);
	$$->push_back(pc2);

	delete $2;
	delete $5;
	delete $7;
}
|
{
	$$ = new vector<PolynomialConstraint>(0);
}
;

linear_constraints: linear_constraints linear_polynomial LEQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval B($4);
	PolynomialConstraint pc(*$2, B);
	$$->push_back(pc);

	delete $2;
}
|
linear_constraints linear_polynomial GEQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval I(-1);
	$2->mul_assign(I);

	Interval B(-$4);
	PolynomialConstraint pc(*$2, B);
	$$->push_back(pc);

	delete $2;
}
|
linear_constraints linear_polynomial EQ NUM
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	$$ = $1;
	Interval B($4);
	PolynomialConstraint pc1(*$2, B);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	Interval mB(-$4);
	PolynomialConstraint pc2(*$2, mB);
	$$->push_back(pc2);

	delete $2;
}
|
linear_constraints linear_polynomial BELONGSTO '[' NUM ',' NUM ']'
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*$2, $7);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	PolynomialConstraint pc2(*$2, -$5);
	$$->push_back(pc2);

	delete $2;
}
|
linear_constraints linear_polynomial LEQ IDENT
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}
	
	int id = continuousProblem.getIDForPar(*$4);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$4);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	PolynomialConstraint pc(*$2, range);
	$$->push_back(pc);

	delete $2;
	delete $4;
}
|
linear_constraints linear_polynomial GEQ IDENT
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*$4);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$4);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	Interval I(-1);
	$2->mul_assign(I);
	range *= I;
	PolynomialConstraint pc(*$2, range);
	$$->push_back(pc);

	delete $2;
	delete $4;
}
|
linear_constraints linear_polynomial EQ IDENT
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*$4);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$4);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	PolynomialConstraint pc1(*$2, range);
	$$->push_back(pc1);

	Interval I(-1);
	$2->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*$2, range);
	$$->push_back(pc2);

	delete $2;
	delete $4;
}
|
linear_constraints linear_polynomial BELONGSTO '[' IDENT ',' IDENT ']'
{
	if($2->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}
	
	int id = continuousProblem.getIDForPar(*$7);
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$7);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$7).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*$2, range);
	$$->push_back(pc1);

	id = continuousProblem.getIDForPar(*$5);

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *$5);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval I(-1);
	$2->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*$2, range);
	$$->push_back(pc2);

	delete $2;
	delete $5;
	delete $7;
}
|
{
	$$ = new vector<PolynomialConstraint>(0);
}
;

continuous: stateVarDecls SETTING '{' settings print '}' POLYODE1 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = ONLY_PICARD;

	delete $9;
	delete $13;
}
|
stateVarDecls SETTING '{' settings print '}' POLYODE2 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LOW_DEGREE;

	delete $9;
	delete $13;
}
|
stateVarDecls SETTING '{' settings print '}' POLYODE3 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = HIGH_DEGREE;

	delete $9;
	delete $13;
}
|
stateVarDecls SETTING '{' settings print '}' NPODE_TAYLOR '{' npode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = NONPOLY_TAYLOR;

	delete $9;
	delete $13;
}
|
stateVarDecls SETTING '{' settings print '}' LINEARODE '{' linear_ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$9, *$13);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LINEAR;

	delete $9;
	delete $13;
}
|
stateVarDecls parDecls SETTING '{' settings print '}' POLYODE1 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$10, *$14);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = ONLY_PICARD;

	delete $10;
	delete $14;
}
|
stateVarDecls parDecls SETTING '{' settings print '}' POLYODE2 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$10, *$14);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LOW_DEGREE;

	delete $10;
	delete $14;
}
|
stateVarDecls parDecls SETTING '{' settings print '}' POLYODE3 '{' ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$10, *$14);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = HIGH_DEGREE;

	delete $10;
	delete $14;
}
|
stateVarDecls parDecls SETTING '{' settings print '}' NPODE_TAYLOR '{' npode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$10, *$14);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = NONPOLY_TAYLOR;

	delete $10;
	delete $14;
}
|
stateVarDecls parDecls SETTING '{' settings print '}' LINEARODE '{' linear_ode '}' INIT '{' init '}'
{
	ContinuousSystem system(*$10, *$14);
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LINEAR;

	delete $10;
	delete $14;
}
;

hybrid: stateVarDecls SETTING '{' settings MAXJMPS NUM print '}' MODES '{' modes '}' JUMPS '{' jumps '}' INIT '{' hybrid_init '}'
{
	if($6 < 0)
	{
		parseError("The maximum jump depth should be a nonnegative integer.", lineNum);
		exit(1);
	}

	hybridProblem.maxJumps = (int)$6;
}
|
stateVarDecls parDecls SETTING '{' settings MAXJMPS NUM print '}' MODES '{' modes '}' JUMPS '{' jumps '}' INIT '{' hybrid_init '}'
{
	if($7 < 0)
	{
		parseError("The maximum jump depth should be a nonnegative integer.", lineNum);
		exit(1);
	}

	hybridProblem.maxJumps = (int)$7;
}
;

hybrid_init: IDENT '{' intervals '}'
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	Flowpipe initialSet(*$3, intZero);
	hybridProblem.initialConfig(id, initialSet);

	int numVars = hybridProblem.stateVarNames.size();

	string tVar("local_t");
	hybridProblem.declareTMVar(tVar);
	continuousProblem.declareTMVar(tVar);

	char name[NAME_SIZE];

	for(int i=0; i<numVars; ++i)
	{
		sprintf(name, "%s%d", local_var_name, i+1);
		string tmVarName(name);
		hybridProblem.declareTMVar(tmVarName);
		continuousProblem.declareTMVar(tmVarName);
	}

	delete $1;
	delete $3;
}
|
IDENT '{' tmVarDecls taylor_model taylor_model_domain '}'
{
	int id = hybridProblem.getIDForMode(*$1);
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Flowpipe initialSet(*$4, *$5);
	hybridProblem.initialConfig(id, initialSet);

	delete $4;
	delete $5;
}
;

modes: modes IDENT local_setting '{' POLYODE1 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$7, *$11, ONLY_PICARD, mode_local_setting);
	mode_local_setting.clear();

	delete $2;
	delete $7;
	delete $11;
}
|
IDENT local_setting '{' POLYODE1 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$6, *$10, ONLY_PICARD, mode_local_setting);
	mode_local_setting.clear();

	delete $1;
	delete $6;
	delete $10;
}
|
modes IDENT local_setting '{' POLYODE2 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$7, *$11, LOW_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete $2;
	delete $7;
	delete $11;
}
|
IDENT local_setting '{' POLYODE2 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$6, *$10, LOW_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete $1;
	delete $6;
	delete $10;
}
|
modes IDENT local_setting '{' POLYODE3 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$7, *$11, HIGH_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete $2;
	delete $7;
	delete $11;
}
|
IDENT local_setting '{' POLYODE3 '{' ode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$6, *$10, HIGH_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete $1;
	delete $6;
	delete $10;
}
|
modes IDENT local_setting '{' NPODE_TAYLOR '{' npode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$7, *$11, NONPOLY_TAYLOR, mode_local_setting);
	mode_local_setting.clear();

	delete $2;
	delete $7;
	delete $11;
}
|
IDENT local_setting '{' NPODE_TAYLOR '{' npode '}' INV '{' polynomial_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$6, *$10, NONPOLY_TAYLOR, mode_local_setting);
	mode_local_setting.clear();

	delete $1;
	delete $6;
	delete $10;
}
|
modes IDENT local_setting '{' LINEARODE '{' linear_ode '}' INV '{' linear_constraints '}' '}'
{
	hybridProblem.declareMode(*$2, *$7, *$11, LINEAR, mode_local_setting);
	mode_local_setting.clear();

	delete $2;
	delete $7;
	delete $11;
}
|
IDENT local_setting '{' LINEARODE '{' linear_ode '}' INV '{' linear_constraints '}' '}'
{
	hybridProblem.declareMode(*$1, *$6, *$10, LINEAR, mode_local_setting);
	mode_local_setting.clear();

	delete $1;
	delete $6;
	delete $10;
}
;



local_setting: SETTING '{' parameters  '}'
{
}
|
{
}
;

parameters: parameters FIXEDST NUM
{
	if($3 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveSteps = false;
	mode_local_setting.step = $3;
}
|
parameters REMEST NUM
{
	if($3 <= 0)
	{
		parseError("Remainder estimation should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$3, $3);

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		mode_local_setting.estimation.push_back(I);
	}
}
|
parameters REMEST '{' remainders '}'
{
	mode_local_setting.estimation = *$4;
}
|
parameters QRPRECOND
{
	mlog1("QR1.\n");
	mode_local_setting.precondition = QR_PRE;
}
|
parameters IDPRECOND
{
	mode_local_setting.precondition = ID_PRE;
}
|
parameters FIXEDORD NUM
{
	int order = (int)$3;

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveOrders = false;
	mode_local_setting.orderType = UNIFORM;
	mode_local_setting.orders.push_back(order);
	mode_local_setting.globalMaxOrder = order;
}
|
parameters CUTOFF NUM
{
	if($3 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$3, $3);
	mode_local_setting.cutoff_threshold = I;
}
|
parameters ADAPTIVEST '{' MIN NUM ',' MAX NUM '}'
{
	if($5 <= 0 || $8 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	if($5 > $8)
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveSteps = true;
	mode_local_setting.step = $8;
	mode_local_setting.miniStep = $5;
	mode_local_setting.bAdaptiveOrders = false;
}
|
parameters ADAPTIVEORD '{' MIN NUM ',' MAX NUM '}'
{
	int minOrder = (int)$5;
	int maxOrder = (int)$8;

	if(minOrder <= 0 || maxOrder <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	if(minOrder > maxOrder)
	{
		parseError("MAX order should be no smaller than MIN order.", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveSteps = false;
	mode_local_setting.bAdaptiveOrders = true;
	mode_local_setting.orderType = UNIFORM;
	mode_local_setting.orders.push_back(minOrder);
	mode_local_setting.maxOrders.push_back(maxOrder);
	mode_local_setting.globalMaxOrder = maxOrder;
}
|
parameters FIXEDORD '{' orders '}'
{
	mode_local_setting.bAdaptiveOrders = false;
	mode_local_setting.orderType = MULTI;
	mode_local_setting.orders = *$4;

	for(int i=0; i<$4->size(); ++i)
	{
		if((*$4)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*$4)[0];
	for(int i=1; i<$4->size(); ++i)
	{
		if(maxOrder < (*$4)[i])
		{
			maxOrder = (*$4)[i];
		}
	}

	mode_local_setting.globalMaxOrder = maxOrder;
}
|
parameters ADAPTIVEORD '{' MIN '{' orders '}' ',' MAX '{' orders '}' '}'
{
	mode_local_setting.bAdaptiveSteps = false;
	mode_local_setting.bAdaptiveOrders = true;
	mode_local_setting.orderType = MULTI;
	mode_local_setting.orders = *$6;
	mode_local_setting.maxOrders = *$11;

	if($6->size() != $11->size())
	{
		parseError("Orders are not properly specified.", lineNum);
		exit(1);
	}

	for(int i=0; i<$11->size(); ++i)
	{
		if((*$6)[i] <= 0 || (*$11)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}

		if((*$6)[i] > (*$11)[i])
		{
			parseError("MAX order should be no smaller than MIN order.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*$11)[0];
	for(int i=1; i<$11->size(); ++i)
	{
		if(maxOrder < (*$11)[i])
		{
			maxOrder = (*$11)[i];
		}
	}

	mode_local_setting.globalMaxOrder = maxOrder;
}
|
{
}
;

jumps: jumps IDENT '-' '>' IDENT GUARD '{' polynomial_constraints '}' RESET '{' reset '}' PARAAGGREG '{' real_valued_vectors '}'
{
	int startID = hybridProblem.getIDForMode(*$2);
	if(startID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int endID = hybridProblem.getIDForMode(*$5);
	if(endID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($16->size() > 0)
	{
		hybridProblem.declareTrans(startID, endID, *$8, *$12, PARA_AGGREG, *$16);
	}
	else
	{
		vector<vector<double> > emptyVec;
		hybridProblem.declareTrans(startID, endID, *$8, *$12, PARA_AGGREG, emptyVec);
	}
}
|
jumps IDENT '-' '>' IDENT GUARD '{' polynomial_constraints '}' RESET '{' reset '}' INTAGGREG
{
	int startID = hybridProblem.getIDForMode(*$2);
	if(startID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int endID = hybridProblem.getIDForMode(*$5);
	if(endID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	vector<vector<double> > empty;
	hybridProblem.declareTrans(startID, endID, *$8, *$12, INTERVAL_AGGREG, empty);
}
|
{
	hybridProblem.declareTrans();
}
;

reset: reset IDENT '\'' ASSIGN reset_polynomial '+' '[' NUM ',' NUM ']'
{
	$$ = $1;

	int id = hybridProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($8 > $10)
	{
		parseError("Invalid remainder interval.", lineNum);
		exit(1);
	}

	Interval I($8, $10);
	TaylorModel tmTemp(*$5, I);
	$$->tmvReset.tms[id] = tmTemp;

	delete $5;
}
|
reset IDENT '\'' ASSIGN reset_polynomial
{
	$$ = $1;

	int id = hybridProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*$5, intZero);
	$$->tmvReset.tms[id] = tmTemp;

	delete $5;
}
|
{
	int numVars = hybridProblem.stateVarNames.size();

	Matrix coefficients_identity_reset(numVars, numVars+1);

	for(int i=0; i<numVars; ++i)
	{
		coefficients_identity_reset.set(1, i, i+1);
	}

	TaylorModelVec tmvReset(coefficients_identity_reset);

	$$ = new ResetMap(tmvReset);
}
;

real_valued_vectors: real_valued_vectors real_valued_vector
{
	$$->push_back(*$2);
	delete $2;
}
|
{
	$$ = new vector<vector<double> >(0);
}
;

real_valued_vector: '[' vector_components ']'
{
	int rangeDim = $2->size();

	if(rangeDim != hybridProblem.stateVarNames.size())
	{
		parseError("The vector dimension should be equivalent to the system dimension.", lineNum);
		exit(1);
	}

	$$ = new vector<double>(0);

	for(int i=0; i<rangeDim; ++i)
	{
		$$->push_back(0);
	}

	bool bZero = true;
	for(int i=0; i<rangeDim; ++i)
	{
		if((*$2)[i] < -THRESHOLD_LOW || (*$2)[i] > THRESHOLD_LOW)
		{
			if(bZero)
			{
				bZero = false;
			}
		}

		(*$$)[i] = (*$2)[i];
	}

	if(bZero)
	{
		parseError("A template vector should not be zero.", lineNum);
		exit(1);
	}

	delete $2;
}
;

vector_components: vector_components ',' IDENT ':' NUM
{
	int id = hybridProblem.getIDForStateVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
	(*$$)[id] = $5;
	delete $3;
}
|
IDENT ':' NUM
{
	int num = hybridProblem.stateVarNames.size();
	$$ = new vector<double>(num);

	int id = hybridProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = $3;
	delete $1;
}
;

stateVarDecls: STATEVAR stateIdDeclList
{
}
;

stateIdDeclList: stateIdDeclList ',' IDENT
{
	//mlog1("stateIdDeclList1");
  logger.inc();
	//mlog1(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	if(!continuousProblem.declareStateVar(*$3))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	//mlog1(sbuilder() << "*$3: " << *$3);
  logger.dec();
	hybridProblem.declareStateVar(*$3);
	delete $3;
}
|
IDENT
{
	//mlog1("stateIdDeclList2");
  logger.inc();
	//mlog1(sbuilder() << "*$1: " << *$1);
	//mlog1(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	if(!continuousProblem.declareStateVar(*$1))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	
	//mlog1(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	hybridProblem.declareStateVar(*$1);
	//mlog1(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
  logger.dec();
	delete $1;
}
;



parDecls: PAR '{' parDeclList '}'
{
}
;

parDeclList: parDeclList IDENT EQ NUM
{
	if(continuousProblem.getIDForStateVar(*$2) != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	else
	{
		Interval range($4);

		if(!continuousProblem.declarePar(*$2, range))
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "Parameter %s has already been declared.", (*$2).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		hybridProblem.declarePar(*$2, range);
	}
}
|
IDENT EQ NUM
{
	if(continuousProblem.getIDForStateVar(*$1) != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	else
	{
		Interval range($3);

		if(!continuousProblem.declarePar(*$1, range))
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "Parameter %s has already been declared.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		hybridProblem.declarePar(*$1, range);
	}
}
;



settings: FIXEDST NUM TIME NUM remainder_estimation precondition plotting FIXEDORD NUM CUTOFF NUM PRECISION NUM OUTPUT IDENT algorithm decomposition
{
	mlog1("settings1");
  mlog1(continuousProblem.algorithm);
	if($2 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	int order = (int)$9;

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}
	mlog1(sbuilder() << "step: " << $2 << ", time: " << time << ", order: " << order);
	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(order);
	continuousProblem.globalMaxOrder = order;

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = $2;
	hybridProblem.time = $4;
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = UNIFORM;
	hybridProblem.global_setting.orders.push_back(order);
	hybridProblem.global_setting.globalMaxOrder = order;

	if($11 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$11,$11);
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)$13;

	strcpy(continuousProblem.outputFileName, (*$15).c_str());
	strcpy(hybridProblem.outputFileName, (*$15).c_str());
	mlog1(sbuilder() << "cutoff_thrs: " << $11 << ", intnumprec: " << $13 << ", outFname: " << *$15);
	

	delete $15;
}
|
FIXEDST NUM TIME NUM remainder_estimation precondition plotting ADAPTIVEORD '{' MIN NUM ',' MAX NUM '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	mlog1("settings2");
	if($2 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	int minOrder = (int)$11;
	int maxOrder = (int)$14;

	if(minOrder <= 0 || maxOrder <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	if(minOrder > maxOrder)
	{
		parseError("MAX order should be no smaller than MIN order.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = true;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(minOrder);
	continuousProblem.maxOrders.push_back(maxOrder);
	continuousProblem.globalMaxOrder = maxOrder;

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = $2;
	hybridProblem.time = $4;
	hybridProblem.global_setting.bAdaptiveOrders = true;
	hybridProblem.global_setting.orderType = UNIFORM;
	hybridProblem.global_setting.orders.push_back(minOrder);
	hybridProblem.global_setting.maxOrders.push_back(maxOrder);
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if($17 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$17,$17);
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)$19;

	strcpy(continuousProblem.outputFileName, (*$21).c_str());
	strcpy(hybridProblem.outputFileName, (*$21).c_str());

	delete $21;
}
|
FIXEDST NUM TIME NUM remainder_estimation precondition plotting FIXEDORD '{' orders '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	mlog1("settings3");
	if($2 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *$10;

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = $2;
	hybridProblem.time = $4;
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = MULTI;
	hybridProblem.global_setting.orders = *$10;

	for(int i=0; i<$10->size(); ++i)
	{
		if((*$10)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*$10)[0];
	for(int i=1; i<$10->size(); ++i)
	{
		if(maxOrder < (*$10)[i])
		{
			maxOrder = (*$10)[i];
		}
	}

	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if($13 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$13,$13);
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)$15;

	strcpy(continuousProblem.outputFileName, (*$17).c_str());
	strcpy(hybridProblem.outputFileName, (*$17).c_str());

	delete $10;
	delete $17;
}
|
FIXEDST NUM TIME NUM remainder_estimation precondition plotting ADAPTIVEORD '{' MIN '{' orders '}' ',' MAX '{' orders '}' '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	mlog1("settings4");
	if($2 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = $2;
	continuousProblem.time = $4;
	continuousProblem.bAdaptiveOrders = true;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *$12;
	continuousProblem.maxOrders = *$17;

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = $2;
	hybridProblem.time = $4;
	hybridProblem.global_setting.bAdaptiveOrders = true;
	hybridProblem.global_setting.orderType = MULTI;
	hybridProblem.global_setting.orders = *$12;
	hybridProblem.global_setting.maxOrders = *$17;

	if($12->size() != $17->size())
	{
		parseError("Orders are not properly specified.", lineNum);
		exit(1);
	}

	for(int i=0; i<$17->size(); ++i)
	{
		if((*$12)[i] <= 0 || (*$17)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}

		if((*$12)[i] > (*$17)[i])
		{
			parseError("MAX order should be no smaller than MIN order.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*$17)[0];
	for(int i=1; i<$17->size(); ++i)
	{
		if(maxOrder < (*$17)[i])
		{
			maxOrder = (*$17)[i];
		}
	}

	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if($21 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$21,$21);
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)$23;

	strcpy(continuousProblem.outputFileName, (*$25).c_str());
	strcpy(hybridProblem.outputFileName, (*$25).c_str());

	delete $12;
	delete $17;
	delete $25;
}
|
ADAPTIVEST '{' MIN NUM ',' MAX NUM '}' TIME NUM remainder_estimation precondition plotting FIXEDORD NUM CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	mlog1("settings5");
	if($4 <= 0 || $7 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	if($4 > $7)
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	int order = (int)$15;

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = true;
	continuousProblem.step = $7;
	continuousProblem.miniStep = $4;
	continuousProblem.time = $10;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(order);
	continuousProblem.globalMaxOrder = order;

	hybridProblem.global_setting.bAdaptiveSteps = true;
	hybridProblem.global_setting.step = $7;
	hybridProblem.global_setting.miniStep = $4;
	hybridProblem.time = $10;
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = UNIFORM;
	hybridProblem.global_setting.orders.push_back(order);
	hybridProblem.global_setting.globalMaxOrder = order;

	if($17 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$17,$17);
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)$19;

	strcpy(continuousProblem.outputFileName, (*$21).c_str());
	strcpy(hybridProblem.outputFileName, (*$21).c_str());

	delete $21;
}
|
ADAPTIVEST '{' MIN NUM ',' MAX NUM '}' TIME NUM remainder_estimation precondition plotting FIXEDORD '{' orders '}' CUTOFF NUM PRECISION NUM OUTPUT IDENT
{
	mlog1("settings6");
	if($4 <= 0 || $7 <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	if($4 > $7)
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	for(int i=0; i<$16->size(); ++i)
	{
		if((*$16)[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	continuousProblem.bAdaptiveSteps = true;
	continuousProblem.step = $7;
	continuousProblem.miniStep = $4;
	continuousProblem.time = $10;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *$16;

	hybridProblem.global_setting.bAdaptiveSteps = true;
	hybridProblem.global_setting.step = $7;
	hybridProblem.global_setting.miniStep = $4;
	hybridProblem.time = $10;
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = MULTI;
	hybridProblem.global_setting.orders = *$16;

	int maxOrder = (*$16)[0];
	for(int i=1; i<$16->size(); ++i)
	{
		if(maxOrder < (*$16)[i])
		{
			maxOrder = (*$16)[i];
		}
	}
	
	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if($19 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$19,$19);
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)$21;

	strcpy(continuousProblem.outputFileName, (*$23).c_str());
	strcpy(hybridProblem.outputFileName, (*$23).c_str());

	delete $16;
	delete $23;
}
;

remainder_estimation: REMEST NUM
{
	//mlog1("remainder_estimation1");
	if($2 <= 0)
	{
		parseError("Remainder estimation should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$2, $2);

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		//mlog1(sbuilder() << "i: " << i << " - " << $2);
		continuousProblem.estimation.push_back(I);
		hybridProblem.global_setting.estimation.push_back(I);
	}
}
|
REMEST '{' remainders '}'
{
	mlog1("remainder_estimation2");
	for(int i=0; i<$3->size(); ++i)
	{
		if((*$3)[i].inf() >= (*$3)[i].sup() - THRESHOLD_LOW)
		{
			parseError("Invalid remainder estimation.", lineNum);
			exit(1);
		}
	}

	continuousProblem.estimation = *$3;
	hybridProblem.global_setting.estimation = *$3;
	delete $3;
}
;

remainders: remainders ',' IDENT ':' '[' NUM ',' NUM ']'
{
	$$ = $1;
	int id = continuousProblem.getIDForStateVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($6 >= $8)
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	Interval I($6,$8);
	(*$$)[id] = I;
	delete $3;
}
|
IDENT ':' '[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<Interval>(numVars);

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($4 >= $6)
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	Interval I($4,$6);
	(*$$)[id] = I;
	delete $1;
}
;

orders: orders ',' IDENT ':' NUM
{
	$$ = $1;
	int id = continuousProblem.getIDForStateVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (int)$5;
	delete $3;
}
|
IDENT ':' NUM
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<int>(numVars);
	for(int i=0; i<numVars; ++i)
	{
		(*$$)[i] = 0;
	}

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (int)$3;
	delete $1;
}
;

precondition: QRPRECOND
{
	mlog1(sbuilder() << "QR_PRE: " << QR_PRE);
	continuousProblem.precondition = QR_PRE;
	hybridProblem.global_setting.precondition = QR_PRE;
  //Transformer *transformer = new QRTransformer();
  //continuousProblem.settings->transformer = transformer;
  Transformer *transformer = new QRTransformerPlain();
  continuousProblem.settings->transformer = transformer;
}
|
QRPRECOND1
{
  mforce("making new");
	mlog1("qr precond1");
	mlog1(sbuilder() << "QR_PRE: " << QR_PRE);
	continuousProblem.precondition = QR_PRE;
	hybridProblem.global_setting.precondition = QR_PRE;
  Transformer *transformer = new QRTransformer1();
  continuousProblem.settings->transformer = transformer;
}
|
QRPRECOND2
{
  mforce("making new");
	mlog1("qr precond1");
	mlog1(sbuilder() << "QR_PRE: " << QR_PRE);
	continuousProblem.precondition = QR_PRE;
	hybridProblem.global_setting.precondition = QR_PRE;
  Transformer *transformer = new QRTransformer2();
  continuousProblem.settings->transformer = transformer;
}
|
QRPRECOND3
{
  mforce("making new");
	mlog1("qr precond1");
	mlog1(sbuilder() << "QR_PRE: " << QR_PRE);
	continuousProblem.precondition = QR_PRE;
	hybridProblem.global_setting.precondition = QR_PRE;
  Transformer *transformer = new QRTransformer3();
  continuousProblem.settings->transformer = transformer;
}
|
IDPRECOND
{
	mlog1("id precond");
	continuousProblem.precondition = ID_PRE;
	hybridProblem.global_setting.precondition = ID_PRE;
  Transformer *transformer = new IdentityTransformer();
  continuousProblem.settings->transformer = transformer;
}|
COMPIDPRECOND
{
	mforce("comp id precond");
	continuousProblem.precondition = ID_PRE;
	hybridProblem.global_setting.precondition = ID_PRE;
  Transformer *transformer = new SingleComponentIdentityTransformer();
  continuousProblem.settings->transformer = transformer;
}
|
NOPRECOND
{
	mlog1("no precond");
  Transformer *transformer = new NullTransformer();
  continuousProblem.settings->transformer = transformer;
}
|
SHRINRWRAPPING NUM
{
	mlog1("shrink num");
	continuousProblem.precondition = SHRINK_WRAPPING;
	ShrinkWrappingCondition *cond = new ShrinkWrappingCondition($2);
  Transformer *transformer;
  transformer = new ShrinkWrapper(cond);
  continuousProblem.settings->transformer = transformer;
}
|
SHRINRWRAPPING REM
{
	mlog1("shrink rem");
	ShrinkWrappingCondition *cond = new ShrinkWrappingCondition();
	continuousProblem.precondition = SHRINK_WRAPPING;
  Transformer *transformer;
  transformer = new ShrinkWrapper(cond);
  continuousProblem.settings->transformer = transformer;
}
;

algorithm: ALG_SIMPLE_IMPL
{
	mlog1("ALG_SIMPLE_IMPL");
	if(continuousProblem.precondition == SHRINK_WRAPPING) {
    parseError(
        "Only small component algorith supports shrink wrapping", lineNum);
    exit(1);
	}
	continuousProblem.algorithm = ALGORITHM_SIMPLE_IMPL;
}
|
ALG_SIMPLE_COMP
{
	mlog1("ALG_SIMPLE_COMP");
	if(continuousProblem.precondition == SHRINK_WRAPPING) {
    parseError(
        "Only small component algorith supports shrink wrapping", lineNum);
    exit(1);
	}
	continuousProblem.algorithm = ALGORITHM_SIMPLE_COMP;
}
|
ALG_SMALL_COMP
{
	mlog1("ALG_SMALL_COMP");
	continuousProblem.algorithm = ALGORITHM_SMALL_COMP;
}
|
ALG_SMALL_COMP FLOW_IMPL
{
	mlog1("ALG_SMALL_COMP");
	continuousProblem.algorithm = ALGORITHM_SMALL_COMP;
	continuousProblem.settings->useFlow = true;
	//exit(0);
}
|
ALG_FLOW
{
	mlog1(sbuilder() << "type: " << typeid(continuousProblem).name());
	mlog1("ALG_DEF");
	continuousProblem.algorithm = ALGORITHM_DEFAULT;
}
| {
  mlog1(sbuilder() << "type: " << typeid(continuousProblem).name());
	mlog1("ALG_DEF");
	continuousProblem.algorithm = ALGORITHM_DEFAULT;
}
;

decomposition: DECOMPOSITION '[' components ']'
{
	mlog1("DECOMPOSITION");
}
| NODECOMPOSITION
{
  logger.force("none");
  int varNumber = continuousProblem.stateVarNames.size();
  vector<int> vs;
  for(int i = 0; i < varNumber; i++) {
    vs.push_back(i);
  }
  continuousProblem.components.push_back(vs);
}
| {
  if(continuousProblem.algorithm == ALGORITHM_SMALL_COMP) {
    parseError("No decompostion given for composition algorithm.", lineNum);
    exit(1);
  }
}

components: components ',' component { }
|
component { }

component: '[' compVarIds ']'
{
  for(vector<int>::iterator varIt = $2->begin(); varIt < $2->end(); varIt++) {
    int varId = *varIt;
    
    //check if variable is in any of the other components
    for(vector< vector<int> >::iterator it = continuousProblem.components.begin(); 
        it < continuousProblem.components.end(); it++) {
        
      //check if variable is in the component
      if(it->end() != find(it->begin(), it->end(), varId)) {
        char errMsg[MSG_SIZE];
        string vName;
        continuousProblem.getStateVarName(vName, varId);
        sprintf(errMsg, "Variable %s is already in the component.", vName.c_str());
        parseError(errMsg, lineNum);
        exit(1);
      }
    }
  }

  //add this component to all components
  continuousProblem.components.push_back(*($2));
  delete $2;
}

compVarIds: compVarIds ',' IDENT
{
	//mlog1("compVarIds1");
	//mlog1(sbuilder() << "*$3: " << *$3);
  
	$$ = $1;
  int varId = continuousProblem.getIDForStateVar(*$3);
  if(varId < 0) {
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
  
  //check if variable is in the current component
  if($$->end() != find($$->begin(), $$->end(), varId)) {
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Variable %s is already in the component.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
  $$->push_back(varId);
	delete $3;
}
|
IDENT
{
	//mlog1("compVarIds2");
	//mlog1(sbuilder() << "*$1: " << *$1);
  int varId = continuousProblem.getIDForStateVar(*$1);
  
  if(varId < 0) {
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
  
	$$ = new vector<int>();
  $$->push_back(varId);
	delete $1;
}
;


plotting: GNUPLOT INTERVAL IDENT ',' IDENT
{
  mlog1("plotting1");
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_INTERVAL;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_INTERVAL;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete $3;
	delete $5;
}
|
GNUPLOT OCTAGON IDENT ',' IDENT
{
  mlog1("plotting2");
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State Variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_OCTAGON;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete $3;
	delete $5;
}
|
GNUPLOT GRID NUM IDENT ',' IDENT
{
  mlog1("plotting3");
	int x = continuousProblem.getIDForStateVar(*$4);
	int y = continuousProblem.getIDForStateVar(*$6);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$6).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_GRID;
	continuousProblem.numSections = (int)$3;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_GRID;
	hybridProblem.numSections = (int)$3;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete $4;
	delete $6;
}
|
MATLAB INTERVAL IDENT ',' IDENT
{
  mlog1("plotting4");
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_INTERVAL;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_INTERVAL;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete $3;
	delete $5;
}
|
MATLAB OCTAGON IDENT ',' IDENT
{
  mlog1("plotting5");
	int x = continuousProblem.getIDForStateVar(*$3);
	int y = continuousProblem.getIDForStateVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State Variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_OCTAGON;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete $3;
	delete $5;
}
|
MATLAB GRID NUM IDENT ',' IDENT
{
  mlog1("plotting6");
	int x = continuousProblem.getIDForStateVar(*$4);
	int y = continuousProblem.getIDForStateVar(*$6);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$6).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_GRID;
	continuousProblem.numSections = (int)$3;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_GRID;
	hybridProblem.numSections = (int)$3;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete $4;
	delete $6;
}
;

init: tmVarDecls taylor_model taylor_model_domain
{
	mlog1("INIT1");
	$$ = new Flowpipe(*$2, *$3);

	delete $2;
	delete $3;
}
|
intervals
{
	Interval intZero;
	$$ = new Flowpipe(*$1, intZero);

	delete $1;
}

tmVarDecls: TMVAR tmIdDeclList
{
}
;

tmIdDeclList: tmIdDeclList ',' IDENT
{
	int id = continuousProblem.getIDForStateVar(*$3);

	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	id = continuousProblem.getIDForPar(*$3);
	
	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a parameter.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(!continuousProblem.declareTMVar(*$3))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareTMVar(*$3);
	delete $3;
}
|
IDENT
{
	string tVar("local_t");
	continuousProblem.declareTMVar(tVar);
	hybridProblem.declareTMVar(tVar);

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	id = continuousProblem.getIDForPar(*$1);
	
	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a parameter.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(!continuousProblem.declareTMVar(*$1))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareTMVar(*$1);
	delete $1;
}
;

taylor_model: taylor_model IDENT EQ polynomial '+' '[' NUM ',' NUM ']'
{
	mlog1("taylor_model");
	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($7 > $9)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I($7,$9);
	TaylorModel tmTemp(*$4, I);
	$$ = $1;
	$$->tms[id] = tmTemp;

	delete $2;
	delete $4;
}
|
{
	TaylorModel tmEmpty;
	$$ = new TaylorModelVec;

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		$$->tms.push_back(tmEmpty);
	}
}
;

taylor_model_domain: taylor_model_domain IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForTMVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($5 > $7)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I($5,$7);
	$$ = $1;
	(*$$)[id] = I;

	delete $2;
}
|
IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	$$ = new vector<Interval>( continuousProblem.tmVarNames.size() );

	Interval intZero;
	(*$$)[0] = intZero;

	int id = continuousProblem.getIDForTMVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($4 > $6)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I($4,$6);
	(*$$)[id] = I;

	delete $1;
}
;

intervals: intervals IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForStateVar(*$2);
	
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($5 > $7)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I($5,$7);
	$$ = $1;
	(*$$)[id] = I;

	delete $2;
}
|
{
	mlog1("intervals2");
	//exit(1);
	int numVars = continuousProblem.stateVarNames.size();
	mlog1(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	$$ = new vector<Interval>(numVars);

	string tVar("local_t");
	continuousProblem.declareTMVar(tVar);

	char name[NAME_SIZE];

	mlog1(local_var_name);
	for(int i=0; i<numVars; ++i)
	{
		sprintf(name, "%s%d", local_var_name, i+1);
		string tmVarName(name);
		mlog1(sbuilder() << "name: " << name);
		continuousProblem.declareTMVar(tmVarName);
	}
}
;

ode: ode IDENT '\'' EQ ODEpolynomial
{
	mlog1("ODE1");
	$$ = $1;
	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*$5, intZero);
	$$->tms[id] = tmTemp;

	delete $2;
	delete $5;
}
|
{	
	mlog1("ODE2");
	int numVars = continuousProblem.stateVarNames.size();
	mlog1(sbuilder() << "numVars: " << numVars);
	$$ = new TaylorModelVec;
	TaylorModel tmTemp;
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		$$->tms.push_back(tmTemp);
	}
}
;

npode: npode IDENT '\'' EQ non_polynomial_rhs_string
{
	$$ = $1;

	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = (*$5);

	delete $2;
	delete $5;
}
|
{
	int numVars = continuousProblem.stateVarNames.size();
	$$ = new vector<string>;

	string empty;
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		$$->push_back(empty);
	}
}
;




linear_ode: linear_ode IDENT '\'' EQ linear_polynomial
{
	$$ = $1;

	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*$5, intZero);
	$$->tms[id] = tmTemp;

	delete $2;
	delete $5;
}
|
{
	int numVars = continuousProblem.stateVarNames.size();

	$$ = new TaylorModelVec;
	TaylorModel tmTemp;
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		$$->tms.push_back(tmTemp);
	}
}
;






// only used in Taylor model declaration
polynomial: polynomial '+' polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
polynomial '-' polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' polynomial ')'
{
	$$ = $2; 
}
|
polynomial '*' polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.tmVarNames.size());
	}
	else
	{
		$$ = new Polynomial;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForTMVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.tmVarNames.size();
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id] = 1;
	Monomial monomial(I, degrees);

	$$ = new Polynomial(monomial);
	delete $1;
}
|
NUM
{
	int numVars = continuousProblem.tmVarNames.size();
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
;











ODEpolynomial: ODEpolynomial '+' ODEpolynomial
{
	mlog1("ODEpolynomial+");
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
ODEpolynomial '-' ODEpolynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' ODEpolynomial ')'
{
	$$ = $2; 
}
|
ODEpolynomial '*' ODEpolynomial
{
	mlog1("ODEpolynomial*");
	$$ = $1;
	(*$$) *= (*$3);
	delete $3;
}
|
ODEpolynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		$$ = new Polynomial;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' ODEpolynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	mlog1("ODEpolynomialIDENT");
	int id = continuousProblem.getIDForPar(*$1);
	mlog1(sbuilder() << "idPar: " << id);
	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *$1);

		int numVars = continuousProblem.stateVarNames.size()+1;
		$$ = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*$1);
		mlog1(sbuilder() << "idSV: " << id);

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		$$ = new Polynomial(monomial);
	}

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	mlog1("ODEpolynomial[,]");
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($2, $4);
	$$ = new Polynomial(I, numVars);
}
|
NUM
{
	mlog1("ODEpolynomialNUM");
	mlog1(sbuilder() << "$1: " << $1);
	int numVars = continuousProblem.stateVarNames.size()+1;
	mlog1(sbuilder() << "numVars+1:" << numVars);
	Interval I($1);
	$$ = new Polynomial(I, numVars);
	mlog1($$);
}
;



constraint_polynomial: constraint_polynomial '+' constraint_polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
constraint_polynomial '-' constraint_polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' constraint_polynomial ')'
{
	$$ = $2; 
}
|
constraint_polynomial '*' constraint_polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
constraint_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		$$ = new Polynomial;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' constraint_polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForPar(*$1);

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *$1);

		int numVars = continuousProblem.stateVarNames.size()+1;
		$$ = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*$1);

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		$$ = new Polynomial(monomial);
	}

	delete $1;
}
|
NUM
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
;




reset_polynomial: reset_polynomial '+' reset_polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
reset_polynomial '-' reset_polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' reset_polynomial ')'
{
	$$ = $2; 
}
|
reset_polynomial '*' reset_polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
reset_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		$$ = new Polynomial;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' reset_polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForPar(*$1);

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *$1);

		int numVars = continuousProblem.stateVarNames.size()+1;
		$$ = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*$1);

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		$$ = new Polynomial(monomial);
	}

	delete $1;
}
|
NUM
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
;





interval_taylor_model: interval_taylor_model IDENT EQ interval_polynomial '+' '[' NUM ',' NUM ']'
{
	int id = continuousProblem.getIDForStateVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($7 > $9)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I($7,$9);
	TaylorModel tmTemp(*$4, I);
	$$ = $1;
	$$->tms[id] = tmTemp;

	delete $2;
	delete $4;
}
|
IDENT EQ interval_polynomial '+' '[' NUM ',' NUM ']'
{
	TaylorModel tmEmpty;
	$$ = new TaylorModelVec;

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		$$->tms.push_back(tmEmpty);
	}

	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($6 > $8)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I($6,$8);
	TaylorModel tmTemp(*$3, I);

	$$->tms[id] = tmTemp;

	delete $1;
	delete $3;
}
;










interval_polynomial: interval_polynomial '+' interval_polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
interval_polynomial '-' interval_polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' interval_polynomial ')'
{
	$$ = $2; 
}
|
interval_polynomial '*' interval_polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
interval_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new Polynomial(I, continuousProblem.tmVarNames.size());
	}
	else
	{
		$$ = new Polynomial;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' interval_polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForTMVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.tmVarNames.size();
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id] = 1;
	Monomial monomial(I, degrees);

	$$ = new Polynomial(monomial);
	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.tmVarNames.size();
	Interval I($2, $4);
	$$ = new Polynomial(I, numVars);
}
|
NUM
{
	int numVars = continuousProblem.tmVarNames.size();
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
;




















non_polynomial_rhs_picard: non_polynomial_rhs_picard '+' non_polynomial_rhs_picard
{
	$$ = $1;
	$1->add_assign(*$3);
	delete $3;
}
|
non_polynomial_rhs_picard '-' non_polynomial_rhs_picard
{
	$$ = $1;
	$1->sub_assign(*$3);
	delete $3;
}
|
non_polynomial_rhs_picard '*' non_polynomial_rhs_picard
{
	$$ = $1;

	Interval intPoly1, intPoly2, intTrunc;

	$3->polyRangeNormal(intPoly2, parseSetting.step_exp_table);
	$1->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, *$3, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);

	parseSetting.ranges.push_back(intPoly1);
	parseSetting.ranges.push_back(intPoly2);
	parseSetting.ranges.push_back(intTrunc);

	delete $3;
}
|
'(' non_polynomial_rhs_picard ')'
{
	$$ = $2;
}
|
non_polynomial_rhs_picard '/' non_polynomial_rhs_picard
{
	TaylorModel tmTemp;
	$3->rec_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	$$ = $1;

	Interval intPoly1, intPoly2, intTrunc;

	tmTemp.polyRangeNormal(intPoly2, parseSetting.step_exp_table);
	$1->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, tmTemp, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);

	parseSetting.ranges.push_back(intPoly1);
	parseSetting.ranges.push_back(intPoly2);
	parseSetting.ranges.push_back(intTrunc);

	delete $3;
}
|
EXP '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->exp_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*$3 = tmTemp;
	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->sin_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*$3 = tmTemp;
	$$ = $3;
}
|
COS '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->cos_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*$3 = tmTemp;
	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->log_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*$3 = tmTemp;
	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_picard ')'
{
	TaylorModel tmTemp;
	$3->sqrt_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*$3 = tmTemp;
	$$ = $3;
}
|
non_polynomial_rhs_picard '^' NUM
{
	int exp = (int)$3;

	if(exp == 0)
	{
		Interval I(1);
		$$ = new TaylorModel(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		TaylorModel res = *$1;
		TaylorModel pow = *$1;
		
		Interval intPoly1, intPoly2, intTrunc;
		
		for(int degree = exp - 1; degree > 0;)
		{
			pow.polyRangeNormal(intPoly2, parseSetting.step_exp_table);
		
			if(degree & 1)
			{
				res.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, pow, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);
				
				parseSetting.ranges.push_back(intPoly1);
				parseSetting.ranges.push_back(intPoly2);
				parseSetting.ranges.push_back(intTrunc);
			}
			
			degree >>= 1;
			
			if(degree > 0)
			{
				pow.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, pow, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);
				
				parseSetting.ranges.push_back(intPoly1);
				parseSetting.ranges.push_back(intPoly2);
				parseSetting.ranges.push_back(intTrunc);
			}
		}
		
		$$ = new TaylorModel(res);
	}

	delete $1;
}
|
'-' non_polynomial_rhs_picard %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new TaylorModel;
	(*$$) = parseSetting.flowpipe.tms[id];

	Interval intTemp;
	(*$$).expansion.ctrunc_normal(intTemp, parseSetting.step_exp_table, parseSetting.order);
	parseSetting.ranges.push_back(intTemp);
	(*$$).remainder += intTemp;

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($2, $4);
	$$ = new TaylorModel(I, numVars);
}
;
















non_polynomial_rhs_remainder: non_polynomial_rhs_remainder '+' non_polynomial_rhs_remainder
{
	$$ = $1;
	(*$$) += (*$3);
	delete $3;
}
|
non_polynomial_rhs_remainder '-' non_polynomial_rhs_remainder
{
	$$ = $1;
	(*$$) -= (*$3);
	delete $3;
}
|
non_polynomial_rhs_remainder '*' non_polynomial_rhs_remainder
{
	$$ = new Interval;

	(*$$) = (*parseSetting.iterRange) * (*$3);
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange) * (*$1);
	(*$$) += (*$1) * (*$3);
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete $1;
	delete $3;
}
|
'(' non_polynomial_rhs_remainder ')'
{
	$$ = $2;
}
|
non_polynomial_rhs_remainder '/' non_polynomial_rhs_remainder
{
	Interval intTemp;
	rec_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	$$ = new Interval;

	(*$$) = (*parseSetting.iterRange) * intTemp;
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange) * (*$1);
	(*$$) += (*$1) * intTemp;
	++parseSetting.iterRange;
	(*$$) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete $1;
	delete $3;
}
|
EXP '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	exp_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	sin_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
COS '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	cos_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	log_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_remainder ')'
{
	Interval intTemp;
	sqrt_taylor_only_remainder(intTemp, *$3, parseSetting.iterRange, parseSetting.order);

	(*$3) = intTemp;
	$$ = $3;
}
|
non_polynomial_rhs_remainder '^' NUM
{
	int exp = (int)$3;

	if(exp == 0)
	{
		Interval intZero;
		(*$1) = intZero;
		$$ = $1;
	}
	else
	{
		Interval res(*$1);
		Interval pow(*$1);
		
		Interval intPoly1, intPoly2, intTrunc;
		
		for(int degree = exp - 1; degree > 0;)
		{
			if(degree & 1)
			{
				Interval intTemp;
				intTemp = (*parseSetting.iterRange) * pow;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange) * res;
				intTemp += pow * res;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange);
				++parseSetting.iterRange;
			
				res = intTemp;
			}
			
			degree >>= 1;
			
			if(degree > 0)
			{
				Interval intTemp;
				intTemp = (*parseSetting.iterRange) * pow;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange) * pow;
				intTemp += pow * pow;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange);
				++parseSetting.iterRange;
			
				pow = intTemp;
			}
		}
		
		$$ = new Interval(res);
	}

	delete $1;
}
|
'-' non_polynomial_rhs_remainder %prec uminus
{
	$$ = $2;
	$$->inv_assign();
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new Interval;
	(*$$) = parseSetting.flowpipe.tms[id].getRemainder();
	
	(*$$) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	$$ = new Interval;
}
;

















non_polynomial_rhs_no_remainder: non_polynomial_rhs_no_remainder '+' non_polynomial_rhs_no_remainder
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
non_polynomial_rhs_no_remainder '-' non_polynomial_rhs_no_remainder
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
non_polynomial_rhs_no_remainder '*' non_polynomial_rhs_no_remainder
{
	$$ = $1;
	(*$$) *= (*$3);
	$$->nctrunc(parseSetting.order);
	$$->cutoff(parseSetting.cutoff_threshold);

	delete $3;
}
|
'(' non_polynomial_rhs_no_remainder ')'
{
	$$ = $2;
}
|
non_polynomial_rhs_no_remainder '/' non_polynomial_rhs_no_remainder
{
	Polynomial polyTemp;
	$3->rec_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*$1) *= polyTemp;
	$1->nctrunc(parseSetting.order);
	$$ = $1;
	$$->cutoff(parseSetting.cutoff_threshold);

	delete $3;
}
|
EXP '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->exp_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*$3) = polyTemp;
	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->sin_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*$3) = polyTemp;
	$$ = $3;
}
|
COS '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->cos_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*$3) = polyTemp;
	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->log_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*$3) = polyTemp;
	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_no_remainder ')'
{
	Polynomial polyTemp;
	$3->sqrt_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*$3) = polyTemp;
	$$ = $3;
}
|
non_polynomial_rhs_no_remainder '^' NUM
{
	int exp = (int) $3;
	
	$$ = new Polynomial;
	
	(*$1).pow(*$$, exp, parseSetting.order);
	$$->cutoff(parseSetting.cutoff_threshold);
}
|
'-' non_polynomial_rhs_no_remainder %prec uminus
{
	$$ = $2;
	$$->inv_assign();
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new Polynomial;
	mlog1($$);
	mlog1(*$1);
	exit(0);
	
	parseSetting.flowpipe.tms[id].getExpansion(*$$);
	
	(*$$).nctrunc(parseSetting.order);
}
|
'[' NUM ',' NUM ']'
{
	Interval I($2, $4);
	$$ = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
}
;














non_polynomial_rhs_string: non_polynomial_rhs_string '+' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '+';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
non_polynomial_rhs_string '-' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '-';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
non_polynomial_rhs_string '*' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '*';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
'(' non_polynomial_rhs_string ')'
{
	string str;
	str += '(';
	str += (*$2);
	str += ')';
	(*$2) = str;

	$$ = $2;
}
|
non_polynomial_rhs_string '/' non_polynomial_rhs_string
{
	(*$1) += ' ';
	(*$1) += '/';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
EXP '(' non_polynomial_rhs_string ')'
{
	string str("exp");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_string ')'
{
	string str("sin");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
COS '(' non_polynomial_rhs_string ')'
{
	string str("cos");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_string ')'
{
	string str("log");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_string ')'
{
	string str("sqrt");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
non_polynomial_rhs_string '^' NUM
{
	(*$1) += '^';

	char strNum[NUM_LENGTH];
	sprintf(strNum, "%d", (int)$3);
	string num(strNum);
	(*$1) += num;

	$$ = $1;
}
|
'-' non_polynomial_rhs_string %prec uminus
{
	string str;
	str += '-';
	str += (*$2);
	(*$2) = str;

	$$ = $2;
}
|
IDENT
{
	int id = continuousProblem.getIDForPar(*$1);

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *$1);

		$$ = new string;
		range.toString(*$$);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*$1);

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		$$ = $1;
	}
}
|
'[' NUM ',' NUM ']'
{
	$$ = new string;
	char strNum_lo[NUM_LENGTH], strNum_up[NUM_LENGTH];
	sprintf(strNum_lo, "%.20e", $2);
	sprintf(strNum_up, "%.20e", $4);

	string num_lo(strNum_lo);
	string num_up(strNum_up);

	(*$$) += '[';
	(*$$) += num_lo;
	(*$$) += ' ';
	(*$$) += ',';
	(*$$) += ' ';
	(*$$) += num_up;
	(*$$) += ']';
}
|
NUM
{
	$$ = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", $1);
	string num(strNum);

	(*$$) += '[';
	(*$$) += num;
	(*$$) += ' ';
	(*$$) += ',';
	(*$$) += ' ';
	(*$$) += num;
	(*$$) += ']';
}
;


















non_polynomial_rhs_center: non_polynomial_rhs_center '+' non_polynomial_rhs_center
{
	(*$1) += ' ';
	(*$1) += '+';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
non_polynomial_rhs_center '-' non_polynomial_rhs_center
{
	(*$1) += ' ';
	(*$1) += '-';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
non_polynomial_rhs_center '*' non_polynomial_rhs_center
{
	(*$1) += ' ';
	(*$1) += '*';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
'(' non_polynomial_rhs_center ')'
{
	string str;
	str += '(';
	str += (*$2);
	str += ')';
	(*$2) = str;

	$$ = $2;
}
|
non_polynomial_rhs_center '/' non_polynomial_rhs_center
{
	(*$1) += ' ';
	(*$1) += '/';
	(*$1) += ' ';
	(*$1) += (*$3);

	$$ = $1;
	delete $3;
}
|
EXP '(' non_polynomial_rhs_center ')'
{
	string str("exp");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
SIN '(' non_polynomial_rhs_center ')'
{
	string str("sin");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
COS '(' non_polynomial_rhs_center ')'
{
	string str("cos");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
LOG '(' non_polynomial_rhs_center ')'
{
	string str("log");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
SQRT '(' non_polynomial_rhs_center ')'
{
	string str("sqrt");
	str += '(';
	str += (*$3);
	str += ')';
	(*$3) = str;

	$$ = $3;
}
|
non_polynomial_rhs_center '^' NUM
{
	(*$1) += '^';

	char strNum[NUM_LENGTH];
	sprintf(strNum, "%d", (int)$3);
	string num(strNum);
	(*$1) += num;

	$$ = $1;
}
|
'-' non_polynomial_rhs_center %prec uminus
{
	string str;
	str += '-';
	str += (*$2);
	(*$2) = str;

	$$ = $2;
}
|
IDENT
{
	int id = continuousProblem.getIDForStateVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = $1;
}
|
'[' NUM ',' NUM ']'
{
	$$ = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", ($2+$4)/2);

	string num(strNum);

	(*$$) += '[';
	(*$$) += num;
	(*$$) += ' ';
	(*$$) += ',';
	(*$$) += ' ';
	(*$$) += num;
	(*$$) += ']';
}
|
NUM
{
	$$ = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", $1);
	string num(strNum);

	(*$$) += '[';
	(*$$) += num;
	(*$$) += ' ';
	(*$$) += ',';
	(*$$) += ' ';
	(*$$) += num;
	(*$$) += ']';
}
;

my_polys: my_poly {
  $$ = new vector<Polynomial>();
  $$->push_back(*$1);
  delete $1;
} | my_poly ',' my_polys {
  $3->insert($3->begin(),*$1);
  $$ = $3;
  delete $1;
};

my_taylor_models: my_taylor_model {
  //mlog1("my taylormodels");
  $$ = new vector<TaylorModel>();
  $$->push_back(TaylorModel(*$1));
  delete $1;
} | my_taylor_model ',' my_taylor_models {
  //mlog1("my taylormodels list");
  //$3->push_back(*$1);
  $3->insert($3->begin(),*$1);
  $$ = $3;
  delete $1;
}
;


mpfr_interval: '[' MPFRNUM ',' MPFRNUM ']' {
  mpfr_t lower, upper;
  parseMpfr($2, lower);
  parseMpfr($4, upper);
  
  //mpfr_printf("lower = %.17Rg\n", lower);
  //mpfr_printf("upper = %.17Rg\n", upper);
  //logger.force("mpfr");
  //logger.force(sbuilder() << *$2);
  //logger.force(sbuilder() << *$4);
  
  Interval* ret = new Interval(lower, upper);
  //mlog1(ret->toString());
  $$ = ret;
}

my_taylor_model: my_poly {
  //mlog1("my taylor model poly");
  
  Interval c;
  $1->constant(c);
  $1->rmConstant();
  
  $$ = new TaylorModel(*$1, c);
  //mlog("tm", *$$);
  delete $1;
} | '<' my_poly ',' mpfr_interval '>' {
  $$ = new TaylorModel(*$2, *$4);
  delete $2;
  delete $4;
};

my_poly: my_poly '+' my_poly{
  //mlog1("+");
  $$ = $1;
  *$$ += *($3);
  delete $3;
} | my_poly '-' my_poly{
  //mlog1("-");
  $$ = $1;
  *$$ -= *($3);
  delete $3;
} | my_poly '^' NUM{
  int exp = (int) $3;
	if(exp == 0) {
		Interval I(1);
		$$ = new Polynomial(I, parseSetting.variables.size());
	} else {
		$$ = new Polynomial;
		(*$1).pow(*$$, exp);
	}
	delete $1;
} | my_poly '*' my_poly{
  //mlog1("*");
  $$ = $1;
  *$$ *= *($3);
  delete $3;
} | IDENT {
  //mlog1(*$1);
  vector<string> vars = parseSetting.variables;
  int dim = parseSetting.variables.size();
  int pos = find(vars.begin(), vars.end(), *$1) - vars.begin();
  if(pos >= dim) {
		char errMsg[MSG_SIZE];
		string s = sbuilder() << "variable '" << *$1 << "' wasn't declared";
		sprintf(errMsg, s.c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	vector<int> degs;
	for(int i = 0; i < dim; i++) {
	  if(i == pos) {
	    degs.push_back(1);
	    continue;
	  }
	  degs.push_back(0);
	}
  Monomial m(Interval(1), degs);
  $$ = new Polynomial(m);
} | '[' NUM ',' NUM ']' {
  int dim = parseSetting.variables.size();
  $$ = new Polynomial(Interval($2, $4), dim);
} | NUM {
  int dim = parseSetting.variables.size();
  Polynomial *p = new Polynomial(Interval($1), dim);
  $$ = new Polynomial(Interval($1), dim);
} | mpfr_interval {
  //mlog1("poly mpfr");
  int dim = parseSetting.variables.size();
  $$ = new Polynomial(Interval(*$1), dim);
  delete $1;
};

my_mono: NUM '<' my_integer_vector '>' {
	Monomial *m = new Monomial(Interval($1), *$3);
  $$ = m;
  delete $3;
} | '<' my_integer_vector '>' {
	Monomial *m = new Monomial(Interval(1), *$2);
  $$ = m;
  delete $2;
};

my_integer_vector: NUM{
  vector<int> *degrees = new vector<int>();
  degrees->push_back($1);
  $$ = degrees;
} | NUM ',' my_integer_vector {
  $3->insert($3->begin(), $1);
  $$ = $3;
};

my_interval_vector: '[' NUM ',' NUM ']' {
  Interval I($2, $4);
  vector<Interval> *ret = new vector<Interval>();
  ret->push_back(I);
  $$ = ret;
} | '[' NUM ',' NUM ']' ',' my_interval_vector {
  Interval I($2, $4);
  $7->insert($7->begin(), I);
  $$ = $7;
}
;








linear_polynomial: linear_polynomial '+' linear_polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
linear_polynomial '-' linear_polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
linear_polynomial '*' linear_polynomial
{
	if($1->degree() + $3->degree() > 1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Polynomial is not linear.");
		parseError(errMsg, lineNum);
		exit(1);
	}
	else
	{
		$$ = $1;
		(*$$) *= (*$3);
		delete $3;
	}
}
|
'(' linear_polynomial ')'
{
	$$ = $2; 
}
|
'-' linear_polynomial %prec uminus
{
	Interval I(-1);
	$$ = $2;
	$$->mul_assign(I);
}
|
IDENT
{
	int id = continuousProblem.getIDForPar(*$1);

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *$1);

		int numVars = continuousProblem.stateVarNames.size()+1;
		$$ = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*$1);

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		$$ = new Polynomial(monomial);
	}

	delete $1;
}
|
NUM
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($1);
	$$ = new Polynomial(I, numVars);
}
|
'[' NUM ',' NUM ']'
{
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I($2, $4);
	$$ = new Polynomial(I, numVars);
}
;





















%%

int yyerror(const char * what)
{
	fprintf(stderr, "Error line %d: %s\n", lineNum,what);
	err=true;
	return 1;
}

int yyerror(string what)
{
	cerr << "Error line "<<lineNum<<" "<<what<<endl;
	err=true;
	return 1;
}
