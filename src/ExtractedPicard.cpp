#include "ExtractedPicard.h"

ExtractedPicard::ExtractedPicard(const vector<Interval> & box, const Interval & I)
: Flowpipe(box, I) {
}

//class TMVSerializer;
//TMVSerializer *pSerializer;
//pSerializer = new TMVSerializer("flow2.txt", 10)

// Taylor model integration by only using Picard iteration
int Flowpipe::advance_picard2(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  if(pSerializer == NULL) {
    //pSerializer = new TMVSerializer("flow.txt", 4 + 7*(10-1), false);
    pSerializer = new TMVSerializer("flow.txt", INT_MAX, true);
  }
  mreset(old);
  mdisable();
  mlog1("Picard1 (ext) <");
  minc();
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();
  
	tstart(fl_eval);
  //mlog("step", step_end_exp_table);
  //mlog("est", estimation);
	mlog1(sbuilder() << ode.size());
	
	// evaluate the the initial set x0
	TaylorModelVec range_of_x0; //construct a empty taylormodel vec
  //mlog1(sbuilder() << "rox0: " << range_of_x0.toString(getVNames(3)));
  
  
  
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);
  
	mlog("left*", range_of_x0);
	
	pSerializer->add(range_of_x0, "leftStar");
	
	/*
	//NOTE BEGIN REMOVE
  TaylorModelVec composedBP;
  vector<Interval> rightBPRange;
  tmv.polyRangeNormal(rightBPRange, step_exp_table);
	range_of_x0.insert_ctrunc_normal(composedBP, tmv, rightBPRange, step_exp_table,
	     domain.size(), order, cutoff_threshold);
  //pSerializer->add(composedBP, "composed_before_precond");
  //END REMOVE
	*/
	
	
	
	
	
	tend(fl_eval);
	
	//pSerializer->add(range_of_x0, "leftStar");
	
  tstart(fl_precond);
  tstart(fl_part_all);
  
  tstart(fl_part1);
  tstart1(fl_pre_start);
  
	//pSerializer->add(range_of_x0, "leftStar"); //0, 4
  // the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	
  //mlog("iVC", intVecCenter);
	
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();
  tend1(fl_pre_start);
  
  tstart1(fl_pre_matrix);

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;
  tend1(fl_pre_matrix);

  tend(fl_part1);
  tstart(fl_part_matrix);
  tstart1(fl_pre_ltr);
  //serializeTMV(range_of_x0, "flow.txt");
	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;

		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);

		break;
	}
	}
  tend1(fl_pre_ltr);
  tend(fl_part_matrix);
  
  mlog("right", tmv);
  tstart(fl_part2);
  
  tstart1(fl_pre_r_range);
	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
  tend1(fl_pre_r_range);
  
  tstart1(fl_pre_insert);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange,
      step_end_exp_table, domain.size(), order, cutoff_threshold);
  tend1(fl_pre_insert);
  
  tstart1(fl_pre_scaling);
	vector<Interval> boundOfr0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0) {
	  //should never happen in our set up
		throw std::runtime_error( "invariant size wasn't zero" );
	} else {
		result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);
	}

	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0,i,i);
			invS.set(1,i,i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}
  tend1(fl_pre_scaling);

  tstart1(fl_pre_lin_trans);
	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table, cutoff_threshold);
  tend1(fl_pre_lin_trans);


  tstart1(fl_pre_left);
	TaylorModelVec x, x0;

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec c0_plus_Ar0;
	Ar0.add(c0_plus_Ar0, c0);
  tend1(fl_pre_left);

  tstart1(fl_pre_end);
	x0 = c0_plus_Ar0;
	x = x0;
	
  tend1(fl_pre_end);
  tend(fl_part2);
  tend(fl_part_all);
  tend(fl_precond);
  
  tstart(fl_int_all);
  
  
  pSerializer->add(x, "left_after_precond");
  pSerializer->add(result.tmv, "right_after_precond");
  
  /*
  //NOTE BEGIN REMOVE
  TaylorModelVec composed;
  vector<Interval> rightAPRange;
  result.tmv.polyRangeNormal(rightAPRange, step_exp_table);
	x.insert_ctrunc_normal(composed, result.tmv, rightAPRange, step_exp_table,
	     domain.size(), order, cutoff_threshold);
  //pSerializer->add(composed, "composed_after_precond");
  //END REMOVE
  */
  
	
  tstart(fl_integrate);
  
  
  
  
  tstart(fl_int_poly);
	for(int i=1; i<=order; ++i)
	{
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, i, cutoff_threshold);
	}	
	x.cutoff(cutoff_threshold);
  tend(fl_int_poly);
  
  tstart(fl_int_rem);
  
	tstart(sc_int_rem_setup);
  //pSerializer->add(x, "no_rem");
	
	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i]; // + step_uncertainties[i];		// apply the remainder estimation
	}
  //pSerializer->add(x, "inspect");

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
  tstart1(fl_ref_first_picard);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode,
      step_exp_table, rangeDim+1, order, cutoff_threshold);
  tend1(fl_ref_first_picard);

  //pSerializer->add(tmvTemp, "inspect");

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += intDifferences[i];
	}
  //pSerializer->add(tmvTemp, "inspect");

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}
	//mforce1(sbuilder() << "bfound: " << bfound);
	//mforce1("initial");

	if(!bfound)
	{
		return 0;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	tend(sc_int_rem_setup);
  
	tstart(fl_int_refine);
	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
	  //mforce1(sbuilder() << "counter: " << (rSteps + 1)); 
		bfinished = true;
    
		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);
    //pSerializer->add(x, "only_rem");

		// add the uncertainties and the cutoff intervals onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];
		}
    //mforce("new rem", newRemainders);
    


		for(int i=0; i<rangeDim; ++i)
		{

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
  	//mforce("x", x);
  	//mforce("fr", x.getRemainders());
	}
  //mforce("x", x);
	//exit(0);
	tend(fl_int_refine);
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	
  tend(fl_int_rem);
	
  tend(fl_int_all);
  
	
	mdec();
	mlog1("Picard1 >");
	mrestore(old);
  tend(fl_integrate);
  
  
  /*
  //pSerializer->add(x, "comp_left");
  //pSerializer->add(result.tmv, "comp_right");
  
  vector<Interval> rightRange;
	result.tmv.polyRangeNormal(rightRange, step_end_exp_table);
	TaylorModelVec composed;
	x.insert_ctrunc_normal(composed, result.tmv, rightRange, step_end_exp_table,
	    domain.size(), order, 0);
  
  
  //pSerializer->add(x, "composed");
  */
  
	return 1;
}



//NOTHING REMOVED FROM THIS VERSION
/*
// Taylor model integration by only using Picard iteration
int Flowpipe::advance_picard2(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  mreset(old);
  //mdisable();
  mlog1("Picard1 (ext) <");
  minc();
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();
	
  //mlog("step", step_end_exp_table);
  //mlog("est", estimation);
	mlog1(sbuilder() << ode.size());
  
	// evaluate the the initial set x0
	TaylorModelVec range_of_x0; //construct a empty taylormodel vec
  //mlog1(sbuilder() << "rox0: " << range_of_x0.toString(getVNames(3)));
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);
	serializeTMV(range_of_x0, "flow.txt");
  //mlog(step_end_exp_table);
  //mlog1(sbuilder() << "rox0: " << range_of_x0.toString(getVNames(3)));
  //mlog1(sbuilder() << "tmvPre: " << tmvPre.toString(getVNames(3)));
	
  // the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	
  mlog("iVC", intVecCenter);
	
	// the center point of the remainder of x0
	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		range_of_x0.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	TaylorModelVec c0(intVecCenter, rangeDim+1);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	range_of_x0.rmConstant();

	// compute the preconditioning matrix
	Matrix A(rangeDim,rangeDim), invA(rangeDim,rangeDim);
	TaylorModelVec range_of_r0;

	switch(precondition)
	{
	case ID_PRE:
	{
		range_of_r0 = range_of_x0;

		break;
	}

	case QR_PRE:
	{
		preconditionQR(A, range_of_x0, rangeDim, rangeDim+1);
		A.transpose(invA);
		range_of_x0.linearTrans(range_of_r0, invA);

		break;
	}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_end_exp_table);
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), order, cutoff_threshold);

	vector<Interval> boundOfr0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		vector<Interval> polyRangeOfr0;
		result.tmv.polyRangeNormal(polyRangeOfr0, step_end_exp_table);

		Matrix matPre_constraints(rangeDim, rangeDim+1);

		switch(precondition)
		{
		case ID_PRE:
		{
			for(int i=0; i<rangeDim; ++i)
			{
				matPre_constraints.set( 1 , i, i+1);
			}

			break;
		}

		case QR_PRE:
		{
			for(int i=0; i<rangeDim; ++i)
			{
				for(int j=0; j<rangeDim; ++j)
				{
					matPre_constraints.set( A.get(i,j) , i, j+1);
				}
			}

			break;
		}
		}

		TaylorModelVec tmvPrecond(matPre_constraints);
		tmvPrecond.add_assign(c0);

		vector<HornerForm> preconditioned_constraints;
		vector<Interval> b;
		for(int i=0; i<invariant.size(); ++i)
		{
			Polynomial polyTemp;
			invariant[i].hf.insert_no_remainder_no_cutoff(polyTemp, tmvPrecond, rangeDim+1);

			HornerForm hf;
			polyTemp.toHornerForm(hf);
			preconditioned_constraints.push_back(hf);

			b.push_back(invariant[i].B);
		}

		vector<Interval> contracted_remainders;
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders.push_back(result.tmv.tms[i].remainder);
		}

		int res = contract_remainder(polyRangeOfr0, contracted_remainders, preconditioned_constraints, b);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			boundOfr0.push_back(polyRangeOfr0[i] + contracted_remainders[i]);
		}
	}
	else
	{
		result.tmv.intEvalNormal(boundOfr0, step_end_exp_table);
	}


	// Compute the scaling matrix S.
	Matrix S(rangeDim, rangeDim);
	Matrix invS(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intSup;
		boundOfr0[i].mag(intSup);

		if(!intSup.subseteq(intZero))
		{
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			boundOfr0[i] = intUnit;
		}
		else
		{
			S.set(0,i,i);
			invS.set(1,i,i);
		}
	}

	// apply the scaling matrix S, i.e., A S S^(-1) A^(-1), then the interval bound of r0 is normalized to [-1,1]^n

	switch(precondition)
	{
	case ID_PRE:
	{
		A = S;
		invA = invS;
		break;
	}

	case QR_PRE:
	{
		A = A * S;
		invA = invS * invA;
		break;
	}
	}

	result.tmv.linearTrans_assign(invS);
	result.tmv.cutoff_normal(step_end_exp_table, cutoff_threshold);


	TaylorModelVec x, x0;

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec c0_plus_Ar0;
	Ar0.add(c0_plus_Ar0, c0);

	x0 = c0_plus_Ar0;
	x = x0;

	for(int i=1; i<=order; ++i)
	{
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, i, cutoff_threshold);
	}

	x.cutoff(cutoff_threshold);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i]; // + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);

		intDifferences.push_back(intTemp);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return 0;
	}
	else
	{
		for(int i=0; i<rangeDim; ++i)
			x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		// add the uncertainties and the cutoff intervals onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];
		}

		for(int i=0; i<rangeDim; ++i)
		{

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}
			}
			else
			{
				bfinished = true;
				break;
			}

			x.tms[i].remainder = newRemainders[i];
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	mdec();
	mlog1("Picard1 >");
	mrestore(old);
	return 1;
}*/

