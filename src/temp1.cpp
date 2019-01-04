int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
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
	return 1;
}


void TaylorModelVec::Picard_ctrunc_normal(TaylorModelVec & result, vector<RangeTree *> & trees, 
    const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, 
    const vector<Interval> & step_exp_table, const int numVars, const int order, 
    const Interval & cutoff_threshold) const {
  //mforce1("my version shouldn't call this1");
  //exit(0);
  mreset(old);
  mdisable();
  minc();
  //mlog1("Picard_ctrunc_normal2");
  //mlog("x0", x0);
  //mlog("ode", ode);
  //mlog("this", *this);
  //mlog("polyRange", polyRange);
  //mlog1(sbuilder() << "order: " << order);
  //mlog1(sbuilder() << "numVars: " << numVars);
	TaylorModelVec tmvTemp;

	trees.clear();
	for(int i=0; i<ode.size(); ++i)
	{
		trees.push_back(NULL);
	}
	if(order <= 1)
	{
		for(int i=0; i<ode.size(); ++i)
		{
			TaylorModel tmTemp;
			//mlog1(ode[i].toString());
			ode[i].insert_ctrunc_normal(tmTemp, trees[i], *this, polyRange, step_exp_table, numVars, 0, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
			//mlog(sbuilder() << "tm[" << i << "]", tmTemp);
		}
	}
	else
	{
		for(int i=0; i<ode.size(); ++i)
		{
		  //order - 1, since the integration makes the degree 1 higher
			TaylorModel tmTemp;
			ode[i].insert_ctrunc_normal(tmTemp, trees[i], *this, polyRange, step_exp_table, numVars, order-1, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
	}
	//mlog("tmvTemp(TM)", tmvTemp);

	TaylorModelVec tmvTemp2;
	tmvTemp.integral(tmvTemp2, step_exp_table[1]);
	//mlog("tmvTemp2(TM)", tmvTemp2);

	x0.add(result, tmvTemp2);
	//mlog("result", result);
	mdec();
	mrestore(old);
  //mforce("result", result);
}

// Taylor model integration by only using Picard iteration
int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, 
    const vector<HornerForm> & ode_centered, const int precondition, 
    vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, 
    const int order, const vector<Interval> & estimation, 
    const vector<PolynomialConstraint> & invariant, 
    const Interval & cutoff_threshold) const
{
  mlog1("Picard1 <");
  minc();
	int rangeDim = ode.size();
	Interval intZero, intOne(1,1), intUnit(-1,1);
	result.clear();
	
	
  vector<Interval>::iterator i1 = step_end_exp_table.begin();
	string s;
	(*i1).toString(s);
	mlog1(sbuilder() << s);
	for (unsigned i=0; i<estimation.size(); i++) {
		string s;
		estimation.at(i).toString(s);
		mlog1(sbuilder() << "estimation[" <<  i << "] = " << s);
	}
	mlog1(sbuilder() << ode.size());

	// evaluate the the initial set x0
	TaylorModelVec range_of_x0;
	tmvPre.evaluate_t(range_of_x0, step_end_exp_table);
	

	// the center point of x0's polynomial part
	vector<Interval> intVecCenter;
	range_of_x0.constant(intVecCenter);
	
	
	for (unsigned i=0; i<intVecCenter.size(); i++) {
		string s;
		intVecCenter.at(i).toString(s);
		//mlog1(sbuilder() << "intVecCenter[" << i << "] = " << s);
	}

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
	return 1;
}