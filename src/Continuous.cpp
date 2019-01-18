/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Continuous.h"
//#include "Utils.h"


void vspode_error(int mode) ;




// vector<LinearConstraint> constraint_template;

Flowpipe::Flowpipe()
{
}

Flowpipe::Flowpipe(const TaylorModelVec & tmvPre_input, const TaylorModelVec & tmv_input, const vector<Interval> & domain_input):
		tmvPre(tmvPre_input), tmv(tmv_input), domain(domain_input)
{
}

Flowpipe::Flowpipe(const vector<Interval> & box, const Interval & I)
{
	int rangeDim = box.size();
	int domainDim = rangeDim + 1;
	Interval intUnit(-1,1), intZero;

	TaylorModelVec tmvCenter;
	vector<double> scalars;

	domain.push_back(I);		// time interval

	// normalize the domain to [-1,1]^n
	for(int i=0; i<rangeDim; ++i)
	{
		double midpoint = box[i].midpoint();
		Interval intMid(midpoint);
		TaylorModel tmTemp(intMid, domainDim);
		tmvCenter.tms.push_back(tmTemp);

		Interval intTemp = box[i];
		intTemp.sub_assign(midpoint);
		scalars.push_back( intTemp.sup() );
		domain.push_back(intUnit);
	}

	Matrix coefficients_of_tmvPre(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients_of_tmvPre.set(scalars[i], i, i+1);
	}

	TaylorModelVec tmvTemp(coefficients_of_tmvPre);
	tmvTemp.add(tmvPre, tmvCenter);

	Matrix coefficients_of_tmv(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients_of_tmv.set(1, i, i+1);
	}

	TaylorModelVec tmvTemp2(coefficients_of_tmv);
	tmv = tmvTemp2;
}

Flowpipe::Flowpipe(const TaylorModelVec & tmv_input, const vector<Interval> & domain_input)
{
	int rangeDim = tmv_input.tms.size();

	tmv = tmv_input;
	domain = domain_input;

	Matrix coefficients_of_tmvPre(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients_of_tmvPre.set(1, i, i+1);
	}

	TaylorModelVec tmvTemp(coefficients_of_tmvPre);
	tmvPre = tmvTemp;

	normalize();
}

Flowpipe::Flowpipe(const Flowpipe & flowpipe):tmvPre(flowpipe.tmvPre), tmv(flowpipe.tmv), domain(flowpipe.domain)
{
}

Flowpipe::~Flowpipe()
{
	clear();
}

void Flowpipe::clear()
{
	tmvPre.clear();
	tmv.clear();
	domain.clear();
}

void Flowpipe::dump(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;

	composition(tmvTemp, cutoff_threshold);

	// dump the Taylor model
	tmvTemp.dump_interval(fp, stateVarNames, tmVarNames);

	//dump the domain
	for(int i=0; i<domain.size(); ++i)
	{
		fprintf(fp, "%s in ", tmVarNames[i].c_str());
		domain[i].dump(fp);
		fprintf(fp, "\n");
	}
}

void Flowpipe::dump_normal(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames, vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;

	composition_normal(tmvTemp, step_exp_table, cutoff_threshold);

	// dump the Taylor model
	tmvTemp.dump_interval(fp, stateVarNames, tmVarNames);

	//dump the domain
	for(int i=0; i<domain.size(); ++i)
	{
		fprintf(fp, "%s in ", tmVarNames[i].c_str());
		domain[i].dump(fp);
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
}

void Flowpipe::composition(TaylorModelVec & result, const Interval & cutoff_threshold) const
{
	vector<int> orders;

	for(int i=0; i<tmv.tms.size(); ++i)
	{
		int d1 = tmv.tms[i].degree();
		int d2 = tmvPre.tms[i].degree();

		if(d1 > d2)
		{
			orders.push_back(d1);
		}
		else
		{
			orders.push_back(d2);
		}
	}

	vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);
	tmvPre.insert_ctrunc(result, tmv, tmvPolyRange, domain, orders, cutoff_threshold);
}

void Flowpipe::composition_normal(TaylorModelVec & result, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const
{
	int domainDim = domain.size();

	vector<int> orders;

	for(int i=0; i<tmv.tms.size(); ++i)
	{
		int d1 = tmv.tms[i].degree();
		int d2 = tmvPre.tms[i].degree();

		if(d1 > d2)
		{
			orders.push_back(d1);
		}
		else
		{
			orders.push_back(d2);
		}
	}
	//mlog1();
	//mlog("tmv", tmv);
	//mlog("tmvPre", tmvPre);
	
	vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_exp_table);
	tmvPre.insert_ctrunc_normal(result, tmv, tmvPolyRange, step_exp_table, 
	    domainDim, orders, cutoff_threshold);
	
	//pSerializer->add(result, "comp_comp");
	//mlog("result", result);
	//mlog1();
}

void Flowpipe::intEval(vector<Interval> & result, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
	composition(tmvTemp, cutoff_threshold);

	tmvTemp.intEval(result, domain);
}

void Flowpipe::intEvalNormal(vector<Interval> & result, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
	composition_normal(tmvTemp, step_exp_table, cutoff_threshold);
	tmvTemp.intEvalNormal(result, step_exp_table);
}

void Flowpipe::normalize()
{
	Interval intZero;

	// we first normalize the Taylor model tmv
	tmv.normalize(domain);

	int rangeDim = tmv.tms.size();

	// compute the center point of tmv
	vector<Interval> intVecCenter;
	tmv.constant(intVecCenter);
	tmv.rmConstant();

	for(int i=0; i<rangeDim; ++i)
	{
		Interval M;
		tmv.tms[i].remainder.remove_midpoint(M);
		intVecCenter[i] += M;
	}

	vector<Interval> tmvRange;
	tmv.intEval(tmvRange, domain);

	vector<vector<Interval> > coefficients;
	vector<Interval> row;

	for(int i=0; i<rangeDim+1; ++i)
	{
		row.push_back(intZero);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		coefficients.push_back(row);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		Interval intScalor;
		tmvRange[i].mag(intScalor);

		if(intScalor.subseteq(intZero))
		{
			coefficients[i][i+1] = intZero;
		}
		else
		{
			coefficients[i][i+1] = intScalor;
			tmv.tms[i].div_assign(intScalor);
		}
	}

	TaylorModelVec newVars(coefficients);
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp(intVecCenter[i], rangeDim+1);
		newVars.tms[i].add_assign(tmTemp);
	}

	for(int i=0; i<tmvPre.tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tmvPre.tms[i].insert_no_remainder_no_cutoff(tmTemp, newVars, rangeDim+1, tmvPre.tms[i].degree());
		tmvPre.tms[i].expansion = tmTemp.expansion;
	}
}











// Taylor model integration by only using Picard iteration
int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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

int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  mlog1("Picard2 <");
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	TaylorModelVec c_plus_Ar0;
	Ar0.add(c_plus_Ar0, c0);

	x0 = c_plus_Ar0;
	x = x0;

	for(int i=1; i<=globalMaxOrder; ++i)
	{
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, i, cutoff_threshold);
	}

	x.nctrunc(orders);

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
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

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
  mlog1("Picard2 >");
	return 1;
}


// adaptive step sizes and fixed orders
int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  mlog1("Picard3 <");
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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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
  mlog1("Picard3 >");
	return 1;
}

int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  mlog1("Picard4 <");
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	TaylorModelVec c_plus_Ar0;
	Ar0.add(c_plus_Ar0, c0);

	x0 = c_plus_Ar0;
	x = x0;

	for(int i=1; i<=globalMaxOrder; ++i)
	{
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, i, cutoff_threshold);
	}

	x.nctrunc(orders);

	x.cutoff(cutoff_threshold);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return -1;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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
  mlog1("Picard4 >");
	return 1;
}

// adaptive orders and fixed step sizes
int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  mlog1("Picard5 <");
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
			S.set(0, i, i);
			invS.set(1, i, i);
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
		x.tms[i].remainder = estimation[i];
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

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return 0;
		}

		// increase the approximation orders by 1
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, newOrder, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i];// + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrder, cutoff_threshold);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		// add the uncertainties onto the result
		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		// add the uncertainties
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

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
  mlog1("Picard5 >");
	return 1;
}

int Flowpipe::advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
  mlog1("Picard5 <");
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	TaylorModelVec c_plus_Ar0;
	Ar0.add(c_plus_Ar0, c0);

	x0 = c_plus_Ar0;
	x = x0;

	for(int i=1; i<=localMaxOrder; ++i)
	{
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, i, cutoff_threshold);
	}

	x.nctrunc(orders);

	x.cutoff(cutoff_threshold);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

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

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return 0;
		}

		// increase the approximation orders
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, newOrders, bIncreased, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = estimation[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
  mlog1("Picard6 >");
	return 1;
}


// for low-degree ODEs
// fixed step sizes and orders

int Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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

	// the center point of the remainder
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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, order, cutoff_threshold);
		x.tms.push_back(tmTemp);
	}

	x.cutoff(cutoff_threshold);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
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

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders.push_back(tmvTemp.tms[i].remainder);
		}

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

int Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, orders[i], cutoff_threshold);
		x.tms.push_back(tmTemp);
	}

	x.cutoff(cutoff_threshold);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

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


// adaptive step sizes and fixed orders
int Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, order, cutoff_threshold);
		x.tms.push_back(tmTemp);
	}

	x.cutoff(cutoff_threshold);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

int Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, orders[i], cutoff_threshold);
		x.tms.push_back(tmTemp);
	}

	x.cutoff(cutoff_threshold);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

// adaptive orders and fixed step sizes
int Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, order, cutoff_threshold);
		x.tms.push_back(tmTemp);
	}

	x.cutoff(cutoff_threshold);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
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

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return 0;
		}

		// increase the approximation orders by 1
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, newOrder, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = estimation[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrder, cutoff_threshold);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

		// add the uncertainties
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

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return 1;
}

int Flowpipe::advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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
			S.set(0, i, i);
			invS.set(1, i, i);
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

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	TaylorModelVec x;

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		TaylorModel tmTemp;
		taylorExpansion[i].insert_no_remainder(tmTemp, x0, rangeDim+1, orders[i], cutoff_threshold);
		x.tms.push_back(tmTemp);
	}

	x.cutoff(cutoff_threshold);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

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

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( !tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return 0;
		}

		// increase the approximation orders
		x.Picard_no_remainder_assign(x0, ode, rangeDim+1, newOrders, bIncreased, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = estimation[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return 1;
}























// for high-degree ODEs
// fixed step sizes and orders

int Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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


	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1,i,i);
			invS.set(1,i,i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode_centered, rangeDim+1, i, cutoff_threshold);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode_centered[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, order - 1, cutoff_threshold);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
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

int Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode_centered, rangeDim+1, i, cutoff_threshold);	// compute c(t)
	}

	c.nctrunc(orders);

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode_centered[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, orders[i] - 1, cutoff_threshold);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

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

// adaptive step sizes and fixed orders
int Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode_centered, rangeDim+1, i, cutoff_threshold);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode_centered[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, order - 1, cutoff_threshold);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, order, cutoff_threshold);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

int Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode_centered, rangeDim+1, i, cutoff_threshold);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode_centered[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, orders[i] - 1, cutoff_threshold);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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


// adaptive orders and fixed step sizes
int Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode_centered, rangeDim+1, i, cutoff_threshold);	// compute c(t)
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode_centered[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, order - 1, cutoff_threshold);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
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

	// add the uncertainties onto the reault
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

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return 0;
		}

		// increase the approximation orders by 1
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, newOrder, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)	// apply the estimation again
		{
			x.tms[i].remainder = estimation[i];// + step_uncertainties[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrder, cutoff_threshold);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return 1;
}

int Flowpipe::advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=localMaxOrder; ++i)
	{
		c.Picard_no_remainder_assign(c0, ode_centered, rangeDim+1, i, cutoff_threshold);	// compute c(t)
	}

	c.nctrunc(orders);

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)
	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		ode_centered[i].insert_no_remainder(tmTemp, c_plus_Ar, rangeDim+1, orders[i] - 1, cutoff_threshold);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];		// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	vector<RangeTree *> trees;

	vector<Interval> xPolyRange;
	x.polyRangeNormal(xPolyRange, step_exp_table);
	x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, orders, cutoff_threshold);

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

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return 0;
		}

		// increase the approximation orders
		x.Picard_no_remainder_assign(x0, ode_centered, rangeDim+1, newOrders, bIncreased, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = estimation[i];
		}

		// compute the Picard operation again
		x.polyRangeNormal(xPolyRange, step_exp_table);
		x.Picard_ctrunc_normal(tmvTemp, trees, x0, xPolyRange, ode, step_exp_table, rangeDim+1, newOrders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_only_remainder(newRemainders, trees, x0, ode, step_exp_table[1]);

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

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	trees.clear();
	return 1;
}





// integration scheme for non-polynomial ODEs (using Taylor approximations)
// fixed step sizes and orders
int Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = strOde.size();
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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde_centered, i, cutoff_threshold);
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.order = order - 1;
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.cutoff_threshold = cutoff_threshold;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde_centered[i] + suffix;

		parseODE();		// call the parser
    //mforce1("1");exit(0);
    
		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, cutoff_threshold);

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
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], order);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders.push_back(tmvTemp.tms[i].remainder);
		}

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

	return 1;
}

int Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = strOde.size();
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde_centered, i, cutoff_threshold);
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.cutoff_threshold = cutoff_threshold;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde_centered[i] + suffix;
		parseSetting.order = orders[i] - 1;

				parseODE();		// call the parser
				//mforce1("2");exit(0);

		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, cutoff_threshold);

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
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], orders);

		for(int i=0; i<rangeDim; ++i)
		{
			newRemainders.push_back(tmvTemp.tms[i].remainder);
		}

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

	return -1;
}


// adaptive step sizes and fixed orders
int Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = strOde.size();
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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde_centered, i, cutoff_threshold);
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.order = order - 1;
	parseSetting.cutoff_threshold = cutoff_threshold;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde_centered[i] + suffix;

		parseODE();		// call the parser
    //mforce1("3");exit(0);
    
		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*order);

		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, cutoff_threshold);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], order);

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

	return 1;
}

int Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = strOde.size();
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=globalMaxOrder; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde_centered, i, cutoff_threshold);
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.cutoff_threshold = cutoff_threshold;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde_centered[i] + suffix;
		parseSetting.order = orders[i] - 1;

		parseODE();		// call the parser
    //mforce1("4");exit(0);
    
		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	if(step > THRESHOLD_HIGH)		// step size is changed
	{
		construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);
	}

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, cutoff_threshold);

	// compute the interval evaluation of the polynomial difference, this part is not able to be reduced by Picard iteration
	vector<Interval> intDifferences;
	vector<Polynomial> polyDifferences;
	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial polyTemp;
		polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTemp);

		Interval intTemp;
		polyTemp.intEvalNormal(intTemp, step_exp_table);
		intDifferences.push_back(intTemp);
	}

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

	for(; !bfound;)
	{
		bfound = true;
		double newStep = step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < miniStep)
		{
			return 0;
		}

		construct_step_exp_table(step_exp_table, step_end_exp_table, newStep, 2*globalMaxOrder);

		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, cutoff_threshold);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], step_exp_table);
		}

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
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], orders);

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

	return 1;
}


// adaptive orders and fixed step sizes
int Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = strOde.size();
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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=order; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde_centered, i, cutoff_threshold);
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.order = order - 1;
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.cutoff_threshold = cutoff_threshold;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde_centered[i] + suffix;
		mlog1(prefix);
		mlog1(suffix);
		mlog1(strOde_centered[i]);

		parseODE();		// call the parser
    //mforce1("5");exit(0);
    
		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, order, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, order, cutoff_threshold);

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

	int newOrder = order;

	for(; !bfound;)
	{
		++newOrder;

		if(newOrder > maxOrder)
		{
			return 0;
		}

		// increase the approximation orders by 1
		x.Picard_non_polynomial_taylor_no_remainder_assign(x0, strOde_centered, newOrder, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = estimation[i];
		}

		// compute the Picard operation again
		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, newOrder, cutoff_threshold);

		// Update the irreducible part
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial polyTemp;
			polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

			Interval intTemp;
			polyTemp.intEvalNormal(intTemp, step_exp_table);
			intDifferences[i] = intTemp;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;
		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], newOrder);

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

	order = newOrder;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return 1;
}

int Flowpipe::advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const
{
	int rangeDim = strOde.size();
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
	range_of_r0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, step_end_exp_table, domain.size(), orders, cutoff_threshold);

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

	// vector storing the zero ranges
	vector<int> zeroIndices;

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
			S.set(1, i, i);
			invS.set(1, i, i);

			zeroIndices.push_back(i+1);
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

	// compute the Taylor expansion of r(t)
	// since r(t) = A^{-1} * (x(t) - c(t)), we have that r'(t) = A^{-1} * (x'(t) - c'(t)) = A^{-1} * (f(c(t) + A*r(t), t) - c'(t))

	TaylorModelVec c = c0;
	for(int i=1; i<=localMaxOrder; ++i)
	{
		c.Picard_non_polynomial_taylor_no_remainder_assign(c0, strOde_centered, i, cutoff_threshold);
	}

	TaylorModelVec dcdt;
	c.derivative(dcdt, 0);	// compute dc/dt

	Matrix matCoefficients_Ar0(rangeDim, rangeDim+1);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			matCoefficients_Ar0.set( A.get(i,j) , i, j+1);
		}
	}

	TaylorModelVec Ar0(matCoefficients_Ar0);
	TaylorModelVec Ar = Ar0;

	TaylorModelVec c_plus_Ar;
	Ar.add(c_plus_Ar, c);

	// compute the Taylor expansion of the ODE of A*r(t)

	parseSetting.clear();
	parseSetting.flowpipe = c_plus_Ar;
	parseSetting.cutoff_threshold = cutoff_threshold;

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	TaylorModelVec Adrdt;
	for(int i=0; i<rangeDim; ++i)
	{
		parseSetting.strODE = prefix + strOde_centered[i] + suffix;
		parseSetting.order = orders[i] - 1;

		parseODE();		// call the parser
    //mforce1("6");exit(0);
    
		TaylorModel tmTemp(parseResult.expansion, intZero);
		Adrdt.tms.push_back(tmTemp);
	}

	Adrdt.sub_assign(dcdt);

	TaylorModelVec drdt;
	Adrdt.linearTrans(drdt, invA);

	TaylorModelVec taylorExp_Ar;

	computeTaylorExpansion(taylorExp_Ar, Adrdt, drdt, orders, cutoff_threshold);

	// remove the zero terms
	taylorExp_Ar.rmZeroTerms(zeroIndices);
	Ar0.rmZeroTerms(zeroIndices);

	taylorExp_Ar.add_assign(Ar0);

	TaylorModelVec x = c;
	x.add_assign(taylorExp_Ar);		// the Taylor expansion of x(t)

	x.cutoff(cutoff_threshold);

	TaylorModelVec x0;
	Ar0.add(x0, c0);

	bool bfound = true;

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = estimation[i];// + step_uncertainties[i];	// apply the remainder estimation
	}

	TaylorModelVec tmvTemp;
	x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, orders, cutoff_threshold);

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

	for(int i=0; i<rangeDim; ++i)
	{
		tmvTemp.tms[i].remainder += intDifferences[i];
	}

	vector<bool> bIncrease;
	for(int i=0; i<rangeDim; ++i)
	{
		if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bIncrease.push_back(true);
			if(bfound)
				bfound = false;
		}
		else
		{
			bIncrease.push_back(false);
		}
	}

	vector<int> newOrders = orders;
	bool bIncreaseOthers = false;
	int numIncrease = 0;

	vector<bool> bIncreased;
	for(int i=0; i<rangeDim; ++i)
	{
		bIncreased.push_back(false);
	}

	for(; !bfound;)
	{
		bool bChanged = false;

		if(bIncreaseOthers)
		{
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(!bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			bIncreaseOthers = false;
		}
		else
		{
			numIncrease = 0;
			for(int i=0; i<bIncrease.size(); ++i)
			{
				if(bIncrease[i] && newOrders[i] < maxOrders[i])
				{
					++newOrders[i];
					++numIncrease;

					bIncreased[i] = true;

					if(!bChanged)
						bChanged = true;
				}
			}

			if(numIncrease < newOrders.size())
				bIncreaseOthers = true;
		}

		if(!bChanged)
		{
			return 0;
		}

		// increase the approximation orders
		x.Picard_non_polynomial_taylor_no_remainder_assign(x0, strOde_centered, newOrders, bIncreased, cutoff_threshold);
		x.cutoff(cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = estimation[i];
		}

		// compute the Picard operation again
		x.Picard_non_polynomial_taylor_ctrunc_normal(tmvTemp, x0, strOde, step_exp_table, newOrders, cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			// update the irreducible part if necessary
			if(bIncreased[i])
			{
				Polynomial polyTemp;
				polyTemp = tmvTemp.tms[i].expansion - x.tms[i].expansion;

				Interval intTemp;
				polyTemp.intEvalNormal(intTemp, step_exp_table);
				intDifferences[i] = intTemp;
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			tmvTemp.tms[i].remainder += intDifferences[i];
		}

		bfound = true;

		for(int i=0; i<rangeDim; ++i)
		{
			if( ! tmvTemp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = true;
					if(bfound)
						bfound = false;
				}
				else
				{
					bfound = false;
					break;
				}
			}
			else
			{
				if(!bIncreaseOthers)
				{
					bIncrease[i] = false;
				}
			}
		}

		for(int i=0; i<rangeDim; ++i)
		{
			bIncreased[i] = false;
		}
	}

	for(int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTemp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		vector<Interval> newRemainders;
		x.Picard_non_polynomial_taylor_only_remainder(newRemainders, x0, strOde, step_exp_table[1], newOrders);

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

	orders = newOrders;
	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = step_exp_table[1];

	return 1;
}

Flowpipe & Flowpipe::operator = (const Flowpipe & flowpipe)
{
	if(this == &flowpipe)
		return *this;

	tmvPre = flowpipe.tmvPre;
	tmv = flowpipe.tmv;
	domain = flowpipe.domain;
	return *this;
}








































// class Continuous_system

ContinuousSystem::ContinuousSystem()
{
}

ContinuousSystem::ContinuousSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input)
{
	int rangeDim = ode_input.tms.size();
	Interval intZero;

	initialSet = initialSet_input;
	tmvOde = ode_input;

	for(int i=0; i<rangeDim; ++i)
	{
		HornerForm hf;
		tmvOde.tms[i].expansion.toHornerForm(hf);
		hfOde.push_back(hf);
	}

	tmvOde_centered = ode_input;
	tmvOde_centered.center_nc();

	for(int i=0; i<rangeDim; ++i)
	{
		HornerForm hf;
		tmvOde_centered.tms[i].expansion.toHornerForm(hf);
		hfOde_centered.push_back(hf);
	}
}

ContinuousSystem::ContinuousSystem(const vector<string> & strOde_input, const Flowpipe & initialSet_input)
{
	strOde = strOde_input;

	string prefix(str_prefix_center);
	string suffix(str_suffix);

	for(int i=0 ;i<strOde.size(); ++i)
	{
		parseSetting.clear();
		parseSetting.strODE = prefix + strOde[i] + suffix;
		parseODE();
    //mforce1("7");exit(0);
    
		strOde_centered.push_back(parseResult.strExpansion);
	}

	initialSet = initialSet_input;
}

ContinuousSystem::ContinuousSystem(const ContinuousSystem & system)
{
	tmvOde				= system.tmvOde;
	hfOde				= system.hfOde;
	initialSet			= system.initialSet;
	strOde				= system.strOde;
	tmvOde_centered		= system.tmvOde_centered;
	hfOde_centered		= system.hfOde_centered;
	strOde_centered		= system.strOde_centered;
}

ContinuousSystem::~ContinuousSystem()
{
	hfOde.clear();
	hfOde_centered.clear();
	strOde.clear();
	strOde_centered.clear();
}

void ContinuousSystem::reach_linear(list<Flowpipe> & results, const double step, const double time, const int order, const bool bPrint, const Interval & cutoff_threshold) const
{
	Interval intZero, intOne(1), intMOne(-1);
	int rangeDim = tmvOde.tms.size();

	vector<Interval> intVecTemp;
	intVecTemp.push_back(intZero);

	for(int i=0; i<rangeDim; ++i)
	{
		intVecTemp.push_back(intZero);
	}

	vector<vector<Interval> > coefficients;
	for(int i=0; i<rangeDim; ++i)
	{
		coefficients.push_back(intVecTemp);
	}

	tmvOde.linearCoefficients(coefficients);

	// compute the matrix A in dx/dt = Ax + u
	for(int i=0; i<rangeDim; ++i)
	{
		coefficients[i].erase( coefficients[i].begin() );
	}

	Interval_matrix A(coefficients);

	Interval_matrix Ar = A;
	Ar.mul_assign(step);

	Interval_matrix identity(rangeDim, rangeDim);

	// identity matrix
	for(int i=0; i<rangeDim; ++i)
	{
		identity.set(intOne, i, i);
	}

	Interval_matrix expAr_ts, expAr_rem;
	Interval_matrix int_expmAs_ts, int_expmAs_rem;

	// compute an interval matrix for exp(A*step)
	exp_int_mat(expAr_ts, expAr_rem, Ar, order);
	int_exp_int_mat(int_expmAs_ts, int_expmAs_rem, Ar, step, order);

	Interval_matrix expAr = expAr_ts + expAr_rem;
	Interval_matrix int_expmAs = int_expmAs_ts + int_expmAs_rem;

	Interval_matrix transConst = expAr*int_expmAs;

	vector<bool> pow_expAr_computed;
	vector<Interval_matrix> pow_expAr;
	Interval_matrix imEmpty;

	int N = time/step + 2;
	for(int i=0; i<N; ++i)
	{
		pow_expAr_computed.push_back(false);
		pow_expAr.push_back(imEmpty);
	}

	pow_expAr_computed[0] = true;
	pow_expAr[0] = identity;

	pow_expAr_computed[1] = true;
	pow_expAr[1] = expAr;

	for(int i=N-1; i>1; --i)
	{
		if(!pow_expAr_computed[i])
		{
			Interval_matrix temp = expAr;
			Interval_matrix result = expAr;
			int pos_temp = 1;
			int pos_result = 1;

			for(int d=i-1; d > 0;)
			{
				if(d & 1)
				{
					pos_result += pos_temp;

					if(pow_expAr_computed[pos_result])
					{
						result = pow_expAr[pos_result];
					}
					else
					{
						result *= temp;
						pow_expAr_computed[pos_result] = true;
						pow_expAr[pos_result] = result;
					}
				}

				d >>= 1;

				if(d > 0)
				{
					pos_temp <<= 1;

					if(pow_expAr_computed[pos_temp])
					{
						temp = pow_expAr[pos_temp];
					}
					else
					{
						temp *= temp;
						pow_expAr_computed[pos_temp] = true;
						pow_expAr[pos_temp] = temp;
					}
				}
			}
		}
	}

	// compute a superset for the constant part
	vector<Interval> constant_part;
	tmvOde.constant(constant_part);

	Interval_matrix imConstant(constant_part);
	Interval_matrix B = transConst * imConstant;


	// compute the first Taylor model flowpipe
	vector<Polynomial> polyOde;
	for(int i=0; i<rangeDim; ++i)
	{
		polyOde.push_back(tmvOde.tms[i].expansion);
	}


	vector<Polynomial> taylorExpansion;
	vector<Polynomial> LieDeriv_n;

	for(int i=0; i<rangeDim; ++i)
	{
		RowVector row(rangeDim+1);
		row.set(1,i+1);
		Polynomial P(row);
		taylorExpansion.push_back(P);
		LieDeriv_n.push_back(P);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=order; ++j)
		{
			Polynomial P1;
			LieDeriv_n[i].LieDerivative(P1, polyOde);
			LieDeriv_n[i] = P1;

			P1.mul_assign(factorial_rec[j]);
			P1.mul_assign(0,j);

			taylorExpansion[i] += P1;
		}
	}

	vector<Interval> rangeOfx0;
	initialSet.intEval(rangeOfx0, cutoff_threshold);
	Interval_matrix imRangeOfx0(rangeOfx0);

	TaylorModelVec x;

	vector<Interval> step_exp_table;
	vector<Interval> step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	vector<Interval> polyRangeX0;
	initialSet.tmvPre.polyRangeNormal(polyRangeX0, step_exp_table);

	// evaluate a remainder for the first TM flowpipe

	double max_Ar = A.max_norm();
	max_Ar *= step;
	Interval intTemp(max_Ar);
	intTemp.exp_assign();

	double maxNorm = intTemp.mag();
	Interval remainder(-maxNorm, maxNorm);

	Interval_matrix enclExpAr(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			enclExpAr.set(remainder, i, j);
		}
	}

	Interval_matrix remainder_x = expAr_rem * imRangeOfx0 + enclExpAr * int_expmAs_rem * imConstant;

	for(int i=0; i<rangeDim; ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		TaylorModel tmTemp;

		hf.insert_ctrunc_normal(tmTemp, initialSet.tmvPre, polyRangeX0, step_exp_table, rangeDim+1, order, cutoff_threshold);

		tmTemp.remainder += remainder_x.get(i, 0);

		x.tms.push_back(tmTemp);
	}

	x.cutoff_normal(step_exp_table, cutoff_threshold);

	vector<Interval> newDomain = initialSet.domain;
	newDomain[0] = step_exp_table[1];

	Flowpipe firstFlowpipe(x, initialSet.tmv, newDomain);

	results.clear();
	results.push_back(initialSet);
	results.push_back(firstFlowpipe);

	if(bPrint)
	{
		printf("time = %f,\t", step);
		printf("step = %f,\t", step);
		printf("order = %d\n", order);
	}


	vector<Polynomial> P0;
	vector<Interval> R0;
	for(int i=0; i<rangeDim; ++i)
	{
		P0.push_back(firstFlowpipe.tmvPre.tms[i].expansion);
		R0.push_back(firstFlowpipe.tmvPre.tms[i].remainder);
	}

	vector<Interval> V0;
	for(int i=0; i<rangeDim; ++i)
	{
		Interval intTemp;
		B.get(intTemp, i, 0);
		V0.push_back(intTemp);
	}

	vector<Interval_matrix> remainder_template;
	vector<Interval> bounds_Vn;

	for(int i=0; i<rangeDim; ++i)
	{
		Interval_matrix imTemp(1, rangeDim);
		imTemp.set(intOne, 0, i);
		remainder_template.push_back(imTemp);

		bounds_Vn.push_back(intZero);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		Interval_matrix imTemp(1, rangeDim);
		imTemp.set(intMOne, 0, i);
		remainder_template.push_back(imTemp);

		bounds_Vn.push_back(intZero);
	}

	// resursively compute the remaining flowpipes
	int k = 1;
	for(double t=step+THRESHOLD_HIGH; t < time; ++k)
	{
		Flowpipe newFlowpipe;

		newFlowpipe.domain = newDomain;
		newFlowpipe.tmv = firstFlowpipe.tmv;

		vector<Polynomial> Pn;
		pow_expAr[k].linear_trans(Pn, P0);

		vector<Interval> bounds_Rn;

		for(int i=0; i<remainder_template.size(); ++i)
		{
			Interval_matrix vec_Rn = remainder_template[i] * pow_expAr[k];

			Interval intSum;
			for(int j=0; j<rangeDim; ++j)
			{
				Interval intTemp;
				vec_Rn.get(intTemp, 0, j);
				intSum += intTemp * R0[j];
			}

			Interval intSup;
			intSum.sup(intSup);
			bounds_Rn.push_back(intSup);

			Interval_matrix vec_Vn = remainder_template[i] * pow_expAr[k-1];
			intSum = intZero;
			for(int j=0; j<rangeDim; ++j)
			{
				Interval intTemp;
				vec_Vn.get(intTemp, 0, j);
				intSum += intTemp * V0[j];
			}

			intSum.sup(intSup);
			bounds_Vn[i] += intSup;
		}

		for(int i=0, j=rangeDim; i<rangeDim; ++i,++j)
		{
			Interval intUp = bounds_Rn[i] + bounds_Vn[i];
			Interval intLo = bounds_Rn[j] + bounds_Vn[j];

			Interval remainder( -intLo.inf() , intUp.sup() );

			Interval M;
			remainder.remove_midpoint(M);
			Polynomial polyMidpoint(M, rangeDim+1);

			TaylorModel tmTemp(Pn[i] + polyMidpoint, remainder);

			newFlowpipe.tmvPre.tms.push_back(tmTemp);
		}

		newFlowpipe.tmvPre.cutoff_normal(step_exp_table, cutoff_threshold);

		results.push_back(newFlowpipe);

		t += step;

		if(bPrint)
		{
			printf("time = %f,\t", t);
			printf("step = %f,\t", step);
			printf("order = %d\n", order);
		}
	}
}

void ContinuousSystem::reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
  mforce1("Reach Pic 1 <");
  minc();
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	//mforce("initial left", initialSet.tmvPre);
	//mforce("initial right", initialSet.tmv);
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;
  double t=THRESHOLD_HIGH;
	for(; t < time;)
	{
		int res = currentFlowpipe.advance_picard2(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, order, estimation, dummy_invariant, cutoff_threshold);
		if(res == 1)
		{ 
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;
			t += step;
      cerr << ".";
			if(bPrint && false)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
			/*
			mlog("pre", newFlowpipe.tmvPre);
			mlog("tmv", newFlowpipe.tmv);
			mdisable();
			*/
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough (time: %f).\n", t);
      settings->writer->info.push_back(sbuilder() << "reason: " << "The remainder estimation is not large enough.");
			break;
			
		}
	}
  settings->writer->info.push_back(sbuilder() << "int progress: " << t);
  //addFlowInfo(settings->writer->info);
  addMyInfo(settings->writer->info);
  
  cout << endl;
  mdec();
  mlog1("Reach Pic 1 >");
}

void ContinuousSystem::reach_picard(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_picard(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders
void ContinuousSystem::reach_picard(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_picard(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_picard(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_picard(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive orders and fixed step sizes
void ContinuousSystem::reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	int newOrder = order;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_picard(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_picard(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	vector<int> newOrders = orders;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		int res = currentFlowpipe.advance_picard(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}


// fixed step sizes and orders
void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde_centered.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde_centered.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde_centered.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde_centered.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, orders);

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, orders, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders
void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde_centered.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde_centered.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order);

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;

			double tDiffer = time - t;

			if(newStep > tDiffer)
			{
				newStep = tDiffer;
			}

			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde_centered.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde_centered.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, orders);

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive orders and fixed step sizes
void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	int newOrder = order;
	int currentMaxOrder = order;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde_centered.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde_centered.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansionHF;
	vector<Polynomial> taylorExpansionMF;
	vector<Polynomial> highestTerms;

	computeTaylorExpansion(taylorExpansionHF, taylorExpansionMF, highestTerms, polyODE, order);

	vector<vector<HornerForm> > expansions;
	expansions.push_back(taylorExpansionHF);

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, hfOde_centered, expansions[newOrder-order], precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;

				if(newOrder > currentMaxOrder)
				{
					for(int i=currentMaxOrder; i<newOrder; ++i)
					{
						vector<HornerForm> newTaylorExpansionHF;
						vector<Polynomial> newTaylorExpansionMF;

						increaseExpansionOrder(newTaylorExpansionHF, newTaylorExpansionMF, highestTerms, taylorExpansionMF, polyODE, i);

						expansions.push_back(newTaylorExpansionHF);
						taylorExpansionMF = newTaylorExpansionMF;
					}

					currentMaxOrder = newOrder;
				}
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<int> newOrders = orders;
	vector<int> localMaxOrders = orders;

	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde_centered.tms.size(); ++i)
	{
		polyODE.push_back(tmvOde_centered.tms[i].expansion);
	}

	vector<HornerForm> taylorExpansionHF;
	vector<Polynomial> taylorExpansionMF;
	vector<Polynomial> highestTerms;

	computeTaylorExpansion(taylorExpansionHF, taylorExpansionMF, highestTerms, polyODE, orders);

	vector<vector<HornerForm> > expansions;
	vector<HornerForm> emptySet;
	for(int i=0; i<taylorExpansionHF.size(); ++i)
	{
		expansions.push_back(emptySet);
		expansions[i].push_back(taylorExpansionHF[i]);
	}

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_low_degree(newFlowpipe, hfOde, hfOde_centered, taylorExpansionHF, precondition, step_exp_table, step_end_exp_table, newOrders, maxOrders, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
				{
					--newOrders[i];

					if(newOrders[i] > localMaxOrders[i])
					{
						for(int j=localMaxOrders[i]; j<newOrders[i]; ++j)
						{
							HornerForm newTaylorExpansionHF;
							Polynomial newTaylorExpansionMF;

							increaseExpansionOrder(newTaylorExpansionHF, newTaylorExpansionMF, highestTerms[i], taylorExpansionMF[i], polyODE, j);

							expansions[i].push_back(newTaylorExpansionHF);
							taylorExpansionMF[i] = newTaylorExpansionMF;
						}
					}

					localMaxOrders[i] = newOrders[i];

					taylorExpansionHF[i] = expansions[i][newOrders[i]-orders[i]];
				}
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}




// for high-degree ODEs
// fixed step sizes and orders

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders
void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive orders and fixed step sizes
void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	int newOrder = order;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	vector<int> newOrders = orders;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		int res = currentFlowpipe.advance_high_degree(newFlowpipe, hfOde, hfOde_centered, precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}













// for non-polynomial ODEs (using Taylor approximations)
// fixed step sizes and orders
void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, strOde_centered, precondition, step_exp_table, step_end_exp_table, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", order);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, strOde_centered, precondition, step_exp_table, step_end_exp_table, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

// adaptive step sizes and fixed orders
void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, strOde_centered, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, order, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("order = %d\n", order);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	double newStep = 0;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, strOde_centered, precondition, step_exp_table, step_end_exp_table, newStep, miniStep, orders, globalMaxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step_exp_table[1].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step_exp_table[1].sup());
				printf("orders:\t");
				int num = orders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), orders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), orders[num]);
			}

			newStep = step_exp_table[1].sup() * LAMBDA_UP;
			if(newStep > step - THRESHOLD_HIGH)
			{
				newStep = 0;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}


// adaptive orders and fixed step sizes
void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder);

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	int newOrder = order;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, strOde_centered, precondition, step_exp_table, step_end_exp_table, newOrder, maxOrder, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", newOrder);
			}

			if(newOrder > order)
			{
				--newOrder;
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

void ContinuousSystem::reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const
{
	vector<Interval> step_exp_table, step_end_exp_table;

	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*globalMaxOrder);

	vector<int> newOrders = orders;

	results.clear();
	results.push_back(initialSet);
	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	vector<PolynomialConstraint> dummy_invariant;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int localMaxOrder = newOrders[0];
		for(int i=1; i<newOrders.size(); ++i)
		{
			if(localMaxOrder < newOrders[i])
				localMaxOrder = newOrders[i];
		}

		int res = currentFlowpipe.advance_non_polynomial_taylor(newFlowpipe, strOde, strOde_centered, precondition, step_exp_table, step_end_exp_table, newOrders, localMaxOrder, maxOrders, estimation, dummy_invariant, cutoff_threshold);

		if(res == 1)
		{
			results.push_back(newFlowpipe);
			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("orders:\t");
				int num = newOrders.size()-1;
				for(int i=0; i<num; ++i)
				{
					printf("%s : %d, ", stateVarNames[i].c_str(), newOrders[i]);
				}
				printf("%s : %d\n", stateVarNames[num].c_str(), newOrders[num]);
			}

			for(int i=0; i<newOrders.size(); ++i)
			{
				if(newOrders[i] > orders[i])
					--newOrders[i];
			}
		}
		else
		{
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
}

ContinuousSystem & ContinuousSystem::operator = (const ContinuousSystem & system)
{
	if(this == &system)
		return *this;

	tmvOde				= system.tmvOde;
	hfOde				= system.hfOde;
	initialSet			= system.initialSet;
	strOde				= system.strOde;
	tmvOde_centered		= system.tmvOde_centered;
	hfOde_centered		= system.hfOde_centered;
	strOde_centered		= system.strOde_centered;

	return *this;
}




































// class ContinuousReachability
ContinuousReachability::ContinuousReachability()
{
  settings = new MySettings();
}

ContinuousReachability::~ContinuousReachability()
{
	outputAxes.clear();
	flowpipes.clear();
	orders.clear();
	maxOrders.clear();
	flowpipesCompo.clear();
	domains.clear();
	unsafeSet.clear();
	stateVarTab.clear();
	stateVarNames.clear();
	tmVarTab.clear();
	tmVarNames.clear();
	parTab.clear();
	parNames.clear();
}

void ContinuousReachability::dump(FILE *fp) const
{
	fprintf(fp,"state var ");
	for(int i=0; i<stateVarNames.size()-1; ++i)
	{
		fprintf(fp, "%s,", stateVarNames[i].c_str());
	}

	fprintf(fp, "%s\n\n", stateVarNames[stateVarNames.size()-1].c_str());

	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		switch(plotSetting)
		{
		case PLOT_INTERVAL:
			fprintf(fp, "gnuplot interval %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_OCTAGON:
			fprintf(fp, "gnuplot octagon %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_GRID:
			fprintf(fp, "gnuplot grid %d %s , %s\n\n", numSections, stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		}
		break;
	case PLOT_MATLAB:
		switch(plotSetting)
		{
		case PLOT_INTERVAL:
			fprintf(fp, "matlab interval %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_OCTAGON:
			fprintf(fp, "matlab octagon %s , %s\n\n", stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		case PLOT_GRID:
			fprintf(fp, "matlab grid %d %s , %s\n\n", numSections, stateVarNames[outputAxes[0]].c_str(), stateVarNames[outputAxes[1]].c_str());
			break;
		}
		break;
	}

	fprintf(fp, "cutoff %e\n\n", cutoff_threshold.sup());

	fprintf(fp, "output %s\n\n", outputFileName);

	if(bSafetyChecking)
	{
		// dump the unsafe set
		fprintf(fp, "unsafe set\n{\n");

		for(int i=0; i<unsafeSet.size(); ++i)
		{
			unsafeSet[i].dump(fp, stateVarNames);
		}

		fprintf(fp, "}\n\n");
	}

	fprintf(fp, "continuous flowpipes\n{\n");

	fprintf(fp, "tm var ");
	for(int i=1; i<tmVarNames.size()-1; ++i)
	{
		fprintf(fp, "%s,", tmVarNames[i].c_str());
	}
	fprintf(fp, "%s\n\n", tmVarNames[tmVarNames.size()-1].c_str());

	list<TaylorModelVec>::const_iterator fpIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	// every Taylor model flowpipe is enclosed by braces

	for(; fpIter != flowpipesCompo.end(); ++fpIter, ++doIter)
	{
		fprintf(fp, "{\n");
		fpIter->dump_interval(fp, stateVarNames, tmVarNames);

		//dump the domain
		for(int i=0; i<doIter->size(); ++i)
		{
			fprintf(fp, "%s in ", tmVarNames[i].c_str());
			(*doIter)[i].dump(fp);
			fprintf(fp, "\n");
		}

		fprintf(fp, "}\n\n");
	}

	fprintf(fp, "}\n");

/*
	// ======================= test begin ========================
	fpIter = flowpipesCompo.end();
	--fpIter;
	doIter = domains.end();
	--doIter;
	TaylorModelVec tmvTemp;
	vector<Interval> step_exp_table;
	Interval intStep((*doIter)[0].sup());
	construct_step_exp_table(step_exp_table, intStep, globalMaxOrder);
	fpIter->evaluate_t(tmvTemp, step_exp_table);
	vector<Interval> enclosure;
	tmvTemp.intEval(enclosure, *doIter);

	for(int i=0; i<enclosure.size(); ++i)
	{
		enclosure[i].dump(stdout);
		printf("\n");
	}
	// ======================= test end ========================
*/
}

void ContinuousReachability::run()
{
  mlog1("ContRun <");
  minc();
  
  compute_factorial_rec(globalMaxOrder+2);
	compute_power_4(globalMaxOrder+2);
	compute_double_factorial(2*globalMaxOrder+4);
  system.settings = settings;

/*
	if(integrationScheme != NONPOLY_TAYLOR)
	{
		bool bLinear = true;
		int rangeDim = system.tmvOde.tms.size();

		for(int i=0; i<rangeDim; ++i)
		{
			if(system.tmvOde.tms[i].expansion.degree() > 1)
			{
				bLinear = false;
				break;
			}
		}

		if(bLinear)
		{
			integrationScheme = LINEAR;
			system.reach_linear(flowpipes, step, time, orders[0], bPrint);
			return;
		}
	}
*/
  mlog1("run switch 1");
	switch(integrationScheme)
	{
	case ONLY_PICARD:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_picard(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_picard(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_picard(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_picard(flowpipes, step, miniStep, time, orders,            globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_picard(flowpipes, step,           time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_picard(flowpipes, step,           time, orders,            globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		}
		break;
	}

	case LOW_DEGREE:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_low_degree(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_low_degree(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_low_degree(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_low_degree(flowpipes, step, miniStep, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_low_degree(flowpipes, step, time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_low_degree(flowpipes, step, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		}
		break;
	}
  mlog1("run switch 2");
	case HIGH_DEGREE:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_high_degree(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_high_degree(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_high_degree(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_high_degree(flowpipes, step, miniStep, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_high_degree(flowpipes, step, time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_high_degree(flowpipes, step, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		}
		break;
	}

	case NONPOLY_TAYLOR:
	{
		switch(orderType)
		{
		case UNIFORM:
			if(bAdaptiveSteps)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, miniStep, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders[0], maxOrders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		case MULTI:
			if(bAdaptiveSteps)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, miniStep, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else if(bAdaptiveOrders)
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders, maxOrders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			else
			{
				system.reach_non_polynomial_taylor(flowpipes, step, time, orders, globalMaxOrder, precondition, estimation, bPrint, stateVarNames, cutoff_threshold);
			}
			break;
		}
		break;
	}

	case LINEAR:
	{
		system.reach_linear(flowpipes, step, time, orders[0], bPrint, cutoff_threshold);
		break;
	}
	}
  mdec();
  mlog1("ContRun >");
}

void ContinuousReachability::composition()
{
	flowpipesCompo.clear();
	domains.clear();

	if(integrationScheme == LINEAR)
	{
		list<Flowpipe>::const_iterator iter;

		for(iter = flowpipes.begin(); iter != flowpipes.end(); ++iter)
		{
			flowpipesCompo.push_back(iter->tmvPre);
			domains.push_back(iter->domain);
		}
	}
	else
	{
		vector<Interval> step_exp_table;
		Interval intStep;
		list<Flowpipe>::const_iterator iter;

		for(iter = flowpipes.begin(); iter != flowpipes.end(); ++iter)
		{
			if(step_exp_table.size() == 0 || intStep != iter->domain[0])
			{
				construct_step_exp_table(step_exp_table, iter->domain[0], globalMaxOrder);
				intStep = iter->domain[0];
			}

			TaylorModelVec tmvTemp;

			iter->composition_normal(tmvTemp, step_exp_table, cutoff_threshold);
			
			
			//mforce3(old3, "composed", tmvTemp);
			//exit(0);
			
			
      //pSerializer->add(iter->tmvPre, "comp_left");
      //pSerializer->add(iter->tmv, "comp_right");
			//pSerializer->add(tmvTemp, "composed");
			
			//mforce3(tt1, "comp_left", iter->tmvPre);
			//mforce3(tt2, "comp_right", iter->tmv);
			//mforce3(tt3, "comp_domain", iter->domain);
			
			/*
			
			vector<Interval> rightRange;
	    //iter->tmv.polyRange(rightRange, iter->domain);
    	iter->tmv.polyRangeNormal(rightRange, step_exp_table);
	
	    TaylorModelVec ret;
	    iter->tmvPre.insert_ctrunc_normal(ret, iter->tmv, rightRange,
	         step_exp_table, iter->domain.size(), 5, cutoff_threshold);
			//pSerializer->add(ret, "composed");
	    */
			
			

			flowpipesCompo.push_back(tmvTemp);
			//pSerializer->add(tmvTemp, "comp_comp");
			//mforce3(tt1, "pipe", tmvTemp.tms[0]);
			//mforce3(tt2, "manu", ret.tms[0]);
			domains.push_back(iter->domain);
			//mlog("domain", iter->domain);
		}
	}
}

int ContinuousReachability::safetyChecking() const
{
	if(unsafeSet.size() == 0)
	{
		return UNSAFE;	// since the whole state space is unsafe, the system is not safe
	}

	bool bDumpCounterexamples = true;

	int mkres = mkdir(counterexampleDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for counterexamples.\n");
		bDumpCounterexamples = false;
	}

	char filename_counterexamples[NAME_SIZE+10];
	FILE *fpDumpCounterexamples;

	if(bDumpCounterexamples)
	{
		sprintf(filename_counterexamples, "%s%s%s", counterexampleDir, outputFileName, str_counterexample_dumping_name_suffix);
		fpDumpCounterexamples = fopen(filename_counterexamples, "w");
	}

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	vector<Interval> step_exp_table;

	int rangeDim = tmvIter->tms.size();
	int domainDim = doIter->size();

	int result = SAFE;

	list<TaylorModelVec> unsafe_flowpipes;
	list<vector<Interval> > unsafe_flowpipe_domains;
	list<Interval> globalTimes;

	int maxOrder = 0;
	Interval globalTime;

	for(; tmvIter!=flowpipesCompo.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || step_exp_table[1] != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, (*doIter)[0], 2*maxOrder);
		}

		bool bsafe = false;

		vector<Interval> tmvPolyRange;
		tmvIter->polyRangeNormal(tmvPolyRange, step_exp_table);

		for(int i=0; i<unsafeSet.size(); ++i)
		{
			TaylorModel tmTemp;

			// interval evaluation on the constraint
			unsafeSet[i].hf.insert_normal(tmTemp, *tmvIter, tmvPolyRange, step_exp_table, domainDim, cutoff_threshold);

			Interval intTemp;
			tmTemp.intEvalNormal(intTemp, step_exp_table);

			if(intTemp > unsafeSet[i].B)
			{
				// no intersection with the unsafe set
				bsafe = true;
				break;
			}
			else
			{
				continue;
			}
		}

		if(!bsafe)
		{
			unsafe_flowpipes.push_back(*tmvIter);
			unsafe_flowpipe_domains.push_back(*doIter);
			globalTimes.push_back(globalTime);

			result = UNKNOWN;
		}

		globalTime += (*doIter)[0];
	}

	if(bDumpCounterexamples && unsafe_flowpipes.size() > 0)
	{
		dump_potential_counterexample(fpDumpCounterexamples, unsafe_flowpipes, unsafe_flowpipe_domains, globalTimes);
	}

	if(bDumpCounterexamples)
	{
		fclose(fpDumpCounterexamples);
	}

	return result;
}

unsigned long ContinuousReachability::numOfFlowpipes() const
{
	return (unsigned long)(flowpipes.size()) - 1;
}

void ContinuousReachability::dump_potential_counterexample(FILE *fp, const list<TaylorModelVec> & flowpipes, const list<vector<Interval> > & domains, const list<Interval> & globalTimes) const
{
	list<TaylorModelVec>::const_iterator fpIter = flowpipes.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	list<Interval>::const_iterator timeIter = globalTimes.begin();

	for(; fpIter!=flowpipes.end(); ++fpIter, ++doIter, ++timeIter)
	{
		fprintf(fp, "starting time %lf\n{\n", timeIter->sup());

		fpIter->dump_interval(fp, stateVarNames, tmVarNames);

		for(int i=0; i<doIter->size(); ++i)
		{
			fprintf(fp, "%s in ", tmVarNames[i].c_str());
			(*doIter)[i].dump(fp);
			fprintf(fp, "\n");
		}

		fprintf(fp, "}\n\n\n");
	}
}

void ContinuousReachability::plot_2D() const
{
	char filename[NAME_SIZE+10];

	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		sprintf(filename, "%s%s.plt", outputDir, outputFileName);
		break;
	case PLOT_MATLAB:
		sprintf(filename, "%s%s.m", outputDir, outputFileName);
		break;
	}

	FILE *fpPlotting = fopen(filename, "w");

	if(fpPlotting == NULL)
	{
		printf("Can not create the plotting file.\n");
		exit(1);
	}

	printf("Generating the plotting file...\n");
	switch(plotFormat)
	{
	case PLOT_GNUPLOT:
		plot_2D_GNUPLOT(fpPlotting);
		break;
	case PLOT_MATLAB:
		plot_2D_MATLAB(fpPlotting);
		break;
	}
	printf("Done.\n");

	fclose(fpPlotting);
}

void ContinuousReachability::plot_2D_GNUPLOT(FILE *fp) const
{
	switch(plotSetting)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_GNUPLOT(fp);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_GNUPLOT(fp);
		break;
	case PLOT_GRID:
		plot_2D_grid_GNUPLOT(fp);
		break;
	}
}

void ContinuousReachability::plot_2D_interval_GNUPLOT(FILE *fp) const
{
	fprintf(fp, "set term png\n");
  //fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.png", imageDir, outputFileName);
	//sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}

		vector<Interval> box;
		tmvIter->intEvalNormal(box, step_exp_table);

		Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

		// output the vertices
		fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
		fprintf(fp, "%lf %lf\n", X.sup(), Y.inf());
		fprintf(fp, "%lf %lf\n", X.sup(), Y.sup());
		fprintf(fp, "%lf %lf\n", X.inf(), Y.sup());
		fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
		fprintf(fp, "\n\n");
	}

	fprintf(fp, "e\n");
}

void ContinuousReachability::plot_2D_octagon_GNUPLOT(FILE *fp) const
{
	int x = outputAxes[0];
	int y = outputAxes[1];

	int rangeDim = stateVarNames.size();
	Matrix output_poly_temp(8, rangeDim);

	output_poly_temp.set(1, 0, x);
	output_poly_temp.set(1, 1, y);
	output_poly_temp.set(-1, 2, x);
	output_poly_temp.set(-1, 3, y);
	output_poly_temp.set(1/sqrt(2), 4, x);
	output_poly_temp.set(1/sqrt(2), 4, y);
	output_poly_temp.set(1/sqrt(2), 5, x);
	output_poly_temp.set(-1/sqrt(2), 5, y);
	output_poly_temp.set(-1/sqrt(2), 6, x);
	output_poly_temp.set(1/sqrt(2), 6, y);
	output_poly_temp.set(-1/sqrt(2), 7, x);
	output_poly_temp.set(-1/sqrt(2), 7, y);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	Matrix sortedTemplate(rows, cols);
	RowVector rowVec(cols);
	vector<RowVector> sortedRows;
	vector<RowVector>::iterator iterp, iterq;

	output_poly_temp.row(rowVec, 0);
	sortedRows.push_back(rowVec);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			double tmp1 = output_poly_temp.get(i,x) * iterp->get(y) - output_poly_temp.get(i,y) * iterp->get(x);
			double tmp2 = output_poly_temp.get(i,x) * iterq->get(y) - output_poly_temp.get(i,y) * iterq->get(x);

			if(tmp1 < 0 && tmp2 > 0)
			{
				output_poly_temp.row(rowVec, i);
				sortedRows.insert(iterq, rowVec);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			output_poly_temp.row(rowVec, i);
			sortedRows.push_back(rowVec);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate.set(iterp->get(j), i, j);
		}
	}

	ColVector b(rows);
	Polyhedron polyTemplate(sortedTemplate, b);

	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}


		templatePolyhedronNormal(polyTemplate, *tmvIter, step_exp_table);

		double f1, f2;

		vector<LinearConstraint>::iterator iterp, iterq;
		iterp = iterq = polyTemplate.constraints.begin();
		++iterq;

		for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
		{
			gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
			gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
			gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
			gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

			gsl_vector_set(d, 0, iterp->B.midpoint());
			gsl_vector_set(d, 1, iterq->B.midpoint());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			if(iterp == polyTemplate.constraints.begin())
			{
				f1 = v1;
				f2 = v2;
			}

			fprintf(fp, "%lf %lf\n", v1, v2);
		}

		iterp = polyTemplate.constraints.begin();
		--iterq;

		gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
		gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
		gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
		gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

		gsl_vector_set(d, 0, iterp->B.midpoint());
		gsl_vector_set(d, 1, iterq->B.midpoint());

		gsl_linalg_HH_solve(C, d, vertex);

		double v1 = gsl_vector_get(vertex, 0);
		double v2 = gsl_vector_get(vertex, 1);

		fprintf(fp, "%lf %lf\n", v1, v2);

		fprintf(fp, "%lf %lf\n", f1, f2);
		fprintf(fp, "\n\n");
	}

	fprintf(fp, "e\n");

	gsl_matrix_free(C);
	gsl_vector_free(d);
	gsl_vector_free(vertex);
}

void ContinuousReachability::plot_2D_grid_GNUPLOT(FILE *fp) const
{
	fprintf(fp, "set terminal postscript\n");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.eps", imageDir, outputFileName);
	fprintf(fp, "set output '%s'\n", filename);

	fprintf(fp, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(fp, "set autoscale\n");
	fprintf(fp, "unset label\n");
	fprintf(fp, "set xtic auto\n");
	fprintf(fp, "set ytic auto\n");
	fprintf(fp, "set xlabel \"%s\"\n", stateVarNames[outputAxes[0]].c_str());
	fprintf(fp, "set ylabel \"%s\"\n", stateVarNames[outputAxes[1]].c_str());
	fprintf(fp, "plot '-' notitle with lines ls 1\n");

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	int domainDim = doIter->size();

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		// decompose the domain
		list<vector<Interval> > grids;

		gridBox(grids, *doIter, numSections);

		// we only consider the output dimensions
		HornerForm hfOutputX;
		Interval remainderX;
		tmvIter->tms[outputAxes[0]].toHornerForm(hfOutputX, remainderX);

		HornerForm hfOutputY;
		Interval remainderY;
		tmvIter->tms[outputAxes[1]].toHornerForm(hfOutputY, remainderY);

		// evaluate the images from all of the grids
		list<vector<Interval> >::const_iterator gIter = grids.begin();
		for(; gIter!=grids.end(); ++gIter)
		{
			Interval X;
			hfOutputX.intEval(X, *gIter);
			X += remainderX;

			Interval Y;
			hfOutputY.intEval(Y, *gIter);
			Y += remainderY;

			// output the vertices
			fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
			fprintf(fp, "%lf %lf\n", X.sup(), Y.inf());
			fprintf(fp, "%lf %lf\n", X.sup(), Y.sup());
			fprintf(fp, "%lf %lf\n", X.inf(), Y.sup());
			fprintf(fp, "%lf %lf\n", X.inf(), Y.inf());
			fprintf(fp, "\n\n");
		}
	}

	fprintf(fp, "e\n");
}

void ContinuousReachability::plot_2D_MATLAB(FILE *fp) const
{
	switch(plotSetting)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_MATLAB(fp);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_MATLAB(fp);
		break;
	case PLOT_GRID:
		plot_2D_grid_MATLAB(fp);
		break;
	}
}

void ContinuousReachability::plot_2D_interval_MATLAB(FILE *fp) const
{
	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}

		vector<Interval> box;
		tmvIter->intEvalNormal(box, step_exp_table);

		Interval X(box[outputAxes[0]]), Y(box[outputAxes[1]]);

		// output all the vertices
		fprintf(fp,"plot( [%lf,%lf,%lf,%lf,%lf] , [%lf,%lf,%lf,%lf,%lf] , 'b');\nhold on;\nclear;\n",
				X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
	}
}

void ContinuousReachability::plot_2D_octagon_MATLAB(FILE *fp) const
{
	int x = outputAxes[0];
	int y = outputAxes[1];

	int rangeDim = stateVarNames.size();
	Matrix output_poly_temp(8, rangeDim);

	output_poly_temp.set(1, 0, x);
	output_poly_temp.set(1, 1, y);
	output_poly_temp.set(-1, 2, x);
	output_poly_temp.set(-1, 3, y);
	output_poly_temp.set(1/sqrt(2), 4, x);
	output_poly_temp.set(1/sqrt(2), 4, y);
	output_poly_temp.set(1/sqrt(2), 5, x);
	output_poly_temp.set(-1/sqrt(2), 5, y);
	output_poly_temp.set(-1/sqrt(2), 6, x);
	output_poly_temp.set(1/sqrt(2), 6, y);
	output_poly_temp.set(-1/sqrt(2), 7, x);
	output_poly_temp.set(-1/sqrt(2), 7, y);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	Matrix sortedTemplate(rows, cols);
	RowVector rowVec(cols);
	vector<RowVector> sortedRows;
	vector<RowVector>::iterator iterp, iterq;

	output_poly_temp.row(rowVec, 0);
	sortedRows.push_back(rowVec);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			double tmp1 = output_poly_temp.get(i,x) * iterp->get(y) - output_poly_temp.get(i,y) * iterp->get(x);
			double tmp2 = output_poly_temp.get(i,x) * iterq->get(y) - output_poly_temp.get(i,y) * iterq->get(x);

			if(tmp1 < 0 && tmp2 > 0)
			{
				output_poly_temp.row(rowVec, i);
				sortedRows.insert(iterq, rowVec);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			output_poly_temp.row(rowVec, i);
			sortedRows.push_back(rowVec);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate.set(iterp->get(j), i, j);
		}
	}

	ColVector b(rows);
	Polyhedron polyTemplate(sortedTemplate, b);

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	vector<Interval> step_exp_table, step_end_exp_table;
	Interval intStep;

	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();

	int maxOrder = 0;

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		int tmp = maxOrder;
		for(int i=0; i<tmvIter->tms.size(); ++i)
		{
			int order = tmvIter->tms[i].expansion.degree();
			if(maxOrder < order)
			{
				maxOrder = order;
			}
		}

		if(step_exp_table.size() == 0 || intStep != (*doIter)[0] || maxOrder > tmp)
		{
			construct_step_exp_table(step_exp_table, step_end_exp_table, (*doIter)[0].sup(), maxOrder);
			intStep = (*doIter)[0];
		}


		templatePolyhedronNormal(polyTemplate, *tmvIter, step_exp_table);

		double f1, f2;

		vector<LinearConstraint>::iterator iterp, iterq;
		iterp = iterq = polyTemplate.constraints.begin();
		++iterq;

		vector<double> vertices_x, vertices_y;

		for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
		{
			gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
			gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
			gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
			gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

			gsl_vector_set(d, 0, iterp->B.midpoint());
			gsl_vector_set(d, 1, iterq->B.midpoint());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			if(iterp == polyTemplate.constraints.begin())
			{
				f1 = v1;
				f2 = v2;
			}

			vertices_x.push_back(v1);
			vertices_y.push_back(v2);
		}

		iterp = polyTemplate.constraints.begin();
		--iterq;

		gsl_matrix_set(C, 0, 0, iterp->A[x].midpoint());
		gsl_matrix_set(C, 0, 1, iterp->A[y].midpoint());
		gsl_matrix_set(C, 1, 0, iterq->A[x].midpoint());
		gsl_matrix_set(C, 1, 1, iterq->A[y].midpoint());

		gsl_vector_set(d, 0, iterp->B.midpoint());
		gsl_vector_set(d, 1, iterq->B.midpoint());

		gsl_linalg_HH_solve(C, d, vertex);

		double v1 = gsl_vector_get(vertex, 0);
		double v2 = gsl_vector_get(vertex, 1);

		vertices_x.push_back(v1);
		vertices_y.push_back(v2);
		vertices_x.push_back(f1);
		vertices_y.push_back(f2);

		fprintf(fp, "plot( ");

		fprintf(fp, "[ ");
		for(int i=0; i<vertices_x.size()-1; ++i)
		{
			fprintf(fp, "%lf , ", vertices_x[i]);
		}
		fprintf(fp, "%lf ] , ", vertices_x.back());

		fprintf(fp, "[ ");
		for(int i=0; i<vertices_y.size()-1; ++i)
		{
			fprintf(fp, "%lf , ", vertices_y[i]);
		}
		fprintf(fp, "%lf ] , ", vertices_y.back());

		fprintf(fp, "'b');\nhold on;\nclear;\n");
	}

	gsl_matrix_free(C);
	gsl_vector_free(d);
	gsl_vector_free(vertex);
}

void ContinuousReachability::plot_2D_grid_MATLAB(FILE *fp) const
{
	list<TaylorModelVec>::const_iterator tmvIter = flowpipesCompo.begin();
	list<vector<Interval> >::const_iterator doIter = domains.begin();
	int domainDim = doIter->size();

	for(; tmvIter != flowpipesCompo.end() && doIter != domains.end(); ++tmvIter, ++doIter)
	{
		// decompose the domain
		list<vector<Interval> > grids;

		gridBox(grids, *doIter, numSections);

		// we only consider the output dimensions
		HornerForm hfOutputX;
		Interval remainderX;
		tmvIter->tms[outputAxes[0]].toHornerForm(hfOutputX, remainderX);

		HornerForm hfOutputY;
		Interval remainderY;
		tmvIter->tms[outputAxes[1]].toHornerForm(hfOutputY, remainderY);

		// evaluate the images from all of the grids
		list<vector<Interval> >::const_iterator gIter = grids.begin();
		for(; gIter!=grids.end(); ++gIter)
		{
			Interval X;
			hfOutputX.intEval(X, *gIter);
			X += remainderX;

			Interval Y;
			hfOutputY.intEval(Y, *gIter);
			Y += remainderY;

			// output the vertices
			fprintf(fp,"plot( [%lf,%lf,%lf,%lf,%lf] , [%lf,%lf,%lf,%lf,%lf] , 'b');\nhold on;\nclear;\n",
					X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
		}
	}
}

bool ContinuousReachability::declareStateVar(const string & vName)
{
	map<string,int>::const_iterator iter;

	if((iter = stateVarTab.find(vName)) == stateVarTab.end())
	{
		stateVarTab[vName] = stateVarNames.size();
		stateVarNames.push_back(vName);
		//settings->varNames.push_back(vName);
		return true;
	}
	else
	{
		return false;
	}
}

int ContinuousReachability::getIDForStateVar(const string & vName) const
{
	map<string,int>::const_iterator iter;
	if((iter = stateVarTab.find(vName)) == stateVarTab.end())
	{
		return -1;
	}

	return iter -> second;
}

bool ContinuousReachability::getStateVarName(string & vName, int id) const
{
	if(id >= 0 && id < stateVarNames.size())
	{
		vName = stateVarNames[id];
		return true;
	}
	else
	{
		return false;
	}
}



bool ContinuousReachability::declareTMVar(const string & vName)
{
	map<string,int>::const_iterator iter;

	if((iter = tmVarTab.find(vName)) == tmVarTab.end())
	{
		tmVarTab[vName] = tmVarNames.size();
		tmVarNames.push_back(vName);
		return true;
	}
	else
	{
		return false;
	}
}

int ContinuousReachability::getIDForTMVar(const string & vName) const
{
	map<string,int>::const_iterator iter;
	if((iter = tmVarTab.find(vName)) == tmVarTab.end())
	{
		return -1;
	}

	return iter->second;
}

bool ContinuousReachability::getTMVarName(string & vName, int id) const
{
	if(id >= 0 && id < tmVarNames.size())
	{
		vName = tmVarNames[id];
		return true;
	}
	else
	{
		return false;
	}
}




bool ContinuousReachability::declarePar(const string & pName, const Interval & range)
{
	map<string,int>::const_iterator iter;

	if((iter = parTab.find(pName)) == parTab.end())
	{
		parTab[pName] = parNames.size();
		parNames.push_back(pName);
		parRanges.push_back(range);
		return true;
	}
	else
	{
		return false;
	}
}

int ContinuousReachability::getIDForPar(const string & pName) const
{
	map<string,int>::const_iterator iter;
	if((iter = parTab.find(pName)) == parTab.end())
	{
		return -1;
	}

	return iter -> second;
}

bool ContinuousReachability::getParName(string & pName, int id) const
{
	if(id >= 0 && id < parNames.size())
	{
		pName = parNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

bool ContinuousReachability::getRangeForPar(Interval & range, const string & pName) const
{
	int id = getIDForPar(pName);

	if(id == -1)
	{
		return false;
	}
	else
	{
		range = parRanges[id];
		return true;
	}
}












void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const int order, const Interval & cutoff_threshold)
{
	vector<Interval> intVecZero;
	Interval intZero, intOne(1,1);

	intVecZero.push_back(intOne);
	intVecZero.push_back(intZero);

	// we compute the Taylor expansion (without the 0-order term)
	TaylorModelVec taylorExpansion;
	first_order_deriv.evaluate_t(taylorExpansion, intVecZero);
//	taylorExpansion.nctrunc(order - 1);
	taylorExpansion.mul_assign(0, 1);

	TaylorModelVec tmvLieDeriv_n = first_order_deriv;

	for(int i=2; i<=order; ++i)
	{
		TaylorModelVec tmvTemp;
		tmvLieDeriv_n.LieDerivative_no_remainder(tmvTemp, ode, order - i, cutoff_threshold);

		TaylorModelVec tmvTemp2;
		tmvTemp.evaluate_t(tmvTemp2, intVecZero);
		tmvTemp2.mul_assign(factorial_rec[i]);
		tmvTemp2.mul_assign(0,i);			// multiplied by t^i
//		tmvTemp2.nctrunc(order);

		taylorExpansion.add_assign(tmvTemp2);

		tmvLieDeriv_n = tmvTemp;
	}
/*
	for(int i=0; i<taylorExpansion.tms.size(); ++i)
	{
		taylorExpansion.tms[i].cutoff(cutoff_threshold);
	}
*/
	result = taylorExpansion;
}

void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const vector<int> & orders, const Interval & cutoff_threshold)
{
	int rangeDim = ode.tms.size();
	vector<Interval> intVecZero;
	Interval intZero, intOne(1,1);

	intVecZero.push_back(intOne);
	intVecZero.push_back(intZero);

	// we compute the Taylor expansion (without the 0-order term)
	TaylorModelVec taylorExpansion;
	first_order_deriv.evaluate_t(taylorExpansion, intVecZero);
/*
	for(int i=0; i<rangeDim; ++i)
	{
		taylorExpansion.tms[i].nctrunc(orders[i] - 1);
	}
*/
	taylorExpansion.mul_assign(0, 1);

	TaylorModelVec tmvLieDeriv_n = first_order_deriv;

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=2; j<=orders[i]; ++j)
		{
			TaylorModel tmTemp;
			tmvLieDeriv_n.tms[i].LieDerivative_no_remainder(tmTemp, ode, orders[i] - j, cutoff_threshold);

			TaylorModel tmTemp2;
			tmTemp.evaluate_t(tmTemp2, intVecZero);
			tmTemp2.mul_assign(factorial_rec[j]);
			tmTemp2.mul_assign(0,j);
//			tmTemp2.nctrunc(orders[i]);

			taylorExpansion.tms[i].add_assign(tmTemp2);

			tmvLieDeriv_n.tms[i] = tmTemp;
		}
	}
/*
	for(int i=0; i<taylorExpansion.tms.size(); ++i)
	{
		taylorExpansion.tms[i].cutoff(cutoff_threshold);
	}
*/
	result = taylorExpansion;
}

void construct_step_exp_table(vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const int order)
{
	step_exp_table.clear();
	step_end_exp_table.clear();

	Interval intProd(1), intStep(0,step);

	for(int i=0; i<=order; ++i)
	{
		step_exp_table.push_back(intProd);

		Interval intTend(intProd.sup());
		step_end_exp_table.push_back(intTend);

		intProd *= intStep;
	}
}

void construct_step_exp_table(vector<Interval> & step_exp_table, const Interval & step, const int order)
{
	step_exp_table.clear();

	Interval intProd(1);
	step_exp_table.push_back(intProd);

	for(int i=1; i<=order; ++i)
	{
		intProd *= step;
		step_exp_table.push_back(intProd);
	}
}

void preconditionQR(Matrix & result, const TaylorModelVec & x0, const int rangeDim, const int domainDim)
{
	Interval intZero;
	vector<vector<Interval> > intCoefficients;

	vector<Interval> intVecTemp;
	for(int i=0; i<domainDim; ++i)
	{
		intVecTemp.push_back(intZero);
	}

	for(int i=0; i<rangeDim; ++i)
	{
		intCoefficients.push_back(intVecTemp);
	}

	x0.linearCoefficients(intCoefficients);
	Matrix matCoefficients(rangeDim, rangeDim);

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=rangeDim; ++j)
		{
			matCoefficients.set(intCoefficients[i][j].midpoint(), i, j-1);
		}
	}
	
	/*
	mlog1("linear coefficients <");
	mlog(matCoefficients);
	mlog1("linear coefficients >");
	//*/

	matCoefficients.sortColumns();
	//mlog(matCoefficients);
	matCoefficients.QRfactor(result);
	//mlog(result);
}

Interval rho(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & domain)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		tmv.tms[i].mul(tmTemp, l[i]);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEval(intRange, domain);

	Interval S;
	intRange.sup(S);

	return S;
}

Interval rhoNormal(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & step_exp_table)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		tmv.tms[i].mul(tmTemp, l[i]);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEvalNormal(intRange, step_exp_table);

	Interval S;
	intRange.sup(S);

	return S;
}

Interval rho(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & domain)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		Interval intTemp(l.get(i));
		tmv.tms[i].mul(tmTemp, intTemp);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEval(intRange, domain);

	Interval S;
	intRange.sup(S);

	return S;
}

Interval rhoNormal(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & step_exp_table)
{
	int d = l.size();
	TaylorModel tmObj;

	for(int i=0; i<d; ++i)
	{
		TaylorModel tmTemp;
		Interval intTemp(l.get(i));
		tmv.tms[i].mul(tmTemp, intTemp);
		tmObj.add_assign(tmTemp);
	}

	Interval intRange;
	tmObj.intEvalNormal(intRange, step_exp_table);

	Interval S;
	intRange.sup(S);

	return S;
}

void templatePolyhedron(Polyhedron & result, const TaylorModelVec & tmv, const vector<Interval> & domain)
{
	for(int i=0; i<result.constraints.size(); ++i)
	{
		result.constraints[i].B = rho(tmv, result.constraints[i].A, domain);
	}
}

void templatePolyhedronNormal(Polyhedron & result, const TaylorModelVec & tmv, vector<Interval> & step_exp_table)
{
	for(int i=0; i<result.constraints.size(); ++i)
	{
		result.constraints[i].B = rhoNormal(tmv, result.constraints[i].A, step_exp_table);
	}
}

/*
int intersection_check_interval_arithmetic(const vector<PolynomialConstraint> & pcs, const vector<HornerForm> & objFuncs, const vector<Interval> & remainders, const vector<Interval> & domain, vector<bool> & bNeeded)
{
	int counter = 0;
	bNeeded.clear();

	for(int i=0; i<pcs.size(); ++i)
	{
		Interval intTemp;
		objFuncs[i].intEval(intTemp, domain);
		intTemp += remainders[i];

		if(intTemp > pcs[i].B)
		{
			// no intersection
			return -1;
		}
		else if(intTemp.smallereq(pcs[i].B))
		{
			// the flowpipe is entirely contained in the guard, domain contraction is not needed
			bNeeded.push_back(false);
			++counter;
		}
		else
		{
			bNeeded.push_back(true);
		}
	}

	return counter;
}

bool boundary_intersected_collection(const vector<PolynomialConstraint> & pcs, const vector<HornerForm> & objFuncs, const vector<Interval> & remainders, const vector<Interval> & domain, vector<bool> & boundary_intersected)
{
	boundary_intersected.clear();

	for(int i=0; i<pcs.size(); ++i)
	{
		Interval intTemp;
		objFuncs[i].intEval(intTemp, domain);
		intTemp += remainders[i];

		if(intTemp > pcs[i].B)
		{
			// no intersection
			return false;
		}
		else if(intTemp < pcs[i].B)
		{
			// do not intersect the boundary
			boundary_intersected.push_back(false);
		}
		else
		{
			boundary_intersected.push_back(true);
		}
	}

	return true;
}
*/

int contract_interval_arithmetic(TaylorModelVec & flowpipe, vector<Interval> & domain, const vector<PolynomialConstraint> & pcs, vector<bool> & boundary_intersected, const Interval & cutoff_threshold)
{
	int rangeDim = flowpipe.tms.size();
	int domainDim = domain.size();

	// contract the remainder firstly
	bool bvalid = true;
	bool bcontinue = true;
	Interval W;
	Interval intZero;

	vector<bool> bNeeded;
	int counter = 0;

	for(int i=0; i<pcs.size(); ++i)
	{
		bNeeded.push_back(true);
	}

	boundary_intersected.clear();

	vector<Interval> remainders;
	for(int i=0; i<rangeDim; ++i)
	{
		remainders.push_back(flowpipe.tms[i].remainder);
	}

	vector<Interval> flowpipePolyRange;
	flowpipe.polyRange(flowpipePolyRange, domain);

	vector<Interval> intVecTemp = flowpipePolyRange;
	intVecTemp.insert(intVecTemp.begin(), intZero);

	// 1: we check the intersection with every constraint
	for(int i=0; i<rangeDim; ++i)
	{
		intVecTemp[i+1] = flowpipePolyRange[i] + remainders[i];
	}

	for(int i=0; i<pcs.size(); ++i)
	{
		Interval intTemp;

		pcs[i].hf.intEval(intTemp, intVecTemp);

		if(intTemp > pcs[i].B)
		{
			// no intersection on the left half
			bvalid = false;
			break;
		}
		else if(intTemp.smallereq(pcs[i].B))
		{
			// do not need to apply domain contraction w.r.t. the current constraint
			boundary_intersected.push_back(false);
			bNeeded[i] = false;
			++counter;
		}
		else
		{
			boundary_intersected.push_back(true);
			bNeeded[i] = true;
			continue;
		}
	}

	if(!bvalid)
	{
		boundary_intersected.clear();
		return -1;	// no intersection is detected
	}
	else if(counter == pcs.size())
	{
		return 0;	// no need to do contraction
	}


	// 2: remainder contraction
	for(; bcontinue; )
	{
		vector<Interval> oldRemainders = remainders;

		for(int i=0; i<rangeDim; ++i)
		{
			Interval newInt = remainders[i];
			vector<bool> localNeeded = bNeeded;
			int localCounter = counter;

			for(int k=0; k<rangeDim; ++k)
			{
				if(k != i)
				{
					intVecTemp[k+1] = flowpipePolyRange[k] + remainders[k];
				}
				else
				{
					intVecTemp[k+1] = flowpipePolyRange[k];
				}
			}

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<pcs.size(); ++j)
				{
					if(localNeeded[j])
					{
						Interval intTemp;
						vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i+1] += intLeft;

						pcs[j].hf.intEval(intTemp, newIntVecTemp);

						if(intTemp > pcs[j].B)
						{
							// no intersection on the left half
							newInt = intRight;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(pcs[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							newInt.width(W);

							continue;
						}
					}
				}

				if(localCounter == pcs.size())
				{
					break;
				}
			}

			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			remainders[i].setInf(Inf);

			newInt = remainders[i];
			newInt.width(W);

			localNeeded = bNeeded;
			localCounter = counter;

			// search an approximation for the upper bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<pcs.size(); ++j)
				{
					if(localNeeded[j])
					{
						Interval intTemp;
						vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i+1] += intRight;

						pcs[j].hf.intEval(intTemp, newIntVecTemp);

						if(intTemp > pcs[j].B)
						{
							// no intersection on the right half
							newInt = intLeft;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(pcs[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							newInt.width(W);
							continue;
						}
					}
				}

				if(localCounter == pcs.size())
				{
					break;
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			remainders[i].setSup(Sup);	// set the upper bound

			if(!remainders[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<rangeDim; ++i)
		{
			if(oldRemainders[i].widthRatio(remainders[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}
	}

	if(!bvalid)
	{
		boundary_intersected.clear();
		return -1;	// no intersection is detected
	}

	for(int i=0; i<rangeDim; ++i)
	{
		flowpipe.tms[i].remainder = remainders[i];
	}


	// the Horner forms of p(T(x))
	vector<HornerForm> objHF;

	vector<Interval> eval_remainders;

	for(int i=0; i<pcs.size(); ++i)
	{
		TaylorModel tmTemp;
		pcs[i].hf.insert(tmTemp, flowpipe, flowpipePolyRange, domain, cutoff_threshold);

		HornerForm hf;
		Interval remainder;
		tmTemp.toHornerForm(hf, remainder);
		objHF.push_back(hf);
		eval_remainders.push_back(remainder);
	}

	Interval intTime = domain[0];

	bvalid = true;
	bcontinue = true;

	// 3: domain contraction
	for(; bcontinue; )
	{
		vector<Interval> oldDomain = domain;

		// contract the domain
		for(int i=0; i<domainDim; ++i)
		{
			Interval newInt = domain[i];
			vector<bool> localNeeded = bNeeded;
			int localCounter = counter;

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<pcs.size(); ++j)
				{
					if(localNeeded[j])
					{
						vector<Interval> newDomain = domain;
						newDomain[i] = intLeft;

						Interval intTemp;
						objHF[j].intEval(intTemp, newDomain);
						intTemp += eval_remainders[j];

						if(intTemp > pcs[j].B)
						{
							// no intersection on the left half
							newInt = intRight;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(pcs[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							newInt.width(W);

							continue;
						}
					}
				}

				if(localCounter == pcs.size())
				{
					break;
				}
			}

			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			domain[i].setInf(Inf);

			newInt = domain[i];

			localNeeded = bNeeded;
			localCounter = counter;

			newInt.width(W);

			// search an approximation for the upper bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<pcs.size(); ++j)
				{
					if(localNeeded[j])
					{
						vector<Interval> newDomain = domain;
						newDomain[i] = intRight;

						Interval intTemp;
						objHF[j].intEval(intTemp, newDomain);
						intTemp += eval_remainders[j];

						if(intTemp > pcs[j].B)
						{
							// no intersection on the right half
							newInt = intLeft;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(pcs[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							newInt.width(W);
							continue;
						}
					}
				}

				if(localCounter == pcs.size())
				{
					break;
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			domain[i].setSup(Sup);	// set the upper bound

			if(!domain[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<domainDim; ++i)
		{
			if(oldDomain[i].widthRatio(domain[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}

		if(bcontinue)
		{
			objHF.clear();
			eval_remainders.clear();

			flowpipe.polyRange(flowpipePolyRange, domain);

			for(int i=0; i<pcs.size(); ++i)
			{
				TaylorModel tmTemp;

				pcs[i].hf.insert(tmTemp, flowpipe, flowpipePolyRange, domain, cutoff_threshold);

				HornerForm hf;
				Interval remainder;
				tmTemp.toHornerForm(hf, remainder);
				objHF.push_back(hf);
				eval_remainders.push_back(remainder);
			}
		}
	}

	if(!bvalid)
	{
		boundary_intersected.clear();
		return -1;
	}

	if(intTime != domain[0])
	{
		return 2;
	}
	else
	{
		return 1;
	}
}


int contract_interval_arithmetic(TaylorModelVec & flowpipe, vector<Interval> & domain, const Polyhedron & inv, vector<bool> & boundary_intersected)
{
	int rangeDim = flowpipe.tms.size();
	int domainDim = domain.size();

	boundary_intersected.clear();

	// contract the remainder firstly
	bool bvalid = true;
	bool bcontinue = true;
	Interval W;
	Interval intZero;

	vector<bool> bNeeded;
	int counter = 0;

	for(int i=0; i<inv.constraints.size(); ++i)
	{
		bNeeded.push_back(true);
	}

	boundary_intersected.clear();

	vector<Interval> remainders;
	for(int i=0; i<rangeDim; ++i)
	{
		remainders.push_back(flowpipe.tms[i].remainder);
	}

	vector<Interval> flowpipePolyRange;
	flowpipe.polyRange(flowpipePolyRange, domain);

	vector<Interval> intVecTemp = flowpipePolyRange;

	// 1: we check the intersection with every constraint
	for(int i=0; i<rangeDim; ++i)
	{
		intVecTemp[i] = flowpipePolyRange[i] + remainders[i];
	}

	for(int i=0; i<inv.constraints.size(); ++i)
	{
		Interval intTemp;

		for(int j=0; j<inv.constraints[i].A.size(); ++j)
		{
			intTemp += inv.constraints[i].A[j] * intVecTemp[j];
		}

		if(intTemp > inv.constraints[i].B)
		{
			// no intersection on the left half
			bvalid = false;
			break;
		}
		else if(intTemp.smallereq(inv.constraints[i].B))
		{
			// do not need to apply domain contraction w.r.t. the current constraint
			boundary_intersected.push_back(false);
			bNeeded[i] = false;
			++counter;
		}
		else
		{
			boundary_intersected.push_back(true);
			bNeeded[i] = true;
			continue;
		}
	}

	if(!bvalid)
	{
		boundary_intersected.clear();
		return -1;	// no intersection is detected
	}
	else if(counter == inv.constraints.size())
	{
		return 0;	// no need to do contraction
	}


	// 2: remainder contraction
	for(; bcontinue; )
	{
		vector<Interval> oldRemainders = remainders;

		for(int i=0; i<rangeDim; ++i)
		{
			Interval newInt = remainders[i];
			vector<bool> localNeeded = bNeeded;
			int localCounter = counter;

			for(int k=0; k<rangeDim; ++k)
			{
				if(k != i)
				{
					intVecTemp[k] = flowpipePolyRange[k] + remainders[k];
				}
				else
				{
					intVecTemp[k] = flowpipePolyRange[k];
				}
			}

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<inv.constraints.size(); ++j)
				{
					if(localNeeded[j])
					{
						Interval intTemp;
						vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i] += intLeft;

						for(int k=0; k<inv.constraints[j].A.size(); ++k)
						{
							intTemp += inv.constraints[j].A[k] * newIntVecTemp[k];
						}

						if(intTemp > inv.constraints[j].B)
						{
							// no intersection on the left half
							newInt = intRight;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(inv.constraints[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							newInt.width(W);

							continue;
						}
					}
				}

				if(localCounter == inv.constraints.size())
				{
					break;
				}
			}

			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			remainders[i].setInf(Inf);

			newInt = remainders[i];
			newInt.width(W);

			localNeeded = bNeeded;
			localCounter = counter;

			// search an approximation for the upper bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<inv.constraints.size(); ++j)
				{
					if(localNeeded[j])
					{
						Interval intTemp;
						vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i] += intRight;

						for(int k=0; k<inv.constraints[j].A.size(); ++k)
						{
							intTemp += inv.constraints[j].A[k] * newIntVecTemp[k];
						}

						if(intTemp > inv.constraints[j].B)
						{
							// no intersection on the right half
							newInt = intLeft;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(inv.constraints[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							newInt.width(W);
							continue;
						}
					}
				}

				if(localCounter == inv.constraints.size())
				{
					break;
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			remainders[i].setSup(Sup);	// set the upper bound

			if(!remainders[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<rangeDim; ++i)
		{
			if(oldRemainders[i].widthRatio(remainders[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}
	}

	if(!bvalid)
	{
		boundary_intersected.clear();
		return -1;	// no intersection is detected
	}

	for(int i=0; i<rangeDim; ++i)
	{
		flowpipe.tms[i].remainder = remainders[i];
	}


	// the Horner forms of p(T(x))
	vector<HornerForm> objHF;
	vector<Interval> eval_remainders;

	for(int i=0; i<inv.constraints.size(); ++i)
	{
		TaylorModel tmTemp;

		for(int j=0; j<inv.constraints[i].A.size(); ++j)
		{
			TaylorModel tmTemp2;
			flowpipe.tms[j].mul(tmTemp2, inv.constraints[i].A[j]);
			tmTemp.add_assign(tmTemp2);
		}

		HornerForm hf;
		Interval remainder;
		tmTemp.toHornerForm(hf, remainder);

		objHF.push_back(hf);
		eval_remainders.push_back(remainder);
	}

	Interval intTime = domain[0];

	bvalid = true;
	bcontinue = true;

	// 3: domain contraction

	for(; bcontinue; )
	{
		vector<Interval> oldDomain = domain;

		// contract the domain
		for(int i=0; i<domainDim; ++i)
		{
			Interval newInt = domain[i];
			vector<bool> localNeeded = bNeeded;
			int localCounter = counter;

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<inv.constraints.size(); ++j)
				{
					if(localNeeded[j])
					{
						vector<Interval> newDomain = domain;
						newDomain[i] = intLeft;

						Interval intTemp;
						objHF[j].intEval(intTemp, newDomain);
						intTemp += eval_remainders[j];

						if(intTemp > inv.constraints[j].B)
						{
							// no intersection on the left half
							newInt = intRight;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(inv.constraints[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							newInt.width(W);

							continue;
						}
					}
				}

				if(localCounter == inv.constraints.size())
				{
					break;
				}
			}

			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			domain[i].setInf(Inf);

			newInt = domain[i];

			localNeeded = bNeeded;
			localCounter = counter;

			newInt.width(W);

			// search an approximation for the upper bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<inv.constraints.size(); ++j)
				{
					if(localNeeded[j])
					{
						vector<Interval> newDomain = domain;
						newDomain[i] = intRight;

						Interval intTemp;
						objHF[j].intEval(intTemp, newDomain);
						intTemp += eval_remainders[j];

						if(intTemp > inv.constraints[j].B)
						{
							// no intersection on the right half
							newInt = intLeft;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(inv.constraints[j].B))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							newInt.width(W);
							continue;
						}
					}
				}

				if(localCounter == inv.constraints.size())
				{
					break;
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			domain[i].setSup(Sup);	// set the upper bound

			if(!domain[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<domainDim; ++i)
		{
			if(oldDomain[i].widthRatio(domain[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}
	}

	if(!bvalid)
	{
		boundary_intersected.clear();
		return -1;
	}

	if(intTime != domain[0])
	{
		return 2;
	}
	else
	{
		return 1;
	}
}

int contract_remainder(const vector<Interval> & polyRange, vector<Interval> & remainders, const vector<HornerForm> & hfs, const vector<Interval> & b)
{
	bool bvalid = true;
	bool bcontinue = true;
	Interval W;
	Interval intZero;

	vector<bool> bNeeded;
	int counter = 0;

	for(int i=0; i<hfs.size(); ++i)
	{
		bNeeded.push_back(true);
	}

	int rangeDim = polyRange.size();

	vector<Interval> intVecTemp = polyRange;
	intVecTemp.insert(intVecTemp.begin(), intZero);		// range of the dummy time variable

	// 1: we check the intersection with every constraint
	for(int i=0; i<rangeDim; ++i)
	{
		intVecTemp[i+1] = polyRange[i] + remainders[i];
	}

	for(int i=0; i<hfs.size(); ++i)
	{
		Interval intTemp;

		hfs[i].intEval(intTemp, intVecTemp);

		if(intTemp > b[i])
		{
			// no intersection on the left half
			bvalid = false;
			break;
		}
		else if(intTemp.smallereq(b[i]))
		{
			// do not need to apply domain contraction w.r.t. the current constraint
			bNeeded[i] = false;
			++counter;
		}
		else
		{
			bNeeded[i] = true;
			continue;
		}
	}

	if(!bvalid)
	{
		return -1;	// no intersection is detected
	}
	else if(counter == hfs.size())
	{
		return 0;	// no need to do contraction
	}

	// 2: contract the remainder
	for(; bcontinue; )
	{
		vector<Interval> oldRemainders = remainders;

		for(int i=0; i<rangeDim; ++i)
		{
			Interval newInt = remainders[i];
			vector<bool> localNeeded = bNeeded;
			int localCounter = counter;

			for(int k=0; k<rangeDim; ++k)
			{
				if(k != i)
				{
					intVecTemp[k+1] = polyRange[k] + remainders[k];
				}
				else
				{
					intVecTemp[k+1] = polyRange[k];
				}
			}

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<hfs.size(); ++j)
				{
					if(localNeeded[j])
					{
						Interval intTemp;
						vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i+1] += intLeft;

						hfs[j].intEval(intTemp, newIntVecTemp);

						if(intTemp > b[j])
						{
							// no intersection on the left half
							newInt = intRight;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(b[j]))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							newInt.width(W);

							continue;
						}
					}
				}

				if(localCounter == hfs.size())
				{
					break;
				}
			}

			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			remainders[i].setInf(Inf);

			newInt = remainders[i];
			newInt.width(W);

			localNeeded = bNeeded;
			localCounter = counter;

			// search an approximation for the upper bound
			for(; W >= DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<hfs.size(); ++j)
				{
					if(localNeeded[j])
					{
						Interval intTemp;
						vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i+1] += intRight;

						hfs[j].intEval(intTemp, newIntVecTemp);

						if(intTemp > b[j])
						{
							// no intersection on the right half
							newInt = intLeft;
							newInt.width(W);
							break;
						}
						else if(intTemp.smallereq(b[j]))
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							newInt.width(W);
							localNeeded[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							newInt.width(W);
							continue;
						}
					}
				}

				if(localCounter == hfs.size())
				{
					break;
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			remainders[i].setSup(Sup);	// set the upper bound

			if(!remainders[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<rangeDim; ++i)
		{
			if(oldRemainders[i].widthRatio(remainders[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}
	}

	if(!bvalid)
	{
		return -1;	// no intersection is detected
	}
	else
	{
		return 1;
	}
}

void gridBox(list<vector<Interval> > & grids, const vector<Interval> & box, const int num)
{
	grids.clear();
	grids.push_back(box);

	for(int i=0; i<box.size(); ++i)
	{
		list<vector<Interval> >::iterator gridIter;
		list<vector<Interval> > newGrids;

		for(; grids.size() > 0;)
		{
			gridIter = grids.begin();

			list<Interval> queue;
			(*gridIter)[i].split(queue, num);

			list<Interval>::iterator iterComponent = queue.begin();
			for(; iterComponent != queue.end(); ++iterComponent)
			{
				vector<Interval> tmpBox = *gridIter;
				tmpBox[i] = *iterComponent;
				newGrids.push_back(tmpBox);
			}

			grids.pop_front();
		}

		grids = newGrids;
	}
}

void exp_int_mat(Interval_matrix & result_ts, Interval_matrix & result_rem, const Interval_matrix & A, const int order)
{
	int dim = A.cols();

	Interval_matrix identity(dim, dim), zero(dim, dim);
	Interval intOne(1);

	for(int i=0; i<dim; ++i)
	{
		identity.set(intOne, i, i);
	}

	if(order == 0)
	{
		result_ts = identity;
		result_rem = zero;
	}
	else
	{
		result_ts = identity;

		for(int k=order; k>0; --k)
		{
			result_ts.div_assign(k);
			result_ts *= A;
			result_ts += identity;
		}

		Interval_matrix matTemp = A;
		matTemp.pow_assign(order+1);
		matTemp.mul_assign(factorial_rec[order+1]);

		double max_A = A.max_norm();
		Interval intTemp(max_A);
		intTemp.exp_assign();

		double maxNorm = intTemp.mag();
		Interval remainder(-maxNorm, maxNorm);

		Interval_matrix matRemainder(dim, dim);
		for(int i=0; i<dim; ++i)
		{
			for(int j=0; j<dim; ++j)
			{
				matRemainder.set(remainder, i, j);
			}
		}

		matTemp *= matRemainder;
		result_rem = matTemp;
	}
}

void int_exp_int_mat(Interval_matrix & result_ts, Interval_matrix & result_rem, const Interval_matrix & A, const double step, const int order)
{
	int dim = A.cols();

	Interval_matrix identity(dim, dim), zero(dim, dim);
	Interval intOne(1);

	for(int i=0; i<dim; ++i)
	{
		identity.set(intOne, i, i);
	}

	if(order == 0)
	{
		result_ts = identity;
		result_rem = zero;
	}
	else
	{
		Interval_matrix R(dim, dim);
		Interval intStep(step);

		for(int i=0; i<dim; ++i)
		{
			R.set(intStep, i, i);
		}

		result_ts = R;
		Interval_matrix mA = A;
		mA.mul_assign(-1.0);

		for(int k=order+1; k>1; --k)
		{
			result_ts.div_assign(k);
			result_ts *= mA;
			result_ts += R;
		}

		Interval_matrix matTemp = mA;
		matTemp.pow_assign(order+1);
		matTemp.mul_assign(factorial_rec[order+2]);
		matTemp.mul_assign(step);

		double max_mA = mA.max_norm();
		Interval intTemp(max_mA);
		intTemp.exp_assign();

		double maxNorm = intTemp.mag();
		Interval remainder(-maxNorm, maxNorm);

		Interval_matrix matRemainder(dim, dim);
		for(int i=0; i<dim; ++i)
		{
			for(int j=0; j<dim; ++j)
			{
				matRemainder.set(remainder, i, j);
			}
		}

		matTemp *= matRemainder;
		result_rem = matTemp;
	}
}



/*
void build_constraint_template(const int d)
{
	vector<Interval> intVecZero;
	Interval intZero, intOne(1), intMOne(-1);

	for(int i=0; i<d; ++i)
	{
		intVecZero.push_back(intZero);
	}

	for(int i=0; i<d; ++i)
	{
		vector<Interval> A = intVecZero;
		A[i] = intOne;

		LinearConstraint lc(A, intZero);
		constraint_template.push_back(lc);
	}

	for(int i=0; i<d; ++i)
	{
		vector<Interval> A = intVecZero;
		A[i] = intMOne;

		LinearConstraint lc(A, intZero);
		constraint_template.push_back(lc);
	}
}
*/









