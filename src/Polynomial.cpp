/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Polynomial.h"
#include "TaylorModel.h"
#include "MyLogger.h"


vector<Interval> factorial_rec;
vector<Interval> power_4;
vector<Interval> double_factorial;

RangeTree::RangeTree()
{
}

RangeTree::RangeTree(const list<Interval> & ranges_input, const list<RangeTree *> & children_input)
{
	ranges = ranges_input;
	children = children_input;
}

RangeTree::RangeTree(const RangeTree & tree)
{
	ranges = tree.ranges;
	children = tree.children;
}

RangeTree::~RangeTree()
{
	list<RangeTree *>::iterator iter = children.begin();

	for(; iter!=children.end(); ++iter)
	{
		delete *iter;
	}

	ranges.clear();
	children.clear();
}

RangeTree & RangeTree::operator = (const RangeTree & tree)
{
	if(this == &tree)
		return *this;

	ranges = tree.ranges;
	children = tree.children;

	return *this;
}





























// class HornerForm

HornerForm::HornerForm()
{
}

HornerForm::HornerForm(const Interval & I):constant(I)
{
}

HornerForm::HornerForm(const Interval & I, const vector<HornerForm> & hfs):constant(I), hornerForms(hfs)
{
}

HornerForm::HornerForm(const HornerForm & hf):constant(hf.constant), hornerForms(hf.hornerForms)
{
}

HornerForm::~HornerForm()
{
	hornerForms.clear();
}

void HornerForm::clear()
{
	constant.set(0,0);
	hornerForms.clear();
}

void HornerForm::intEval(Interval & result, const vector<Interval> & domain) const
{
  //minc();
  //mlog1(sbuilder() << "h_result0: " << result.toString(60));
	result = constant;
  //mlog1(sbuilder() << "constant: " << constant.toString(60));

	for(int i=0; i<hornerForms.size(); ++i)
	{
	  //mlog1(i);
		Interval intHF;
		hornerForms[i].intEval(intHF, domain);
		//mlog1(sbuilder() << "intHF0: " << intHF.toString(60));
		intHF *= domain[i];
		//mlog1(sbuilder() << "domain: " << domain[i].toString(60));
		//mlog1(sbuilder() << "intHF1: " << intHF.toString(60));
		result += intHF;
	}
  //mlog1(sbuilder() << "h_result1: " << result.toString(60));
	//mdec();
}

void HornerForm::insert(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const
{
  //minc();
  //mlog1(sbuilder() << "ode: " << this->toString());
  //mlog("vars", vars);
	Interval intZero;
	int numVars = domain.size();
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert(tmTemp, vars, varsPolyRange, domain, cutoff_threshold);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= domain[0];
    //mlog1(sbuilder() << "rem0: " << result.remainder.toString());
		result.add_assign(tmTemp);
    //mlog1(sbuilder() << "rem1: " << result.remainder.toString());

		for(int i=1; i<hornerForms.size(); ++i) {
      //mlog1(sbuilder() << "rem2(" << i << "): " << 
      //    result.remainder.toString());
			hornerForms[i].insert(tmTemp, vars, varsPolyRange, domain,
			    cutoff_threshold);	// recursive call
      //mlog1(sbuilder() << "rem3: " << result.remainder.toString());
      //mlog("tm(insert)", vars.tms[i-1]);
			tmTemp.mul_insert_assign(vars.tms[i-1], varsPolyRange[i-1], domain,
			    cutoff_threshold);
			result.add_assign(tmTemp);
      //mlog1(sbuilder() << "rem4: " << result.remainder.toString());
		}
	}
	//mlog("result", result);
	//mlog1(sbuilder() << "rem: " << result.remainder.toString());
	//mdec();
}

void HornerForm::insert_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const Interval & cutoff_threshold) const
{
	Interval intZero;

	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, cutoff_threshold);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= step_exp_table[1];
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, cutoff_threshold);	// recursive call
			tmTemp.mul_insert_normal_assign(vars.tms[i-1], varsPolyRange[i-1], step_exp_table, cutoff_threshold);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_ctrunc(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const
{
	Interval intZero;
	int numVars = domain.size();

	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_ctrunc(tmTemp, vars, varsPolyRange, domain, order, cutoff_threshold);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= domain[0];

		tmTemp.ctrunc(domain, order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_ctrunc(tmTemp, vars, varsPolyRange, domain, order, cutoff_threshold);	// recursive call
			tmTemp.mul_insert_ctrunc_assign(vars.tms[i-1], varsPolyRange[i-1], domain, order, cutoff_threshold);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_no_remainder(TaylorModel & result, 
    const TaylorModelVec & vars, const int numVars, const int order, 
    const Interval & cutoff_threshold) const {
	Interval intZero;
	result.clear();
  minc();
	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}
	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_no_remainder(tmTemp, vars, numVars, order, cutoff_threshold);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.nctrunc(order);
    //mforce("temp0", tmTemp);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); i++) {
      //mforce1(sbuilder() << "i: " << i);
      //mforce("hf[i]", hornerForms[i]);
			hornerForms[i].insert_no_remainder(tmTemp, vars, numVars, order,
          cutoff_threshold);	// recursive call
      //mforce("vars[i-1]", vars.tms[i-1]);
      //mforce1(sbuilder() << "order: " << order);
      //mforce1(sbuilder() << "vars: " << numVars);
      //mforce1(sbuilder() << vars.tms[i-1].getParamCount());
      //mforce1(sbuilder() << tmTemp.getParamCount());
			tmTemp.mul_no_remainder_assign(vars.tms[i-1], order, cutoff_threshold);
      //mforce("tempi", tmTemp);
			result.add_assign(tmTemp);
		}
	}
  mdec();
}

void HornerForm::insert_no_remainder_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_no_remainder_no_cutoff(tmTemp, vars, numVars, order);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.nctrunc(order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_no_remainder_no_cutoff(tmTemp, vars, numVars, order);	// recursive call
			tmTemp.mul_no_remainder_no_cutoff_assign(vars.tms[i-1], order);
			result.add_assign(tmTemp);
		}
	}
}

void HornerForm::insert_no_remainder_no_cutoff(Polynomial & result, const TaylorModelVec & vars, const int numVars) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		Polynomial polyConstant(constant, numVars);
		result = polyConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		Polynomial polyTemp;
		hornerForms[0].insert_no_remainder_no_cutoff(polyTemp, vars, numVars);

		polyTemp.mul_assign(0,1);			// multiplied by t
		result += polyTemp;

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_no_remainder_no_cutoff(polyTemp, vars, numVars);	// recursive call
			polyTemp *= vars.tms[i-1].expansion;
			result += polyTemp;
		}
	}
}

void HornerForm::insert_no_remainder_no_cutoff(Polynomial & result, const vector<Polynomial> & vars, const int numVars) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		Polynomial polyConstant(constant, numVars);
		result = polyConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		Polynomial polyTemp;
		hornerForms[0].insert_no_remainder_no_cutoff(polyTemp, vars, numVars);

		polyTemp.mul_assign(0,1);			// multiplied by t
		result += polyTemp;

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_no_remainder_no_cutoff(polyTemp, vars, numVars);	// recursive call
			polyTemp *= vars[i-1];
			result += polyTemp;
		}
	}
}

void HornerForm::insert_ctrunc_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval intZero;
	result.clear();

	if(!constant.subseteq(intZero))
	{
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}

	if(hornerForms.size() > 0)						// the first variable is t
	{
		TaylorModel tmTemp;
		hornerForms[0].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);

		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		tmTemp.remainder *= step_exp_table[1];

		tmTemp.ctrunc_normal(step_exp_table, order);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
			hornerForms[i].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);	// recursive call

			tmTemp.mul_insert_ctrunc_normal_assign(vars.tms[i-1], varsPolyRange[i-1], step_exp_table, order, cutoff_threshold);

			result.add_assign(tmTemp);
		}
	}
}

//substitutes vars in to HF (need to have var for each of the HFs (besides time)
//numVars is the number of parameters in result and numVars
void HornerForm::insert_ctrunc_normal(TaylorModel & result, RangeTree * & tree, 
    const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, 
    const vector<Interval> & step_exp_table, const int numVars, const int order,
    const Interval & cutoff_threshold) const {
  mreset(old);
  mdisable();
  mlog1("hf_i_c_n");
  minc();
  mlog1(sbuilder() << "hf: " << this->toString() << ", order: " << order << 
      " hfs: " << hornerForms.size());
	Interval intZero;
	result.clear();
  

	if(!constant.subseteq(intZero)) {
		TaylorModel tmConstant(constant, numVars);
		result = tmConstant;
	}
	RangeTree *pnode = new RangeTree;

	if(hornerForms.size() > 0) {
    //mforce1(sbuilder() << "numVars: " << numVars);
    //mforce1(sbuilder() << "hornerForms.size(): " << hornerForms.size());
    //mforce1(sbuilder() << "vars.tms.size(): " << vars.tms.size());
    
    #ifdef do_checks
    if(vars.tms.size() != hornerForms.size() - 1) {
	    cout << "vars.tms.size(): " << vars.tms.size() << endl;
	    cout << "hornerForms.size() - 1: " << hornerForms.size() - 1 << endl;
	    throw std::runtime_error("variables in hf different than given");
      mforce1("not equal");
      exit(0);
    }
    #endif
    
    // the first variable is t
		TaylorModel tmTemp;
		RangeTree *child;
		
		
	  mlog1(sbuilder() << "i: " << 0);
	  mlog1(sbuilder() << "hfi: " << hornerForms[0].toString());
    
    //tmTmp is the tmv inserted in hornerforms
		hornerForms[0].insert_ctrunc_normal(tmTemp, child, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);

    //multiply the polynomial part by t (t's index is 0)
		tmTemp.expansion.mul_assign(0,1);			// multiplied by t
		//multiply the remainder by the [0,t_step]
		tmTemp.remainder *= step_exp_table[1];

    //truncate the too high terms (add them to remainder)
		Interval intTrunc;
		tmTemp.expansion.ctrunc_normal(intTrunc, step_exp_table, order);
		tmTemp.remainder += intTrunc;
		
		//too high terms in the hf for t (after multiplying with t)
		mlog1(sbuilder() << "pushing0: " << intTrunc.toString());

		pnode->ranges.push_back(intTrunc); //truncation for hf_t after insert
		pnode->children.push_back(child);

    mlog("tmTemp", tmTemp);
		result.add_assign(tmTemp);

		for(int i=1; i<hornerForms.size(); ++i)
		{
		  mlog1(sbuilder() << "i: " << i);
		  mlog1(sbuilder() << "hfi: " << hornerForms[i].toString());
			TaylorModel tmTemp;
			RangeTree *child;

      //insert initial set for variables in ode
			hornerForms[i].insert_ctrunc_normal(tmTemp, child, vars, varsPolyRange, 
			    step_exp_table, numVars, order, cutoff_threshold);	// recursive call
			mlog("hf_insert", tmTemp);
      mlog("vars[i-1]", vars.tms[i-1]);

      //multiply inserted TM by the variable the hf corresponds too
			Interval tm1, intTrunc2;
      // here coefficient_range = tm1
			tmTemp.mul_insert_ctrunc_normal_assign(tm1, intTrunc2, vars.tms[i-1], 
          varsPolyRange[i-1], step_exp_table, order, cutoff_threshold);
      
			
			mlog1(sbuilder() << "pushing[1]: " << tm1.toString());
			mlog1(sbuilder() << "pushing[2]: " << varsPolyRange[i-1].toString());
			mlog1(sbuilder() << "pushing[3]: " << intTrunc2.toString());
			
      //range of the poly in inserted variables for hf[i] (before multiplying 
      //with x1) after truncation (but only if vars[i] has nonzero polynomial <- removed this condition)
			pnode->ranges.push_back(tm1); 
			pnode->ranges.push_back(varsPolyRange[i-1]); //range of picard (k-1)
			pnode->ranges.push_back(intTrunc2); //higher terms gotten from after multiplying (purely in polynom)
			pnode->children.push_back(child);
      mlog(sbuilder() << "x[" << i << "] * hf_insert", tmTemp);

			result.add_assign(tmTemp);
		}
	}

	tree = pnode;
	mdec();
	mrestore(old);
}

void HornerForm::insert_only_remainder(Interval & result, RangeTree *tree, const TaylorModelVec & vars, const Interval & timeStep) const
{
  mreset(old);
  mdisable();
  mlog1("insert_only_remainder");
  //mlog("vars", vars);
  minc();
	Interval intZero;

	result = intZero;
	list<Interval>::const_iterator iter = tree->ranges.begin();
	list<RangeTree *>::const_iterator child = tree->children.begin();

	if(hornerForms.size() > 0)						// the first variable is t
	{
		Interval intTemp;
		hornerForms[0].insert_only_remainder(intTemp, *child, vars, timeStep);
		mlog1(sbuilder() << "intTemp: " << intTemp.toString());
		intTemp *= timeStep;
  	mlog1(sbuilder() << "hf[0] iter[0]: " << iter->toString());

		intTemp += (*iter);
		result += intTemp;

		++iter;
		++child;

		for(int i=1; i<hornerForms.size(); ++i,++child)
		{
			Interval intTemp2;
			hornerForms[i].insert_only_remainder(intTemp2, *child, vars, timeStep);
			
  		mlog1(sbuilder() << "hf[" << i << "] iter[1]: " << iter->toString());
  		//+ range of polynomial of hf[i] * remainder of vars[i-1]
			Interval newRemainder = (*iter) * vars.tms[i-1].remainder;
			++iter;
  		mlog1(sbuilder() << "hf[" << i << "] iter[2]: " << iter->toString());
  		//range of polynomial of vars[i-1] * remainder of hf[i]
			newRemainder += (*iter) * intTemp2;
			//remainder of vars[i-1].remainder * remainder of hf[i]
			newRemainder += vars.tms[i-1].remainder * intTemp2;
			++iter;
  		mlog1(sbuilder() << "hf[" << i << "] iter[3]: " << iter->toString());
  		//interval evaluation of higher terms
			newRemainder += (*iter);

			result += newRemainder;
			++iter;
		}
	}
	mdec();
	mrestore(old);
}

void HornerForm::dump(FILE *fp, const vector<string> & varNames) const
{
	int numVars = hornerForms.size();

	Interval intZero;
	bool bPlus = false;

	fprintf(fp, " ( ");
	if(!constant.subseteq(intZero))
	{
		bPlus = true;
		constant.dump(fp);
	}

	if(numVars == 0)
	{
		fprintf(fp, " ) ");
		return;
	}

	for(int i=0; i<numVars; ++i)
	{
		if(hornerForms[i].hornerForms.size() != 0 || !hornerForms[i].constant.subseteq(intZero))
		{
			if(bPlus)		// only used to print the "+" symbol
				fprintf(fp, " + ");
			else
				bPlus = true;

			hornerForms[i].dump(fp, varNames);
			fprintf(fp, "* %s", varNames[i].c_str());
		}
	}

	fprintf(fp, " ) ");
}

string HornerForm::toString() const {
  minc();
  std::stringstream ss;
  //don't add "[0] + " in the beginning
  if(string("[0]").compare(constant.toString()) != 0)
    ss << constant.toString() << " + ";
  for (unsigned i=0; i<hornerForms.size(); i++) {
		string s = hornerForms.at(i).toString();
    //if(string("[0.]").compare(s) == 0)
    //  continue;
    if(s.length() == 0)
      continue;
    ss << "x" << i << "(" << s << ") + ";
	}
  //remove last " + "
  
  string ret = ss.str().substr(0, ss.str().length() - 3);
  //if(ret.length() != 0)
  //  mlog1(ret);
  mdec();
  return ret;
}

HornerForm & HornerForm::operator = (const HornerForm & hf)
{
	if(this == &hf)
		return *this;

	constant = hf.constant;
	hornerForms = hf.hornerForms;
	return *this;
}



































// class Polynomial

Polynomial::Polynomial()
{
}

Polynomial::Polynomial(const Interval & constant, const int numVars)
{
	Interval intZero;

	if(!constant.subseteq(intZero))
	{
		Monomial monomial(constant, numVars);
		monomials.push_back(monomial);
	}
}

Polynomial::Polynomial(const RowVector & coefficients)
{
	int numVars = coefficients.size();

	for(int i=0; i<numVars; ++i)
	{
		double dTemp = coefficients.get(i);
		if(dTemp <= THRESHOLD_LOW && dTemp >= -THRESHOLD_LOW)		// dTemp is zero
			continue;

		Interval intTemp(dTemp);
		Monomial monoTemp(intTemp, numVars);
		monoTemp.degrees[i] = 1;
		monoTemp.d = 1;
		monomials.push_back(monoTemp);
	}

	reorder();
}

Polynomial::Polynomial(const vector<Interval> & coefficients)
{
	int numVars = coefficients.size();
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		if(coefficients[i].subseteq(intZero))		// the coefficient is zero
			continue;

		Monomial monoTemp(coefficients[i], numVars);
		monoTemp.degrees[i] = 1;
		monoTemp.d = 1;
		monomials.push_back(monoTemp);
	}

	reorder();
}

Polynomial::Polynomial(const Monomial & monomial)
{
	monomials.push_back(monomial);
}

Polynomial::Polynomial(const list<Monomial> & monos):monomials(monos)
{
	reorder();
}

Polynomial::Polynomial(const Polynomial & polynomial):monomials(polynomial.monomials)
{
}

Polynomial::~Polynomial()
{
	monomials.clear();
}

void Polynomial::reorder()
{
	monomials.sort();
}

void Polynomial::clear()
{
	monomials.clear();
}

void Polynomial::dump_interval(FILE *fp, const vector<string> & varNames) const
{
	if(monomials.size() == 0)
	{
		fprintf(fp, "[0,0]");
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		iter->dump_interval(fp, varNames);
		fprintf(fp, " + ");
	}

	monomials.back().dump_interval(fp, varNames);
}


void Polynomial::serialize(FILE *fp, const vector<string> & tmParams) const {
  
	if(monomials.size() == 0)
	{ 
    ZERO_INTERVAL.serialize(fp);
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		iter->serialize(fp, tmParams);
		fprintf(fp, " + ");
	}

	monomials.back().serialize(fp, tmParams);
}

void Polynomial::dump_constant(FILE *fp, const vector<string> & varNames) const
{
	if(monomials.size() == 0)
	{
		fprintf(fp, "[0,0]");
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		iter->dump_constant(fp, varNames);
		fprintf(fp, " + ");
	}

	monomials.back().dump_constant(fp, varNames);
}

void Polynomial::constant(Interval & result) const
{
	Interval intZero;

	if(monomials.size() > 0 && (monomials.begin())->d == 0)
	{
    //mlog1(sbuilder() << "result: " << monomials.begin());
		//mlog1(sbuilder() << (monomials.back())->coefficient.toString());
		//mlog1(sbuilder() << "size: " << monomials.size());
    //mlog1(sbuilder() << "begin: " << (monomials.begin())->coefficient.toString());
		result = (monomials.begin())->coefficient;
	}
	else
	{
		result = intZero;
	}
}

void Polynomial::intEval(Interval & result, const vector<Interval> & domain) const
{
  //mlog1("here");
  //mlog1(sbuilder() << "poly: " << toString(getVNames(10)));
	HornerForm hf;
	toHornerForm(hf);
  //mlog1(sbuilder() << "hf: " << hf.toString());
  //mlog1(sbuilder() << "p_result0: " << result.toString(60));
	hf.intEval(result, domain);
  //mlog1(sbuilder() << "p_result1: " << result.toString(60));
}

void Polynomial::intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	Interval intZero;
	result = intZero;

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		Interval intTemp;
		iter->intEvalNormal(intTemp, step_exp_table);

		result += intTemp;
	}
}

void Polynomial::inv(Polynomial & result) const
{
	result = *this;
	result.inv_assign();
}

void Polynomial::inv_assign()
{
	list<Monomial>::iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		iter->coefficient.inv_assign();
	}
}

void Polynomial::pow(Polynomial & result, const int degree) const
{
	Polynomial temp = *this;
	result = *this;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= temp;
		}

		d >>= 1;

		if(d > 0)
		{
			temp *= temp;
		}
	}
}

void Polynomial::pow_assign(const int degree)
{
	Polynomial temp = *this;
	Polynomial result = *this;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= temp;
		}

		d >>= 1;

		if(d > 0)
		{
			temp *= temp;
		}
	}

	*this = result;
}

void Polynomial::pow(Polynomial & result, const int degree, const int order) const
{
	Polynomial p = *this;
	p.nctrunc(order);

	Polynomial temp = p;
	result = p;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= temp;
			result.nctrunc(order);
		}

		d >>= 1;

		if(d > 0)
		{
			temp *= temp;
			temp.nctrunc(order);
		}
	}
}

void Polynomial::pow_assign(const int degree, const int order)
{
	Polynomial p = *this;
	p.nctrunc(order);

	Polynomial temp = p;
	Polynomial result = p;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= temp;
			result.nctrunc(order);
		}

		d >>= 1;

		if(d > 0)
		{
			temp *= temp;
			temp.nctrunc(order);
		}
	}

	*this = result;
}

void Polynomial::center()
{
	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		bool bvalid = iter->center();

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}

void Polynomial::add_assign(const Monomial & monomial)
{
	bool bAdded = false;

	list<Monomial>::iterator iter;

	Interval intZero;

	if(monomial.coefficient.subseteq(intZero))
	{
		return;
	}

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		if(monomial < *iter)
		{
			monomials.insert(iter, monomial);
			bAdded = true;
			break;
		}
		else if(monomial == *iter)
		{
			(*iter) += monomial;
/*
			if(iter->coefficient.subseteq(intZero))
			{
				monomials.erase(iter);
			}
*/
			bAdded = true;
			break;
		}
	}

	if(!bAdded)
	{
		monomials.push_back(monomial);
	}
}

void Polynomial::sub_assign(const Monomial & monomial)
{
	Monomial monoTemp;
	monomial.inv(monoTemp);
	add_assign(monoTemp);
}

void Polynomial::mul_assign(const Monomial & monomial)
{
	Interval intZero;

	if(monomial.coefficient.subseteq(intZero))	// the monomial is zero
	{
		clear();
	}
	else
	{
		list<Monomial>::iterator iter;
		for(iter = monomials.begin(); iter != monomials.end(); )
		{
			(*iter) *= monomial;
/*
			if(iter->coefficient.subseteq(intZero))
			{
				iter = monomials.erase(iter);
			}
			else
			{
*/
				++iter;
//			}
		}
	}
}

void Polynomial::mul_assign(const Interval & I)
{
	Interval intZero;

	if(I.subseteq(intZero))	// the interval is zero
	{
		clear();
	}
	else
	{
		list<Monomial>::iterator iter;
		for(iter = monomials.begin(); iter != monomials.end(); )
		{
			iter->coefficient *= I;
			if(iter->coefficient.subseteq(intZero))
			{
				iter = monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
}

void Polynomial::div_assign(const Interval & I)
{
//	Interval intZero;

	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		iter->coefficient /= I;
/*
		if(iter->coefficient.subseteq(intZero))
		{
			iter = monomials.erase(iter);
		}
		else
		{
*/
			++iter;
//		}
	}
}

void Polynomial::mul(Polynomial & result, const Interval & I) const
{
	result = *this;
	result.mul_assign(I);
}

void Polynomial::div(Polynomial & result, const Interval & I) const
{
	result = *this;
	result.div_assign(I);
}

void Polynomial::mul_assign(const int varIndex, const int degree)
{
	list<Monomial>::iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		iter->degrees[varIndex] += degree;
		iter->d += degree;
	}
}

void Polynomial::mul(Polynomial result, const int varIndex, const int degree) const
{
	result = *this;
	result.mul_assign(varIndex, degree);
}

Polynomial & Polynomial::operator = (const Polynomial & polynomial)
{
	if(this == &polynomial)
		return *this;

	monomials = polynomial.monomials;
	return *this;
}

Polynomial & Polynomial::operator += (const Polynomial & polynomial)
{
	Polynomial result;

//	Interval intZero;

	list<Monomial>::const_iterator iterA;	// polynomial A
	list<Monomial>::const_iterator iterB;	// polynomial B

	for(iterA = monomials.begin(), iterB = polynomial.monomials.begin(); ; )
	{
		if(iterA == monomials.end() || iterB == polynomial.monomials.end())
			break;

		if((*iterA) < (*iterB))
	    {
			result.monomials.push_back(*iterA);
			++iterA;
	    }
		else if((*iterB) < (*iterA))
	    {
			result.monomials.push_back(*iterB);
			++iterB;
	    }
		else
		{
			Interval intTemp;
			intTemp = iterA->coefficient + iterB->coefficient;

//			if(!intTemp.subseteq(intZero))
//			{
				Monomial monoTemp(*iterA);
				monoTemp.coefficient = intTemp;
				result.monomials.push_back(monoTemp);
//			}

			++iterA;
			++iterB;
		}
	}

	if(iterA == monomials.end() && iterB != polynomial.monomials.end())
	{
		for(; iterB != polynomial.monomials.end(); ++iterB)
			result.monomials.push_back(*iterB);
	}
	else if(iterA != monomials.end() && iterB == polynomial.monomials.end())
	{
		for(; iterA != monomials.end(); ++iterA)
			result.monomials.push_back(*iterA);
	}

	*this = result;
	return *this;
}

Polynomial & Polynomial::operator -= (const Polynomial & polynomial)
{
	Polynomial polyTemp = polynomial;
	polyTemp.inv_assign();
	*this += polyTemp;

	return *this;
}

Polynomial & Polynomial::operator *= (const Polynomial & polynomial)
{
	Polynomial result;

	if((monomials.size() == 0) || (polynomial.monomials.size() == 0))
	{
		this->clear();
		return *this;
	}

	list<Monomial>::const_iterator iterB;	// polynomial B

	for(iterB = polynomial.monomials.begin(); iterB != polynomial.monomials.end(); ++iterB)
	{
		Polynomial polyTemp = *this;
    //mforce("polyTemp1", polyTemp);
    //mforce("m", *iterB);
		polyTemp.mul_assign(*iterB);
    //mforce("polyTemp2", polyTemp);
		result += polyTemp;
	}

	*this = result;
	return *this;
}

const Polynomial Polynomial::operator + (const Polynomial & polynomial) const
{
	Polynomial result = *this;
	result += polynomial;
	return result;
}

const Polynomial Polynomial::operator - (const Polynomial & polynomial) const
{
	Polynomial result = *this;
	result -= polynomial;
	return result;
}

const Polynomial Polynomial::operator * (const Polynomial & polynomial) const
{
	Polynomial result = *this;
	result *= polynomial;
	return result;
}

void Polynomial::ctrunc(Interval & remainder, const vector<Interval> & domain, const int order)
{
	Polynomial polyTemp;
	Monomial monoTemp;

	for(; monomials.size() > 0;)
	{
		monoTemp = monomials.back();

		if(monoTemp.d > order)
		{
			polyTemp.monomials.insert(polyTemp.monomials.begin(), monoTemp);
			monomials.pop_back();
		}
		else
		{
			break;
		}
	}

	polyTemp.intEval(remainder, domain);
}

void Polynomial::nctrunc(const int order)
{
	Monomial monoTemp;

	for(; monomials.size() > 0;)
	{
		monoTemp = monomials.back();

		if(monoTemp.d > order)
		{
			monomials.pop_back();
		}
		else
		{
			break;
		}
	}
}

//moves terms higher than order to remainder
void Polynomial::ctrunc_normal(Interval & remainder, const vector<Interval> & step_exp_table, const int order)
{
	//mforce("HERE");
	//mlog1(order);
	//mlog1(this);
	Polynomial polyTemp;
	Monomial monoTemp;

	for(; monomials.size() > 0;)
	{
		monoTemp = monomials.back();

		if(monoTemp.d > order)
		{
			polyTemp.monomials.insert(polyTemp.monomials.begin(), monoTemp);
			monomials.pop_back();
		}
		else
		{
			break;
		}
	}

	polyTemp.intEvalNormal(remainder, step_exp_table);
	//mlog1(this);
	//mlog1(&polyTemp);
}

void Polynomial::linearCoefficients(vector<Interval> & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			result[i] = iter->coefficient;
		}
	}
}

void Polynomial::linearCoefficients(RowVector & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			result.set(iter->coefficient.sup(), i);
		}
	}
}

void Polynomial::constraintCoefficients(RowVector & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			if(i > 0)
			{
				result.set(iter->coefficient.sup(), i-1);
			}
		}
	}
}

void Polynomial::constraintCoefficients(vector<Interval> & result) const
{
	// initially, the result should be filled with 0

	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			if(i > 0)
			{
				result[i-1] = iter->coefficient;
			}
		}
	}
}

void Polynomial::toHornerForm(HornerForm & result) const
{
	result.clear();

	if(monomials.size() == 0)
		return;

	int numVars = (monomials.begin())->degrees.size();

	list<Monomial> lstMono = monomials;
	list<Monomial>::iterator iter = lstMono.begin();

	if(iter->d == 0)
	{
		result.constant = iter->coefficient;
		iter = lstMono.erase(iter);

		if(lstMono.size() == 0)
			return;
	}

	vector<list<Monomial> > vlMono;

	for(int i=0; i<numVars; ++i)
	{
		list<Monomial> lst_ith;

		for(iter = lstMono.begin(); iter != lstMono.end();)
		{
			if(iter->degrees[i] > 0)
			{
				iter->degrees[i] -= 1;
				iter->d -= 1;
				lst_ith.push_back(*iter);
				iter = lstMono.erase(iter);
			}
			else
			{
				++iter;
			}
		}

		vlMono.push_back(lst_ith);
	}

	for(int i=0; i<numVars; ++i)
	{
		Polynomial polyTemp(vlMono[i]);
		HornerForm hf;
		polyTemp.toHornerForm(hf);
		result.hornerForms.push_back(hf);
	}
}

void Polynomial::rmConstant()
{
	if(monomials.size() > 0 && (monomials.begin())->d == 0)
	{
    if(monomials.size() == 1) {
      monomials.push_back(
          Monomial(ZERO_INTERVAL, monomials.begin()->degrees.size()));
    }
		monomials.erase( monomials.begin() );
	}
}

int Polynomial::degree() const
{
	if(monomials.size() > 0)
	{
		list<Monomial>::const_iterator iter = monomials.end();
		--iter;
		return iter->d;
	}
	else
	{
		return 0;
	}
}

bool Polynomial::isZero() const
{
	if(monomials.size() == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Polynomial::rmZeroTerms(const vector<int> & indices)
{
	if(indices.size() == 0)
	{
		return;
	}

	list<Monomial>::iterator iter = monomials.begin();

	for(; iter != monomials.end();)
	{
		bool bDeleted = false;

		for(int i=0; i<indices.size(); ++i)
		{
			if(iter->degrees[indices[i]] > 0)
			{
				iter = monomials.erase(iter);
				bDeleted = true;
				break;
			}
		}

		if(bDeleted == false)
		{
			++iter;
		}
	}
}

void Polynomial::cutoff_normal(Interval & intRem, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	Polynomial polyTemp;

	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		Monomial monoTemp;
		bool bvalid = iter->cutoff(monoTemp, cutoff_threshold);

		polyTemp.monomials.push_back(monoTemp);

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}

	polyTemp.intEvalNormal(intRem, step_exp_table);
}

void Polynomial::cutoff(Interval & intRem, const vector<Interval> & domain, const Interval & cutoff_threshold)
{
	Polynomial polyTemp;

	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		Monomial monoTemp;
		bool bvalid = iter->cutoff(monoTemp, cutoff_threshold);

		polyTemp.monomials.push_back(monoTemp);

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}

	polyTemp.intEval(intRem, domain);
}

void Polynomial::cutoff(const Interval & cutoff_threshold)
{
	list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); )
	{
		bool bvalid = iter->cutoff(cutoff_threshold);

		if(!bvalid)
		{
			iter = monomials.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}

void Polynomial::derivative(Polynomial & result, const int varIndex) const
{
	result = *this;

	list<Monomial>::iterator iter;

	for(iter = result.monomials.begin(); iter != result.monomials.end(); )
	{
		if(iter->degrees[varIndex] > 0)
		{
			double tmp = iter->degrees[varIndex];
			iter->degrees[varIndex] -= 1;
			iter->d -= 1;
			iter->coefficient.mul_assign(tmp);
			++iter;
		}
		else
		{
			iter = result.monomials.erase(iter);
		}
	}
}

void Polynomial::LieDerivative(Polynomial & result, const vector<Polynomial> & f) const
{
	derivative(result, 0);

	int rangeDim = f.size();

	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial P;
		derivative(P, i+1);
		P *= f[i];
		result += P;
	}
}

void Polynomial::sub(Polynomial & result, const Polynomial & P, const int order) const
{
	list<Monomial> monomials1, monomials2;
	list<Monomial>::const_iterator iter;

	for(iter = monomials.begin(); iter != monomials.end(); ++iter)
	{
		if(iter->d == order)
		{
			monomials1.push_back(*iter);
		}
	}

	for(iter = P.monomials.begin(); iter != P.monomials.end(); ++iter)
	{
		if(iter->d == order)
		{
			monomials2.push_back(*iter);
		}
	}

	Polynomial P1(monomials1), P2(monomials2);
	result = P1 - P2;
}

void Polynomial::exp_taylor(Polynomial & result, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	const_part.exp_assign();	// exp(c)

	if(F.isZero())				// tm = c
	{
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval I(1);
	Polynomial polyOne(I, numVars);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result = polyOne;

	for(int i=order; i>0; --i)
	{
		Interval intFactor(1);
		intFactor.div_assign((double)i);

		result.mul_assign(intFactor);

		result *= F;
		result.nctrunc(order);
		result.cutoff(cutoff_threshold);

		result += polyOne;
	}

	result.mul_assign(const_part);
}

void Polynomial::rec_taylor(Polynomial & result, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	const_part.rec_assign();	// 1/c

	if(F.isZero())				// tm = c
	{
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval I(1);
	Polynomial polyOne(I, numVars);
	Polynomial F_c;
	F.mul(F_c, const_part);

	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result = polyOne;

	for(int i=order; i>0; --i)
	{
		result.inv_assign();

		result *= F_c;
		result.nctrunc(order);
		result.cutoff(cutoff_threshold);

		result += polyOne;
	}

	result.mul_assign(const_part);
}

void Polynomial::sin_taylor(Polynomial & result, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	if(F.isZero())			// tm = c
	{
		const_part.sin_assign();
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();
	sinc.inv(msinc);
	cosc.inv(mcosc);

	Polynomial polyTemp(sinc, numVars);
	result = polyTemp;

	int k=1;
	Interval I(1);

	Polynomial polyPowerF(I, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		I.div_assign((double)i);

		switch(k)
		{
		case 0:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * sinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 1:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * cosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 2:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * msinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 3:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * mcosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		}
	}

	result.cutoff(cutoff_threshold);
}

void Polynomial::cos_taylor(Polynomial & result, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	if(F.isZero())			// tm = c
	{
		const_part.cos_assign();
		Polynomial polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Interval sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();
	sinc.inv(msinc);
	cosc.inv(mcosc);

	Polynomial polyTemp(cosc, numVars);
	result = polyTemp;

	int k=1;
	Interval I(1);

	Polynomial polyPowerF(I, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		I.div_assign((double)i);

		switch(k)
		{
		case 0:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * cosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 1:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * msinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 2:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * mcosc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		case 3:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);

			polyTemp = polyPowerF;
			Interval intTemp = I * sinc;
			polyTemp.mul_assign(intTemp);

			result += polyTemp;

			break;
		}
		}
	}

	result.cutoff(cutoff_threshold);
}

void Polynomial::log_taylor(Polynomial & result, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	Interval C = const_part;

	const_part.log_assign();	// log(c)

	if(F.isZero())			// tm = c
	{
		Polynomial polyLog(const_part, numVars);
		result = polyLog;

		return;
	}

	Polynomial F_c;
	F.div(F_c, C);

	result = F_c;

	Interval I((double)order);
	result.div_assign(I);			// F/c * (1/order)

	for(int i=order; i>=2; --i)
	{
		Interval J(1);
		J.div_assign((double)(i-1));
		Polynomial polyJ(J, numVars);

		result -= polyJ;
		result.inv_assign();

		result *= F_c;
		result.nctrunc(order);
		result.cutoff(cutoff_threshold);
	}

	Polynomial const_part_poly(const_part, numVars);
	result += const_part_poly;
}

void Polynomial::sqrt_taylor(Polynomial & result, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	Interval C = const_part;
	const_part.sqrt_assign();	// log(c)

	if(F.isZero())			// tm = c
	{
		Polynomial polySqrt(const_part, numVars);
		result = polySqrt;

		return;
	}

	Polynomial F_2c;
	F.div(F_2c, C);

	Interval intTwo(2);
	F_2c.div_assign(intTwo);	// F/2c

	Interval intOne(1);
	Polynomial polyOne(intOne, numVars);

	result = F_2c;

	Interval K(1), J(1);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		// i
		Interval K((double)i);

		// j = 2*i-3
		Interval J((double)j);

		result.inv_assign();
		result.mul_assign( J / K );

		result += polyOne;
		result *= F_2c;
		result.nctrunc(order);
		result.cutoff(cutoff_threshold);
	}

	result += polyOne;

	result.mul_assign(const_part);
}

void Polynomial::toString(string & result, const vector<string> & varNames) const
{
	string strPoly;

	if(monomials.size() == 0)
	{
		strPoly = "(0)";
		return;
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;

	strPoly += '(';

	for(iter = monomials.begin(); iter != iter_last; ++iter)
	{
		string strTemp;
		iter->toString(strTemp, varNames);

		strPoly += strTemp;
		strPoly += ' ';
		strPoly += '+';
		strPoly += ' ';
	}

	string strTemp2;
	monomials.back().toString(strTemp2, varNames);
	strPoly += strTemp2;
	strPoly += ')';

	result = strPoly;
}








































//1, 1/2!, 1/3!, 1/4!, ...
void compute_factorial_rec(const int order)
{
	Interval I(1);
	//mlog1("compute_factorial <");
	//minc();
	factorial_rec.push_back(I);

	for(int i=1; i<=order; ++i)
	{
		I.div_assign((double)i);
		factorial_rec.push_back(I);
		//string s;
		//I.toString(s);
		//mlog1(sbuilder() << s);
	}
	//mdec();
	//mlog1("compute_factorial >");
}

//4, 4^2, 4^3, 4^4, ...
void compute_power_4(const int order)
{
	Interval I(1);

	power_4.push_back(I);

	for(int i=1; i<=order; ++i)
	{
		I.mul_assign(4.0);
		//string s;
		//I.toString(s);
		//mlog1(sbuilder() << s);
		power_4.push_back(I);
	}
}
//1, 2, 2*4, 2*4*6, ...
//1, 3, 3*5, 3*5*7, ...
void compute_double_factorial(const int order)
{
	Interval odd(1), even(1);

	double_factorial.push_back(even);
	double_factorial.push_back(odd);

	for(int i=2; i<=order; ++i)
	{
		if(i%2 == 0)
		{
			even.mul_assign((double)i);
			double_factorial.push_back(even);
			//string s;
			//even.toString(s);
			//mlog1(sbuilder() << s);
		}
		else
		{
			odd.mul_assign((double)i);
			double_factorial.push_back(odd);
			//string s;
			//odd.toString(s);
			//mlog1(sbuilder() << s);
		}
		
	}
}

void computeTaylorExpansion(vector<HornerForm> & result, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

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
			LieDeriv_n[i].LieDerivative(P1, ode);
			LieDeriv_n[i] = P1;

			P1.mul_assign(factorial_rec[j]);
			P1.mul_assign(0,j);

			taylorExpansion[i] += P1;
		}
	}

	result.clear();

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		result.push_back(hf);
	}
}

void computeTaylorExpansion(vector<HornerForm> & result, const vector<Polynomial> & ode, const vector<int> & orders)
{
	int rangeDim = ode.size();

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
		for(int j=1; j<=orders[i]; ++j)
		{
			Polynomial P1;
			LieDeriv_n[i].LieDerivative(P1, ode);
			LieDeriv_n[i] = P1;

			P1.mul_assign(factorial_rec[j]);
			P1.mul_assign(0,j);

			taylorExpansion[i] += P1;
		}
	}

	result.clear();

	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		result.push_back(hf);
	}
}

void computeTaylorExpansion(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

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

	highest.clear();

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=order; ++j)
		{
			Polynomial P;
			LieDeriv_n[i].LieDerivative(P, ode);
			LieDeriv_n[i] = P;

			if(j == order)
			{
				highest.push_back(P);
			}

			P.mul_assign(factorial_rec[j]);
			P.mul_assign(0,j);

			taylorExpansion[i] += P;
		}
	}

	resultMF = taylorExpansion;

	resultHF.clear();
	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		resultHF.push_back(hf);
	}
}

void computeTaylorExpansion(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & ode, const vector<int> & orders)
{
	int rangeDim = ode.size();

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

	highest.clear();

	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=1; j<=orders[i]; ++j)
		{
			Polynomial P;
			LieDeriv_n[i].LieDerivative(P, ode);
			LieDeriv_n[i] = P;

			if(j == orders[i])
			{
				highest.push_back(P);
			}

			P.mul_assign(factorial_rec[j]);
			P.mul_assign(0,j);

			taylorExpansion[i] += P;
		}
	}

	resultMF = taylorExpansion;

	resultHF.clear();
	for(int i=0; i<taylorExpansion.size(); ++i)
	{
		HornerForm hf;
		taylorExpansion[i].toHornerForm(hf);
		resultHF.push_back(hf);
	}
}

void increaseExpansionOrder(vector<HornerForm> & resultHF, vector<Polynomial> & resultMF, vector<Polynomial> & highest, const vector<Polynomial> & taylorExpansion, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

	vector<Polynomial> expansion = taylorExpansion;
	vector<Polynomial> LieDeriv_n = highest;

	for(int i=0; i<rangeDim; ++i)
	{
		Polynomial P1;
		LieDeriv_n[i].LieDerivative(P1, ode);

		highest[i] = P1;

		P1.mul_assign(factorial_rec[order+1]);
		P1.mul_assign(0, order+1);

		expansion[i] += P1;
	}

	resultMF = expansion;

	resultHF.clear();
	for(int i=0; i<expansion.size(); ++i)
	{
		HornerForm hf;
		expansion[i].toHornerForm(hf);
		resultHF.push_back(hf);
	}
}

void increaseExpansionOrder(HornerForm & resultHF, Polynomial & resultMF, Polynomial & highest, const Polynomial & taylorExpansion, const vector<Polynomial> & ode, const int order)
{
	int rangeDim = ode.size();

	Polynomial expansion = taylorExpansion;
	Polynomial LieDeriv_n = highest;

	Polynomial P1;
	LieDeriv_n.LieDerivative(P1, ode);

	highest = P1;

	P1.mul_assign(factorial_rec[order+1]);
	P1.mul_assign(0, order+1);

	expansion += P1;

	resultMF = expansion;
	expansion.toHornerForm(resultHF);
}

string Polynomial::toString(const vector<string> & varNames) const {
  string strPoly;
	if(monomials.size() == 0) {
		return "(0)";
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;
	strPoly += '(';
	for(iter = monomials.begin(); iter != iter_last; ++iter) {
		string strTemp = iter->toString(varNames);
		strPoly += strTemp;
		strPoly += ' ';
		strPoly += '+';
		strPoly += ' ';
	}

	string strTemp2 = monomials.back().toString(varNames);
	strPoly += strTemp2;
	strPoly += ')';

	return strPoly;
}

string Polynomial::toMathematicaString() const {
  string strPoly;
	if(monomials.size() == 0) {
		return "0";
	}

	list<Monomial>::const_iterator iter, iter_last;
	iter_last = monomials.end();
	--iter_last;
	for(iter = monomials.begin(); iter != iter_last; ++iter) {
		string strTemp = iter->toMathematicaString();
		strPoly += strTemp;
		strPoly += ' ';
		strPoly += '+';
		strPoly += ' ';
	}

	string strTemp2 = monomials.back().toMathematicaString();
	strPoly += strTemp2;

	return strPoly;
  return "math poly";
}

HornerForm HornerForm::transform(map<int, int> lookup, int size) const{
  //mlog1("h transform");
  if(hornerForms.size() == 0) {
    return HornerForm(constant);
  }
  mlog1(hornerForms.size());
  Interval intZero;
  for(int i = 0; i < hornerForms.size(); i++) {
    //mlog1(sbuilder() << i << ": " << hornerForms.at(i).toString());
    //mlog1(hornerForms.at(i).constant.subseteq(intZero));
    //check that all non zero parts are retained
    if(!hornerForms.at(i).constant.subseteq(intZero)) {
      if(lookup.find(i) == lookup.end()) {
        throw std::invalid_argument("error: non zero hf in transforming");
      }
    }
  }
  vector<HornerForm> hfs;
  for(int i = 0; i < size; i++) {
    hfs.push_back(HornerForm(intZero));
  }
  
  //remap the new hfs
  for(map<int,int>::iterator it=lookup.begin(); it!=lookup.end(); ++it) {
    //mlog1(sbuilder() << (it->second) << "->" << (it->first));
    hfs.at(it->second) = hornerForms.at(it->first).transform(lookup, size);
  }
  
  HornerForm ret(constant, hfs);
  return ret;
}


HornerForm HornerForm::transform(vector<int> indexes) const {
  if(hornerForms.size() == 0) {
    return HornerForm(constant);
  }
  Interval intZero;
  mlog("indexes", indexes);
  mlog1(sbuilder() << "hf size: " << hornerForms.size());
  for(int i = 0; i < hornerForms.size(); i++) {
    //mlog1(sbuilder() << i << ": " << hornerForms.at(i).toString());
    //mlog1(hornerForms.at(i).constant.subseteq(intZero));
    //check that all non zero parts are retained
    if(!hornerForms.at(i).constant.subseteq(intZero)) {
      //find nonzero variable in component indexes (-1 since time is the first one)
      if(find(indexes.begin(), indexes.end(), i - 1) != indexes.end()) {
        //mlog1(sbuilder() << "if, i: " << i);
      } else {
        throw std::invalid_argument("error: non zero hf in transforming");
      }
    }
  }
  vector<HornerForm> hfs;
  
  for(int i = 0; i < indexes.size() + 1; i++) {
    hfs.push_back(HornerForm(intZero));
  }
  
  hfs.at(0) = hornerForms.at(0).transform(indexes);
  
  for(int i = 0; i < indexes.size(); i++) {
    int varIndex = indexes.at(i) + 1; // +1 since time is 0
    hfs.at(i+1) = hornerForms.at(varIndex).transform(indexes);
  }
  
  
  HornerForm ret(constant, hfs);
  return ret;
}

int Polynomial::equals(const Polynomial & p2) const {
  return 0;
}


int Polynomial::getVariableCount() const {
  if(monomials.size() == 0)
    return -1;
  Monomial m = monomials.front();
  return m.getVariableCount();
}

Polynomial Polynomial::addNVariables(int n) const {
  list<Monomial>::const_iterator it, it_last;
	it_last = monomials.end();
	
	list<Monomial> retMonos;
  for(it = monomials.begin(); it != it_last; it++) {
    Monomial m = it->addNVariables(n);
    retMonos.push_back(m);
	}
	Polynomial ret(retMonos);
	return ret;
}


//returns true if the polynomial is not zero
bool HornerForm::getVars(int vars[]) const {
  minc();
  //mlog1(sbuilder() << "this: " << toString());
  //mlog1(sbuilder() << constant.toString() << ", " << hornerForms.size());
  
  //return true if constant is non zero
  bool ret = !constant.subseteq(Interval());
  int varIndex = 0;
  for(vector<HornerForm>::const_iterator it = hornerForms.begin(); 
      it < hornerForms.end(); it++, varIndex++) {    
    bool nonZero = it->getVars(vars);
    
    //update the variable count for all non zero sub polynomials
    if(nonZero) {
      //-1 because of time
      vars[varIndex - 1]++;
    }
    //return true if any of the sub polynomials are non zero
    ret |= nonZero;
    //mlog1(sbuilder() << "b:" << nonZero);
  }
  //mlog1(sbuilder() << "ret: " << ret);
  mdec();
  return ret;
}


bool HornerForm::isClose(const HornerForm & hf, double d) const {
//  mlog1("checking close");
//  mlog1(sbuilder() << constant.toString() << " and " << 
//      hf.constant.toString());
  if(constant.isClose(hf.constant, d) == false)
    return false;
//  mlog1(sbuilder() << "s1: " << hornerForms.size() << ", s2: " 
//      << hf.hornerForms.size());
  if(hornerForms.size() != hf.hornerForms.size())
    return false;
  for(int i = 0; i < hornerForms.size(); i++) {
    if(hornerForms[i].isClose(hf.hornerForms[i], d) == false)
      return false;
  }
  return true;
}

void Polynomial::filter(vector<int> powers) {
  mlog1("fitlering");
  
  list<Monomial> retMonos;
  
  list<Monomial>::iterator iter;
	for(iter = monomials.begin(); iter != monomials.end(); iter++) {
	  vector<int> monoPowers = iter->degrees;
	  bool good = true;
	  for(int i = 0; i < powers.size(); i++) {
	    if(powers[i] != monoPowers[i]) {
        good = false;
	    }
	  }
	  if(good)
	    retMonos.push_back(*iter);
	}
	monomials = retMonos;
}

vector<int> HornerForm::getVars() const {
  minc();
  mreset(old);
  mlog("hf", *this);
  //mlog1(sbuilder() << "this: " << toString());
  //mlog1(sbuilder() << constant.toString() << ", " << hornerForms.size());
  
  //return true if constant is non zero
  vector<int> ret;
  
  mrestore(old);
  mdec();
  return ret;
  /*
  bool ret = !constant.subseteq(Interval());
  
  int varIndex = 0;
  for(vector<HornerForm>::const_iterator it = hornerForms.begin(); 
      it < hornerForms.end(); it++, varIndex++) {    
    bool nonZero = it->getVars(vars);
    
    //update the variable count for all non zero sub polynomials
    if(nonZero) {
      //-1 because of time
      vars[varIndex - 1]++;
    }
    //return true if any of the sub polynomials are non zero
    ret |= nonZero;
    //mlog1(sbuilder() << "b:" << nonZero);
  }
  
  //mlog1(sbuilder() << "ret: " << ret);
  mdec();
  return ret;*/
}

int HornerForm::getVarCount() const {
	return hornerForms.size();
}