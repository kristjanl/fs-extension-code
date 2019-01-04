/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "TaylorModel.h"
#include "MyLogger.h"
#include "MyComponent.h"
#include "Utils.h"
#include "Exceptions.h"

ParseSetting parseSetting;
ParseResult parseResult;

// class TaylorModel

TaylorModel::TaylorModel()
{
}

TaylorModel::TaylorModel(const Interval & I, const int numVars)
{
	Interval intZero;
	Polynomial polyTemp(I, numVars);
	expansion = polyTemp;
	remainder = intZero;
}

TaylorModel::TaylorModel(const Polynomial & polyExp, const Interval & I):expansion(polyExp), remainder(I)
{
}

TaylorModel::TaylorModel(const RowVector & coefficients)
{
	Interval intZero;
	Polynomial polyTemp(coefficients);
	expansion = polyTemp;
	remainder = intZero;
}

TaylorModel::TaylorModel(const RowVector & coefficients, const Interval & I)
{
	Polynomial polyTemp(coefficients);
	expansion = polyTemp;
	remainder = I;
}

TaylorModel::TaylorModel(const vector<Interval> & coefficients)
{
	Interval intZero;
	Polynomial polyTemp(coefficients);
	expansion = polyTemp;
	remainder = intZero;
}

TaylorModel::TaylorModel(const vector<Interval> & coefficients, const Interval & I)
{
	Polynomial polyTemp(coefficients);
	expansion = polyTemp;
	remainder = I;
}

TaylorModel::TaylorModel(const TaylorModel & tm):expansion(tm.expansion), remainder(tm.remainder)
{
}

TaylorModel::~TaylorModel()
{
	expansion.clear();
}

void TaylorModel::clear()
{
	Interval intZero;
	expansion.clear();
	remainder = intZero;
}

void TaylorModel::dump_interval(FILE *fp, vector<string> const & varNames) const
{
	expansion.dump_interval(fp, varNames);
	fprintf(fp, " + ");
	remainder.dump(fp);
	fprintf(fp, "\n");
}

void TaylorModel::serialize(FILE *fp, vector<string> const & tmParams) const
{
  fprintf(fp, "<");
	expansion.serialize(fp, tmParams);
	fprintf(fp, " , ");
	remainder.serialize(fp);
	fprintf(fp, ">");
	fprintf(fp, "\n");
}

void TaylorModel::dump_constant(FILE *fp, vector<string> const & varNames) const
{
	expansion.dump_constant(fp, varNames);
	fprintf(fp, " + ");
	remainder.dump(fp);
	fprintf(fp, "\n");
}

void TaylorModel::constant(Interval & result) const
{
	expansion.constant(result);
}

void TaylorModel::intEval(Interval & result, const vector<Interval> & domain) const
{
	expansion.intEval(result, domain);
	//mlog1(sbuilder() << "res: " << result.toString(70));
	//mlog1(sbuilder() << "rem: " << remainder.toString(70));
	result += remainder;
}

void TaylorModel::intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	expansion.intEvalNormal(result, step_exp_table);
	result += remainder;
}

void TaylorModel::ctrunc(const vector<Interval> & domain, const int order)
{
	Interval intTemp;
	expansion.ctrunc(intTemp, domain, order);
	remainder += intTemp;
}

void TaylorModel::nctrunc(const int order)
{
	expansion.nctrunc(order);
}

void TaylorModel::ctrunc_normal(const vector<Interval> & step_exp_table, const int order)
{
	Interval intTemp;
	expansion.ctrunc_normal(intTemp, step_exp_table, order);
	remainder += intTemp;
}

void TaylorModel::inv(TaylorModel & result) const
{
	expansion.inv(result.expansion);
	remainder.inv(result.remainder);
}

void TaylorModel::inv_assign()
{
	expansion.inv_assign();
	remainder.inv_assign();
}

void TaylorModel::add(TaylorModel & result, const TaylorModel & tm) const
{
	result.expansion = expansion + tm.expansion;
	result.remainder = remainder + tm.remainder;
}

void TaylorModel::sub(TaylorModel & result, const TaylorModel & tm) const
{
	result.expansion = expansion - tm.expansion;
	result.remainder = remainder - tm.remainder;
}

void TaylorModel::add_assign(const TaylorModel & tm)
{
	expansion += tm.expansion;
	remainder += tm.remainder;
}

void TaylorModel::sub_assign(const TaylorModel & tm)
{
	expansion -= tm.expansion;
	remainder -= tm.remainder;
}

void TaylorModel::mul_ctrunc(TaylorModel & result, const TaylorModel & tm, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const
{
	Interval intZero;
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(!tm.remainder.subseteq(intZero))
	{
		expansion.intEval(P1xI2, domain);
		P1xI2 *= tm.remainder;
	}

	if(!remainder.subseteq(intZero))
	{
		tm.expansion.intEval(P2xI1, domain);
		P2xI1 *= remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc(domain, order);

	result.cutoff(domain, cutoff_threshold);
}

void TaylorModel::mul_ctrunc_normal(TaylorModel & result, const TaylorModel & tm, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const
{
	Interval intZero;
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(!tm.remainder.subseteq(intZero))
	{
		expansion.intEvalNormal(P1xI2, step_exp_table);
		P1xI2 *= tm.remainder;
	}

	if(!remainder.subseteq(intZero))
	{
		tm.expansion.intEvalNormal(P2xI1, step_exp_table);
		P2xI1 *= remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc_normal(step_exp_table, order);

	result.cutoff_normal(step_exp_table, cutoff_threshold);
}

void TaylorModel::mul_no_remainder(TaylorModel & result, const TaylorModel & tm, const int order, const Interval & cutoff_threshold) const
{
  //mforce("res1", *this);
  //mforce("tm", tm);
  //mforce1(sbuilder() << this->getParamCount());
  //mforce1(sbuilder() << tm.getParamCount());
	result.expansion = expansion * tm.expansion;
  //mforce("res2", *this);
	result.expansion.nctrunc(order);

	result.expansion.cutoff(cutoff_threshold);
}

void TaylorModel::mul_no_remainder_no_cutoff(TaylorModel & result, const TaylorModel & tm, const int order) const
{
	result.expansion = expansion * tm.expansion;
	result.expansion.nctrunc(order);
}

void TaylorModel::mul(TaylorModel & result, const Interval & I) const
{
	expansion.mul(result.expansion, I);
	result.remainder = remainder * I;
}

void TaylorModel::mul_ctrunc_assign(const TaylorModel & tm, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold)
{
	TaylorModel result;
	mul_ctrunc(result, tm, domain, order, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_ctrunc_normal_assign(const TaylorModel & tm, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold)
{
	TaylorModel result;
	mul_ctrunc_normal(result, tm, step_exp_table, order, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_no_remainder_assign(const TaylorModel & tm, const int order, const Interval & cutoff_threshold)
{
	TaylorModel result;
	mul_no_remainder(result, tm, order, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_no_remainder_no_cutoff_assign(const TaylorModel & tm, const int order)
{
	TaylorModel result;
	mul_no_remainder_no_cutoff(result, tm, order);
	*this = result;
}

void TaylorModel::mul_assign(const Interval & I)
{
	TaylorModel result;
	mul(result, I);
	*this = result;
}

void TaylorModel::mul_insert(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const
{
  //mlog("tm(mul_insert)", tm);
	Interval intZero;
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(!tm.remainder.subseteq(intZero))
	{
		expansion.intEval(P1xI2, domain);
		P1xI2 *= tm.remainder;
		//mlog1(sbuilder() << "P1xI2: " << P1xI2.toString());
	}

	if(!remainder.subseteq(intZero))
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;
	//mlog("r1", result);

	result.remainder = I1xI2;
	//mlog("r2", result);
	result.remainder += P2xI1;
	//mlog("r3", result);
	result.remainder += P1xI2;
	//mlog1(P1xI2.toString());
	//mlog("r4", result);

	result.cutoff(domain, cutoff_threshold);
}

void TaylorModel::mul_insert_normal(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const
{
	Interval intZero;
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(!tm.remainder.subseteq(intZero))
	{
		expansion.intEvalNormal(P1xI2, step_exp_table);
		P1xI2 *= tm.remainder;
	}

	if(!remainder.subseteq(intZero))
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.cutoff_normal(step_exp_table, cutoff_threshold);
}


void TaylorModel::mul_insert_ctrunc(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const
{
	Interval intZero;
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(!tm.remainder.subseteq(intZero))
	{
		expansion.intEval(P1xI2, domain);
		P1xI2 *= tm.remainder;
	}

	if(!remainder.subseteq(intZero))
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc(domain, order);
	result.cutoff(domain, cutoff_threshold);
}

void TaylorModel::mul_insert_ctrunc_normal(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const
{
	Interval intZero;
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(!tm.remainder.subseteq(intZero))
	{
		expansion.intEvalNormal(P1xI2, step_exp_table);
		P1xI2 *= tm.remainder;
	}

	if(!remainder.subseteq(intZero))
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc_normal(step_exp_table, order);
	result.cutoff_normal(step_exp_table, cutoff_threshold);
}

void TaylorModel::mul_insert_ctrunc_normal(TaylorModel & result, Interval & tm1, Interval & intTrunc, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const
{
  mreset(old);
  mdisable();
  //mlog1("mul_insert_ctrunc_normal2");
  //mlog("this", *this);
  //mlog("tm", tm);
  //mlog1(sbuilder() << "order: " << order);
	Polynomial P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	Interval intZero;
	tm1 = intZero;
	intTrunc = intZero;

  //mforce1(tm.remainder.toString());
	if(!tm.remainder.subseteq(intZero) || true)
	{
	  //mlog1("evaling");
		expansion.intEvalNormal(P1xI2, step_exp_table);
		tm1 = P1xI2;
		//mlog1(tm1.toString());
		P1xI2 *= tm.remainder;
	}

	if(!remainder.subseteq(intZero))
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

  //mlog("result (before ctrunc)", result);
	result.expansion.ctrunc_normal(intTrunc, step_exp_table, order);
  //mlog("result (before cutoff)", result);
	//mlog1(sbuilder() << "intTrunc: " << intTrunc.toString());

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	intTrunc += intRound;
	//mlog1(sbuilder() << "intRound: " << intRound.toString());

	result.remainder += intTrunc;
  //mlog("result (after cutoff)", result);
  mrestore(old);
}

void TaylorModel::mul_insert_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold)
{
	TaylorModel result;
	//mlog("tm(mul_insert_assign)", tm);
	mul_insert(result, tm, tmPolyRange, domain, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_insert_normal_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	TaylorModel result;
	mul_insert_normal(result, tm, tmPolyRange, step_exp_table, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_insert_ctrunc_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold)
{
	TaylorModel result;
	mul_insert_ctrunc(result, tm, tmPolyRange, domain, order, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_insert_ctrunc_normal_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold)
{
	TaylorModel result;
	mul_insert_ctrunc_normal(result, tm, tmPolyRange, step_exp_table, order, cutoff_threshold);
	*this = result;
}

void TaylorModel::mul_insert_ctrunc_normal_assign(Interval & tm1, Interval & intTrunc, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold)
{
  mreset(old);
  mdisable();
  //mlog1("mul_insert_ctrunc_normal_assign2");
  minc();
	TaylorModel result;
	mul_insert_ctrunc_normal(result, tm1, intTrunc, tm, tmPolyRange, step_exp_table, order, cutoff_threshold);
	//mlog("result", result);
	//mlog1(sbuilder() << "tm1: " << tm1.toString());
	//mlog1(sbuilder() << "intTrunc: " << intTrunc.toString());
	*this = result;
	mdec();
	mrestore(old);
}

void TaylorModel::div(TaylorModel & result, const Interval & I) const
{
	expansion.div(result.expansion, I);
	result.remainder = remainder / I;
}

void TaylorModel::div_assign(const Interval & I)
{
	expansion.div_assign(I);
	remainder /= I;
}

void TaylorModel::derivative(TaylorModel & result, const int varIndex) const
{
	Interval intZero;
	result = *this;

	list<Monomial>::iterator iter;

	for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end();)
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
			iter = result.expansion.monomials.erase(iter);
		}
	}

	result.remainder = intZero;
}

void TaylorModel::LieDerivative_no_remainder(TaylorModel & result, const TaylorModelVec & f, const int order, const Interval & cutoff_threshold) const
{
	derivative(result, 0);

	int rangeDim = f.tms.size();

	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp;
		derivative(tmTemp, i+1);
		tmTemp.mul_no_remainder_assign(f.tms[i], order, cutoff_threshold);
		result.add_assign(tmTemp);
	}
}

void TaylorModel::integral(TaylorModel & result, const Interval & I) const
{
	result = *this;

	list<Monomial>::iterator iter;

	for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end(); ++iter)
	{
		if(iter->degrees[0] > 0)
		{
			iter->degrees[0] += 1;
			iter->d += 1;
			double tmp = iter->degrees[0];
			iter->coefficient.div_assign(tmp);
		}
		else
		{
			iter->degrees[0] += 1;
			iter->d += 1;
		}
	}

	result.remainder *= I;
}

void TaylorModel::integral_no_remainder(TaylorModel & result) const
{
	result = *this;

	list<Monomial>::iterator iter;

	for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end(); ++iter)
	{
		if(iter->degrees[0] > 0)
		{
			iter->degrees[0] += 1;
			iter->d += 1;
			double tmp = iter->degrees[0];
			iter->coefficient.div_assign(tmp);
		}
		else
		{
			iter->degrees[0] += 1;
			iter->d += 1;
		}
	}
}

void TaylorModel::linearCoefficients(vector<Interval> & result) const
{
	expansion.linearCoefficients(result);
}

void TaylorModel::toHornerForm(HornerForm & result, Interval & I) const
{
	expansion.toHornerForm(result);
	I = remainder;
}

void TaylorModel::insert(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		list<Monomial>::iterator iter;

		for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end(); )
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm hf;
		expansion.toHornerForm(hf);

		hf.insert(result, vars, varsPolyRange, domain, cutoff_threshold);
		result.remainder += remainder;
	}
}

void TaylorModel::insert_ctrunc(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		list<Monomial>::iterator iter;

		for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end(); )
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm hf;
		expansion.toHornerForm(hf);

		hf.insert_ctrunc(result, vars, varsPolyRange, domain, order, cutoff_threshold);
		result.remainder += remainder;
	}
}

void TaylorModel::insert_no_remainder(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		list<Monomial>::iterator iter;

		for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end();)
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm hf;
		expansion.toHornerForm(hf);

		hf.insert_no_remainder(result, vars, numVars, order, cutoff_threshold);
	}
}

void TaylorModel::insert_no_remainder_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		list<Monomial>::iterator iter;

		for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end();)
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm hf;
		expansion.toHornerForm(hf);

		hf.insert_no_remainder_no_cutoff(result, vars, numVars, order);
	}
}

//numVars is the number of parameters on the vars TMV
void TaylorModel::insert_ctrunc_normal(TaylorModel & result, const TaylorModelVec & vars, 
    const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, 
    const int numVars, const int order, const Interval & cutoff_threshold) const 
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		list<Monomial>::iterator iter;

		for(iter = result.expansion.monomials.begin(); iter != result.expansion.monomials.end(); )
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.monomials.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm hf;
		expansion.toHornerForm(hf);

		hf.insert_ctrunc_normal(result, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);
		result.remainder += remainder;
	}
}

void TaylorModel::evaluate_t(TaylorModel & result, const vector<Interval> & step_exp_table) const
{
	result.expansion.clear();
	result.remainder = remainder;

	if(expansion.monomials.size() == 0)
		return;

	list<Monomial>::const_iterator iter;
	Interval intZero;

	if(step_exp_table[1].subseteq(intZero))		// t = 0
	{
		for(iter = expansion.monomials.begin(); iter != expansion.monomials.end(); ++iter)
		{
			if(iter->degrees[0] == 0)
			{
				result.expansion.add_assign(*iter);
			}
		}
	}
	else
	{
		for(iter = expansion.monomials.begin(); iter != expansion.monomials.end(); ++iter)
		{
			Monomial monoTemp = *iter;
			int tmp = monoTemp.degrees[0];

			if(tmp > 0)
			{
				monoTemp.coefficient *= step_exp_table[tmp];
				monoTemp.d -= tmp;
				monoTemp.degrees[0] = 0;
			}

			result.expansion.add_assign(monoTemp);
		}
	}
}

void TaylorModel::mul(TaylorModel & result, const int varIndex, const int degree) const
{
	result = *this;
	result.mul_assign(varIndex, degree);
}

void TaylorModel::mul_assign(const int varIndex, const int degree)
{
	expansion.mul_assign(varIndex, degree);
}

void TaylorModel::rmConstant()
{
	expansion.rmConstant();
}

void TaylorModel::cutoff_normal(const vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	Interval intRem;
	expansion.cutoff_normal(intRem, step_exp_table, cutoff_threshold);
	remainder += intRem;
}

void TaylorModel::cutoff(const vector<Interval> & domain, const Interval & cutoff_threshold)
{
	Interval intRem;
	expansion.cutoff(intRem, domain, cutoff_threshold);
	remainder += intRem;
}

void TaylorModel::cutoff(const Interval & cutoff_threshold)
{
	expansion.cutoff(cutoff_threshold);
}

int TaylorModel::degree() const
{
	return expansion.degree();
}

bool TaylorModel::isZero() const
{
	Interval intZero;

	if(expansion.isZero())
	{
		if(remainder.subseteq(intZero))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

void TaylorModel::center_nc()
{
	expansion.center();
}

void TaylorModel::rmZeroTerms(const vector<int> & indices)
{
	expansion.rmZeroTerms(indices);
}

void TaylorModel::normalize(vector<Interval> & domain)
{
	int domainDim = domain.size();

	// compute the center of the original domain and make it origin-centered
	vector<Interval> intVecCenter;
	for(int i=1; i<domainDim; ++i)		// we omit the time dimension
	{
		Interval M;
		domain[i].remove_midpoint(M);
		intVecCenter.push_back(M);
	}

	// compute the scalars
	Interval intZero;
	vector<vector<Interval> > coefficients;
	vector<Interval> row;

	for(int i=0; i<domainDim; ++i)
	{
		row.push_back(intZero);
	}

	for(int i=0; i<domainDim-1; ++i)
	{
		coefficients.push_back(row);
	}

	for(int i=1; i<domainDim; ++i)
	{
		Interval M;
		domain[i].mag(M);
		coefficients[i-1][i] = M;
	}

	TaylorModelVec newVars(coefficients);
	for(int i=0; i<domainDim-1; ++i)
	{
		TaylorModel tmTemp(intVecCenter[i], domainDim);
		newVars.tms[i].add_assign(tmTemp);
	}

	Interval intUnit(-1,1);
	for(int i=1; i<domainDim; ++i)
	{
		domain[i] = intUnit;
	}

	TaylorModel tmTemp;
	insert_no_remainder_no_cutoff(tmTemp, newVars, domainDim, degree());
	expansion = tmTemp.expansion;
}

void TaylorModel::polyRange(Interval & result, const vector<Interval> & domain) const
{
	expansion.intEval(result, domain);
}

void TaylorModel::polyRangeNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	expansion.intEvalNormal(result, step_exp_table);
}

void TaylorModel::exp_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	TaylorModel tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	const_part.exp_assign();	// exp(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel tmExp(const_part, numVars);
		result = tmExp;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	ranges.push_back(const_part);			// keep the unchanged part

	Interval I(1);
	TaylorModel tmOne(I, numVars);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result = tmOne;
	Interval tmFPolyRange;
	tmF.polyRangeNormal(tmFPolyRange, step_exp_table);

	for(int i=order; i>0; --i)
	{
		Interval intFactor(1);
		intFactor.div_assign((double)i);

		result.mul_assign(intFactor);

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmFPolyRange);		// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part

		result.add_assign(tmOne);
	}

	result.mul_assign(const_part);

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	ranges.push_back(intRound);				// keep the unchanged part
	result.remainder += intRound;

	Interval rem, tmRange;
	ranges.push_back(tmFPolyRange);			// keep the unchanged part
	tmRange = tmFPolyRange + tmF.remainder;
	exp_taylor_remainder(rem, tmRange, order+1);

	result.remainder += const_part * rem;
}

void TaylorModel::rec_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	TaylorModel tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	const_part.rec_assign();	// 1/c

	if(tmF.isZero())			// tm = c
	{
		TaylorModel tmRec(const_part, numVars);
		result = tmRec;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	Interval I(1);
	TaylorModel tmOne(I, numVars);
	TaylorModel tmF_c;

	ranges.push_back(const_part);			// keep the unchanged part
	tmF.mul(tmF_c, const_part);

	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result = tmOne;
	Interval tmF_cPolyRange;
	tmF_c.polyRangeNormal(tmF_cPolyRange, step_exp_table);

	for(int i=order; i>0; --i)
	{
		result.inv_assign();

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF_c, tmF_cPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmF_cPolyRange);	// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part

		result.add_assign(tmOne);
	}

	result.mul_assign(const_part);

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	ranges.push_back(intRound);				// keep the unchanged part
	result.remainder += intRound;

	Interval rem, tmF_cRange;
	ranges.push_back(tmF_cPolyRange);		// keep the unchanged part
	tmF_cRange = tmF_cPolyRange + tmF_c.remainder;

	rec_taylor_remainder(rem, tmF_cRange, order+1);

	result.remainder += rem * const_part;
}

void TaylorModel::sin_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	TaylorModel tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	if(tmF.isZero())			// tm = c
	{
		const_part.sin_assign();
		TaylorModel tmSin(const_part, numVars);
		result = tmSin;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	ranges.push_back(const_part);				// keep the unchanged part

	Interval sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();
	sinc.inv(msinc);
	cosc.inv(mcosc);

	TaylorModel tmTemp(sinc, numVars);
	result = tmTemp;

	Interval tmFPolyRange;
	tmF.polyRangeNormal(tmFPolyRange, step_exp_table);

	int k=1;
	Interval I(1);

	TaylorModel tmPowerTmF(I, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		I.div_assign((double)i);

		switch(k)
		{
		case 0:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * sinc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		case 1:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * cosc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		case 2:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * msinc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		case 3:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * mcosc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		}
	}

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	ranges.push_back(intRound);					// keep the unchanged part
	result.remainder += intRound;

	// evaluate the remainder
	Interval tmRange, rem;
	tmRange = tmFPolyRange + tmF.remainder;

	ranges.push_back(tmFPolyRange);				// keep the unchanged part
	sin_taylor_remainder(rem, const_part, tmRange, order+1);

	result.remainder += rem;
}

void TaylorModel::cos_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	TaylorModel tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	if(tmF.isZero())			// tm = c
	{
		const_part.cos_assign();
		TaylorModel tmCos(const_part, numVars);
		result = tmCos;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	ranges.push_back(const_part);				// keep the unchanged part

	Interval sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();
	sinc.inv(msinc);
	cosc.inv(mcosc);

	TaylorModel tmTemp(cosc, numVars);
	result = tmTemp;

	Interval tmFPolyRange;
	tmF.polyRangeNormal(tmFPolyRange, step_exp_table);

	int k=1;
	Interval I(1);

	TaylorModel tmPowerTmF(I, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		I.div_assign((double)i);

		switch(k)
		{
		case 0:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * cosc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		case 1:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * msinc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		case 2:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * mcosc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		case 3:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			tmTemp = tmPowerTmF;

			Interval intTemp = I * sinc;
			ranges.push_back(intTemp);			// keep the unchanged part
			tmTemp.mul_assign(intTemp);

			result.add_assign(tmTemp);

			break;
		}
		}
	}

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	ranges.push_back(intRound);					// keep the unchanged part
	result.remainder += intRound;

	// evaluate the remainder
	Interval tmRange, rem;
	tmRange = tmFPolyRange + tmF.remainder;

	ranges.push_back(tmFPolyRange);				// keep the unchanged part
	cos_taylor_remainder(rem, const_part, tmRange, order+1);

	result.remainder += rem;
}

void TaylorModel::log_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	TaylorModel tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	Interval C = const_part;
	ranges.push_back(const_part);			// keep the unchanged part

	const_part.log_assign();	// log(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel tmLog(const_part, numVars);
		result = tmLog;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	TaylorModel tmF_c;
	tmF.div(tmF_c, C);

	result = tmF_c;

	Interval I((double)order);
	result.div_assign(I);			// F/c * (1/order)

	Interval tmF_cPolyRange;
	tmF_c.polyRangeNormal(tmF_cPolyRange, step_exp_table);

	for(int i=order; i>=2; --i)
	{
		Interval J(1);
		J.div_assign((double)(i-1));
		TaylorModel tmJ(J, numVars);

		result.sub_assign(tmJ);
		result.inv_assign();

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF_c, tmF_cPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmF_cPolyRange);	// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part
	}

	TaylorModel const_part_tm(const_part, numVars);
	result.add_assign(const_part_tm);

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	ranges.push_back(intRound);				// keep the unchanged part
	result.remainder += intRound;

	Interval rem, tmF_cRange;
	ranges.push_back(tmF_cPolyRange);		// keep the unchanged part
	tmF_cRange = tmF_cPolyRange + tmF_c.remainder;

	log_taylor_remainder(rem, tmF_cRange, order+1);

	result.remainder += rem;
}

void TaylorModel::sqrt_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	Interval const_part;

	TaylorModel tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	Interval C = const_part;
	ranges.push_back(const_part);			// keep the unchanged part

	const_part.sqrt_assign();	// log(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel tmSqrt(const_part, numVars);
		result = tmSqrt;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	TaylorModel tmF_2c;
	tmF.div(tmF_2c, C);

	Interval intTwo(2);
	tmF_2c.div_assign(intTwo);	// F/2c

	Interval intOne(1);
	TaylorModel tmOne(intOne, numVars);

	result = tmF_2c;

	Interval K(1), J(1);

	Interval tmF_2cPolyRange;
	tmF_2c.polyRangeNormal(tmF_2cPolyRange, step_exp_table);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		// i
		Interval K((double)i);

		// j = 2*i-3
		Interval J((double)j);

		result.inv_assign();
		result.mul_assign( J / K );

		result.add_assign(tmOne);

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF_2c, tmF_2cPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmF_2cPolyRange);	// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part
	}

	result.add_assign(tmOne);

	result.mul_assign(const_part);

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);
	ranges.push_back(intRound);				// keep the unchanged part
	result.remainder += intRound;

	Interval rem, tmF_cRange = tmF_2cPolyRange;
	tmF_cRange.mul_assign(2.0);

	ranges.push_back(tmF_cRange);			// keep the unchanged part

	Interval tmF_c_remainder = tmF_2c.remainder;
	tmF_c_remainder.mul_assign(2.0);
	tmF_cRange += tmF_c_remainder;

	sqrt_taylor_remainder(rem, tmF_cRange, order+1);

	result.remainder += rem * const_part;
}

Interval TaylorModel::getRemainder() const
{
	return remainder;
}

void TaylorModel::getExpansion(Polynomial & P) const
{
	P = expansion;
}

TaylorModel & TaylorModel::operator = (const TaylorModel & tm)
{
	if(this == &tm)
		return *this;

	expansion = tm.expansion;
	remainder = tm.remainder;
	return *this;
}











































// class TaylorModelVec

TaylorModelVec::TaylorModelVec()
{
}

TaylorModelVec::TaylorModelVec(int n) {
  tms.reserve(n);
}

TaylorModelVec::TaylorModelVec(const vector<TaylorModel> & tms_input):tms(tms_input)
{
}

TaylorModelVec::TaylorModelVec(const Matrix & coefficients)
{
	int cols = coefficients.cols();
	RowVector rowVec(cols);

	int rows = coefficients.rows();
  tms.reserve(rows);

	for(int i=0; i<rows; ++i)
	{
		coefficients.row(rowVec, i);
		TaylorModel tmTemp(rowVec);
		tms.push_back(tmTemp);
	}
}

TaylorModelVec::TaylorModelVec(const Matrix & coefficients, const vector<Interval> & remainders)
{
	int cols = coefficients.cols();
	RowVector rowVec(cols);

	int rows = coefficients.rows();
  tms.reserve(rows);

	for(int i=0; i<rows; ++i)
	{
		coefficients.row(rowVec, i);
		TaylorModel tmTemp(rowVec, remainders[i]);
		tms.push_back(tmTemp);
	}
}

TaylorModelVec::TaylorModelVec(const vector<Interval> & constants, const int numVars)
{
  tms.reserve(constants.size());
	for(int i=0; i<constants.size(); ++i)
	{
		TaylorModel tmTemp(constants[i], numVars);
		tms.push_back(tmTemp);
	}
}

TaylorModelVec::TaylorModelVec(const vector<vector<Interval> > & coefficients)
{
  tms.reserve(coefficients.size());
	for(int i=0; i<coefficients.size(); ++i)
	{
		TaylorModel tmTemp(coefficients[i]);
		tms.push_back(tmTemp);
	}
}

TaylorModelVec::TaylorModelVec(const vector<vector<Interval> > & coefficients, const vector<Interval> & remainders)
{
  tms.reserve(coefficients.size());
	for(int i=0; i<coefficients.size(); ++i)
	{
		TaylorModel tmTemp(coefficients[i], remainders[i]);
		tms.push_back(tmTemp);
	}
}

TaylorModelVec::TaylorModelVec(const vector<Interval> & intVec, vector<Interval> & domain)
{
	int rangeDim = intVec.size();
	domain = intVec;
	Interval intZero;
	domain.insert(domain.begin(), intZero);

	vector<Interval> intVecCenter;
	for(int i=0; i<rangeDim; ++i)		// we omit the time dimension
	{
		double center = intVec[i].midpoint();
		Interval intTemp(center);
		intVecCenter.push_back(intTemp);
		domain[i+1].sub_assign(center);
	}

	// compute the scalars
	Matrix coefficients(rangeDim, rangeDim+1);
	for(int i=0; i<rangeDim; ++i)
	{
		coefficients.set( domain[i+1].sup() , i, i+1);
	}

	TaylorModelVec tmvTemp(coefficients);
	tms = tmvTemp.tms;

	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel tmTemp(intVecCenter[i], rangeDim+1);
		tms[i].add_assign(tmTemp);
	}

	Interval intUnit(-1,1);
	for(int i=0; i<rangeDim; ++i)
	{
		domain[i+1] = intUnit;
	}
}

TaylorModelVec::TaylorModelVec(const TaylorModelVec & tmv):tms(tmv.tms)
{
}

TaylorModelVec::~TaylorModelVec()
{
	tms.clear();
}

void TaylorModelVec::clear()
{
	tms.clear();
}

void TaylorModelVec::dump_interval(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames) const
{
	for(int i=0; i<tms.size(); ++i)
	{
		fprintf(fp, "%s = ", stateVarNames[i].c_str());
		tms[i].dump_interval(fp, tmVarNames);
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
}

void TaylorModelVec::serialize_old(FILE *fp, const vector<string> & tmParams) const {
	fprintf(fp, "my models{\n");
  bool first = true;
	for(int i=0; i<tms.size(); ++i) {
	  if(first == false)
	    fprintf(fp, ",");
	  first = false;
		tms[i].dump_interval(fp, tmParams);
	}
	fprintf(fp, "}\n");
}
void TaylorModelVec::serialize(FILE *fp, const vector<string> & tmParams) const {
	fprintf(fp, "my models{\n");
  bool first = true;
	for(int i=0; i<tms.size(); ++i) {
	  if(first == false)
	    fprintf(fp, ",");
	  first = false;
		tms[i].serialize(fp, tmParams);
	}
	fprintf(fp, "}\n");
}

void TaylorModelVec::dump_constant(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames) const
{
	for(int i=0; i<tms.size(); ++i)
	{
		fprintf(fp, "%s = ", stateVarNames[i].c_str());
		tms[i].dump_constant(fp, tmVarNames);
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
}

void TaylorModelVec::constant(vector<Interval> & result) const
{
	result.clear();
  result.reserve(tms.size());

	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].constant(I);
		result.push_back(I);
	}
}

void TaylorModelVec::intEval(vector<Interval> & result, const vector<Interval> & domain) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].intEval(I, domain);
		result.push_back(I);
	}
}

void TaylorModelVec::intEvalNormal(vector<Interval> & result, const vector<Interval> & step_exp_table) const
{
	result.clear();
  result.reserve(tms.size());
	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].intEvalNormal(I, step_exp_table);
		result.push_back(I);
	}
}

void TaylorModelVec::ctrunc(const vector<Interval> & domain, const int order)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].ctrunc(domain, order);
}

void TaylorModelVec::nctrunc(const int order)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].nctrunc(order);
}

void TaylorModelVec::ctrunc_normal(const vector<Interval> & step_exp_table, const int order)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].ctrunc_normal(step_exp_table, order);
}

void TaylorModelVec::ctrunc(const vector<Interval> & domain, const vector<int> & orders)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].ctrunc(domain, orders[i]);
}

void TaylorModelVec::nctrunc(const vector<int> & orders)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].nctrunc(orders[i]);
}

void TaylorModelVec::ctrunc_normal(const vector<Interval> & step_exp_table, const vector<int> & orders)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].ctrunc_normal(step_exp_table, orders[i]);
}

void TaylorModelVec::inv(TaylorModelVec & result) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].inv(tmTemp);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::add(TaylorModelVec & result, const TaylorModelVec & tmv) const
{
	result.clear();
  result.tms.reserve(tmv.tms.size());

	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].add(tmTemp, tmv.tms[i]);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::sub(TaylorModelVec & result, const TaylorModelVec & tmv) const
{
	result.clear();
  result.tms.reserve(tmv.tms.size());

	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].sub(tmTemp, tmv.tms[i]);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::add_assign(const TaylorModelVec & tmv)
{
	TaylorModelVec result;
	add(result, tmv);
	*this = result;
}

void TaylorModelVec::sub_assign(const TaylorModelVec & tmv)
{
	TaylorModelVec result;
	sub(result, tmv);
	*this = result;
}

void TaylorModelVec::mul(TaylorModelVec & result, const Interval & I) const
{
	result.clear();
  result.tms.reserve(tms.size());

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].mul(tmTemp, I);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::mul_assign(const Interval & I)
{
	TaylorModelVec result;
	mul(result, I);
	*this = result;
}

void TaylorModelVec::div(TaylorModelVec & result, const Interval & I) const
{
	result.clear();
  result.tms.reserve(tms.size());

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].div(tmTemp, I);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::div_assign(const Interval & I)
{
	TaylorModelVec result;
	div(result, I);
	*this = result;
}

void TaylorModelVec::derivative(TaylorModelVec & result, const int varIndex) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].derivative(tmTemp, varIndex);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::LieDerivative_no_remainder(TaylorModelVec & result, const TaylorModelVec & f, const int order, const Interval & cutoff_threshold) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].LieDerivative_no_remainder(tmTemp, f, order, cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::LieDerivative_no_remainder(TaylorModelVec & result, const TaylorModelVec & f, const vector<int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].LieDerivative_no_remainder(tmTemp, f, orders[i], cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::integral(TaylorModelVec & result, const Interval & I) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].integral(tmTemp, I);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::integral_no_remainder(TaylorModelVec & result) const
{
	result.clear();
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].integral_no_remainder(tmTemp);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::linearCoefficients(vector<vector<Interval> > & result) const
{
	for(int i=0; i<tms.size(); ++i)
	{
		tms[i].linearCoefficients(result[i]);
	}
}

void TaylorModelVec::rmZeroTerms(const vector<int> & indices)
{
	if(indices.size() == 0)
	{
		return;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		tms[i].rmZeroTerms(indices);
	}
}

void TaylorModelVec::insert(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert(tmTemp, vars, varsPolyRange, domain, cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_ctrunc(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_ctrunc(tmTemp, vars, varsPolyRange, domain, order, cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_no_remainder(TaylorModelVec & result, const TaylorModelVec & vars, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_no_remainder(tmTemp, vars, numVars, order, cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_ctrunc(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const vector<int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_ctrunc(tmTemp, vars, varsPolyRange, domain, orders[i], cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_no_remainder(TaylorModelVec & result, const TaylorModelVec & vars, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_no_remainder(tmTemp, vars, numVars, orders[i], cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, step_exp_table, numVars, orders[i], cutoff_threshold);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::insert_ctrunc_normal(TaylorModelVec & result, 
		const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, 
		const MySettings *settings, MyComponent *component, const int numVars) const {
	//settings.step_end_exp_table, rightParamCount, settings.order, *settings.cutoff

	//int var = comp2->solveIndexes[j]; //local index inside of component
	//int globalIndex = comp2->compVars[var]; //global index out of all vars
	//int varOrder = settings->orderLookup[globalIndex];
	//mlog1(sbuilder() << "var: " << var);

	result.clear();
	for(int i=0; i<tms.size(); i++) {
		//tms is ordered the same as allVars, which references the global index
		int globalIndex = component->allVars[i]; 
		int varOrder = settings->orderLookup.at(globalIndex);
		//mforce1(sbuilder() << "globalIndex: " << globalIndex << ", order: " << varOrder);
		TaylorModel tmTemp;
		tms[i].insert_ctrunc_normal(tmTemp, vars, varsPolyRange, 
				settings->step_exp_table, numVars, varOrder, *settings->cutoff);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::evaluate_t(TaylorModelVec & result, const vector<Interval> & step_exp_table) const
{
	result.clear();
  result.tms.reserve(tms.size());
	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].evaluate_t(tmTemp, step_exp_table);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::mul(TaylorModelVec & result, const int varIndex, const int degree) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].mul(tmTemp, varIndex, degree);
		result.tms.push_back(tmTemp);
	}
}

void TaylorModelVec::mul_assign(const int varIndex, const int degree)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].mul_assign(varIndex, degree);
}

void TaylorModelVec::linearTrans(TaylorModelVec & result, const Matrix & A) const
{
	result.clear();
	if(tms.size() != A.cols())
	{
		printf("Dimensions do not coincide.\n");
		return;
	}

	int rows = A.rows();
  result.tms.reserve(rows);
	for(int i=0; i<rows; ++i)
	{
		TaylorModel tm1;

		for(int j=0; j<A.cols(); ++j)
		{
			TaylorModel tm2;
			Interval I( A.get(i,j) );
			tms[j].mul(tm2, I);
			tm1.add_assign(tm2);
		}

		result.tms.push_back(tm1);
	}
}

void TaylorModelVec::linearTrans_assign(const Matrix & A)
{
	TaylorModelVec result;
	linearTrans(result, A);
	*this = result;
}

void TaylorModelVec::rmConstant()
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].rmConstant();
}

void TaylorModelVec::cutoff_normal(const vector<Interval> & step_exp_table, const Interval & cutoff_threshold)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].cutoff_normal(step_exp_table, cutoff_threshold);
}

void TaylorModelVec::cutoff(const vector<Interval> & domain, const Interval & cutoff_threshold)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].cutoff(domain, cutoff_threshold);
}


void TaylorModelVec::cutoff(const Interval & cutoff_threshold)
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].cutoff(cutoff_threshold);
}

void TaylorModelVec::center_nc()
{
	for(int i=0; i<tms.size(); ++i)
		tms[i].center_nc();
}

void TaylorModelVec::Picard_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
  /*
  mreset(old);
  mlog1(ode[0].toString());
  mrestore(old);
  */
	if(order <= 1)
	{
		for(int i=0; i<ode.size(); ++i)
		{
			TaylorModel tmTemp;
			ode[i].insert_no_remainder(tmTemp, *this, numVars, 0, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
	}
	else
	{
		for(int i=0; i<ode.size(); ++i)
		{
			TaylorModel tmTemp;
			ode[i].insert_no_remainder(tmTemp, *this, numVars, order-1, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral_no_remainder(tmvTemp2);

	x0.add(result, tmvTemp2);
}


//find next iteration of with Picard operator (no remainder version)
//if variables is not to be solved in component, then it's just copied over
//otherwise it's substituted into odes, integrated and added it's intial 
//condition (from component.initSet)
void TaylorModelVec::Picard_no_remainder(TaylorModelVec & result, 
    MyComponent *comp, const int numVars, const int order, 
    const Interval & cutoff_threshold) const {
  mreset(old);
  mdisable();
  //mlog1(sbuilder() << "order: " << order);
  
	TaylorModelVec tmvTemp;
	vector<int> & sIndexes = comp->solveIndexes;
	int lowerOrder = max(0, order - 1);
	
	//find picard iteration for variables that need to be solved
  TaylorModelVec compPicard(sIndexes.size());
	for(int i=0; i<sIndexes.size(); ++i) {
	  TaylorModel tmTemp;
		comp->odes[sIndexes[i]].insert_no_remainder(tmTemp, *this, numVars, order-1,
		    cutoff_threshold);
		compPicard.tms.push_back(tmTemp);
	}
	//mlog("compPic", compPicard);
	
	TaylorModelVec integrated(tms.size());
	compPicard.integral_no_remainder(integrated);
	
	for(int i = 0; i < tms.size(); i++) {
    //mlog1(sbuilder() << "i: " << i);
    vector<int>::iterator it = find(comp->solveIndexes.begin(), 
        comp->solveIndexes.end(), i);
    //old variable (alrady has solution, just copy it)
    if(it == comp->solveIndexes.end()) {
      //mlog("old", tms[i]);
      result.tms.push_back(tms[i]);
    }else {
      //new variable
      int indexInComp = it - comp->solveIndexes.begin();
      TaylorModel tmTemp;
      //add initial conditon to integration result
  		comp->initSet.tms[i].add(tmTemp, integrated.tms[indexInComp]);
  		//mlog("new", tmTemp);
      result.tms.push_back(tmTemp);
    }
	}
  mrestore(old);
}

void TaylorModelVec::Picard_no_remainder_assign(const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const int order, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_no_remainder(result, x0, ode, numVars, order, cutoff_threshold);
	*this = result;
}


void TaylorModelVec::Picard_no_remainder_assign(MyComponent *component, 
    const int numVars, const int order, const Interval & cutoff_threshold) {
	TaylorModelVec result(tms.size());
	Picard_no_remainder(result, component, numVars, order, cutoff_threshold);
	*this = result;
}

void TaylorModelVec::Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
  mforce1("Picard_ctrunc_normal");

	if(order <= 1)
	{
		for(int i=0; i<ode.size(); ++i)
		{
			TaylorModel tmTemp;
			ode[i].insert_ctrunc_normal(tmTemp, *this, polyRange, step_exp_table, numVars, 0, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
	}
	else
	{
		for(int i=0; i<ode.size(); ++i)
		{
			TaylorModel tmTemp;
			ode[i].insert_ctrunc_normal(tmTemp, *this, polyRange, step_exp_table, numVars, order-1, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral(tmvTemp2, step_exp_table[1]);

	x0.add(result, tmvTemp2);
}

void TaylorModelVec::Picard_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_ctrunc_normal(result, x0, polyRange, ode, step_exp_table, numVars, order, cutoff_threshold);
	*this = result;
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

void TaylorModelVec::Picard_ctrunc_normal(TaylorModelVec & result, 
    vector<RangeTree *> & trees, MyComponent *comp2, MySettings *settings,
    const vector<Interval> & polyRange, const vector<Interval> & step_exp_table, 
    const int numVars, const Interval & cutoff_threshold) const {
  mreset(old);
  mdisable();
  minc();
  mlog1("Picard_ctrunc_normal (2, i vec)");
  //mlog("x0", comp2->initSet);
  //mlog("ode", comp2->odes);
  //mlog("this", *this);
  //mlog("polyRange", polyRange);
  //mlog1(sbuilder() << "order: " << order);
  //mlog1(sbuilder() << "numVars: " << numVars);
  //mlog("comp", comp2->solveIndexes);
	TaylorModelVec inserted(comp2->solveIndexes.size());

	trees.clear();
	for(int i=0; i<comp2->odes.size(); ++i) {
		trees.push_back(NULL);
	}

	for(int j = 0; j < comp2->solveIndexes.size(); j++) {
		int var = comp2->solveIndexes[j]; //local index inside of component
		int globalIndex = comp2->compVars[var]; //global index out of all vars
		int varOrder = settings->orderLookup.at(globalIndex);
		//mforce1(sbuilder() << "var: " << var << ", order: " << varOrder);
		//TODO two branches are the same one if varOrder = max(1, varOrder)
		if(varOrder <= 1) {
			TaylorModel tmTemp;
			//mlog1(comp2->odes[var].toString());
			comp2->odes[var].insert_ctrunc_normal(tmTemp, trees[var], *this, 
          polyRange, step_exp_table, numVars, 0, cutoff_threshold);
			inserted.tms.push_back(tmTemp);
			//mlog(sbuilder() << "tm[" << var << "]", tmTemp);
		} else {
	    mlog1(sbuilder() << "var: " << var);
			mlog1(comp2->odes[var].toString());
      mlog1(sbuilder() << numVars);
		  //varOrder - 1, since the integration makes the degree 1 higher
			TaylorModel tmTemp;
			comp2->odes[var].insert_ctrunc_normal(tmTemp, trees[var], *this, 
          polyRange, step_exp_table, numVars, varOrder-1, cutoff_threshold);
			inserted.tms.push_back(tmTemp);
		}
	}
	//mlog("inserted(TM)", inserted);

	TaylorModelVec integrated(inserted.tms.size());
	for(int i = 0; i < inserted.tms.size(); i++) {
		TaylorModel tmTemp;
		inserted.tms[i].integral(tmTemp, step_exp_table[1]);
		integrated.tms.push_back(tmTemp);
	}
	//mlog("integrated(TM)", integrated);
	TaylorModelVec added(comp2->solveIndexes.size());
	for(int i = 0; i < comp2->solveIndexes.size(); i++) {
	  //i is index in integrated, var is index in original
	  int var = comp2->solveIndexes[i];
		TaylorModel tmTemp;
		comp2->initSet.tms[var].add(tmTemp, integrated.tms[i]);
		added.tms.push_back(tmTemp);
	}
	
	result.clear();
	result.tms.reserve(tms.size());
	for(int i = 0; i < tms.size(); i++) {
    //mlog1(sbuilder() << "i: " << i);
    vector<int>::iterator it = find(comp2->solveIndexes.begin(), comp2->solveIndexes.end(), i);
    if(it == comp2->solveIndexes.end()) {
      result.tms.push_back(tms[i]);
    }else {
      int indexInComp = it - comp2->solveIndexes.begin();
      result.tms.push_back(added.tms[indexInComp]);
    }
	}
	//mlog("result", result);
	//mforce("result", result);
	mdec();
	mrestore(old);
}

void TaylorModelVec::Picard_ctrunc_normal(TaylorModelVec & result, vector<RangeTree *> & trees, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const
{
  //mforce1("Picard_ctrunc_normal3");
	TaylorModelVec tmvTemp;

	trees.clear();
	for(int i=0; i<ode.size(); ++i)
	{
		trees.push_back(NULL);
	}
	for(int i=0; i<ode.size(); ++i)
	{
		if(orders[i] <= 1)
		{
			TaylorModel tmTemp;
			ode[i].insert_ctrunc_normal(tmTemp, trees[i], *this, polyRange, step_exp_table, numVars, 0, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
		else
		{
			TaylorModel tmTemp;
			ode[i].insert_ctrunc_normal(tmTemp, trees[i], *this, polyRange, step_exp_table, numVars, orders[i]-1, cutoff_threshold);
			tmvTemp.tms.push_back(tmTemp);
		}
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral(tmvTemp2, step_exp_table[1]);

	x0.add(result, tmvTemp2);
}

void TaylorModelVec::Picard_only_remainder(vector<Interval> & result, vector<RangeTree *> & trees, const TaylorModelVec & x0, const vector<HornerForm> & ode, const Interval & timeStep) const
{
  //mforce1("my version shoudln't call this (remainder)3");
  //exit(0);
	result.clear();

	for(int i=0; i<ode.size(); ++i)
	{
		Interval intTemp;
		ode[i].insert_only_remainder(intTemp, trees[i], *this, timeStep);
		intTemp *= timeStep;
		result.push_back(intTemp);
	}
}
void TaylorModelVec::Picard_only_remainder(vector<Interval> & result, 
    vector<RangeTree *> & trees, MyComponent *comp, 
    const Interval & timeStep) const {
	result.clear();
	mreset(old);
	mdisable();
	//mlog("this", *this);
	
	for(int i = 0; i < tms.size(); i++) {
    //mlog1(sbuilder() << "i: " << i);
    vector<int>::iterator it = find(comp->solveIndexes.begin(), comp->solveIndexes.end(), i);
    if(it == comp->solveIndexes.end()) {
      //mlog1("copying");
      result.push_back(tms[i].remainder);
    }else {
      //mlog1("computing");
      Interval intTemp;
		  comp->odes[i].insert_only_remainder(intTemp, trees[i], *this, timeStep);
		  intTemp *= timeStep; //multiply by time interval (integration)
      intTemp += comp->initSet.tms[0].remainder; //add the initial condition rem
		  result.push_back(intTemp);
    }
	}
	//mlog("ret", result);
	mrestore(old);
}

void TaylorModelVec::Picard_update_remainder(vector<RangeTree *> & trees, 
  	  MyComponent *comp, const Interval & timeStep) {
  vector<Interval> result;
  Picard_only_remainder(result, trees, comp, timeStep);
  for(int i = 0; i < tms.size(); i++) {
    tms[i].remainder = result[i];
  }
}



void TaylorModelVec::Picard_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold) const
{
	result = *this;

	for(int i=0; i<ode.size(); ++i)
	{
		if(bIncreased[i])
		{
			TaylorModel tmTemp;
			if(orders[i] <= 1)
			{
				ode[i].insert_no_remainder(tmTemp, *this, numVars, 0, cutoff_threshold);
			}
			else
			{
				ode[i].insert_no_remainder(tmTemp, *this, numVars, orders[i]-1, cutoff_threshold);
			}

			TaylorModel tmTemp2;
			tmTemp.integral_no_remainder(tmTemp2);
			x0.tms[i].add(result.tms[i], tmTemp2);
		}
	}
}

void TaylorModelVec::Picard_no_remainder_assign(const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_no_remainder(result, x0, ode, numVars, orders, bIncreased, cutoff_threshold);
	*this = result;
}

void TaylorModelVec::Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;

	for(int i=0; i<ode.size(); ++i)
	{
		TaylorModel tmTemp;
		if(orders[i] <= 1)
		{
			ode[i].insert_ctrunc_normal(tmTemp, *this, polyRange, step_exp_table, numVars, 0, cutoff_threshold);
		}
		else
		{
			ode[i].insert_ctrunc_normal(tmTemp, *this, polyRange, step_exp_table, numVars, orders[i]-1, cutoff_threshold);
		}
		tmvTemp.tms.push_back(tmTemp);
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral(tmvTemp2, step_exp_table[1]);

	x0.add(result, tmvTemp2);
}

void TaylorModelVec::Picard_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_ctrunc_normal(result, x0, polyRange, ode, step_exp_table, numVars, orders, cutoff_threshold);
	*this = result;
}





// using Taylor approximation
void TaylorModelVec::Picard_non_polynomial_taylor_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
	Interval intZero;
	int rangeDim = strOde.size();

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	parseSetting.clear();

	parseSetting.cutoff_threshold = cutoff_threshold;
	parseSetting.flowpipe = *this;

	if(order <= 1)
	{
		parseSetting.order = 0;
	}
	else
	{
		parseSetting.order = order-1;
	}

	for(int i=0; i<strOde.size(); ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;

		//parseODE();		// call the parser //put back for flowstar

		TaylorModel tmTemp(parseResult.expansion, intZero);

		tmvTemp.tms.push_back(tmTemp);
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral_no_remainder(tmvTemp2);

	x0.add(result, tmvTemp2);
}

void TaylorModelVec::Picard_non_polynomial_taylor_no_remainder_assign(const TaylorModelVec & x0, const vector<string> & strOde, const int order, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_non_polynomial_taylor_no_remainder(result, x0, strOde, order, cutoff_threshold);
	*this = result;
}

void TaylorModelVec::Picard_non_polynomial_taylor_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold) const
{
	result = *this;

	Interval intZero;
	int rangeDim = strOde.size();

	string prefix(str_prefix_taylor_polynomial);
	string suffix(str_suffix);

	parseSetting.clear();

	parseSetting.cutoff_threshold = cutoff_threshold;
	parseSetting.flowpipe = *this;

	for(int i=0; i<strOde.size(); ++i)
	{
		if(bIncreased[i])
		{
			if(orders[i] <= 1)
			{
				parseSetting.order = 0;
			}
			else
			{
				parseSetting.order = orders[i]-1;
			}

			parseSetting.strODE = prefix + strOde[i] + suffix;

			//parseODE();		// call the parser

			TaylorModel tmTemp(parseResult.expansion, intZero);

			TaylorModel tmTemp2;
			tmTemp.integral_no_remainder(tmTemp2);
			x0.tms[i].add(result.tms[i], tmTemp2);
		}
	}
}

void TaylorModelVec::Picard_non_polynomial_taylor_no_remainder_assign(const TaylorModelVec & x0, const vector<string> & strOde, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_non_polynomial_taylor_no_remainder(result, x0, strOde, orders, bIncreased, cutoff_threshold);
	*this = result;
}

void TaylorModelVec::Picard_non_polynomial_taylor_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
	Interval intZero;
	int rangeDim = strOde.size();

	string prefix(str_prefix_taylor_picard);
	string suffix(str_suffix);

	parseSetting.clear();

	parseSetting.flowpipe = *this;
	parseSetting.step_exp_table = step_exp_table;
	parseSetting.cutoff_threshold = cutoff_threshold;

	if(order <= 1)
	{
		parseSetting.order = 0;
	}
	else
	{
		parseSetting.order = order-1;
	}

	for(int i=0; i<strOde.size(); ++i)
	{
		parseSetting.strODE = prefix + strOde[i] + suffix;

		//parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, parseResult.remainder);

		tmvTemp.tms.push_back(tmTemp);
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral(tmvTemp2, step_exp_table[1]);

	x0.add(result, tmvTemp2);
}

void TaylorModelVec::Picard_non_polynomial_taylor_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_non_polynomial_taylor_ctrunc_normal(result, x0, strOde, step_exp_table, order, cutoff_threshold);
	*this = result;
}

void TaylorModelVec::Picard_non_polynomial_taylor_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const vector<int> & orders, const Interval & cutoff_threshold) const
{
	TaylorModelVec tmvTemp;
	Interval intZero;
	int rangeDim = strOde.size();

	string prefix(str_prefix_taylor_picard);
	string suffix(str_suffix);

	parseSetting.clear();

	parseSetting.flowpipe = *this;
	parseSetting.step_exp_table = step_exp_table;
	parseSetting.cutoff_threshold = cutoff_threshold;

	for(int i=0; i<strOde.size(); ++i)
	{
		if(orders[i] <= 1)
		{
			parseSetting.order = 0;
		}
		else
		{
			parseSetting.order = orders[i]-1;
		}

		parseSetting.strODE = prefix + strOde[i] + suffix;

		//parseODE();		// call the parser

		TaylorModel tmTemp(parseResult.expansion, parseResult.remainder);

		tmvTemp.tms.push_back(tmTemp);
	}

	TaylorModelVec tmvTemp2;
	tmvTemp.integral(tmvTemp2, step_exp_table[1]);

	x0.add(result, tmvTemp2);
}

void TaylorModelVec::Picard_non_polynomial_taylor_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const vector<int> & orders, const Interval & cutoff_threshold)
{
	TaylorModelVec result;
	Picard_non_polynomial_taylor_ctrunc_normal(result, x0, strOde, step_exp_table, orders, cutoff_threshold);
	*this = result;
}

void TaylorModelVec::Picard_non_polynomial_taylor_only_remainder(vector<Interval> & result, const TaylorModelVec & x0, const vector<string> & strOde, const Interval & timeStep, const int order) const
{
	result.clear();

	string prefix(str_prefix_taylor_remainder);
	string suffix(str_suffix);

	parseSetting.flowpipe = *this;
	parseSetting.iterRange = parseSetting.ranges.begin();

	if(order <= 1)
	{
		parseSetting.order = 0;
	}
	else
	{
		parseSetting.order = order-1;
	}

	for(int i=0; i<strOde.size(); ++i)
	{
		Interval intTemp;

		parseSetting.strODE = prefix + strOde[i] + suffix;

		//parseODE();		// call the parser

		intTemp = parseResult.remainder * timeStep;
		result.push_back(intTemp);
	}
}

void TaylorModelVec::Picard_non_polynomial_taylor_only_remainder(vector<Interval> & result, const TaylorModelVec & x0, const vector<string> & strOde, const Interval & timeStep, const vector<int> & orders) const
{
	result.clear();

	string prefix(str_prefix_taylor_remainder);
	string suffix(str_suffix);

	parseSetting.flowpipe = *this;
	parseSetting.iterRange = parseSetting.ranges.begin();

	for(int i=0; i<strOde.size(); ++i)
	{
		if(orders[i] <= 1)
		{
			parseSetting.order = 0;
		}
		else
		{
			parseSetting.order = orders[i]-1;
		}

		Interval intTemp;

		parseSetting.strODE = prefix + strOde[i] + suffix;

		//parseODE();		// call the parser

		intTemp = parseResult.remainder * timeStep;
		result.push_back(intTemp);
	}
}

void TaylorModelVec::normalize(vector<Interval> & domain)
{
	int domainDim = domain.size();
	int rangeDim = tms.size();

	// compute the center of the original domain and make it origin-centered
	vector<Interval> intVecCenter;
	for(int i=1; i<domainDim; ++i)		// we omit the time dimension
	{
		Interval intTemp;
		domain[i].remove_midpoint(intTemp);
		intVecCenter.push_back(intTemp);
	}

	// compute the scalars
	Interval intZero;
	vector<vector<Interval> > coefficients;
	vector<Interval> row;

	for(int i=0; i<domainDim; ++i)
	{
		row.push_back(intZero);
	}

	for(int i=0; i<domainDim-1; ++i)
	{
		coefficients.push_back(row);
	}

	for(int i=1; i<domainDim; ++i)
	{
		Interval M;
		domain[i].mag(M);
		coefficients[i-1][i] = M;
	}

	TaylorModelVec newVars(coefficients);
	for(int i=0; i<domainDim-1; ++i)
	{
		TaylorModel tmTemp(intVecCenter[i], domainDim);
		newVars.tms[i].add_assign(tmTemp);
	}

	Interval intUnit(-1,1);
	for(int i=1; i<domainDim; ++i)
	{
		domain[i] = intUnit;
	}

	for(int i=0; i<tms.size(); ++i)
	{
		TaylorModel tmTemp;
		tms[i].insert_no_remainder_no_cutoff(tmTemp, newVars, domainDim, tms[i].degree());
		tms[i].expansion = tmTemp.expansion;
	}
}

void TaylorModelVec::polyRange(vector<Interval> & result, const vector<Interval> & domain) const
{
	result.clear();

	for(int i=0; i<tms.size(); ++i)
	{
		Interval intTemp;
		tms[i].polyRange(intTemp, domain);
		result.push_back(intTemp);
	}
}

void TaylorModelVec::polyRangeNormal(vector<Interval> & result, const vector<Interval> & step_exp_table) const
{
	result.clear();
  result.reserve(tms.size());

	for(int i=0; i<tms.size(); ++i)
	{
		Interval intTemp;
		tms[i].polyRangeNormal(intTemp, step_exp_table);
		result.push_back(intTemp);
	}
}

TaylorModelVec & TaylorModelVec::operator = (const TaylorModelVec & tmv)
{
	if(this == &tmv)
		return *this;

	tms = tmv.tms;
	return *this;
}




















ParseSetting::ParseSetting()
{
	iterRange = ranges.begin();
	order = 0;
}

ParseSetting::ParseSetting(const ParseSetting & setting)
{
	strODE = setting.strODE;

	ranges = setting.ranges;
	iterRange = setting.iterRange;

	cutoff_threshold = setting.cutoff_threshold;

	step_exp_table = setting.step_exp_table;
	flowpipe = setting.flowpipe;
	order = setting.order;
}

ParseSetting::~ParseSetting()
{
	ranges.clear();
	step_exp_table.clear();
}

void ParseSetting::clear()
{
	ranges.clear();
	iterRange = ranges.begin();

	step_exp_table.clear();

	Interval intZero;
	cutoff_threshold = intZero;
	variables = vector<string>();
}

void ParseSetting::addVar(string v) {
  variables.push_back(v);
}

ParseSetting & ParseSetting::operator = (const ParseSetting & setting)
{
	if(this == &setting)
		return *this;

	strODE = setting.strODE;

	ranges = setting.ranges;
	iterRange = setting.iterRange;

	cutoff_threshold = setting.cutoff_threshold;

	step_exp_table = setting.step_exp_table;
	flowpipe = setting.flowpipe;
	order = setting.order;

	return *this;
}































ParseResult::ParseResult()
{
}

ParseResult::ParseResult(const ParseResult & result)
{
	expansion		= result.expansion;
	remainder		= result.remainder;
	strExpansion	= result.strExpansion;
}

ParseResult::~ParseResult()
{
}

ParseResult & ParseResult::operator = (const ParseResult & result)
{
	if(this == &result)
		return *this;

	expansion		= result.expansion;
	remainder		= result.remainder;
	strExpansion	= result.strExpansion;

	return *this;
}
















// class Interval_matrix

Interval_matrix::Interval_matrix()
{
}

Interval_matrix::Interval_matrix(const vector<vector<Interval> > & elements_input)
{
	elements = elements_input;
}

Interval_matrix::Interval_matrix(const vector<Interval> & intVector)
{
	for(int i=0; i<intVector.size(); ++i)
	{
		vector<Interval> intVecTemp;

		intVecTemp.push_back(intVector[i]);
		elements.push_back(intVecTemp);
	}
}

Interval_matrix::Interval_matrix(const int m, const int n)
{
	Interval intZero;

	vector<Interval> intVecTemp;
	for(int i=0; i<n; ++i)
	{
		intVecTemp.push_back(intZero);
	}

	for(int i=0; i<m; ++i)
	{
		elements.push_back(intVecTemp);
	}
}

Interval_matrix::Interval_matrix(const Interval_matrix & A)
{
	elements = A.elements;
}

Interval_matrix::~Interval_matrix()
{
}

int Interval_matrix::rows() const
{
	return (int)(elements.size());
}

int Interval_matrix::cols() const
{
	if(elements.size() > 0)
	{
		return (int)(elements[0].size());
	}
	else
	{
		return 0;
	}
}

void Interval_matrix::get(Interval & I, const int i, const int j) const
{
	I = elements[i][j];
}

Interval Interval_matrix::get(const int i, const int j) const
{
	return elements[i][j];
}

void Interval_matrix::set(const Interval & I, const int i, const int j)
{
	elements[i][j] = I;
}

void Interval_matrix::mul(Interval_matrix & result, const double v) const
{
	result = *this;

	for(int i=0; i<result.elements.size(); ++i)
	{
		for(int j=0; j<result.elements[i].size(); ++j)
		{
			result.elements[i][j].mul_assign(v);
		}
	}
}

void Interval_matrix::mul_assign(const double v)
{
	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			elements[i][j].mul_assign(v);
		}
	}
}

void Interval_matrix::mul(Interval_matrix & result, const Interval & I) const
{
	result = *this;

	for(int i=0; i<result.elements.size(); ++i)
	{
		for(int j=0; j<result.elements[i].size(); ++j)
		{
			result.elements[i][j] *= I;
		}
	}
}

void Interval_matrix::mul_assign(const Interval & I)
{
	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			elements[i][j] *= I;
		}
	}
}

void Interval_matrix::div(Interval_matrix & result, const double v) const
{
	result = *this;

	for(int i=0; i<result.elements.size(); ++i)
	{
		for(int j=0; j<result.elements[i].size(); ++j)
		{
			result.elements[i][j].div_assign(v);
		}
	}
}

void Interval_matrix::div_assign(const double v)
{
	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			elements[i][j].div_assign(v);
		}
	}
}

void Interval_matrix::decompose(Interval_matrix & det, Interval_matrix & nondet) const
{
	det = nondet = *this;

	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			nondet.elements[i][j].remove_midpoint(det.elements[i][j]);
		}
	}
}

void Interval_matrix::pow(Interval_matrix & result, const int degree) const
{
	Interval_matrix temp = *this;
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

void Interval_matrix::pow_assign(const int degree)
{
	Interval_matrix temp = *this;
	Interval_matrix result = *this;

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

double Interval_matrix::max_norm() const
{
	double max = 0;

	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			double tmp = elements[i][j].mag();

			if(tmp > max)
			{
				max = tmp;
			}
		}
	}

	return max;
}

void Interval_matrix::linear_trans(vector<Polynomial> & result, const vector<Polynomial> & polynomial) const
{
	if(polynomial.size() != elements[0].size())
	{
		printf("Dimensions do not match.\n");
		exit(1);
	}

	result.clear();

	for(int i=0; i<elements.size(); ++i)
	{
		Polynomial polyTemp1;

		for(int j=0; j<elements[i].size(); ++j)
		{
			Polynomial polyTemp2;

			polynomial[j].mul(polyTemp2, elements[i][j]);

			polyTemp1 += polyTemp2;
		}

		result.push_back(polyTemp1);
	}
}

void Interval_matrix::transpose(Interval_matrix & result) const
{
	int rows = elements.size();
	int cols = elements[0].size();

	Interval_matrix imTemp(cols, rows);

	for(int i=0; i<rows; ++i)
	{
		for(int j=0; j<cols; ++j)
		{
			imTemp.elements[j][i] = elements[i][j];
		}
	}

	result = imTemp;
}

void Interval_matrix::output(FILE *fp) const
{
	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			elements[i][j].dump(fp);
			fprintf(fp, "\t");
		}

		fprintf(fp, "\n");
	}
}

Interval_matrix & Interval_matrix::operator += (const Interval_matrix & A)
{
	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			elements[i][j] += A.elements[i][j];
		}
	}

	return *this;
}

Interval_matrix & Interval_matrix::operator -= (const Interval_matrix & A)
{
	for(int i=0; i<elements.size(); ++i)
	{
		for(int j=0; j<elements[i].size(); ++j)
		{
			elements[i][j] -= A.elements[i][j];
		}
	}

	return *this;
}

Interval_matrix & Interval_matrix::operator *= (const Interval_matrix & A)
{
	int m = elements.size();
	int n = A.elements[0].size();
	int k = A.elements.size();

	Interval_matrix result(m,n);

	for(int i=0; i<m; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			Interval tmp;

			for(int p=0; p<k; ++p)
			{
				tmp += elements[i][p] * A.elements[p][j];
			}

			result.elements[i][j] = tmp;
		}
	}

	*this = result;
	return *this;
}

Interval_matrix Interval_matrix::operator + (const Interval_matrix & A) const
{
	Interval_matrix result = *this;

	result += A;

	return result;
}

Interval_matrix Interval_matrix::operator - (const Interval_matrix & A) const
{
	Interval_matrix result = *this;

	result -= A;

	return result;
}

Interval_matrix Interval_matrix::operator * (const Interval_matrix & A) const
{
	int m = elements.size();
	int n = A.elements[0].size();
	int k = A.elements.size();

	Interval_matrix result(m,n);

	for(int i=0; i<m; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			Interval tmp;

			for(int p=0; p<k; ++p)
			{
				tmp += elements[i][p] * A.elements[p][j];
			}

			result.elements[i][j] = tmp;
		}
	}

	return result;
}

Interval_matrix & Interval_matrix::operator = (const Interval_matrix & A)
{
	if(this == &A)
		return *this;

	elements = A.elements;

	return *this;
}

























void exp_taylor_remainder(Interval & result, const Interval & tmRange, const int order)
{
	Interval intProd = tmRange.pow(order);

	Interval J(0,1);
	J *= tmRange;
	J.exp_assign();

	result = factorial_rec[order] * intProd * J;
}

void rec_taylor_remainder(Interval & result, const Interval & tmRange, const int order)
{
	Interval J(0,1), intOne(1), intMOne(-1);
	J *= tmRange;
	J += intOne;
	J.rec_assign();

	Interval intProd = J;
	intProd *= tmRange;
	intProd *= intMOne;

	result = intProd.pow(order);
	result *= J;
}

void sin_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const int order)
{
	Interval intProd = tmRange.pow(order);

	Interval J(0,1);
	J *= tmRange;
	J += C;

	int k = order % 4;

	switch(k)
	{
	case 0:
		J.sin_assign();
		break;
	case 1:
		J.cos_assign();
		break;
	case 2:
		J.sin_assign();
		J.inv_assign();
		break;
	case 3:
		J.cos_assign();
		J.inv_assign();
		break;
	}

	result = factorial_rec[order] * intProd * J;
}

void cos_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const int order)
{
	Interval intProd = tmRange.pow(order);

	Interval J(0,1);
	J *= tmRange;
	J += C;

	int k = order % 4;

	switch(k)
	{
	case 0:
		J.cos_assign();
		break;
	case 1:
		J.sin_assign();
		J.inv_assign();
		break;
	case 2:
		J.cos_assign();
		J.inv_assign();
		break;
	case 3:
		J.sin_assign();
		break;
	}

	result = factorial_rec[order] * intProd * J;
}

void log_taylor_remainder(Interval & result, const Interval & tmRange, const int order)
{
	Interval J(0,1);
	J *= tmRange;
	J.add_assign(1.0);
	J.rec_assign();

	Interval I = tmRange;
	I *= J;

	result = I.pow(order);

	result.div_assign((double)order);

	if((order+1)%2 == 1)		// order+1 is odd
	{
		result.inv_assign();
	}
}

void sqrt_taylor_remainder(Interval & result, const Interval & tmRange, const int order)
{
	Interval I(0,1);
	I *= tmRange;
	I.add_assign(1.0);
	I.rec_assign();

	Interval intTemp;
	I.sqrt(intTemp);

	I *= tmRange;
	I.div_assign(2.0);

	Interval intProd = I.pow(order-1);

	intProd /= intTemp;
	intProd *= tmRange;
	intProd.div_assign(2.0);

	result = double_factorial[2*order-3] * factorial_rec[order] * intProd;

	if(order % 2 == 0)
	{
		result.inv_assign();
	}
}



















void exp_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order)
{
	Interval intZero;
	result = intZero;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	for(int i=order; i>0; --i)
	{
		Interval intFactor(1);
		intFactor.div_assign((double)i);

		result *= intFactor;

		Interval intTemp;
		intTemp = (*iterRange) * remainder;		// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;		// P2 x I1
		intTemp += remainder * result;			// I2 x I1
		++iterRange;
		intTemp += (*iterRange);				// truncation
		++iterRange;

		result = intTemp;
	}

	result *= const_part;

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval tmRange = (*iterRange) + remainder;
	++iterRange;

	Interval rem;
	exp_taylor_remainder(rem, tmRange, order+1);
	result += const_part * rem;
}

void rec_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order)
{
	Interval intZero;
	result = intZero;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	Interval tmF_c_remainder = remainder * const_part;

	for(int i=order; i>0; --i)
	{
		result.inv_assign();

		Interval intTemp;
		intTemp = (*iterRange) * tmF_c_remainder;	// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;			// P2 x I1
		intTemp += tmF_c_remainder * result;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);					// truncation
		++iterRange;

		result = intTemp;
	}

	result *= const_part;

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval rem, tmF_cRange;
	tmF_cRange = (*iterRange) + tmF_c_remainder;
	++iterRange;

	rec_taylor_remainder(rem, tmF_cRange, order+1);

	result += rem * const_part;
}

void sin_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order)
{
	Interval intZero;
	result = intZero;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	Interval tmPowerTmF_remainder;

	for(int i=1; i<=order; ++i)
	{
		Interval intTemp;
		intTemp = (*iterRange) * remainder;					// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * tmPowerTmF_remainder;		// P2 x I1
		intTemp += remainder * tmPowerTmF_remainder;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);							// truncation
		++iterRange;

		tmPowerTmF_remainder = intTemp;

		Interval intTemp2 = tmPowerTmF_remainder;

		intTemp2 *= (*iterRange);
		++iterRange;

		result += intTemp2;
	}

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval tmRange, rem;
	tmRange = (*iterRange) + remainder;
	++iterRange;

	sin_taylor_remainder(rem, const_part, tmRange, order+1);

	result += rem;
}

void cos_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order)
{
	Interval intZero;
	result = intZero;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	int k=1;

	Interval tmPowerTmF_remainder;

	for(int i=1; i<=order; ++i, ++k)
	{
		Interval intTemp;
		intTemp = (*iterRange) * remainder;					// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * tmPowerTmF_remainder;		// P2 x I1
		intTemp += remainder * tmPowerTmF_remainder;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);							// truncation
		++iterRange;

		tmPowerTmF_remainder = intTemp;

		Interval intTemp2 = tmPowerTmF_remainder;

		intTemp2 *= (*iterRange);
		++iterRange;

		result += intTemp2;
	}

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval tmRange, rem;
	tmRange = (*iterRange) + remainder;
	++iterRange;

	cos_taylor_remainder(rem, const_part, tmRange, order+1);

	result += rem;
}

void log_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order)
{
	Interval intZero;
	result = intZero;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval C = *iterRange;
	++iterRange;

	Interval const_part = C;
	const_part.log_assign();

	Interval tmF_c_remainder = remainder / C;

	result = tmF_c_remainder;
	result.div_assign((double)order);

	for(int i=order; i>=2; --i)
	{
		result.inv_assign();

		Interval intTemp;
		intTemp = (*iterRange) * tmF_c_remainder;	// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;			// P2 x I1
		intTemp += tmF_c_remainder * result;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);					// truncation
		++iterRange;

		result = intTemp;
	}

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval rem, tmF_cRange;
	tmF_cRange = (*iterRange) + tmF_c_remainder;
	++iterRange;

	log_taylor_remainder(rem, tmF_cRange, order+1);

	result += rem;
}

void sqrt_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order)
{
	Interval intZero;
	result = intZero;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval C = *iterRange;
	++iterRange;

	Interval const_part = C;
	const_part.sqrt_assign();

	Interval intTwo(2);
	Interval tmF_2c_remainder = (remainder / C) / intTwo;

	result = tmF_2c_remainder;

	Interval K(1), J(1);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		// i
		Interval K((double)i);

		// j = 2*i-3
		Interval J((double)j);

		result.inv_assign();
		result *= J / K;

		Interval intTemp;
		intTemp = (*iterRange) * tmF_2c_remainder;	// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;			// P2 x I1
		intTemp += tmF_2c_remainder * result;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);					// truncation
		++iterRange;

		result = intTemp;
	}

	result *= const_part;

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval rem, tmF_cRange;
	tmF_cRange = (*iterRange);
	++iterRange;

	tmF_cRange += tmF_2c_remainder * intTwo;

	sqrt_taylor_remainder(rem, tmF_cRange, order+1);

	result += rem * const_part;
}

string TaylorModel::toString(const vector<string> & varNames) const {
  //Polynomial expansion;	// Taylor expansion
	//Interval remainder;
  return sbuilder() << expansion.toString(varNames) << " + " << remainder.toString();
}

string TaylorModelVec::toString(const vector<string> & varNames) const {
  sbuilder sb;
  for (unsigned i=0; i<tms.size(); i++) {
		sb << tms.at(i).toString(varNames) << "\n";
	}
  return sb;
}

string TaylorModel::toMathematicaString() const {
  stringstream ss;
  ss << "{" << expansion.toMathematicaString() << 
      "," << remainder.toMathematicaString() << "}";
  
  string ret = ss.str();
  
  //replace scientific notation
  boost::replace_all(ret, "e", " * 10^");
  
  
  return ret;
}
string TaylorModelVec::toMathematicaString() const {
  sbuilder sb;
  sb << "{";
  
  for (int i=0; i < tms.size() - 1; i++) {
		sb << tms[i].toMathematicaString() << ",";
	}
  sb << tms[tms.size()-1].toMathematicaString();
  sb << "}";
  return sb;
}

TaylorModel TaylorModel::transform(map<int, int> lookup, int size) {
  //mlog1("transform");
  list<Monomial> l = expansion.monomials;
  list<Monomial> ms;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    Monomial m = (*it).transform(lookup, size);
    ms.push_back(m);
    //mlog1(m->toString(getVNames(10)));
    //mlog1((*it).toString(getVNames(10)));
  }
  Polynomial p(ms);
  TaylorModel* tm = new TaylorModel(p, remainder);
  //mlog("tm2", *tm);
  return *tm;
}

TaylorModelVec TaylorModelVec::transform(vector<int> indexes) {
  vector<TaylorModel> temp;
  for(vector<TaylorModel>::iterator it = tms.begin(); it < tms.end(); it++) {
    minc();
    temp.push_back((*it).transform(indexes));
    mdec();
  }
  return TaylorModelVec(temp);
}

TaylorModelVec TaylorModelVec::addNParams(int n) const {
//  mlog1("adding tmv");
  vector<TaylorModel> retTms;
  for(vector<TaylorModel>::const_iterator it = tms.begin(); 
      it < tms.end(); it++) {
    retTms.push_back(it->addNParams(n));
  }
  return TaylorModelVec(retTms);
}

TaylorModel TaylorModel::addNParams(int n) const{
//  mlog1("adding tm");
  Polynomial retExpansion = expansion.addNVariables(n);
  return TaylorModel(retExpansion, remainder);
}

bool TaylorModelVec::isClose(const TaylorModelVec & tmv, double d) const {
  //mlog1("comparing tmv");
  if(tms.size() != tmv.tms.size())
    return false;
  for(int i = 0; i < tms.size(); i++) {
    if(false == tms.at(i).isClose(tmv.tms.at(i), d))
      return false;
  }
  return true;
}

bool TaylorModel::isClose(const TaylorModel & tm, double d) const {
  mreset(old);
  mdisable();
  if(getParamCount() != tm.getParamCount()) {
    mlog1("not same param count");
    mrestore(old);
    return false;
  }
  TaylorModel subbed;
	sub(subbed, tm);
	vector<Interval> domain;
	for(int i = 0; i < getParamCount(); i++) {
	  domain.push_back(Interval(-1,1));
	}
	Interval result;
	subbed.polyRange(result, domain);
	//mlog1(result.mag() > d);
	if(result.mag() > d) {
    mlog1("poly diff larger");
    mlog1(result.toString());
    mrestore(old);
	  return false;
	}
	
	bool ret = remainder.isClose(tm.remainder, d);
	mrestore(old);
	return ret;
}

TaylorModel TaylorModel::transform(vector<int> indexes) {
  mreset(old);
  mdisable();
  //mlog1("transform");
  //mlog("tm", *this);
  minc();
  list<Monomial> l = expansion.monomials;
  list<Monomial> ms;
  //mlog1(sbuilder() << "l: " <<  l.size());
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); it++) {
    //mlog1("before");
    //mlog1(it->toString());
    //mlog1("mid1");
    //mlog1(it->getVariableCount());
    //mlog1("mid2");
    Monomial m = it->transform(indexes);
    //mlog1("mid3");
    //mlog1(sbuilder() << "res:" << m.toString());
    //Monomial rm = (*it).rtransform(indexes);
    //mlog1(sbuilder() << "rm: " << &rm);
    //mlog1("after");
    ms.push_back(m);
  }
  Polynomial p(ms);
  TaylorModel* tm = new TaylorModel(p, remainder);
  mdec();
  //mlog1("done");
  mrestore(old);
  return *tm;
}


TaylorModel* TaylorModel::ptransform(vector<int> indexes) {
  //mlog1("transform");
  list<Monomial> l = expansion.monomials;
  list<Monomial> ms;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    //Monomial m = (*it).transform(indexes);
    //ms.push_back(m);
    
    Monomial *pm = (*it).ptransform(indexes);
    //mlog1(sbuilder() << pm);
    ms.push_back(*pm);
  }
  Polynomial *p = new Polynomial(ms);
  //mlog1(sbuilder() << "&p1: " << p);
  TaylorModel* tm = new TaylorModel(*p, remainder);
  //mlog1(sbuilder() << "&p2: " << &(tm->expansion));
  //mlog1(sbuilder() << "&tm1: " << tm);
  return tm;
}

TaylorModelVec TaylorModelVec::prepareSecondary(int prefix) const {
  //mlog1("preparing");
  TaylorModelVec *ret = new TaylorModelVec();
  for(int i = 0; i < prefix; i++) {
    TaylorModel tm;
    ret->tms.push_back(tm);
  }
  for(int i = 0; i < tms.size(); i++) {
    TaylorModel tm = tms.at(i).prepareSecondary(prefix);
    //mlog1(prefix);
    //mlog("1", tms.at(i));
    //mlog("2", tm);
    ret->tms.push_back(tm);
  }
  return *ret;
}
TaylorModel TaylorModel::prepareSecondary(int prefix) const {
  //mlog1("preparing tm");
  
  list<Monomial> l = expansion.monomials;
  list<Monomial> ms;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    Monomial m = (*it).prepareSecondary(prefix);
    ms.push_back(m);
    //mlog1((*it).toString(getVNames(10)));
    //mlog1(m.toString(getVNames(10)));
  }
  Polynomial p(ms);
  TaylorModel* tm = new TaylorModel(p, remainder);
  return *tm;
}

int TaylorModel::getParamCount() const {
  int ret = expansion.getVariableCount();
  if(ret == -1) {
    throw ArgumentException("asked paramCount from zero TM");
  }
  return ret;
}
int TaylorModel::getIgnoringParamCount() const {
  int ret = expansion.getVariableCount();
  if(ret == -1) {
    return 0;
  }
  return ret;
}
int TaylorModelVec::getParamCount() const {
  if(tms.size() == 0) {
    throw ArgumentException("no tms to get paramCount from");
  }
  for(int i = 0; i < tms.size() - 1; i++) {
    int temp = tms[i].getIgnoringParamCount();
    if(temp > 0) {
      return temp;
    }
  }
  //also throws exection if last possiblity doesn't have params
  return tms[tms.size() - 1].getParamCount();
}
int TaylorModelVec::getIgnoringParamCount() const {
  if(tms.size() == 0) {
    throw ArgumentException("no tms to get paramCount from");
  }
  return tms[0].getIgnoringParamCount();
}

vector<int> TaylorModel::getParams() const {
  vector<int> ret;
  
  list<Monomial> l = expansion.monomials;
  
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    //mlog1(*it);
    //mlog("mparams", it->getParams());
    vector<int> mparams = it->getParams();
    //mlog("mparams", mparams);
    ret.insert(ret.end(), mparams.begin(), mparams.end());
  }
  
  sort(ret.begin(), ret.end());
  ret.erase( unique( ret.begin(), ret.end() ), ret.end() );
  
  return ret;
}


void TaylorModel::getLinearPart(Interval ret[], vector<int> variables) const {
  //mlog1("lin part");
  //expansion.getLinearPart();
  //mlog1(expansion.toString());
  int dim = variables.size();
  
  list<Monomial> l = expansion.monomials;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    //mlog1("--------------------------------");
    //mlog1(it->coefficient.toString());
    int linearVar = it->linearVariable(); //variable in linear monomial
    if(linearVar > -1) {
      //position of variable in component
      int pos = (find(variables.begin(), variables.end(), linearVar-1) - variables.begin());
      //mlog1(sbuilder() << "var: " << linearVar << ", pos: " << pos);
      ret[pos] = it->coefficient;
      //ret[linearVar - 1] = it->coefficient; //-1 cause time is on index 0
      
    }
  }
}

void TaylorModel::getLinearPart2(vector<int> variables) const {
  //mlog1("lin part2");
  //expansion.getLinearPart();
  //mlog1(expansion.toString());
  int dim = variables.size();
  
  list<Monomial> l = expansion.monomials;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    //mlog1("--------------------------------");
    //mlog1(it->coefficient.toString());
    int linearVar = it->linearVariable(); //variable in linear monomial
    if(linearVar > -1) {
      //position of variable in component
      int pos = (find(variables.begin(), variables.end(), linearVar-1) - variables.begin());
      //mlog1(sbuilder() << "var: " << linearVar << ", pos: " << pos);
      //ret[pos] = it->coefficient;
      //ret[linearVar - 1] = it->coefficient; //-1 cause time is on index 0
      
    }
  }
}

void TaylorModel::getParts(TaylorModel & constant, TaylorModel & linear, 
      TaylorModel & nonLinear, TaylorModel & remainder) const {
  //mlog1("tm parts");
  
  //mlog("tm", *this);
  
  list<Monomial> l = expansion.monomials;
  list<Monomial> cMonos;
  list<Monomial> linMonos;
  list<Monomial> nlMonos;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); ++it) {
    //mlog1("--------------------------------");
    int degree = it->degree();
    
    if(degree == 0) {
      cMonos.push_back(*it);
    }
    if(degree == 1) {
      linMonos.push_back(*it);
    }
    if(degree > 1) {
      nlMonos.push_back(*it);
    }
  }
  Polynomial linP(linMonos);
  
  constant.expansion = Polynomial(cMonos);
  constant.remainder = ZERO_INTERVAL;
  
  linear.expansion = Polynomial(linMonos);
  constant.remainder = ZERO_INTERVAL;
  
  nonLinear.expansion = Polynomial(nlMonos);
  constant.remainder = ZERO_INTERVAL;
  
  remainder.expansion = Polynomial();
  remainder.remainder = this->remainder;
}
void TaylorModelVec::getParts(TaylorModelVec & constant, TaylorModelVec & linear, 
      TaylorModelVec & nonLinear, TaylorModelVec & remainder) const {
  constant.tms.clear();
  linear.tms.clear();
  nonLinear.tms.clear();
  remainder.tms.clear();
  
  for(int i = 0; i < tms.size(); i++) {
    TaylorModel c, lin, nl, rem;
    tms[i].getParts(c, lin, nl, rem);
    
    constant.tms.push_back(c);
    linear.tms.push_back(lin);
    nonLinear.tms.push_back(nl);
    remainder.tms.push_back(rem);
  }
}


bool TaylorModelVec::compare(const TaylorModelVec & tmv, 
      const vector<Interval> & domain) const {
  mreset(old);
  mdisable();
  //mlog1("comparing");
  
  vector<Interval> first;
  vector<Interval> second;
  intEval(first, domain);
	tmv.intEval(second, domain);
	//mlog("first", first);
	//mlog("second", second);
	return subseteq(first, second);
  
  mrestore(old);
}


void TaylorModelVec::removeHighTerms(int order) {
  for(int i = 0; i < tms.size(); i++)
    tms[i].removeHighTerms(order);
}
void TaylorModel::removeHighTerms(int order) {
  //mlog1("removing tm");
  list<Monomial> l = expansion.monomials;
  list<Monomial> retained;
  for (list<Monomial>::iterator it = l.begin(); it != l.end(); it++) {
    if(it->degree() <= order)
      retained.push_back(*it);
  }
  expansion.monomials = retained;
}


vector<Interval> TaylorModelVec::getRemainders() {
  vector<Interval> ret;
  for(int i = 0; i < tms.size(); i++) {
    ret.push_back(tms[i].remainder);
  }
  return ret;
}

void TaylorModel::pushConstantToRemainder() {
  //mlog1("pushing");
  
  //mlog("pre", *this);
  
  Interval c;
  constant(c);
  rmConstant();
  
  //mlog1(c.toString());
  remainder += c;
  //mlog("aft", *this);
}

void TaylorModelVec::pushConstantToRemainder() {
  for(int i = 0; i < tms.size(); i++) {
    tms[i].pushConstantToRemainder();
  }
}


TaylorModel TaylorModel::distance(const TaylorModel & tmv) const {
  list<Monomial> l1;
  list<Monomial> l2;
  {
    const list<Monomial> & l1r = expansion.monomials;
    const list<Monomial> & l2r = tmv.expansion.monomials;
    copy(l1r.begin(), l1r.end(), back_inserter(l1));
    copy(l2r.begin(), l2r.end(), back_inserter(l2));
  }
  l1.sort();
  l2.sort();
  
  
  list<Monomial>::iterator it1 = l1.begin();
  list<Monomial>::iterator it2 = l2.begin();
  
  list<Monomial> ret;
  
  //stop when one list is empty
  while(it1 != l1.end() && it2 != l2.end()) {
    //both iterators are at the same monomial
    if(*it1 == *it2) {
      //mlog1("equal");
      
      //measure coefficent distance
      Interval dis = it1->coefficient.distance(it2->coefficient);
      //mlog1(dis.toString());
      vector<int> & degs = it1->degrees;
      it1++;
      it2++;
      
      //skip if coefficents were the same
      if(dis.subseteq(ZERO_INTERVAL))
        continue;
      
      //store the distance with the same monomial
      ret.push_back(*memCreateMonomial(dis, degs));
      continue;
    }
    
    //store *it1 magnitude and advance iterator
    if(*it1 < *it2) {
      //mlog1("it1 < it2");
      Interval dis = ZERO_INTERVAL.distance(it1->coefficient);
      ret.push_back(*memCreateMonomial(dis, it1->degrees));
      it1++;
      continue;
    }
    
    //store *it2 magnitude and advance iterator
    //mlog1("it1 > it2");
    Interval dis = ZERO_INTERVAL.distance(it2->coefficient);
    ret.push_back(*memCreateMonomial(dis, it2->degrees));
    it2++;
  }
  
  //add all the remainding monomials
  while(it1 != l1.end()) {
    Interval dis = ZERO_INTERVAL.distance(it1->coefficient);
    ret.push_back(*memCreateMonomial(dis, it1->degrees));
    it1++;
  }
  while(it2 != l2.end()) {
    Interval dis = ZERO_INTERVAL.distance(it2->coefficient);
    ret.push_back(*memCreateMonomial(dis, it2->degrees));
    it2++;
  }
  
  //measure remainder term coefficient
  Interval remDis = remainder.distance(tmv.remainder);
  
  return TaylorModel(Polynomial(ret), remDis);
}

TaylorModelVec TaylorModelVec::distance(const TaylorModelVec & tmv) const {
  //mlog1("tmv dist");
  if(tms.size() != tmv.tms.size()) {
    mforce1("different tmv sizes");
    exit(0);
  }
  TaylorModelVec ret;
  for(int i = 0; i < tms.size(); i++) {
    ret.tms.push_back(tms[i].distance(tmv.tms[i]));
  }
  
  return ret;
}

void TaylorModel::centerRemainder() {
  mreset(old);
  mdisable();
  Interval mid;
  //mlog("th", *this);
  remainder.remove_midpoint(mid);
  int numOfParams = getParamCount(); //TODO maybe catch -1 error, and return
  TaylorModel centeredModel(mid, numOfParams);
  this->add_assign(centeredModel);
  //mlog1(mid.toString());
  //mlog("cm", centeredModel);
  //mlog("th", *this);
  //mlog1("centering");
  mrestore(old);
}
void TaylorModelVec::centerRemainder() {
  for(int i = 0; i < tms.size(); i++) {
    tms[i].centerRemainder();
  }
}
void TaylorModelVec::centerRemainder(int paramCount) {
  for(int i = 0; i < tms.size(); i++) {
    tms[i].centerRemainder(paramCount);
  }
}
void TaylorModel::centerRemainder(int paramCount) {
  mreset(old);
  mdisable();
  Interval mid;
  //mlog("th", *this);
  remainder.remove_midpoint(mid);
  TaylorModel centeredModel(mid, paramCount);
  this->add_assign(centeredModel);
  //mlog1(mid.toString());
  //mlog("cm", centeredModel);
  //mlog("th", *this);
  //mlog1("centering");
  mrestore(old);
}

//TODO remove and use TaylorModelVec::centerRemainder(int n)
void TaylorModelVec::preconditionCenterRemainder() {
  for(int i = 0; i < tms.size(); i++) {
    tms[i].centerRemainder(tms.size() + 1);
  }
}
