/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Monomial.h"

Monomial::Monomial()
{
}

Monomial::Monomial(const Interval & I, const vector<int> & degs):coefficient(I), degrees(degs), d(0)
{
	for(int i=0; i<degs.size(); ++i)
	{
		d += degs[i];
	}
}

Monomial::Monomial(const Monomial & monomial): coefficient(monomial.coefficient), degrees(monomial.degrees), d(monomial.d)
{
}

Monomial::Monomial(const Interval & I, const int numVars):d(0)
{
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	coefficient = I;
}

Monomial::~Monomial()
{
	degrees.clear();
}

int Monomial::degree() const
{
	return d;
}

int Monomial::dimension() const
{
	return degrees.size();
}

void Monomial::intEval(Interval & result, const vector<Interval> & domain) const
{
	result = coefficient;

	for(int i=0; i<degrees.size(); ++i)
	{
		Interval tmpI(1,1);
		for(int j=0; j<degrees[i]; ++j)
		{
			tmpI *= domain[i];
		}
		result *= tmpI;
	}
}

void Monomial::intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const
{
	Interval intZero;
	result = intZero;

	if(degrees.size() == 0)
		return;

	result = coefficient;
	result *= step_exp_table[degrees[0]];

	Interval evenInt(0,1), oddInt(-1,1);
	Interval intFactor(1);
	bool bSet = false;

	for(int i=1; i<degrees.size(); ++i)
	{
		if(degrees[i] == 0)			// degree is zero
		{
			continue;
		}
		else if(degrees[i]%2 == 0)	// degree is an even number
		{
			if(!bSet)
			{
				intFactor = evenInt;
				bSet = true;
			}
		}
		else						// degree is an odd number
		{
			intFactor = oddInt;
			break;
		}
	}

	result *= intFactor;
}

void Monomial::inv(Monomial & result) const
{
	result = *this;
	coefficient.inv(result.coefficient);
}

Monomial & Monomial::operator = (const Monomial & monomial)
{
	if(this == &monomial)
		return *this;

	coefficient = monomial.coefficient;
	degrees = monomial.degrees;
	d = monomial.d;

	return *this;
}

Monomial & Monomial::operator += (const Monomial & monomial)
{
	coefficient += monomial.coefficient;
	return *this;
}

Monomial & Monomial::operator *= (const Monomial & monomial)
{
	coefficient *= monomial.coefficient;
	
	#ifdef do_checks
	  //cout << "doing check" << endl;
	  if(degrees.size() != monomial.degrees.size()) {
	    cout << "degrees1.size(): " << degrees.size() << endl;
	    cout << "degrees2.size(): " << monomial.degrees.size() << endl;
	    throw std::runtime_error("multiplying monomials with different var count");
	    exit(0);
	  }
  #endif

	for(int i=0; i<degrees.size(); ++i)
	{
		degrees[i] += monomial.degrees[i];
	}

	d += monomial.d;
	return *this;
}

const Monomial Monomial::operator + (const Monomial & monomial) const
{
	Monomial result = *this;
	result += monomial;
	return result;
}

const Monomial Monomial::operator * (const Monomial & monomial) const
{
	Monomial result = *this;
	result *= monomial;
	return result;
}

bool Monomial::isLinear(int & index) const
{
	if(d == 1)
	{
		for(int i=0; i<degrees.size(); ++i)
		{
			if(degrees[i] == 1)
			{
				index = i;
				return true;
			}
		}
	}

	return false;
}

void Monomial::dump_interval(FILE *fp, const vector<string> & varNames) const
{
	coefficient.dump(fp);

	for(int i=0; i<degrees.size()-1; i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
				fprintf(fp, " * %s", varNames[i].c_str());
			else
				fprintf(fp, " * %s^%d", varNames[i].c_str(), degrees[i]);
		}
	}

	if(degrees[degrees.size()-1] != 0)
	{
		if(degrees[degrees.size()-1] == 1)
			fprintf(fp, " * %s", varNames[degrees.size()-1].c_str());
		else
			fprintf(fp, " * %s^%d", varNames[degrees.size()-1].c_str(), degrees[degrees.size()-1]);
	}
}

void Monomial::serialize(FILE *fp, const vector<string> & varNames) const {
	coefficient.serialize(fp);

	for(int i=0; i<degrees.size()-1; i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
				fprintf(fp, " * %s", varNames[i].c_str());
			else
				fprintf(fp, " * %s^%d", varNames[i].c_str(), degrees[i]);
		}
	}

	if(degrees[degrees.size()-1] != 0)
	{
		if(degrees[degrees.size()-1] == 1)
			fprintf(fp, " * %s", varNames[degrees.size()-1].c_str());
		else
			fprintf(fp, " * %s^%d", varNames[degrees.size()-1].c_str(), degrees[degrees.size()-1]);
	}
}

void Monomial::dump_constant(FILE *fp, const vector<string> & varNames) const
{
	double c = coefficient.sup();
	fprintf(fp, "(%lf)", c);

	for(int i=0; i<degrees.size()-1; i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
				fprintf(fp, " * %s", varNames[i].c_str());
			else
				fprintf(fp, " * %s^%d", varNames[i].c_str(), degrees[i]);
		}
	}

	if(degrees[degrees.size()-1] != 0)
	{
		if(degrees[degrees.size()-1] == 1)
			fprintf(fp, " * %s", varNames[degrees.size()-1].c_str());
		else
			fprintf(fp, " * %s^%d", varNames[degrees.size()-1].c_str(), degrees[degrees.size()-1]);
	}
}

void Monomial::toString(string & result, const vector<string> & varNames) const
{
	string strMono;

	strMono += '(';

  /* old flowstar version
	string strInt;
	coefficient.toString(strInt);
	strMono += strInt; */
  strMono += coefficient.toString();
  //cout << degrees.size() << endl;

	for(int i=0; i<degrees.size(); i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
			{
				strMono += ' ';
				strMono += '*';
				strMono += ' ';
				strMono += varNames[i];
			}
			else
			{
				strMono += ' ';
				strMono += '*';
				strMono += ' ';
				strMono += varNames[i];
				strMono += '^';

				char strNum[NUM_LENGTH];
				sprintf(strNum, "%d", degrees[i]);
				string num(strNum);
				strMono += num;
			}
		}
	}

	strMono += ')';

	result = strMono;
}

bool Monomial::classInvariantOK() const
{
	int sum = 0;

	for(int i = 0; i<degrees.size(); ++i)
		sum += degrees[i];
	return (sum == d);
}

bool Monomial::cutoff(Monomial & monoRem, const Interval & cutoff_threshold)
{
	Interval M;
	coefficient.midpoint(M);

	monoRem = *this;

	if(M.subseteq(cutoff_threshold))
	{
		return false;
	}
	else
	{
		monoRem.coefficient -= M;
		coefficient = M;
		return true;
	}
}

bool Monomial::cutoff(const Interval & cutoff_threshold)
{
	Interval M;
	coefficient.midpoint(M);

	if(M.subseteq(cutoff_threshold))
	{
		return false;
	}
	else
	{
		coefficient = M;
		return true;
	}
}

bool Monomial::center()
{
	Interval M, intZero;

	coefficient.midpoint(M);
	if(M.subseteq(intZero))
	{
		return false;
	}
	else
	{
		coefficient = M;
		return true;
	}
}

bool operator == (const Monomial & a, const Monomial & b)
{
	if (a.d == b.d)
	{
		for(int i=0; i<a.degrees.size(); i++)
		{
			if(a.degrees[i] != b.degrees[i])
				return false;
		}
		return true;	// The two monomials are identical without considering the coefficients.
	}
	else
		return false;
}

bool operator < (const Monomial & a, const Monomial & b)
{
	if(a.d < b.d)
		return true;
	else if(a.d > b.d)
		return false;
	else	// a.d == b.d
	{
		for(int i=0; i<a.degrees.size(); ++i)
		{
			if(a.degrees[i] < b.degrees[i])
				return true;
			else if(a.degrees[i] > b.degrees[i])
				return false;
		}
	}
	return false;	// a == b
}

string Monomial::toString(const vector<string> & varNames) const
{
	string ret;
	toString(ret, varNames);
	return ret;
}

string Monomial::toString() const {
	vector<string> variables;
	for(int i = 0; i < degrees.size(); i++) {
	  char buffer[10];
    sprintf (buffer, "a%d", i);
	  variables.push_back(buffer);
	}
  //cout << "deg size: " << degrees.size();
	
	return toString(variables);
}

string Monomial::toMathematicaString() const {
  stringstream ss;
  ss.precision(10);
  
  //ss << coefficient.toString(); //TODO modify
  //ss << coefficient.midpoint();
  ss << coefficient.midToString();
  
	for(int i=0; i < degrees.size(); i++) {
		if(degrees[i] == 0) {
      continue;
    }
    ss << " * a" << i;
    if(degrees[i] != 1) {
      ss << "^" << degrees[i];
    }
	}
	return ss.str();
}

Monomial Monomial::transform(map<int, int> lookup, int size) {
  //cout << "m transform" << endl;
  int i = 0;
  for(vector<int>::iterator it=degrees.begin(); it < degrees.end(); it++) {
    //cout << i << ": " << *it << endl;
    
    //check that all nonzero degrees are retained
    if(*it != 0) {
      if(lookup.find(i) == lookup.end()) {
        cout << "error: non zero degreee in transforming" << endl;
        exit(12);
      }
    }
    i++;
  }
  vector<int> v;
  for(int j = 0; j < size; j++) {
    v.push_back(0);
  }
  for(map<int,int>::iterator it=lookup.begin(); it!=lookup.end(); ++it) {
    v.at(it->second) = degrees.at(it->first);
  }
  Monomial* ret = new Monomial(coefficient, v);
  return *ret;
}

//transforms monomials variables into a new one where, variables are
//t, indexes(0), indexes(1), ... indexes(last)
Monomial Monomial::transform(vector<int> indexes) {
  //cout << "m transform" << endl;
  //cout << "before" << endl;
  int i = 0;
  for(vector<int>::iterator it=degrees.begin(); it < degrees.end(); it++) {
    if(i == 0) { //skip time
      i++;
      continue;
    }
    //cout << i << ": " << *it << endl;
    
    //check that all nonzero degrees are retained
    if(*it != 0) {
      if(find(indexes.begin(), indexes.end(), i - 1) != indexes.end()) {
        //mlog1(sbuilder() << "if, i: " << i);
      } else {
        throw std::invalid_argument(
            "no index for a non zero degree in transforming monomial");
      }
    }
    i++;
  }
  vector<int> v;
  v.push_back(degrees.at(0)); // copy time
  //cout << "middle" << endl;
  for(vector<int>::iterator it=indexes.begin(); it < indexes.end(); it++) {
    //cout << "before2" << endl;
    //add zero degree for padding variable
    if(*it == PADDING_VARIABLE) {
      v.push_back(0);
      continue;
    }
    //cout << "middle2" << endl;
    v.push_back(degrees.at(*it + 1));
    //cout << "after2" << endl;
  }
  //cout << "after" << endl;
  Monomial* ret = new Monomial(coefficient, v);
  return *ret;
}

//transforms monomials variables into a new one where, variables are
//t, indexes(0), indexes(1), ... indexes(last)
Monomial* Monomial::ptransform(vector<int> indexes) {
  //cout << "m transform" << endl;
  
  int i = 0;
  for(vector<int>::iterator it=degrees.begin(); it < degrees.end(); it++) {
    if(i == 0) { //skip time
      i++;
      continue;
    }
    //cout << i << ": " << *it << endl;
    
    //check that all nonzero degrees are retained
    if(*it != 0) {
      if(find(indexes.begin(), indexes.end(), i - 1) != indexes.end()) {
        //mlog1(sbuilder() << "if, i: " << i);
      } else {
        throw std::invalid_argument(
            "no index for a non zero degree in transforming monomial");
      }
    }
    i++;
  }
  vector<int> v;
  v.push_back(degrees.at(0)); // copy time
  for(vector<int>::iterator it=indexes.begin(); it < indexes.end(); it++) {
    //add zero degree for padding variable
    if(*it == PADDING_VARIABLE) {
      v.push_back(0);
      continue;
    }
    v.push_back(degrees.at(*it + 1));
  }
  Monomial* ret = new Monomial(coefficient, v);
  cout << ret << endl;
  cout << "here" << endl;
  return ret;
}

Monomial & Monomial::rtransform(vector<int> indexes) {
  //cout << "m transform" << endl;
  
  int i = 0;
  for(vector<int>::iterator it=degrees.begin(); it < degrees.end(); it++) {
    if(i == 0) { //skip time
      i++;
      continue;
    }
    //cout << i << ": " << *it << endl;
    
    //check that all nonzero degrees are retained
    if(*it != 0) {
      if(find(indexes.begin(), indexes.end(), i - 1) != indexes.end()) {
        //mlog1(sbuilder() << "if, i: " << i);
      } else {
        throw std::invalid_argument(
            "no index for a non zero degree in transforming monomial");
      }
    }
    i++;
  }
  vector<int> v;
  v.push_back(degrees.at(0)); // copy time
  for(vector<int>::iterator it=indexes.begin(); it < indexes.end(); it++) {
    //add zero degree for padding variable
    if(*it == PADDING_VARIABLE) {
      v.push_back(0);
      continue;
    }
    v.push_back(degrees.at(*it + 1));
  }
  Monomial* ret = new Monomial(coefficient, v);
  cout << ret << endl;
  cout << "rhere" << endl;
  return *ret;
}

Monomial Monomial::prepareSecondary(int prefix) {
  vector<int> v;
  v.push_back(degrees.at(0)); // copy time
  for(int i = 0; i < prefix; i++) {
    v.push_back(0);
  }
  for(vector<int>::iterator it=degrees.begin(); it < degrees.end(); it++) {
    if(it == degrees.begin()) {
      continue; //skip time cause already copied
    }
    v.push_back(*it);
  }
  Monomial* ret = new Monomial(coefficient, v);
  return *ret;
}

int Monomial::equals(const Monomial & m2) const {
  if(degrees.size() != m2.degrees.size())
    return 1;
  for(int i = 0; i < degrees.size(); i++) {
    if(degrees.at(i) != m2.degrees.at(i))
      return 1;
  }
  //TODO check coefficient (find a fix of equal intervals not being equal);
  return 0;
}

int Monomial::








getVariableCount() const {
  return degrees.size();
}

vector<int> Monomial::getParams() const {
  vector<int> ret;
  int index = 0;
  for(vector<int>::const_iterator it = degrees.begin(); it < degrees.end();
      it++, index++) {
    if((*it) != 0) {
      ret.push_back(index-1); // -1 because of time
    }
  }
  return ret;
}

//return the index of the linear variable of the monomial
//if the monomial is not linear then returns -1
//index 0 is for time
int Monomial::linearVariable() const {
  int variable = -1;
  int var = 0;
  for(vector<int>::const_iterator it = degrees.begin(); it < degrees.end();
      it++, var++) {
    if((*it) == 1) {
      if(variable == -1) //needs to be previously unset
        variable = var;
      else //if multiple linear variables return -1
        return -1;
    }
    if((*it) > 1) { //if nonlinear variables return -1
      return -1;
    }
  }
  return variable;
}


Monomial Monomial::addNVariables(int n) const {
  vector<int> deg = degrees;
  for(int i = 0; i < n; i++)
    deg.push_back(0);
  
	Monomial ret(coefficient, deg);
	return ret;
}

Monomial* memCreateMonomial(Interval & coef, vector<int> & powers) {
  //cout << "creating" << endl;
  Interval *coefCopy = new Interval(coef);
  vector<int>* degCopy = new vector<int>();
  copy(powers.begin(), powers.end(), back_inserter(*degCopy));
  
  Monomial* m = new Monomial(*coefCopy, *degCopy);
  return m;
}


