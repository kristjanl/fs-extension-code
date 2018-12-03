/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

/*
 * MPFR_RNDU: round toward plus infinity (roundTowardPositive in IEEE 754-2008),
 * MPFR_RNDD: round toward minus infinity (roundTowardNegative in IEEE 754-2008),
 */

#include "DoubleInterval.h"

#include <fenv.h>
#pragma STDC FENV_ACCESS ON
#include <math.h>
#include <iomanip> //setprecision
#include <sstream> //stringstream

using namespace std;

//TODO remove this later likely
//because of modelparser.y
mpfr_prec_t intervalNumPrecision = normal_precision;

Interval::Interval()
{
  lo = 0.0;
  up = 0.0;
}

Interval::Interval(const double c)
{
	//std::fesetround(FE_DOWNWARD);
  lo = c;
	//std::fesetround(FE_UPWARD);
  up = c;
}

Interval::Interval(const double l, const double u)
{
  lo = l;
  up = u;
}

Interval::Interval(const char *strLo, const char *strUp)
{
	std::fesetround(FE_DOWNWARD);
  lo = stod(strLo);
	std::fesetround(FE_UPWARD);
  up = stod(strUp);
}

Interval::Interval(const mpfr_t lower, const mpfr_t upper) {
  throw runtime_error("Shouldn't construct double Intervals from mpfr numbers");
}

Interval::Interval(const Interval & I)
{
	lo = I.lo;
	up = I.up;
}

//TODO is this necessary?
Interval::~Interval() {
	
}

void Interval::set(const double l, const double u)
{
	lo = l;
  up = u;
}

void Interval::set(const double c)
{
	lo = c;
  up = c;
}

void Interval::setInf(const double l)
{
	lo = l;
}

void Interval::setInf(const Interval & I)
{
	lo = I.lo;
}

void Interval::setSup(const double u)
{
	up = u;
}

void Interval::setSup(const Interval & S)
{
	up = S.up;
}

void Interval::split(Interval & left, Interval & right) const
{
	left.lo = lo;
	std::fesetround(FE_UPWARD);
  left.up = (lo + up)/2.0;

	std::fesetround(FE_DOWNWARD);
  right.lo = (lo + up)/2.0;
  right.up = up;
}

void Interval::split(list<Interval> & result, const int n) const
{
  throw runtime_error("'Interval::split(list<Interval> & result, const int n) const' not ported!");
}

void Interval::set_inf()
{
	lo = -1;
  up = 1;
}

double Interval::sup() const
{
	return up;
}

double Interval::inf() const
{
	return lo;
}

void Interval::sup(Interval & S) const
{
	S.lo = up;
  S.up = up;
}

void Interval::inf(Interval & I) const
{
	I.lo = up;
  I.up = up;
}

double Interval::midpoint() const
{
	std::fesetround(FE_TONEAREST);
	return (lo + up) / 2.0;
}

void Interval::midpoint(Interval & M) const
{
	std::fesetround(FE_DOWNWARD);
  M.lo = (lo + up)/2.0;
	std::fesetround(FE_UPWARD);
  M.up = (lo + up)/2.0;

}

void Interval::remove_midpoint(Interval & M)
{
	std::fesetround(FE_DOWNWARD);
  M.lo = (lo + up)/2.0;
	std::fesetround(FE_UPWARD);
  M.up = (lo + up)/2.0;
	//std::fesetround(FE_UPWARD);
	up -= M.lo;
	std::fesetround(FE_DOWNWARD);
	lo -= M.up;
}

void Interval::remove_midpoint()
{
	double tmp1, tmp2;

	std::fesetround(FE_DOWNWARD);
  tmp1 = (lo + up)/2.0;
	std::fesetround(FE_UPWARD);
  tmp2 = (lo + up)/2.0;
	//std::fesetround(FE_UPWARD);
  up -= tmp2;
	std::fesetround(FE_DOWNWARD);
  lo -= tmp1;
}

Interval Interval::intersect(const Interval & I) const
{
	Interval result;

	if(lo > I.lo) {
    result.lo = lo;
	} else {
    result.lo = I.lo;
	}

	if(up > I.up) {
    result.up = I.up;
	} else {
    result.up = up;
	}

	return result;
}

void Interval::intersect_assign(const Interval & I)
{
	if(lo < I.lo) {
		lo = I.lo;
	}

	if(up > I.up) {
		up = I.up;
	}
}

void Interval::bloat(const double e)
{
	std::fesetround(FE_DOWNWARD);
  lo -= e;
	std::fesetround(FE_UPWARD);
  up += e;
}

bool Interval::within(const Interval & I, const double e) const
{
	double d;

  //TODO don't this is necessary, but leaving it here now

	std::fesetround(FE_UPWARD);
	if(up >= I.up) {
    //std::fesetround(FE_UPWARD);
    d = up - I.up;
	} else {
	  //std::fesetround(FE_DOWNWARD);
    //d = up - I.up;
    d = I.up - up;
	}

  //not needed now
	//std::fesetround(FE_UPWARD);
  //d = fabs(d);

	if(d > e) {
		return false;
	}

	if(lo >= I.lo) {
	  //std::fesetround(FE_UPWARD);
		d = lo - I.lo;
	} else {
	  //std::fesetround(FE_DOWNWARD);
		d = I.lo - lo;
	}

  //not needed now
	std::fesetround(FE_DOWNWARD);
  d = fabs(d);

	if(d > e) {
		return false;
	}

  return true;
}

double Interval::width() const
{
	std::fesetround(FE_UPWARD);
	return up - lo;
}

void Interval::width(Interval & W) const
{
	std::fesetround(FE_UPWARD);
	W.lo = up - lo;
	W.up = W.up;
}

double Interval::mag() const
{
  //I assume that fabs will only change sign bit (no rounding needed)
	double inf = fabs(lo);
	double sup = fabs(up);
	return inf < sup ? sup : inf;
}

void Interval::mag(Interval & M) const
{
  //I assume that fabs will only change sign bit (no rounding needed)
	double tmp1 = fabs(lo);
  double tmp2 = fabs(up);;

	if(tmp1 > tmp2) {
		M.lo = tmp1;
    M.up = tmp1;
	} else {
		M.lo = tmp2;
    M.up = tmp2;
	}
}

void Interval::abs(Interval & result) const
{
  //I assume that fabs will only change sign bit (no rounding needed)
	double tmp1 = fabs(lo);
  double tmp2 = fabs(up);

	if(tmp1 > tmp2) {
		result.lo = tmp2;
		result.up = tmp1;
	} else {
		result.lo = tmp1;
		result.up = tmp2;
	}
}

void Interval::abs_assign()
{
	double tmp1 = fabs(lo);
  double tmp2 = fabs(up);

	if(tmp1 > tmp2) {
		lo = tmp2;
		up = tmp1;
	} else {
		lo = tmp1;
		up = tmp2;
	}
}

bool Interval::subseteq(const Interval & I) const
{
	if( (I.lo <= lo) && (I.up >= up) )
		return true;
	return false;
}

bool Interval::supseteq(const Interval & I) const
{
	if( (lo <= I.lo) && (up >= I.up) )
		return true;
	return false;
}

bool Interval::valid() const
{
	if(up >= lo) {
		return true;
	}
  return false;
}

bool Interval::operator == (const Interval & I) const
{
	return (lo == I.lo) && (up == I.up);
}

bool Interval::operator != (const Interval & I) const
{
	return (lo != I.lo) || (up != I.up);
}

bool Interval::operator > (const Interval & I) const
{
	return lo > I.up;
}

bool Interval::operator < (const Interval & I) const
{
	return up < I.lo;
}

bool Interval::operator <= (const Interval & I) const
{
  //why not <=?
	//return ( (mpfr_cmp(lo, I.lo) < 0) );
  return lo < I.lo;
}

bool Interval::operator >= (const Interval & I) const
{
	return up > I.up;
}

bool Interval::smallereq(const Interval & I) const
{
	return up <= I.lo;
}

Interval & Interval::operator = (const Interval & I)
{
	if(this == &I)
		return *this;	// check for self assignment

  lo = I.lo;
	up = I.up;

	return *this;
}

Interval & Interval::operator += (const Interval & I)
{
	std::fesetround(FE_DOWNWARD);
	lo += I.lo;
	std::fesetround(FE_UPWARD);
	up += I.up;
	return *this;
}

Interval & Interval::operator -= (const Interval & I)
{
	std::fesetround(FE_DOWNWARD);
	lo -= I.up;
	std::fesetround(FE_UPWARD);
	up -= I.lo;
	return *this;
}

Interval & Interval::operator *= (const Interval & I)
{
  double lolo, loup, uplo, upup, min, max;

	// compute the lower bound
	std::fesetround(FE_DOWNWARD);
	lolo = lo * I.lo;
	loup = lo * I.up;
	uplo = up * I.lo;
	upup = up * I.up;

	min = lolo;
	if(min > loup) {
		min = loup;
	}
	if(min > uplo) {
		min = uplo;
	}
	if(min > upup) {
		min = upup;
	}

	std::fesetround(FE_UPWARD);
	lolo = lo * I.lo;
	loup = lo * I.up;
	uplo = up * I.lo;
	upup = up * I.up;

	max = lolo;
	if(max < loup) {
		max = loup;
	}
	if(max < uplo) {
		max = uplo;
	}
	if(max < upup) {
		max = upup;
	}

	lo = min;
	up = max;

	return *this;
}

Interval & Interval::operator /= (const Interval & I)
{
	Interval tmp;

	I.rec(tmp);
	*this *= tmp;

	return *this;
}

Interval & Interval::operator ++ ()
{
	std::fesetround(FE_DOWNWARD);
  lo++;
	std::fesetround(FE_UPWARD);
  up++;

	return *this;
}

Interval & Interval::operator -- ()
{
	std::fesetround(FE_DOWNWARD);
  lo--;
	std::fesetround(FE_UPWARD);
  up--;

	return *this;
}

const Interval Interval::operator + (const Interval & I) const
{
	Interval result = *this;
	result += I;
	return result;
}

const Interval Interval::operator - (const Interval & I) const
{
	Interval result = *this;
	result -= I;
	return result;
}

const Interval Interval::operator * (const Interval & I) const
{
	Interval result = *this;
	result *= I;
	return result;
}

const Interval Interval::operator / (const Interval & I) const
{
	Interval result = *this;
	result /= I;
	return result;
}

void Interval::sqrt(Interval & result) const
{
	if(lo < 0) {
		printf("Exception: Square root of a negative number.\n");
		exit(1);
	}
	std::fesetround(FE_DOWNWARD);
	result.lo = std::sqrt(lo);
	std::fesetround(FE_UPWARD);
	result.up = std::sqrt(up);
}

void Interval::inv(Interval & result) const
{
  //I assume that '-' will only change sign bit (no rounding needed)
	//std::fesetround(FE_DOWNWARD);
	result.lo = -up;
	//std::fesetround(FE_UPWARD);
	result.up = -lo;
}

void Interval::rec(Interval & result) const
{
	if (lo <= 0 && up >= 0) {
		printf("Exception: Divided by 0.\n");
		exit(1);
	} else {
	  std::fesetround(FE_DOWNWARD);
    result.lo = 1.0/up;
	  std::fesetround(FE_UPWARD);
    result.up = 1.0/lo;
	}
}

void Interval::sqrt_assign()
{
	if(lo < 0) {
		printf("Exception: Square root of a negative number.\n");
		exit(1);
	}
	std::fesetround(FE_DOWNWARD);
	lo = std::sqrt(lo);
	std::fesetround(FE_UPWARD);
	up = std::sqrt(up);
}

void Interval::inv_assign()
{
	Interval result;
	this->inv(result);
	*this = result;
}

void Interval::rec_assign()
{
	Interval result;
	this->rec(result);
	*this = result;
}

void Interval::add_assign(const double c)
{
	std::fesetround(FE_DOWNWARD);
	lo += c;
	std::fesetround(FE_UPWARD);
	up += c;
}

void Interval::sub_assign(const double c)
{
	std::fesetround(FE_DOWNWARD);
	lo -= c;
	std::fesetround(FE_UPWARD);
	up -= c;
}

void Interval::mul_assign(const double c)
{
	Interval result;

	if(c > 0) {
    std::fesetround(FE_DOWNWARD);
    result.lo = lo * c;
    std::fesetround(FE_UPWARD);
    result.up = up * c;
	} else {
    std::fesetround(FE_DOWNWARD);
    result.lo = up * c;
    std::fesetround(FE_UPWARD);
    result.up = lo * c;
	}

	*this = result;
}

void Interval::div_assign(const double c)
{
	Interval result;

	if(c > 0) {
    std::fesetround(FE_DOWNWARD);
    result.lo = lo / c;
    std::fesetround(FE_UPWARD);
    result.up = up / c;
	} else {
    std::fesetround(FE_DOWNWARD);
    result.lo = up / c;
    std::fesetround(FE_UPWARD);
    result.up = lo / c;
	}

	*this = result;
}

Interval Interval::pow(const int n) const
{
	Interval result;

	if(n % 2 == 1) {	// n is odd
    std::fesetround(FE_DOWNWARD);
    result.lo = std::pow(lo, n);
    std::fesetround(FE_UPWARD);
    result.up = std::pow(up, n);
		return result;
	}
	// n is even
	if(lo >= 0) {			// 0 <= lo <= up
		std::fesetround(FE_DOWNWARD);
		result.lo = std::pow(lo, n);
		std::fesetround(FE_UPWARD);
		result.up = std::pow(up, n);
	} else if(up <= 0)	{	// lo <= up <= 0
		std::fesetround(FE_DOWNWARD);
		result.lo = std::pow(up, n);
		std::fesetround(FE_UPWARD);
		result.up = std::pow(lo, n);
	} else {									// lo < 0 < up
		double tmp1, tmp2;
		std::fesetround(FE_UPWARD);
		tmp1 = std::pow(lo, n);
		tmp2 = std::pow(up, n);
		if(tmp1 >= tmp2) {
			result.up = tmp1;
		} else {
			result.up = tmp2;
		}
		//wasn't present in Flow* (but putting it here for being explicit)
		result.lo = 0;
	}

	// return [a,b]
	return result;
}

Interval Interval::exp() const
{
	Interval result;
	
  throw runtime_error("std::exp function is not precise");
  std::fesetround(FE_DOWNWARD);
  result.lo = std::exp(lo);
  std::fesetround(FE_UPWARD);
  result.up = std::exp(up);

	return result;
}

Interval Interval::sin() const
{
  throw runtime_error("'Interval::sin() not ported!");
}

Interval Interval::cos() const
{
  throw runtime_error("'Interval::cos() not ported!");
}

Interval Interval::log() const
{
	if(lo <= 0) {
		printf("Exception: Logarithm of a non-positive number.\n");
		exit(1);
	} else {
		Interval result;

  	throw runtime_error("std::exp function is not precise (after 17 digits");
  	std::fesetround(FE_DOWNWARD);
  	result.lo = std::log(lo);
  	std::fesetround(FE_UPWARD);
  	result.up = std::log(up);
		return result;
	}
}

void Interval::pow_assign(const int n)
{
	if(n % 2 == 1)		// n is odd
	{		
    std::fesetround(FE_DOWNWARD);
    lo = std::pow(lo, n);
    std::fesetround(FE_UPWARD);
    up = std::pow(up, n);
		return;
	}
	// n is even
	if(lo >= 0) {		// 0 <= lo <= up
		std::fesetround(FE_DOWNWARD);
		lo = std::pow(lo, n);
		std::fesetround(FE_UPWARD);
		up = std::pow(up, n);
	} else if(up <= 0) {		// lo <= up <= 0
		std::fesetround(FE_DOWNWARD);
		double tmp = std::pow(up, n);
		std::fesetround(FE_UPWARD);
		up = std::pow(lo, n);
		lo = tmp;
	} else {							// lo < 0 < up
		double tmp1, tmp2;
		std::fesetround(FE_UPWARD);
		tmp1 = std::pow(lo, n);
		tmp2 = std::pow(up, n);
		if(tmp1 >= tmp2) {
			up = tmp1;
		} else {
			up = tmp2;
		}
		lo = 0;
	}
}

void Interval::exp_assign()
{
	throw runtime_error("std::exp function is not precise");
  std::fesetround(FE_DOWNWARD);
  lo = std::exp(lo);
  std::fesetround(FE_UPWARD);
  up = std::exp(up);
}

void Interval::sin_assign()
{
  throw runtime_error("'Interval::sin_assign() not ported!");
}

void Interval::cos_assign()
{
  throw runtime_error("'Interval::cos_assign() not ported!");
}

void Interval::log_assign()
{
	if(lo <= 0) {
		printf("Exception: Logarithm of a non-positive number.\n");
		exit(1);
	} else {
		throw runtime_error("std::exp function is not precise (after 17 digits");
  	std::fesetround(FE_DOWNWARD);
  	lo = std::log(lo);
  	std::fesetround(FE_UPWARD);
  	up = std::log(up);
	}
}

double Interval::widthRatio(const Interval & I) const
{
	double width1, width2, ratio;

  std::fesetround(FE_UPWARD);
	width1 = up - lo;
	width2 = I.up - I.lo;

	ratio = width2 / width1;		// we assume that width1 >= width2

	return ratio;
}

void Interval::toString(string & result) const
{
  //mlog1("(string) version");
	char strTemp[30];

	string strInt;
	strInt += '[';

  std::fesetround(FE_DOWNWARD);
  sprintf(strTemp, "%.10f", lo);
	string strLo(strTemp);

	strInt += strLo;
	strInt += ' ';
	strInt += ',';
	strInt += ' ';

  std::fesetround(FE_UPWARD);
  sprintf(strTemp, "%.10f", up);
	string strUp(strTemp);
	strInt += strUp;
	strInt += ']';

	result = strInt;
}

bool Interval::isClose(const Interval & I, double d) const {
  double infDiff = inf() - I.inf();
  double supDiff = sup() - I.sup();
  if(infDiff < -d)
    return false;
  if(infDiff > d)
    return false;
  if(supDiff < -d)
    return false;
  if(supDiff > d)
    return false;
  return true;
}

void Interval::dump(FILE *fp) const
{
	fprintf (fp, "[");

	//15 is PN in Flow* (but it's arbitrary anyway)
  std::fesetround(FE_DOWNWARD);
  fprintf (fp, "%.15f", lo);
	fprintf(fp, " , ");
  std::fesetround(FE_UPWARD);
  fprintf (fp, "%.15f", up);

	fprintf(fp, "]");
}


void serializeMpfr(FILE *fp, const double number) {
	throw runtime_error("serialize not ported yet");
	/*
	fprintf (fp, "0m");
  
  //get the mpfr radix into string
  char* str = NULL;
  mpfr_exp_t e;
  str = mpfr_get_str (NULL, &e, MPFR_SERIALIZATION_BASE, 0, number, MPFR_RNDN);
    
  //add radix point and exponent
  char buffer[64];
  sprintf (buffer, ".%s@%ld", str, (long) e);
  
  //swap sign and radix point for negative numbers
  if(str[0] == '-') {
    buffer[0] = buffer[1];
    buffer[1] = '.';
    //cout << "buffer: " << buffer << endl;
  }
	fprintf (fp, buffer); 
	*/
}

void Interval::serialize(FILE *fp) const
{
	fprintf (fp, "[");
  serializeMpfr(fp, lo);
	fprintf(fp, " , ");
  serializeMpfr(fp, up);
	fprintf(fp, "]");
}

void Interval::output(FILE * fp, const char * msg, const char * msg2) const
{
	fprintf (fp, "%s [ ", msg);

	//15 is PN in Flow* (but it's arbitrary anyway)
  std::fesetround(FE_DOWNWARD);
  fprintf (fp, "%.15f", lo);
	fprintf(fp, " , ");
  std::fesetround(FE_UPWARD);
  fprintf (fp, "%.15f", up);

	fprintf(fp, " ] %s", msg2);
}

string toStringHelper(const double data, int prec, char roundingMode) {
	if(roundingMode == 'N') { //round to nearest
    std::fesetround(FE_TONEAREST);
  } else if (roundingMode == 'Y') { // round away from zero
    if(data < 0) {
      std::fesetround(FE_DOWNWARD);
    } else {
      std::fesetround(FE_UPWARD);
    }
  } else {
    runtime_error("illegal rounding mode in toStringHelper");
  }

  
  stringstream ss;
  ss << setprecision(prec);
  ss << data;
  return ss.str();
}


string toStringHelper(const double data, int prec) {
  //round away from zero
  return toStringHelper(data, prec, 'Y');
}

string Interval::toString() const {
  //mlog1("() version");
  return toString(15);
}

string Interval::toString(int prec) const {
  //don't bother with intervals if width is small
	if(width() < ::pow(10, -prec)) {
		return "[" + toStringHelper(lo, prec) + "]";
  }
  string los = toStringHelper(lo, prec);
  string ups = toStringHelper(up, prec);
  return "[" + los + "," + ups + "]";
}

string Interval::toMathematicaString() const {
  string los = toStringHelper(lo, 20);
  string ups = toStringHelper(up, 20);
  return los + "," + ups;
}

string Interval::midToString() const {
	return toStringHelper(midpoint(), 10);
}

void Interval::printFull() const {
	//TODO probably remove this function
	throw runtime_error("printfull not ported");
  /*
  //mpfr_printf ("%.55Rg\n", lo);
  mpfr_out_str (stdout, 2, 0, lo, MPFR_RNDU);
  cout << endl;
  mpfr_out_str (stdout, 2, 0, lo, MPFR_RNDD);
  cout << endl;
	*/
}


string Interval::getLower(int prec) const {
  //round to nearest
  return toStringHelper(lo, prec, 'N');
}


string Interval::getHigher(int prec) const {
  //round to nearest
  return toStringHelper(up, prec, 'N');
}

string Interval::getLower() const {
  return getLower(15);
}


string Interval::getHigher() const {
  return getHigher(15);
}

bool subseteq(const vector<Interval> & v1, const vector<Interval> & v2) {
  if(v1.size() != v2.size())
    throw invalid_argument("Different domains of vectors");
  for(int i = 0; i < v1.size(); i++) {
    if(v1[i].subseteq(v2[i]) == false) {
      return false;
    }
  }
  return true;
}


Interval Interval::distance(const Interval & I) const {
  double lower, upper;
	//this is not rounding away, but shouldn't matter that much
  lower = fabs(lo - I.lo);
	upper = fabs(up - I.up);
  
  Interval ret;
  ret.up = max(lower, upper);
  return ret;
}





