/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <math.h>
#include <boost/archive/text_oarchive.hpp>

#include "include.h"
#include "Exceptions.h"

extern mpfr_prec_t intervalNumPrecision;


class Interval
{
private:
	mpfr_t lo;		// the lower bound
	mpfr_t up;		// the upper bound
public:
	Interval();
	Interval(const double c);
	Interval(const double l, const double u);
	Interval(const char *strLo, const char *strUp);
	Interval(const Interval & I);
	Interval(const mpfr_t lower, const mpfr_t upper);
	~Interval();

	void set(const double l, const double u);
	void set(const double c);

	void setInf(const double l);
	void setInf(const Interval & I);

	void setSup(const double u);
	void setSup(const Interval & S);

	void split(Interval & left, Interval & right) const;		// split the interval at the midpoint
	void split(list<Interval> & result, const int n) const;		// split the interval uniformly by n parts

	void set_inf();

	double sup() const;
	double inf() const;

	void sup(Interval & S) const;
	void inf(Interval & I) const;

	double midpoint() const;
	void midpoint(Interval & M) const;
	void remove_midpoint(Interval & M);
	void remove_midpoint();

	Interval intersect(const Interval & I) const;
	void intersect_assign(const Interval & I);

	void bloat(const double e);
	bool within(const Interval & I, const double e) const;

	double width() const;
	void width(Interval & W) const;

	double mag() const;		// max{|lo|,|up|}
	void mag(Interval & M) const;

	void abs(Interval & result) const;
	void abs_assign();		// absolute value

	bool subseteq(const Interval & I) const;	// returns true if the interval is a subset of I
	bool supseteq(const Interval & I) const;	// returns true if the interval is a superset of I
	bool valid() const;
	bool operator == (const Interval & I) const;
	bool operator != (const Interval & I) const;
	bool operator > (const Interval & I) const;		// lo > up
	bool operator < (const Interval & I) const;		// up < lo
	bool operator <= (const Interval & I) const;	// lo < lo
	bool operator >= (const Interval & I) const;	// up > up

	bool smallereq(const Interval & I) const; 		// up <= lo

	Interval & operator = (const Interval & I);
	Interval & operator += (const Interval & I);
	Interval & operator -= (const Interval & I);
	Interval & operator *= (const Interval & I);
	Interval & operator /= (const Interval & I);
	Interval & operator ++ ();
	Interval & operator -- ();

	const Interval operator + (const Interval & I) const;
	const Interval operator - (const Interval & I) const;
	const Interval operator * (const Interval & I) const;
	const Interval operator / (const Interval & I) const;

	void sqrt(Interval & result) const;		// square root
	void inv(Interval & result) const;		// additive inverse
	void rec(Interval & result) const;		// reciprocal
	void sqrt_assign();
	void inv_assign();
	void rec_assign();

	void add_assign(const double c);
	void sub_assign(const double c);
	void mul_assign(const double c);
	void div_assign(const double c);

	Interval pow(const int n) const;
	Interval exp() const;
	Interval sin() const;
	Interval cos() const;
	Interval log() const;

	void pow_assign(const int n);
	void exp_assign();
	void sin_assign();
	void cos_assign();
	void log_assign();

	double widthRatio(const Interval & I) const;

	void toString(string & result) const;
	void printFull() const;
	string toString() const;
	string toString(int prec) const;
	string getLower(int prec) const;
	string getHigher(int prec) const;
	string getLower() const;
	string getHigher() const;
	bool isClose(const Interval & I, double d) const;
	Interval distance(const Interval & I) const;
	
	void compare(const Interval & I) const;
	
	void dump(FILE *fp) const;
	void serialize(FILE *fp) const;
	void output(FILE *fp, const char * msg, const char * msg2) const;
};

// returns true if all the intervals in I1 are subseteq to intervals in I2
bool subseteq(const vector<Interval> & I1, const vector<Interval> & I2);

const static Interval ZERO_INTERVAL;
const static Interval ONE_INTERVAL(1,1);
const static Interval UNIT_INTERVAL(-1,1);

#endif /* INTERVAL_H_ */
