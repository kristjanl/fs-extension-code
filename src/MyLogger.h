#ifndef MYLOGGER_H_
#define MYLOGGER_H_

#include <iostream>
#include <fstream>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <vector>
#include <stdexcept>


#include "Interval.h"
#include "Polynomial.h"
#include "TaylorModel.h"

using namespace std;

//class TaylorModelVec;
//class TaylorModel;
//class Interval;
//class Polynomial;

struct sbuilder
{
	std::stringstream ss;
	template<typename T>
	sbuilder & operator << (const T &data) {
		ss << data;
		return *this;
	}
	operator std::string() { return ss.str(); }
};


class mylogger2 {
	public:
		void logVI(vector<Interval>);
		void logVI(string, vector<Interval>);
		void logVHF(string, vector<HornerForm>);
		void logVi(string, vector<int>);
		void listVi(string, vector<int>);
		void logTMV(string, TaylorModelVec);
		void logTMVRem(string, TaylorModelVec);
		void logTM(string, TaylorModel);
    void logM(Monomial m);
    void logMatrix(Matrix m);
    void force(string s);
		void log(string s);
    void log(int i);
    void logd(double d);
    void printLevel();
		void logPoly(Polynomial *);
		void inc();
		void dec();
		mylogger2();
		mylogger2(bool disabled);
		void disable();
		void enable();
		int reset();
		void restore(int old);
	private:
		int ltab;
		int disabled;
	
};

vector<string> getVNames(int i);


extern mylogger2 logger;
extern mylogger2 logger2;


#endif /* MYLOGGER_H_ */
