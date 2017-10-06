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


#define mreset(name) int name = logger.reset()
#define mrestore(name) logger.restore(name)

#define minc() logger.inc()
#define mdec() logger.dec()
#define mdisable() logger.disable()

#define mlog(name, o) logger.log(name, o)
#define mlog1(o) logger.log(o)
#define mforce(o) logger.force(o)
#define mlist(name, o) logger.log(name, o) 

#ifdef enablelog
  #define LOGGING_STATUS "logging enabled"
#else
  #define LOGGING_STATUS "logging disabled"
  #undef mreset
  #define mreset(name)
  #undef mrestore
  #define mrestore(name)
  
  #undef minc
  #define minc()
  #undef mdec
  #define mdec()
  #undef mdisable
  #define mdisable()
  
  #undef mlog
  #define mlog(name, o)
  #undef mlog1
  #define mlog1(o)
  #undef mlist
  #define mlist(name, o)
#endif



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

/*
sbuilder & operator << (sbuilder & sb, const T &data) {
	ss << data;
	return *this;
}
*/

class mylogger2 {
	public:
		mylogger2();
		mylogger2(bool disabled);
		
		void log(string, TaylorModelVec);
		void log(string, TaylorModel);
		void log(string, vector<Interval>);
		void log(string, vector<int>);
		void log(string, vector<string>);
		void log(string, vector<HornerForm>);
		void log(string, vector<RangeTree *>);
		void log(string, RangeTree *);
    void log(string, Matrix m);
		
		void log(Polynomial *);
    void log(Monomial m);
    
    void force(string s);
		void log(string s);
    void log(int i);
		
		void inc();
		void dec();
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
