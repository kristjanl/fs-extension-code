#ifndef MYLOGGER_H_
#define MYLOGGER_H_

#include <fstream>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <vector>
#include <stdexcept>
#include <set>

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
#define mforce1(o) logger.force(o)
#define mforce(n,o) {mreset(old); mlog(n,o); mrestore(old);}
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
  #undef mforce1
  #define mforce1(o)
  #undef mforce2
  #define mforce2(b,c)
  #undef mlist
  #define mlist(name, o)
#endif



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
		
		void log(std::string, TaylorModelVec);
		void log(std::string, TaylorModel);
		void log(std::string, Polynomial);
		void log(std::string, Monomial);
		void log(std::string, std::vector<Interval>);
		void log(std::string, std::vector<int>);
		void log(std::string, std::vector<std::string>);
		void log(std::string, std::vector<HornerForm>);
		void log(std::string, HornerForm);
		void log(std::string, std::vector<RangeTree *>);
		void log(std::string, RangeTree *);
    void log(std::string, Matrix m);
    void log(std::string, std::set<int> myset);
		
		void log(Polynomial *);
    void log(Monomial m);
    
    void force(std::string s);
		void log(std::string s);
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

std::vector<std::string> getVNames(int i);


extern mylogger2 logger;
extern mylogger2 logger2;


#endif /* MYLOGGER_H_ */
