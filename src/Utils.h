#ifndef UTILS_H_
#define UTILS_H_

#include "include.h"
#include "TaylorModel.h"
#include "MyLogger.h"
#include "MyComponent.h"
#include "Interval.h"
#include "OutputWriter.h"
#include "Continuous.h"

using namespace std;

class ShrinkWrappingCondition {
	public:
		ShrinkWrappingCondition(int steps);
		ShrinkWrappingCondition(); //constructor for using remainder
		bool checkApplicability(const vector<MyComponent *> comps, 
		    const vector<Interval> & estimation);
		void log() const;
		int getCount() const;
  private:
    bool useRemainder, useSteps;
    int steps, cycleSteps;
    int count;
};

class MySettings {
  public:
    OutputWriter & writer;
    const int order;
    const double step; 
    const double time;
    const vector<Interval> & estimation;
    const vector<Interval> step_end_exp_table;
    const vector<Interval> domain;
    MySettings(OutputWriter & writer, const int order, const double step, 
        const double time, const vector<Interval> & estimation, 
        const vector<Interval> & step_end_exp_table, 
        const vector<Interval> & domain);
  
};

#endif /* UTILS_H_ */
