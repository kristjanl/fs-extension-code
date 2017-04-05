#ifndef UTILS_H_
#define UTILS_H_

#include "include.h"
#include "TaylorModel.h"
#include "MyLogger.h"
#include "MyComponent.h"
#include "Interval.h"

using namespace std;

void foo();

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



#endif /* UTILS_H_ */
