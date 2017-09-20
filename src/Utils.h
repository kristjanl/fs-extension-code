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

/*
class MySettings {
  public:
    OutputWriter *writer;
    const int order;
    const double step; 
    const double time;
    const vector<Interval> & estimation;
    const vector<Interval> step_exp_table;
    const vector<Interval> step_end_exp_table;
    const vector<Interval> domain;
    const Interval & cutoff;
    MySettings(OutputWriter *writer, const int order, const double step, 
        const double time, const vector<Interval> & estimation, 
        const vector<Interval> step_exp_table, 
        const vector<Interval> & step_end_exp_table, 
        const vector<Interval> & domain, const Interval & cutoff);
  
};*/

class MySettings;

class PrecondModel {
  public:
    PrecondModel(TaylorModelVec left, TaylorModelVec right);
    
    TaylorModelVec left;
    TaylorModelVec right;
    
    TaylorModelVec composed(MySettings *settings);
};

void serializeTMV(TaylorModelVec & tmv, string filename);
void serializeFlows(MyComponent *comp, string filename);
vector<TaylorModelVec> & deserializeFlows(string filename);
vector<TaylorModelVec *> pDeserializeFlows(string filename);

void compareFlows(vector<TaylorModelVec *> & first, 
    vector<TaylorModelVec *> & second);
void compareFlows(vector<TaylorModelVec> & first, 
    vector<TaylorModelVec> & second);
double compareIntervalVecs(vector<Interval> & f, vector<Interval> & s);

void printTMVFiles(string file1, string file2, int index1, int index2);

vector<Interval> getUnitBox(int n);

class TMVSerializer {
  public:
    TMVSerializer(string filename);
    TMVSerializer(string filename, int maxSize);
    void add(const TaylorModelVec & tmv);
    void serialize();
  private:
    string filename;
    vector<TaylorModelVec> tmvs;
    int maxSize;
};


#endif /* UTILS_H_ */
