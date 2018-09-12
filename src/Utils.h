#ifndef UTILS_H_
#define UTILS_H_

#include <map>

#include "include.h"
#include "TaylorModel.h"
#include "MyLogger.h"
#include "MyComponent.h"
#include "Interval.h"
#include "OutputWriter.h"
#include "Continuous.h"

//start the clock with variable name <name>Start
#define tstart(name) clock_t name##Start = clock();
//start the clock with variable <name>End, subtract <name>
// Start from it and store it in timeLookup with key <name>
#define tend(name) clock_t name##End = clock(); \
  double name##Dur = double(name##End - name##Start) / CLOCKS_PER_SEC; \
  timeLookup[#name] += name##Dur;
  

#define tstart1(name) clock_t name##Start = clock();
#define tend1(name) clock_t name##End = clock(); \
  double name##Dur = double(name##End - name##Start) / CLOCKS_PER_SEC; \
  timeLookup[#name] += name##Dur;  
  
#define treset(name) timeLookup[#name] = 0;

#define tprint(prefix) printTimes(prefix);

#define taddToInfo(infoName, clockName, infos) addTimeToInfo(infoName, #clockName , infos)

extern map<string, double> timeLookup;

void printTimes(string prefix);

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

class MySettings;

class PrecondModel {
  public:
    PrecondModel(TaylorModelVec left, TaylorModelVec right);
    
    TaylorModelVec left;
    TaylorModelVec right;
    
    TaylorModelVec composed(MySettings *settings);
};


class NamedTMV {
  public:
    string name;
    TaylorModelVec tmv;
    NamedTMV(string name, TaylorModelVec tmv);
};

void serializeTMV(TaylorModelVec & tmv, string filename);
void serializeFlows(MyComponent *comp, string filename);
vector<TaylorModelVec> & deserializeFlows(string filename);
vector<TaylorModelVec *> pDeserializeFlows(string filename);
vector<NamedTMV> pDeserializeNamedFlows(string filename);

void compareFlows(vector<TaylorModelVec *> & first, 
    vector<TaylorModelVec *> & second);
void compareFlows(vector<TaylorModelVec> & first, 
    vector<TaylorModelVec> & second);
double compareIntervalVecs(vector<Interval> & f, vector<Interval> & s);

void printTMVFiles(string file1, string file2, string name, 
    int index1, int index2);
    
void toMathematica(string file);

vector<Interval> getUnitBox(int n);

class TMVSerializer {
  public:
    TMVSerializer(string filename);
    TMVSerializer(string filename, int maxSize);
    TMVSerializer(string filename, int maxSize, bool active);
    void add(const TaylorModelVec & tmv);
    void add(const TaylorModelVec & tmv, string name);
    void serialize();
    void activate();
  private:
    string filename;
    vector<TaylorModelVec> tmvs;
    vector<string> names;
    int maxSize;
    bool active;
};

extern TMVSerializer *pSerializer;

void addTimeToInfo(string name, string clockName, vector<string> & infos);

int findPos(int value, vector<int> *v);

int isIn(int value, vector<int> *v);


void addMyInfo(vector<string> & info);
void addFlowInfo(vector<string> & info);

void printComponents(MySettings *settings);

#endif /* UTILS_H_ */

