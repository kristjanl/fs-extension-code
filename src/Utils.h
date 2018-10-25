#ifndef UTILS_H_
#define UTILS_H_

#include <map>
#include <string>
#include <vector>

class MySettings;
class MyComponent;
class TaylorModelVec;
class Interval;
class Transformer;
class OutputWriter;

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

extern std::map<std::string, double> timeLookup;

void printTimes(std::string prefix);


class MySettings {
  public:
    OutputWriter *writer;
    std::vector< std::vector<int> > intComponents;
    bool autoComponents;
    int order;
    double step; 
    double time;
    std::vector<Interval> estimation;
    std::vector<Interval> step_exp_table;
    std::vector<Interval> step_end_exp_table;
    std::vector<Interval> domain;
    const Interval *cutoff;
    bool useFlow;//TODO remove
    bool discardEmptyParams;
    std::vector<std::string> varNames;
    Transformer *transformer; // determines how are initials sets transformed for each timestep
    MySettings();
    MySettings(OutputWriter *writer, int order, double step, 
        double time, std::vector<Interval> estimation, 
        std::vector<Interval> step_exp_table, 
        std::vector<Interval> step_end_exp_table, 
        std::vector<Interval> domain, const Interval *cutoff);
    void log();
    //MySettings* toOld();
};



class ShrinkWrappingCondition {
	public:
		ShrinkWrappingCondition(int steps);
		ShrinkWrappingCondition(); //constructor for using remainder
		bool checkApplicability(const std::vector<MyComponent *> comps, 
		    const std::vector<Interval> & estimation);
		void log() const;
		int getCount() const;
  private:
    bool useRemainder, useSteps;
    int steps, cycleSteps;
    int count;
};




void serializeTMV(TaylorModelVec & tmv, std::string filename);
void serializeFlows(MyComponent *comp, std::string filename);
std::vector<TaylorModelVec> & deserializeFlows(std::string filename);
std::vector<TaylorModelVec *> pDeserializeFlows(std::string filename);

void compareFlows(std::vector<TaylorModelVec *> & first, 
    std::vector<TaylorModelVec *> & second);
void compareFlows(std::vector<TaylorModelVec> & first, 
    std::vector<TaylorModelVec> & second);
double compareIntervalVecs(std::vector<Interval> & f, std::vector<Interval> & s);

void printTMVFiles(std::string file1, std::string file2, std::string name, 
    int index1, int index2);
    
void toMathematica(std::string file);

std::vector<Interval> getUnitBox(int n);

class TMVSerializer {
  public:
    TMVSerializer(std::string filename);
    TMVSerializer(std::string filename, int maxSize);
    TMVSerializer(std::string filename, int maxSize, bool active);
    void add(const TaylorModelVec & tmv);
    void add(const TaylorModelVec & tmv, std::string name);
    void serialize();
    void activate();
  public:
    std::string filename;
    std::vector<TaylorModelVec> tmvs;
    std::vector<std::string> names;
    int maxSize;
    bool active;
};

extern TMVSerializer *pSerializer;

void addTimeToInfo(std::string name, std::string clockName, std::vector<std::string> & infos);

int findPos(int value, const std::vector<int> *v);

int isIn(int value, const std::vector<int> *v);


void addMyInfo(std::vector<std::string> & info);
void addFlowInfo(std::vector<std::string> & info);

void printComponents(MySettings *settings);

TaylorModelVec getUnitTmv(int varCount);
TaylorModelVec getNVarMParam(int varCount, int paramCount);
TaylorModelVec getNVarMParam(int varCount, std::vector<int> params);

void createOutput(std::vector<MyComponent *> comps, MyComponent & all, 
    Transformer *transformer, MySettings *settings);


#endif /* UTILS_H_ */

