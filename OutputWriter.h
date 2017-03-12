#ifndef OUTPUTWRITER_H_
#define OUTPUTWRITER_H_

#include <iostream>
#include <fstream>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream


#include "MyLogger.h"
#include "MyComponent.h"

class TaylorModelVec;
class Interval;

class OutputWriter {
	public:
		OutputWriter(string name, int var1, int var2);
    void init();
    void writeFlowpipe(const Interval & timeInt, const TaylorModelVec & tmv,
        vector<Interval> & domain) const;
    void writeFlowpipe(const vector<int> comp, const Interval & timeInt, 
        const TaylorModelVec & tmv, vector<Interval> domain) const;
    void addCompomentData(MyComponent & comp, vector<Interval> & domain);
    void addComponents(vector<MyComponent *> comps, vector<Interval> & domain);
    void writeCSV();
    void writeInfo();
		void finish();
		vector<string> info;
		double swTime; //separate to info file
	private:
    string name;
    ofstream *outfile;
    ofstream *csvfile;
    int var1Index;
    int var2Index;
    vector<vector <Interval> > data;
    vector<vector <string> > data2;
};

Interval evalVarToInterval(const Interval & timeInt, const TaylorModelVec & tmv, vector<Interval> & domain, int index);

#endif /* OUTPUTWRITER_H_ */
