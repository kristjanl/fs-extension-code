#ifndef UTILS2_H_
#define UTILS2_H_

//move to Utils file after validating that in can be included everywhere where needed

#include <vector>
#include <string>

class Interval;
class Transformer;
class OutputWriter;

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
#endif /* UTILS_H_ */