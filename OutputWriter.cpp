#include "OutputWriter.h"

OutputWriter::OutputWriter(string name2, int var1, int var2): 
    name(name2), var1Index(var1),var2Index(var2) {
  outfile = new ofstream();
  csvfile = new ofstream();
}
void OutputWriter::init() {
  string outfname = "outputs/" + name + ".plt";
  outfile->open(outfname.c_str());
  
  string var1, var2;
  if(var1Index == -1) {
    var1 = "t";
  } else {
    var1 = sbuilder() << "x" << (var1Index + 1);
  }
  if(var2Index == -1) {
    var2 = "t";
  } else {
    var2 = sbuilder() << "x" << (var2Index + 1);
  }
  //*outfile << "set terminal postscript\n";
  //*outfile << "set output './images/" << name << ".eps'\n";
  *outfile << "set term png\n";
  *outfile << "set output './images/" << name << ".png'\n";
  *outfile << "set style line 1 linecolor rgb \"blue\"\n";
  *outfile << "set autoscale\n";
  *outfile << "unset label\n";
  *outfile << "set xtic auto\n";
  *outfile << "set ytic auto\n";
  *outfile << "set xlabel \"" << var1 << "\"\n";
  *outfile << "set ylabel \"" << var2 << "\"\n";
  *outfile << "plot '-' notitle with lines ls 1\n";
}

Interval evalVarToInterval(const TaylorModelVec & tmv, vector<Interval> & domain, int index) {
  Interval temp;
  tmv.tms[index].intEval(temp, domain);
  return temp;
}

Interval evalVarToInterval(const Interval & timeInt, const TaylorModelVec & tmv, vector<Interval> & domain, int index) {
  if(index == -1) {
    return timeInt;
  }
  return evalVarToInterval(tmv, domain, index);
}

void OutputWriter::writeFlowpipe(const Interval & timeInt, 
    const TaylorModelVec & tmv, vector<Interval> & domain) const {
  logger.enable();
  //intervals for the timestep
  Interval xInterval = evalVarToInterval(timeInt, tmv, domain, var1Index);
  Interval yInterval = evalVarToInterval(timeInt, tmv, domain, var2Index);
  
  //corners of 2-d interval box
  *outfile << xInterval.getLower() << " " << yInterval.getLower() << "\n";
  *outfile << xInterval.getHigher() << " " << yInterval.getLower() << "\n";
  *outfile << xInterval.getHigher() << " " << yInterval.getHigher() << "\n";
  *outfile << xInterval.getLower() << " " << yInterval.getHigher() << "\n";
  *outfile << xInterval.getLower() << " " << yInterval.getLower() << "\n";
  *outfile << "\n\n";
  logger.disable();
}

bool checkIndex(const vector<int> comp, int index) {
  for(int i = 0; i < comp.size(); i++) {
    if(comp.at(i) == index)
      return true;
  }
  return false;
}

void OutputWriter::writeFlowpipe(const vector<int> comp, const Interval & timeInt, 
    const TaylorModelVec & tmv, vector<Interval> domain) const {
  if(checkIndex(comp, var1Index) == false & checkIndex(comp, var2Index) == false) {
    return;
  }
  this->writeFlowpipe(timeInt, tmv, domain);
}

void OutputWriter::addComponents(vector<MyComponent *> comps, 
    vector<Interval> & domain) {
  logger.log("adding components");
  
  //for time
  data.push_back(vector<Interval>());
  //for variables
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    for(int i = 0; i < (*it)->varIndexes.size(); i++) {
      data.push_back(vector<Interval>());
    }
  }
  
  for(int i = 0; i < comps.at(0)->pipes.size(); i++) {
    data.at(0).push_back(domain.at(0) + Interval(i*domain.at(0).sup()));
  }
  
  
  
  
  
  
  logger.inc();
  int dim = 0;
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    dim += (*it)->varIndexes.size();
  }
  logger.log(sbuilder() << "dim: " << dim);
  //2 for time + (per dim) 2 for interval + remainder + range
  int csvSize = 2 + 4*dim;
  for(int i = 0; i < csvSize; i++) {
    data2.push_back(vector<string>());
  }
  
  for(int i = 0; i < comps.at(0)->pipes.size(); i++) {
    Interval timeInt = domain.at(0) + Interval(i*domain.at(0).sup());
    data2.at(0).push_back(timeInt.getLower());
    data2.at(1).push_back(timeInt.getHigher());
  }
  
  logger.dec();
  
  
  
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    addCompomentData(**it, domain);
  }
}

void OutputWriter::addCompomentData(MyComponent & comp, 
    vector<Interval> & domain) {
  /*
  int index = 0;
  for(vector<int>::iterator it = comp.varIndexes.begin(); 
      it < comp.varIndexes.end(); it++, index++) {
    int varIndex = *it;
    int compIndex = comp.solveIndexes.at(index);
    
    for(int i = 0; i < comp.pipes.size(); i++) {
      
      Interval pipe = evalVarToInterval(comp.pipes.at(i), domain, compIndex);
      data.at(varIndex+1).push_back(pipe);
    }
  }*/
  
  
  int dim = domain.size() - 1;
  
  int index = 0;
  for(vector<int>::iterator it = comp.varIndexes.begin(); 
      it < comp.varIndexes.end(); it++, index++) {
    int varIndex = *it;
    int compIndex = comp.solveIndexes.at(index);
    
    for(int i = 0; i < comp.pipes.size(); i++) {
      TaylorModelVec tmv = comp.pipes.at(i);
      Interval pipe = evalVarToInterval(tmv, domain, compIndex);
      data.at(varIndex+1).push_back(pipe);
      
      //shift by all previous variables and time
      data2.at((varIndex +1)*2).push_back(pipe.getLower());
      data2.at((varIndex +1)*2+1).push_back(pipe.getHigher());
      data2.at((dim + 1)*2 + varIndex).push_back(sbuilder() << pipe.width());
      data2.at((dim + 1)*2 + dim + varIndex).push_back(
          sbuilder() << tmv.tms[compIndex].remainder.width());
    }
  }
}

void createDir(string pathname) {
  struct stat info;
  if( stat( pathname.c_str(), &info ) != 0 ) {
    logger.force(sbuilder() << "creating directory '" << pathname << "'");
    string cmd = sbuilder() << "mkdir " << pathname;
    system(cmd.c_str());
  } else if( info.st_mode & S_IFDIR ) {
//    printf( "%s is a directory\n", pathname.c_str() );
  } else
    throw std::invalid_argument(sbuilder() << pathname << " is not directory");
}

void OutputWriter::writeCSV() {
  string csvfname = "csvs/" + name + ".csv";
  
  createDir("csvs");
  
  logger.log(sbuilder() << "writingCSV (" << csvfname << ")");
  csvfile->open(csvfname.c_str());
  
  /*
  int vars = data.size();
  int steps = data.at(0).size();
  for(int step = 0; step < steps; step++) {
    bool doBreak = false;
    for(int var = 0; var < vars; var++) {
      if(data.at(var).size() <= step) {
        doBreak = true;
      }
    }
    if(doBreak)
      break;
      
    for(int var = 0; var < vars; var++) {
      Interval varInt = data.at(var).at(step);
      //logger.log(
      //    sbuilder() << varInt.getLower() << "," <<  varInt.getHigher());
      *csvfile << varInt.getLower() << "," << varInt.getHigher() << ",";
    }
    *csvfile << endl;
    //logger.log("");
    
  }
  */
  int values = data2.size();
  int steps = data2.at(0).size();
  for(int step = 0; step < steps; step++) {
    for(int value = 0; value < values; value++) {
      if(data2.at(value).size() <= step) {
        *csvfile << ",";
        continue;
      }
      *csvfile << data2.at(value).at(step) << ",";
    }
    *csvfile << endl;
  }
  
  csvfile->close();
}
void OutputWriter::writeInfo() {
  string infoName = "infos/" + name + ".txt";
  
  createDir("infos");
  logger.log(sbuilder() << "writingInfo (" << infoName << ")");
  ofstream infoFile;
  infoFile.open(infoName.c_str());
  
  info.push_back(sbuilder() << "shrink wrapping time: " << swTime);
  
  for(vector<string>::iterator it = info.begin(); it < info.end(); it++) {
    //logger.log(*it);  
    infoFile << *it << endl;
  }
  infoFile.close();
}


void OutputWriter::finish() {
  //logger.log(sbuilder() << "finishing writing " << name);
  *outfile << "e\n";
  outfile->close();
}
