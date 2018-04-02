#include "OutputWriter.h"

OutputWriter::OutputWriter(string name2, int var1, int var2): 
    name(name2), var1Index(var1),var2Index(var2) {
  outfile = new ofstream();
  csvfile = new ofstream();
}
OutputWriter::OutputWriter(string name): 
    name(name) {
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


void OutputWriter::addPreconditioned(vector<MyComponent *> comps, 
    vector<Interval> & domain, MyComponent & all) {
  if(all.output.size() == 0)
    return;
  int dim = all.output[0].tms.size();
  //2 for time + (per dim) 2 for interval + remainder + range
  int csvSize = 2 + 4*dim;
  for(int i = 0; i < csvSize; i++) {
    data2.push_back(vector<string>());
  }
  
  //if it's preconditioned then the first timeinterval is [0,0]
  data2.at(0).push_back("0");
  data2.at(1).push_back("0");
  //size - 1, since first one was added manually
  for(int i = 0; i < all.output.size() - 1; i++) {
    Interval timeInt = domain.at(0) + Interval(i*domain.at(0).sup());
    data2.at(0).push_back(timeInt.getLower(5));
    data2.at(1).push_back(timeInt.getHigher(5));
  }
  
  
  for(int i = 0; i < all.output.size(); i++) {
    TaylorModelVec tmv = all.output[i];
    for(int var = 0; var < tmv.tms.size(); var++) {
      Interval pipe = evalVarToInterval(tmv, domain, var);
      
      data2.at((var +1)*2).push_back(pipe.getLower());
      data2.at((var +1)*2+1).push_back(pipe.getHigher());
      data2.at((dim + 1)*2 + var).push_back(sbuilder() << pipe.width());
      data2.at((dim + 1)*2 + dim + var).push_back(
          sbuilder() << tmv.tms[var].remainder.width());
    }
  }
}

void OutputWriter::addComponents(vector<MyComponent *> comps, 
    vector<Interval> & domain, MyComponent & all, bool isPreconditioned) {
  mlog1("adding components");
  
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
  
  if(isPreconditioned) {
    addPreconditioned(comps, domain, all);
    return;
  }
  
  minc();
  int dim = 0;
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    dim += (*it)->varIndexes.size();
  }
  mlog1(sbuilder() << "dim: " << dim);
  //2 for time + (per dim) 2 for interval + remainder + range
  int csvSize = 2 + 4*dim;
  for(int i = 0; i < csvSize; i++) {
    data2.push_back(vector<string>());
  }
  
  for(int i = 0; i < comps.at(0)->pipes.size(); i++) {
    Interval timeInt = domain.at(0) + Interval(i*domain.at(0).sup());
    data2.at(0).push_back(timeInt.getLower(5));
    data2.at(1).push_back(timeInt.getHigher(5));
  }
  
  mdec();
  
  
  
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    addCompomentData(**it, domain);
  }
}

void OutputWriter::addCompomentData(MyComponent & comp, 
    vector<Interval> & domain) {
  
  int dim = domain.size() - 1;
  int index = 0;
  for(vector<int>::iterator it = comp.varIndexes.begin(); 
      it < comp.varIndexes.end(); it++, index++) {
    int varIndex = *it;
    int compIndex = comp.solveIndexes.at(index);
    mforce(sbuilder() << comp.pipes.size());
    
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
    mforce(sbuilder() << "creating directory '" << pathname << "'");
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
  mlog1(sbuilder() << "writingCSV (" << csvfname << ")");
  csvfile->open(csvfname.c_str());
  
  int values = data2.size();
  int steps = data2.at(0).size();
  int dim = (values - 2) / 4;
  
  //for(int i = 0; i < values; i++) {
  //  mforce(sbuilder() << i << ": " << data2[i][0]);
  //}
  
  int samplePoint = data2[0].size() / 8;
  
  mforce(sbuilder() << "steps: " << data2[0].size());
  mforce(sbuilder() << "sample point: " << (samplePoint));
  double sampleWidth = 0;
  for(int i = 0; i < dim; i++) {
    int varWidthIndex = 2 + dim * 2 + i;
    //mforce(sbuilder() << i << ": " << data2[varWidthIndex][0]);
    sampleWidth += atof(data2[varWidthIndex][samplePoint].c_str());
  }
  mforce(sbuilder() << "sampleWidth: " << sampleWidth);
  
  for(int step = 0; step < steps; step++) {
    //mforce(sbuilder() << "time: " << data2[0][step]);
    //mforce(sbuilder() << "size: " << data2.size());
    
    double stepWidth = 0;
    for(int i = 0; i < dim; i++) {
      int varWidthIndex = 2 + dim * 2 + i;
      //mforce(sbuilder() << i << ": " << data2[varWidthIndex][step]);
      double value = atof(data2[varWidthIndex][step].c_str());
      stepWidth += value;
    }
    //mforce(sbuilder() << "stepWidth: " << stepWidth);
    //mforce(sbuilder() << "ratio: "<< (stepWidth/sampleWidth));
    if(stepWidth/sampleWidth > 10) {
      //mforce(sbuilder() << "breaking at step #"<< step);
      break;
    }
    
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
  mlog1(sbuilder() << "writingInfo (" << infoName << ")");
  ofstream infoFile;
  infoFile.open(infoName.c_str());
  
  //info.push_back(sbuilder() << "shrink wrapping time: " << swTime);
  
  for(vector<string>::iterator it = info.begin(); it < info.end(); it++) {
    //mlog1(*it);  
    infoFile << *it << endl;
  }
  infoFile.close();
}


void OutputWriter::finish() {
  //mlog1(sbuilder() << "finishing writing " << name);
  *outfile << "e\n";
  outfile->close();
}


void OutputWriter::fromFlowstar(list<TaylorModelVec> & flowpipesCompo, 
      list<vector<Interval> > & domains) {
  mforce("from flowstar");
  list<TaylorModelVec>::const_iterator fIt = flowpipesCompo.begin();
  list< vector<Interval> >::const_iterator dIt = domains.begin();
  if(fIt == flowpipesCompo.end())
    return;
  
  int i = 1;
  
  int dim = fIt->tms.size();
  //2 for time + (per dim) 2 for interval + remainder + range
  int csvSize = 2 + 4*dim;
  for(int i = 0; i < csvSize; i++) {
    data2.push_back(vector<string>());
  }
  
  Interval shift = Interval(0);
  
  while(fIt != flowpipesCompo.end()) {
    //mlog("tmv", *fIt);
    //mlog("d", *dIt);
    Interval time = (*dIt)[0];
    time += shift;
    data2.at(0).push_back(time.getLower(5));
    data2.at(1).push_back(time.getHigher(5));
    
    shift += (*dIt)[0].sup();
    //mlog1(time.toString());
    
    const TaylorModelVec & tmv = *fIt;
    vector<Interval> domain = *dIt;
    for(int var = 0; var < tmv.tms.size(); var++) {
      Interval pipe = evalVarToInterval(tmv, domain, var);
      
      data2.at((var +1)*2).push_back(pipe.getLower());
      data2.at((var +1)*2+1).push_back(pipe.getHigher());
      data2.at((dim + 1)*2 + var).push_back(sbuilder() << pipe.width());
      data2.at((dim + 1)*2 + dim + var).push_back(
          sbuilder() << tmv.tms[var].remainder.width());
    }
    
    fIt++;
    dIt++;
    i++;
  }
}

