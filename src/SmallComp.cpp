#include "SmallComp.h"


namespace smallComp {
  struct MyException : public exception {
    const char * what () const throw () {
      return "C++ Exception";
    }
  };


  void picardIter(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, TaylorModelVec init, MySettings & settings) {
    mreset(old);
    mdisable();
    mlog1("picardIter <");
    minc();
    TaylorModelVec integratedTMV;
    
    vector<TaylorModel> newModels;
    mlog("in pipe", pipe);
    
    int paramCount = -1;
    for(vector<TaylorModel>::iterator it = pipe.tms.begin(); 
        it < pipe.tms.end(); it++) {
      paramCount = max(paramCount, (*it).getParamCount());
      mlog1(sbuilder() << "paramCount: " << paramCount);
    }
    
    for(int i = 0; i < comp.size(); i++) {
      //mforce(sbuilder() << "comp[" << i << "] = " << comp[i]);
      //mforce(ode.at(comp[i]).toString());
      
      
      
      TaylorModel insertedTM;
      //substitute variables with taylormodels
      //computes just the polynomial part 
      ode.at(comp[i]).insert_no_remainder_no_cutoff(insertedTM, pipe,
          paramCount, settings.order);
      mlog("inserted", insertedTM);
      
      TaylorModel tmInt;
      //integrate with respect to time
      insertedTM.integral_no_remainder(tmInt);
      
      //mlog1(sbuilder() << "integratedTMV[" << comp[i] <<"]=" << 
      //    tmInt.toString(getVNames(tmInt.getParamCount())));
      
      TaylorModel currentInit = init.tms.at(comp[i]);
      
      
      TaylorModel added;
      tmInt.add(added, currentInit);
      mlog("added", added);
      newModels.push_back(added);
    }
    for(int i = 0; i < comp.size(); i++) {
      pipe.tms.at(comp[i]) = newModels.at(i);
    }
    mlog("picard pipe", pipe);
    
    mdec();
    mlog1("picardIter >");
    mrestore(old);
  }


  void computeNewRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    mreset(old);
    mdisable();
    mlog1("picardIterRem <");
    minc();
    TaylorModelVec integratedTMV;
    //mlog("tmv: ", tmv);
    
    mlog("pipe", pipe);
    mlog("domain", domain);
    vector<Interval> remainders;
    for(int i = 0; i < pipe.tms.size(); i++) {
      remainders.push_back(pipe.tms.at(i).remainder);    
    }
    //mlog("rems", remainders);
    for(int i = 0; i < comp.size(); i++) {
      TaylorModel insertedTM;
      
      //mlog("tmv: ", tmv);
      vector<Interval> polyRange;
      pipe.polyRange(polyRange, domain);
      //mlog("range", polyRange);
      
      Interval cutoff(-1e-55,1e-55); //don't want cutoff atm
      
      //substitute variables in TM, includes remainder
      //mlog("inserted", insertedTM);
      //mlog("---pipe", pipe);
      ode.at(comp[i]).insert(insertedTM, pipe, polyRange, domain, cutoff);
      mlog1(sbuilder() << ode.at(comp[i]).toString());
      mlog("inserted", insertedTM);
      
      
      TaylorModel tmInt;
      //integrate with respect to time
      insertedTM.integral(tmInt, domain.at(0));
      
      mlog("integrated", tmInt);
      
      TaylorModel added;
      
      //add the initial conditions
      tmInt.add(added, init.tms[comp[i]]);
      
      mlog("added", added);
      
      Polynomial higherTerms;
      higherTerms = added.expansion - pipe.tms[comp[i]].expansion;
      
      mlog("higher", TaylorModel(higherTerms, Interval(0)));
      Interval termsBound;
      higherTerms.intEval(termsBound, domain);
      mlog1(sbuilder() << "higherBounds: " << termsBound.toString());
      //mlog1(sbuilder() << "added rem: " << added.remainder.toString());
      Interval newRemainder = termsBound + added.remainder;
      //pipe.tms[comp[i]].remainder = newRemainder;
      remainders.at(comp[i]) = newRemainder;
      //mlog1(newRemainder.toString());
    }
    for(int i = 0; i < pipe.tms.size(); i++) {
      pipe.tms.at(i).remainder = remainders.at(i);    
    }
    mlog("rems", remainders);
    //mlog("rems2", remainders);
    //mlog("end", pipe);
    mdec();
    mlog1("picardIterRem >");
    mrestore(old);
  }

  void findDecreasingRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain, const vector<Interval> & estimation) {
    mreset(old);
    mdisable();
    mlog1("decRem <");
    minc();

    vector<Interval> guess;
    for(int i=0; i<comp.size(); i++) {
      guess.push_back(estimation[0]);
    }
    for(int i=0; i<comp.size(); i++) {
      pipe.tms.at(comp[i]).remainder = guess.at(i);
    }
    
    // flowstar doesn't try to increase, so they will have it as 1
    bool redo = false;
    for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
      computeNewRemainder(comp, pipe, ode, init, domain);
      //mlog("guess", guess);
      //mlog("new", pipe.getRemainders());
      redo = false;
      for(int i=0; i<comp.size(); i++) {
        //mlog1(tmv.tms[i].remainder.toString());
        //mlog1(guess.at(i).toString());
        
        //new remainder is not a subset of the old one - so it's bad
        if(pipe.tms[comp[i]].remainder.subseteq(guess.at(i)) == false) {
          redo = true;
          guess.at(i) *= 2; //TODO think/find a better way to increase it
        }
        //reset the remainder
        pipe.tms[comp[i]].remainder = guess.at(i);
      }
      if(redo == false) {
        break;
      }
    }
    if(redo) {
      throw IntegrationException("max increase couldn't find a remainder");
    }
    mdec();
    mlog1("decRem >");
    mrestore(old);
  }

  void refineRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    mreset(old);
    mdisable();
    mlog1("refRem <");
    minc();
    
    
    //store the original remainders
    vector<Interval> prevRems;
    for(int i=0; i<comp.size(); i++) {
      prevRems.push_back(pipe.tms[comp[i]].remainder);
    }
    mlog("prevRems", prevRems);
    
    for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
      bool redo = false;
      //apply picard operator to get new remainders
      //mlog("prenew", pipe);
      mlog("rem", prevRems);
      computeNewRemainder(comp, pipe, ode, init, domain);
      //mlog("aftnew", pipe);
      
      for(int i=0; i<comp.size(); i++) {
        //mlog1(sbuilder() << prevRems[i].width() << ", " << pipe.tms[comp[i]].remainder.width());
        //mlog1(sbuilder() << "ratio[" << i << "]: " << prevRems[i].widthRatio(pipe.tms[comp[i]].remainder));
        //STOP_RATIO = 0.99 (from flowstar)
        //stop if all of the remainders have almost equal widths after picard
        if(prevRems[i].widthRatio(pipe.tms[comp[i]].remainder) <=
            STOP_RATIO) {
          redo = true;
        }
        prevRems.at(i) = pipe.tms[comp[i]].remainder;
        //mlog("r", prevRems);
      }
      if(redo == false) {
        break;
      }
      if(j == MAX_REFINEMENT_STEPS-1) {
        //throw std::invalid_argument("max refinement steps");
        throw IntegrationException(sbuilder() << 
          "max refinement steps");
      }
    }
    mdec();
    mlog1("refRem >");
    mrestore(old);
  }

  void advance_step(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain, MySettings & settings) {
    mreset(old);
    mdisable();
    mlog1("advancing <");
    minc();
    //mlog("pipe", pipe);
    //mlog("init", init);
    
    //comy the initial condtions of the variables to be solved to pipe
    for (unsigned i=0; i<comp.size(); i++) {
      pipe.tms.at(comp.at(i)) = init.tms.at(comp.at(i));
    }
    //mlog("pipe", pipe);
    //mlog("init", init);
    //mlog1("");
    //apply picard operator (order times)
    for(int i = 0; i < settings.order; i++) {
      picardIter(comp, pipe, ode, init, settings);
      //mlog("picard pipe", pipe);
    }
    //remove too high terms
    pipe.removeHighTerms(settings.order);
    mlog("pipe", pipe);
    
    //inflate the remainders until picard operator is contracting
    findDecreasingRemainder(comp, pipe, ode, init, domain, settings.estimation);
    mlog("dec", pipe);
    //contract the remainder with picard operator
    refineRemainder(comp, pipe, ode, init, domain);
    mlog("ref", pipe);
    mdec();
    mlog1("advancing >");
    mrestore(old);
  }
  
  //TODO duplicated in MyComponent
  void addEmptyTM(vector<TaylorModelVec> & pipes, int dim) {  
    TaylorModelVec temp;
    TaylorModel t;
    for(int i = 0; i < dim; i++) {
      temp.tms.push_back(t);
    }
    pipes.push_back(temp);
  }
  
  void findDecreasingRemainderFlow(TaylorModelVec & p, vector<Interval> & pPolyRange, 
	    vector<RangeTree *> & trees, MyComponent & comp, MySettings & settings, 
      vector<Interval> & cutoffInt) {
	  mreset(old);
    mdisable();
	  mlog1("findDecreasingRemainderFlow <");
	  minc();
	  mlog("init", comp.initSet);
    
    const vector<Interval> & step_exp_table = settings.step_exp_table;
    int order = settings.order;
    const Interval & cutoff_threshold = settings.cutoff;

	  //number of taylor model parameters
    int paramCount = comp.initSet.tms[0].getParamCount();
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
    mlog1(sbuilder() << "varCount: " << varCount);
	  
	  //need to find remainders more efficiently
	  p.polyRangeNormal(pPolyRange, step_exp_table);
	  
	  //initial guess for the remainder
	  vector<Interval> guess;
	  for(int i = 0; i < varCount; i++) {
	    guess.push_back(settings.estimation[0]); //TODO maybe support each variable
	  }
	  
	  //set the remainder to be initial guess
	  for(int i = 0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
		  p.tms[i].remainder = guess[i];
	  }
	  //evaluate this one seperately to get the cutoff measures
	  TaylorModelVec compTemp;
    p.Picard_ctrunc_normal(compTemp, trees, &comp, pPolyRange, 
      step_exp_table, paramCount, order, cutoff_threshold);
    
    mlog("comTemp", compTemp);
    mlog("tmv", compTemp);
    mlog("tm0", compTemp.tms[0]);
	  mlog("solve", comp.solveIndexes);
	  
	  //should be because of the turncated parts and uncertainties (?)
	  for(int i=0; i < varCount; i++) {
		  Polynomial polyTemp;
		  polyTemp = compTemp.tms[i].expansion - p.tms[i].expansion;

		  Interval intTemp;
		  polyTemp.intEvalNormal(intTemp, step_exp_table);
		  
		  cutoffInt.push_back(intTemp);
		  
      compTemp.tms[i].remainder += intTemp;
      p.tms[i].remainder = compTemp.tms[i].remainder;
	  }
	  
	  bool notSubset = false;
	  for(int i=0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
	    //mlog1(guess[i].toString());
	    //mlog1(compTemp.tms[i].remainder.toString());
		  if(compTemp.tms[i].remainder.subseteq(guess[i]) == false ) {
		    mlog1("not subset");
		    notSubset = true;
		    guess[i] *= 2;
		  }
	  }
	  mlog("guess", guess);
	  mlog("compTemp", compTemp);
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
		
	  for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
	    //mlog1(sbuilder() << "j: " << j);
	    
		  for(int i = 0; i < varCount; i++) {        
        if(comp.isSolveVar(i) == false)
          continue;
		    p.tms[i].remainder = guess[i];
		  }
	    //mlog("pb", p);
  		p.Picard_only_remainder(newRemainders, trees, &comp, step_exp_table[1]);
	    mlog("nrb", newRemainders);
	    
	    mlog("guess", guess);
	    mlog("newre", newRemainders);
	    
	    
	    notSubset = false;
		  for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
		    newRemainders[i] += cutoffInt[i];
		    if(newRemainders[i].subseteq(guess[i]) == false ) {
		      //mlog1("not subset");
		      notSubset = true;
		      guess[i] *= 2;
		    }
		  }
		  if(notSubset == false)
		    break;
		}
    if(notSubset) {
      throw IntegrationException("max increase couldn't find a remainder");
    }
		
	  for(int i = 0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
	    p.tms[i].remainder = newRemainders[i];
	  }
	  mlog("last", p);
	  
	  mdec();
	  mlog1("findDecreasingRemainderFlow >");
	  mrestore(old);
	}
	
	
	void refineRemainderFlow(TaylorModelVec & p, vector<Interval> & pPolyRange, 
	    vector<RangeTree *> & trees, MyComponent & comp, MySettings & settings, 
      vector<Interval> & cutoffInt) {
	  mreset(old);
    mdisable();
	  mlog1("refineRemainderFlow <");
	  minc();
	  
    const vector<Interval> & step_exp_table = settings.step_exp_table;
    int order = settings.order;
    const Interval & cutoff_threshold = settings.cutoff;
    
	  //number of taylor model parameters
    int paramCount = comp.initSet.tms[0].getParamCount();
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
	  
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
		
	  bool stop = true;
	  for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
      p.Picard_only_remainder(newRemainders, trees, &comp, step_exp_table[1]);
      
	    for(int i = 0; i < varCount; i++)
	      newRemainders[i] += cutoffInt[i];
	    
	    stop = true;
	    for(int i = 0; i < varCount; i++) {
	      mlog1(sbuilder() << "ratio: " << p.tms[i].remainder.widthRatio(newRemainders[i]));
	      if(p.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO) {
          //mlog1("redoing");
          stop = false;
          break;
        }
	    }
	    //mlog("new", newRemainders);
	    mlog1(newRemainders[0].toString());
	    
	    for(int i = 0; i < varCount; i++) {
	      p.tms[i].remainder = newRemainders[i];
	    }
	    if(stop)
	      break;
	  }
	  if(stop == false)
	    throw IntegrationException(sbuilder() << "max refinement steps");
	  	  
	  mdec();
	  mlog1("refineRemainderFlow >");
	  mrestore(old);
	}
  
  void advanceFlow(MyComponent & component, MySettings & settings) {
    mreset(old);
    mdisable();
    mlog("init", component.initSet);
    
    //variable when picard approximation is stored
    TaylorModelVec & p = component.timeStepPipe;
    mlog("tsp", component.timeStepPipe);
    
    vector<int> & sIndexes = component.solveIndexes;
    
    for (unsigned i=0; i<component.solveIndexes.size(); i++) {
      p.tms.at(sIndexes[i]) = component.initSet.tms.at(sIndexes[i]);
    }
    
    mlog("after init", p);
    
    int paramCount = p.tms[0].getParamCount();
    int varCount = p.tms.size();
    
    
    //find the picard polynomial
    for(int i = 1; i <= settings.order; i++) {
      //p.Picard_no_remainder_assign(component.initSet, component.odes, paramCount, i, settings.cutoff);
      p.Picard_no_remainder_assign(&component, paramCount, i, settings.cutoff);
    }
    
    mlog("after", p);
    
    
    p.cutoff(settings.cutoff);
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
    
    mlog("before decr", p);
	  findDecreasingRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
    

    mlog("dec", p);
    refineRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
	  
    mlog("ref", p);
    
    
    mlog1(component.pipes.size());
    
    mlog("comp.tsp", component.timeStepPipe);
    mlog("tsp", p);
        
    mrestore(old);
    
  }
  
  void singleStepIntegrate(MyComponent & component, MySettings & settings) {
    TaylorModelVec nextInit = component.initSet;
    vector<TaylorModelVec> & pipes = component.pipes;
    
    
    double t=THRESHOLD_HIGH;
    //mlog1(sbuilder() << "t: " << t);
    //mlog("init", nextInit);
    Interval stepTime = Interval(t, t + settings.step);
    //mlog1(sbuilder() << "counter: " << counter);
    //add empty Taylor model in first component
    //mlog("nextInit", nextInit);

    //mlog("start", pipe);
    //mlog("nextInit", nextInit);
    
    
    if(settings.useFlow == false) {
      //mlog("initp", nextInit);
      //mforce("plain");
      TaylorModelVec & pipe = component.timeStepPipe;
      smallComp::advance_step(component.solveIndexes, pipe, component.odes, 
          nextInit, component.dom, settings);
    } else {
      //mlog("initf", component.initSet);
      //mforce("flow");
      smallComp::advanceFlow(component, settings);
    }
    //mlog("last2", component.lastPipe());
    
    //mlog("step pipe", pipe);
    
    //mlog("end", pipe);
    //evaluate TM at the end of the timestep
    //pipe.evaluate_t(nextInit, settings.step_end_exp_table);
    //mlog("next", nextInit);
    
    //output the flowpipe for plotting
    //settings.writer.writeFlowpipe(component.varIndexes, stepTime, pipe, component.dom);
    
    //advance the time
    t += settings.step;
    
    if(false) {
      fprintf(stdout, 
          "Terminated -- The remainder estimation is not large enough.\n");
    }
    //exit(0);
  }
  
  void singleStepPrepareIntegrate(MyComponent & component, MySettings & settings) {
    //if component has already been solved return
    if(component.isSolved) {
      //mlog1("solved already");
      //mlog("component vars: ", component.compVars);
      return;
    }
    mreset(old);
    mdisable();
    mlog1("singleStepPrepareIntegrate <");
    minc();
    for(vector<CompDependency *>::iterator it = component.dependencies.begin(); 
        it < component.dependencies.end(); it++) {
      //mlog1(sbuilder() << "link: " << (*it)->linkVar);
      MyComponent *pComp = (*it)->pComp;
      //solve all dependencies
      smallComp::singleStepPrepareIntegrate(*pComp, settings);
    }
    mlog1("solving");
    mlog("component vars: ", component.compVars);
    
    
    //remaps previous components flowpipes
    component.remapLastFlowpipe();
    mlog1("remapping done");
    
    //component.log();
    
    singleStepIntegrate(component, settings);
    mlog1("single done");
    component.isSolved = true;
    
    mdec();
    mlog1("singleStepPrepareIntegrate >");
    mrestore(old);
  }
  
  
  void foo(MyComponent & comp, int order, const Interval cutoff_threshold, 
      vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table) {
    mlog1("foo");
    
    //variable when picard approximation is stored
    TaylorModelVec p = TaylorModelVec(comp.initSet);
    
    int paramCount = comp.initSet.tms[0].getParamCount();
    int varCount = comp.initSet.tms.size();
    
    //find the picard polynomial
    for(int i = 1; i <= order; i++)
      p.Picard_no_remainder_assign(comp.initSet, comp.odes, paramCount, i, cutoff_threshold);
    p.cutoff(cutoff_threshold);
    
    
    mlog("p", p);
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
	  
    /*
	  findDecreasingRemainder(p, pPolyRange, trees, comp, step_exp_table, order, 
	      cutoff_threshold, cutoffInt);
	  
	  refineRemainder(p, pPolyRange, trees, comp, step_exp_table, order, 
	      cutoff_threshold, cutoffInt);
	  */
    mlog("p", p);
    
    //TaylorModelVec temp;
    //p.evaluate_t(temp, step_end_exp_table);
    
    
    //precond(temp, step_end_exp_table);
	  
    exit(0);
  }
}

SmallCompReachability::SmallCompReachability()
: ContinuousReachability() {
  mlog1("simple comp reach constructor");
}



SmallCompSystem::SmallCompSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input)
      : ContinuousSystem(ode_input, initialSet_input) {
  mlog1("simple comp system cons (ode, pipe)");
}
SmallCompSystem::SmallCompSystem(const ContinuousSystem & system, vector< vector<int> > components)
      : ContinuousSystem(system), components(components) {
  mlog1("simple comp system cons (sys)");
}

void foo() {
  mreset(old);
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	int order = 1;
	Interval cutoff_threshold = Interval(0.0,0.0);
	
  //TaylorModelVec pipe = parseTMV("my models{a - a*t  + a*t^2 * 0.5 - a*t^3 * 0.5 * 0.333 +  [-1,1]}");
  TaylorModelVec pipe = parseTMV("my models{a + [-0.001,0.001]}");
  TaylorModelVec init = parseTMV("my models{a}");
  vector<HornerForm> ode = parseHFFromPoly("my hfs {-1*a}");
  vector<int> comp = parseiVec("my iv <0>");
  vector<Interval> domain = parseIVec("my Iv <[0,0.1],[-1,1]>");
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, domain[0].sup(), 2*order);
	int paramCount = pipe.tms[0].getParamCount();
  
  
  exit(0);
  TaylorModelVec oldPipe = TaylorModelVec(pipe);
  smallComp::computeNewRemainder(comp, oldPipe, ode, init, domain);
  mlog("oldPipe", oldPipe);
	mlog1("-------");
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	MyComponent comp1;
	TaylorModelVec comp1Pipe = parseTMV("my models{a + [-0.001,0.001]}");
  comp1.initSet = parseTMV("my models{a + [-1.1,11.1]}");
  comp1.odes = parseHFFromPoly("my hfs {-1*a}");
  comp1.solveIndexes = parseiVec("my iv <0>");
  vector<Interval> comp1Domain = parseIVec("my Iv <[0,0.1],[-1,1]>");
  
	vector<Interval> step_exp_table1, step_end_exp_table1;
	construct_step_exp_table(step_exp_table1, step_end_exp_table1, comp1Domain[0].sup(), 2*order);
  
  
  vector<Interval> comp1Range;
  comp1Pipe.polyRangeNormal(comp1Range, step_exp_table1);
  
  MySettings settings;
  
  settings.step_exp_table = step_exp_table1;
  settings.order = order;
  settings.cutoff = cutoff_threshold;
	settings.estimation.push_back(Interval(-1,1));
  
	vector<Interval> pPolyRange;
	vector<RangeTree *> trees;
	vector<Interval> cutoffInt;
  
  mlog("bef", comp1Pipe);
  smallComp::findDecreasingRemainderFlow(comp1Pipe, pPolyRange, trees, comp1, settings, cutoffInt);
  smallComp::refineRemainderFlow(comp1Pipe, pPolyRange, trees, comp1, settings, cutoffInt);
  
  mlog("dec", comp1Pipe);
  exit(0);
  
  TaylorModelVec comp1Temp;
	vector<RangeTree *> comp1Tree;
  comp1Pipe.Picard_ctrunc_normal(comp1Temp, comp1Tree, &comp1, comp1Range, 
      step_exp_table1, paramCount, order, cutoff_threshold);
  //mlog("comp1_", comp1Tree);
  mlog("comp1", comp1Temp);
  for(int i = 0; i < 1; i++)
    comp1Temp.Picard_update_remainder(comp1Tree, &comp1, comp1Domain[0]);
  mlog("comp1", comp1Temp);
  
  
  
  
  return;
  MyComponent comp2;
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	TaylorModelVec comp2Pipe = parseTMV("my models{a - a*t + a*t^2*0.5 +[-0.001,0.001], b + t*a*b}");
  comp2Pipe.tms[0] = comp1Temp.tms[0];
  comp2.initSet = parseTMV("my models{a,b}");
  comp2.odes = parseHFFromPoly("my hfs {0, a*b}");
  comp2.solveIndexes = parseiVec("my iv <1>");
  vector<Interval> comp2Domain = parseIVec("my Iv <[0,0.1],[-1,1],[-1,1]>");
	
  vector<Interval> comp2Range;
  comp2Pipe.polyRangeNormal(comp2Range, step_exp_table);
  
  TaylorModelVec comp2Temp;
	vector<RangeTree *> comp2Tree;
  /*comp2Pipe.Picard_ctrunc_normal(comp2Temp, comp2Tree, comp2Vars, comp2Init, comp2Range, comp2Ode, 
      step_exp_table, paramCount, order, cutoff_threshold);*/
  comp2Pipe.Picard_ctrunc_normal(comp2Temp, comp2Tree, &comp2, comp2Range, step_exp_table, paramCount, order, cutoff_threshold);
      
  mlog("comp2", comp2Temp);
  
  
  for(int i = 0; i < 10; i++)
    comp2Temp.Picard_update_remainder(comp2Tree, &comp2, domain[0]);
  mlog("comp2", comp2Temp);
  
  
  
  
  
  /*
  vector<Interval> rems;
  comp2Temp.Picard_only_remainder(rems, comp2Tree, &comp2, domain[0]);
  mlog("rems", rems);
  */
  
  
  
  
  
	exit(0);
}



void SmallCompReachability::myRun() {
  mreset(old);
  mlog1("Simple Comp Run <");
  minc();
  
  clock_t begin, end;
	begin = clock();
	
	
  //copy-paste from flowstar
	compute_factorial_rec(globalMaxOrder+2);
	compute_power_4(globalMaxOrder+2);
	compute_double_factorial(2*globalMaxOrder+4);
  
  
  OutputWriter writer = OutputWriter(sbuilder() << outputFileName, -1, 0);
  writer.init();
  
  //analog of flowstar function
  pSystem->my_reach_picard(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold, writer);
  
	end = clock();
	//printf("simple comp time cost: %lf\n", 
	//    (double)(end - begin) / CLOCKS_PER_SEC);
		
  mdec();
  mlog1("Simple Comp Run >");
  mrestore(old);
}


//TODO move to OutputWriter
void createOutput(vector<MyComponent *> comps, MyComponent & all, 
      Transformer *transformer, MySettings *settings) {
      
  mlog1("here");
  if(transformer->isPreconditioned == false)
    return;
  
  for(int i = 0; i < all.pipePairs.size(); i++) {
    cout << ".";
    TaylorModelVec composed = all.pipePairs[i]->composed(settings);
    //mlog("com", com);
	  all.output.push_back(composed);
  }
  cout << endl;
}

void SmallCompSystem::my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const {
  mlog1("sc reach <");
  minc();
  mlog1(sbuilder() << "# of components: " <<components.size());
  
  if(pSerializer == NULL) {
    //pSerializer = new TMVSerializer("comp.txt", 8*11);
    pSerializer = new TMVSerializer("comp.txt");
  }
  
  vector<MyComponent *> comps = createComponents(components, hfOde);
  
  //copy-paste from flowstar 
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
  //end of copy-paste
  
  //create output writer object
  //arguments are name and 2 indexes of variables (-1 is time)
  

  TaylorModelVec currentTMV = initialSet.tmvPre;
  //currentTMV.pushConstantToRemainder();
  
  //domain of the TM variable
  vector<Interval> domain = initialSet.domain;
  domain.at(0) = step_exp_table[1]; //set domain[0] to timestep
  
  settings->writer = &writer;
  settings->order = order;
  settings->step_exp_table = step_exp_table;
  settings->step_end_exp_table = step_end_exp_table;
  settings->domain = domain;
  
  vector<TaylorModelVec> pipes;
  
  if(precondition == SHRINK_WRAPPING || 
      settings->transformer->isPreconditioned) {
    for(int i = 0; i < comps.size(); i++) {
      comps.at(i)->retainEmptyParams = true;
    }
  }
  
  //mlog1("before preparing");
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(currentTMV, hfOde, domain);
  }
  MyComponent all = getSystemComponent(comps, currentTMV, hfOde, domain);
  
  settings->transformer->transform(all, comps, *settings);
  clock_t integrClock = clock();
  double t;
  for(t = 0; (t + THRESHOLD_HIGH) < time; t+= step) {
    //mlog1(sbuilder() << "t: " << t);
    cerr << ".";
    
    try{
      tstart(sc_integrate);
      //solve components
      for(vector<MyComponent *>::iterator it = comps.begin(); 
          it < comps.end(); it++) {
        //mlog("init", (*it)->initSet);
        smallComp::singleStepPrepareIntegrate(**it, *settings);
        //mlog("last", (*it)->lastPipe());
      }
      tend(sc_integrate);
      //don't intdicate that components are solved anymore
      for(int i = 0; i < comps.size(); i++) {
        comps[i]->isSolved = false;
      }
      //make apppropriate next initial set
      tstart(sc_transfrom);
      settings->transformer->transform(all, comps, *settings);
      tend(sc_transfrom);
      //cout << "integrate: " << timeLookup["sc_integrate"] << endl;
      //cout << "transform: " << timeLookup["sc_transfrom"] << endl;
    }catch(IntegrationException& e) {
      mreset(old);
      mlog1("IntegrationException caught");
      mlog1(e.what());
      writer.info.push_back(sbuilder() << "reason: " << e.what());
      break;
    }
    //break; //only the first step, REMOVE!
  }
  cerr << endl;
  
  tprint("tr_part");
  
  settings->transformer->addInfo(writer.info);
  
  writer.info.push_back(sbuilder() << "int progress: " << t);
  clock_t end = clock();
  double integrTime = double(end - integrClock) / CLOCKS_PER_SEC;
  cout << "computation time: " << integrTime << endl;
  writer.info.push_back(sbuilder() << "computation time: " << integrTime);
  
  
  #ifdef no_output
    cout << "not creating my output" << endl;
  #else
    cout << "creating my output" << endl;
    createOutput(comps, all, settings->transformer, settings);
    writer.addComponents(comps, domain, all, 
        settings->transformer->isPreconditioned);
    writer.writeCSV();
    writer.writeInfo();
    writer.finish();
  #endif
  
  
  /*
  if(settings->transformer->isWrapper) {
    //TODO
    //writer.info.push_back(sbuilder() << "shrink wraps: " << swChecker->getCount());
  } else
    writer.info.push_back(sbuilder() << "shrink wraps: 0");
  */
  
  

  mdec();
  mlog1("sc reach >");
  
  if(pSerializer != NULL)
    pSerializer->serialize();
  //comps[0]->serializeFlows();
  
  
  //serializeFlows(comps[0], "plain.txt");
  //vector<TaylorModelVec> & parsed = deserializeFlows("fcomp.txt");
  //compareFlows(comps[0]->pipes, parsed);
  exit(0);
  
}
