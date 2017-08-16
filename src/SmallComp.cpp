#include "SmallComp.h"


namespace smallComp {
  struct MyException : public exception {
    const char * what () const throw () {
      return "C++ Exception";
    }
  };


  void picardIter(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, TaylorModelVec init, MySettings & settings) {
    int old = logger.reset();
    logger.disable();
    logger.log("picardIter <");
    logger.inc();
    TaylorModelVec integratedTMV;
    
    vector<TaylorModel> newModels;
    logger.logTMV("in pipe", pipe);
    
    int paramCount = -1;
    for(vector<TaylorModel>::iterator it = pipe.tms.begin(); 
        it < pipe.tms.end(); it++) {
      paramCount = max(paramCount, (*it).getParamCount());
      logger.log(sbuilder() << "paramCount: " << paramCount);
    }
    
    for(int i = 0; i < comp.size(); i++) {
      //logger.force(sbuilder() << "comp[" << i << "] = " << comp[i]);
      //logger.force(ode.at(comp[i]).toString());
      
      
      
      TaylorModel insertedTM;
      //substitute variables with taylormodels
      //computes just the polynomial part 
      ode.at(comp[i]).insert_no_remainder_no_cutoff(insertedTM, pipe,
          paramCount, settings.order);
      logger.logTM("inserted", insertedTM);
      
      TaylorModel tmInt;
      //integrate with respect to time
      insertedTM.integral_no_remainder(tmInt);
      
      //logger.log(sbuilder() << "integratedTMV[" << comp[i] <<"]=" << 
      //    tmInt.toString(getVNames(tmInt.getParamCount())));
      
      TaylorModel currentInit = init.tms.at(comp[i]);
      
      
      TaylorModel added;
      tmInt.add(added, currentInit);
      logger.logTM("added", added);
      newModels.push_back(added);
    }
    for(int i = 0; i < comp.size(); i++) {
      pipe.tms.at(comp[i]) = newModels.at(i);
    }
    logger.logTMV("picard pipe", pipe);
    
    logger.dec();
    logger.log("picardIter >");
    logger.restore(old);
  }


  void computeNewRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    int old = logger.reset();
    logger.disable();
    logger.log("picardIterRem <");
    logger.inc();
    TaylorModelVec integratedTMV;
    //logger.logTMV("tmv: ", tmv);
    
    logger.logTMV("pipe", pipe);
    logger.logVI("domain", domain);
    vector<Interval> remainders;
    for(int i = 0; i < pipe.tms.size(); i++) {
      remainders.push_back(pipe.tms.at(i).remainder);    
    }
    //logger.logVI("rems", remainders);
    for(int i = 0; i < comp.size(); i++) {
      TaylorModel insertedTM;
      
      //logger.logTMV("tmv: ", tmv);
      vector<Interval> polyRange;
      pipe.polyRange(polyRange, domain);
      //logger.logVI("range", polyRange);
      
      Interval cutoff(-1e-55,1e-55); //don't want cutoff atm
      
      //substitute variables in TM, includes remainder
      //logger.logTM("inserted", insertedTM);
      //logger.logTMV("---pipe", pipe);
      ode.at(comp[i]).insert(insertedTM, pipe, polyRange, domain, cutoff);
      logger.log(sbuilder() << ode.at(comp[i]).toString());
      logger.logTM("inserted", insertedTM);
      
      
      TaylorModel tmInt;
      //integrate with respect to time
      insertedTM.integral(tmInt, domain.at(0));
      
      logger.logTM("integrated", tmInt);
      
      TaylorModel added;
      
      //add the initial conditions
      tmInt.add(added, init.tms[comp[i]]);
      
      logger.logTM("added", added);
      
      Polynomial higherTerms;
      higherTerms = added.expansion - pipe.tms[comp[i]].expansion;
      
      logger.logTM("higher", TaylorModel(higherTerms, Interval(0)));
      Interval termsBound;
      higherTerms.intEval(termsBound, domain);
      logger.log(sbuilder() << "higherBounds: " << termsBound.toString());
      //logger.log(sbuilder() << "added rem: " << added.remainder.toString());
      Interval newRemainder = termsBound + added.remainder;
      //pipe.tms[comp[i]].remainder = newRemainder;
      remainders.at(comp[i]) = newRemainder;
      //logger.log(newRemainder.toString());
    }
    for(int i = 0; i < pipe.tms.size(); i++) {
      pipe.tms.at(i).remainder = remainders.at(i);    
    }
    logger.logVI("rems", remainders);
    //logger.logVI("rems2", remainders);
    //logger.logTMV("end", pipe);
    logger.dec();
    logger.log("picardIterRem >");
    logger.restore(old);
  }

  void findDecreasingRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain, const vector<Interval> & estimation) {
    int old = logger.reset();
    logger.disable();
    logger.log("decRem <");
    logger.inc();

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
      //logger.logVI("guess", guess);
      //logger.logVI("new", pipe.getRemainders());
      redo = false;
      for(int i=0; i<comp.size(); i++) {
        //logger.log(tmv.tms[i].remainder.toString());
        //logger.log(guess.at(i).toString());
        
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
    logger.dec();
    logger.log("decRem >");
    logger.restore(old);
  }

  void refineRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    int old = logger.reset();
    logger.disable();
    logger.log("refRem <");
    logger.inc();
    
    
    //store the original remainders
    vector<Interval> prevRems;
    for(int i=0; i<comp.size(); i++) {
      prevRems.push_back(pipe.tms[comp[i]].remainder);
    }
    logger.logVI("prevRems", prevRems);
    
    for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
      bool redo = false;
      //apply picard operator to get new remainders
      //logger.logTMV("prenew", pipe);
      logger.logVI("rem", prevRems);
      computeNewRemainder(comp, pipe, ode, init, domain);
      //logger.logTMV("aftnew", pipe);
      
      for(int i=0; i<comp.size(); i++) {
        //logger.log(sbuilder() << prevRems[i].width() << ", " << pipe.tms[comp[i]].remainder.width());
        //logger.log(sbuilder() << "ratio[" << i << "]: " << prevRems[i].widthRatio(pipe.tms[comp[i]].remainder));
        //STOP_RATIO = 0.99 (from flowstar)
        //stop if all of the remainders have almost equal widths after picard
        if(prevRems[i].widthRatio(pipe.tms[comp[i]].remainder) <=
            STOP_RATIO) {
          redo = true;
        }
        prevRems.at(i) = pipe.tms[comp[i]].remainder;
        //logger.logVI("r", prevRems);
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
    logger.dec();
    logger.log("refRem >");
    logger.restore(old);
  }

  void advance_step(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain, MySettings & settings) {
    int old = logger.reset();
    logger.disable();
    logger.log("advancing <");
    logger.inc();
    //logger.logTMV("pipe", pipe);
    //logger.logTMV("init", init);
    
    //comy the initial condtions of the variables to be solved to pipe
    for (unsigned i=0; i<comp.size(); i++) {
      pipe.tms.at(comp.at(i)) = init.tms.at(comp.at(i));
    }
    //logger.logTMV("pipe", pipe);
    //logger.logTMV("init", init);
    //logger.log("");
    //apply picard operator (order times)
    for(int i = 0; i < settings.order; i++) {
      picardIter(comp, pipe, ode, init, settings);
      //logger.logTMV("picard pipe", pipe);
    }
    //remove too high terms
    pipe.removeHighTerms(settings.order);
    logger.logTMV("pipe", pipe);
    
    //inflate the remainders until picard operator is contracting
    findDecreasingRemainder(comp, pipe, ode, init, domain, settings.estimation);
    logger.logTMV("dec", pipe);
    //contract the remainder with picard operator
    refineRemainder(comp, pipe, ode, init, domain);
    logger.logTMV("ref", pipe);
    logger.dec();
    logger.log("advancing >");
    logger.restore(old);
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
	  int old = logger.reset();
    logger.disable();
	  logger.log("findDecreasingRemainderFlow <");
	  logger.inc();
	  logger.logTMV("init", comp.initSet);
    
    const vector<Interval> & step_exp_table = settings.step_exp_table;
    int order = settings.order;
    const Interval & cutoff_threshold = settings.cutoff;

	  //number of taylor model parameters
    int paramCount = comp.initSet.tms[0].getParamCount();
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
    logger.log(sbuilder() << "varCount: " << varCount);
	  
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
    
    logger.logTMV("comTemp", compTemp);
    
	  logger.logVi("solve", comp.solveIndexes);
	  
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
	    //logger.log(guess[i].toString());
	    //logger.log(compTemp.tms[i].remainder.toString());
		  if(compTemp.tms[i].remainder.subseteq(guess[i]) == false ) {
		    logger.log("not subset");
		    notSubset = true;
		    guess[i] *= 2;
		  }
	  }
	  logger.logVI("guess", guess);
	  logger.logTMV("compTemp", compTemp);
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
		
	  for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
	    //logger.log(sbuilder() << "j: " << j);
	    
		  for(int i = 0; i < varCount; i++) {        
        if(comp.isSolveVar(i) == false)
          continue;
		    p.tms[i].remainder = guess[i];
		  }
	    //logger.logTMV("pb", p);
  		p.Picard_only_remainder(newRemainders, trees, &comp, step_exp_table[1]);
	    logger.logVI("nrb", newRemainders);
	    
	    logger.logVI("guess", guess);
	    logger.logVI("newre", newRemainders);
	    
	    
	    notSubset = false;
		  for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
		    newRemainders[i] += cutoffInt[i];
		    if(newRemainders[i].subseteq(guess[i]) == false ) {
		      //logger.log("not subset");
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
	  logger.logTMV("last", p);
	  
	  logger.dec();
	  logger.log("findDecreasingRemainderFlow >");
	  logger.restore(old);
	}
	
	
	void refineRemainderFlow(TaylorModelVec & p, vector<Interval> & pPolyRange, 
	    vector<RangeTree *> & trees, MyComponent & comp, MySettings & settings, 
      vector<Interval> & cutoffInt) {
	  int old = logger.reset();
    logger.disable();
	  logger.log("refineRemainderFlow <");
	  logger.inc();
	  
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
	      logger.log(sbuilder() << "ratio: " << p.tms[i].remainder.widthRatio(newRemainders[i]));
	      if(p.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO) {
          //logger.log("redoing");
          stop = false;
          break;
        }
	    }
	    //logger.logTMVRem("pre", p);
	    //logger.logVI("new", newRemainders);
	    logger.log(newRemainders[0].toString());
	    
	    for(int i = 0; i < varCount; i++) {
	      p.tms[i].remainder = newRemainders[i];
	    }
	    if(stop)
	      break;
	  }
	  if(stop == false)
	    throw IntegrationException(sbuilder() << "max refinement steps");
	  	  
	  logger.dec();
	  logger.log("refineRemainderFlow >");
	  logger.restore(old);
	}
  
  void advanceFlow(MyComponent & component, MySettings & settings) {
    int old = logger.reset();
    logger.disable();
    //variable when picard approximation is stored
    //TaylorModelVec p = TaylorModelVec(component.initSet);
    TaylorModelVec p = TaylorModelVec(component.lastPipe());
    logger.logTMV("lastpipe", component.lastPipe());
    
    vector<int> & sIndexes = component.solveIndexes;
    
    
    for (unsigned i=0; i<component.solveIndexes.size(); i++) {
      p.tms.at(sIndexes[i]) = component.initSet.tms.at(sIndexes[i]);
    }
    
    logger.logTMV("after init", p);
    
    int paramCount = p.tms[0].getParamCount();
    int varCount = p.tms.size();
    
    //find the picard polynomial
    for(int i = 1; i <= settings.order; i++) {
      //p.Picard_no_remainder_assign(component.initSet, component.odes, paramCount, i, settings.cutoff);
      p.Picard_no_remainder_assign(&component, paramCount, i, settings.cutoff);
    }
    
    p.cutoff(settings.cutoff);
    
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
    
    logger.logTMV("before decr", p);
	  findDecreasingRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
    

    logger.logTMV("dec", p);
    refineRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
	  
    logger.logTMV("ref", p);
    
    logger.log(component.pipes.size());
    
    component.pipes[component.pipes.size() - 1] = p;
        
    logger.restore(old);
    
  }
  int tempp = 0;
  void singleStepIntegrate(MyComponent & component, MySettings & settings) {
    TaylorModelVec nextInit = component.initSet;
    vector<TaylorModelVec> & pipes = component.pipes;
    
    
    double t=THRESHOLD_HIGH;
    //logger.log(sbuilder() << "t: " << t);
    //logger.logTMV("init", nextInit);
    Interval stepTime = Interval(t, t + settings.step);
    //logger.log(sbuilder() << "counter: " << counter);
    //add empty Taylor model in first component
    
    //logger.logTMV("nextInit", nextInit);
    
    TaylorModelVec & pipe = pipes.at(pipes.size() - 1);
    //logger.logTMV("start", pipe);
    //logger.logTMV("nextInit", nextInit);
    
    //int old2 = logger.reset();
    if(settings.useFlow == false) {
      //logger.logTMV("initp", nextInit);
      //logger.force("plain");
      smallComp::advance_step(component.solveIndexes, pipe, component.odes, nextInit, component.dom, settings);
    } else {
      //logger.logTMV("initf", component.initSet);
      //logger.force("flow");
      smallComp::advanceFlow(component, settings);
    }
    tempp++;
    /*logger.log(tempp);
    logger.log(sbuilder() << "int: " << component.lastPipe().tms[0].remainder.toString());
    logger.restore(old2);
    if(tempp == 3)
      exit(0);
    */
    //logger.logTMV("last2", component.lastPipe());
    
    //logger.logTMV("step pipe", pipe);
    
    //logger.logTMV("end", pipe);
    //evaluate TM at the end of the timestep
    //pipe.evaluate_t(nextInit, settings.step_end_exp_table);
    //logger.logTMV("next", nextInit);
    
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
      //logger.log("solved already");
      //logger.listVi("component vars: ", component.compVars);
      return;
    }
    int old = logger.reset();
    logger.disable();
    logger.log("singleStepPrepareIntegrate <");
    logger.inc();
    for(vector<CompDependency *>::iterator it = component.dependencies.begin(); 
        it < component.dependencies.end(); it++) {
      //logger.log(sbuilder() << "link: " << (*it)->linkVar);
      MyComponent *pComp = (*it)->pComp;
      //solve all dependencies
      smallComp::singleStepPrepareIntegrate(*pComp, settings);
    }
    //logger.log("solving");
    logger.listVi("component vars: ", component.compVars);
    
    
    //remaps previous components flowpipes
    component.remapLastFlowpipe();
    
    //component.log();
    
    singleStepIntegrate(component, settings);
    component.isSolved = true;
    
    logger.dec();
    logger.log("singleStepPrepareIntegrate >");
    logger.restore(old);
  }
  
  
  void foo(MyComponent & comp, int order, const Interval cutoff_threshold, 
      vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table) {
    logger.log("foo");
    
    //variable when picard approximation is stored
    TaylorModelVec p = TaylorModelVec(comp.initSet);
    
    int paramCount = comp.initSet.tms[0].getParamCount();
    int varCount = comp.initSet.tms.size();
    
    //find the picard polynomial
    for(int i = 1; i <= order; i++)
      p.Picard_no_remainder_assign(comp.initSet, comp.odes, paramCount, i, cutoff_threshold);
    p.cutoff(cutoff_threshold);
    
    
    logger.logTMV("p", p);
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
	  
    /*
	  findDecreasingRemainder(p, pPolyRange, trees, comp, step_exp_table, order, 
	      cutoff_threshold, cutoffInt);
	  
	  refineRemainder(p, pPolyRange, trees, comp, step_exp_table, order, 
	      cutoff_threshold, cutoffInt);
	  */
    logger.logTMV("p", p);
    
    //TaylorModelVec temp;
    //p.evaluate_t(temp, step_end_exp_table);
    
    
    //precond(temp, step_end_exp_table);
	  
    exit(0);
  }
}

SmallCompReachability::SmallCompReachability()
: ContinuousReachability() {
  logger.log("simple comp reach constructor");
}



SmallCompSystem::SmallCompSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input)
      : ContinuousSystem(ode_input, initialSet_input) {
  logger.log("simple comp system cons (ode, pipe)");
}
SmallCompSystem::SmallCompSystem(const ContinuousSystem & system, vector< vector<int> > components)
      : ContinuousSystem(system), components(components) {
  logger.log("simple comp system cons (sys)");
}

void foo() {
  logger.reset();
  
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
  logger.logTMV("oldPipe", oldPipe);
	logger.log("-------");
	
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
  
  logger.logTMV("bef", comp1Pipe);
  smallComp::findDecreasingRemainderFlow(comp1Pipe, pPolyRange, trees, comp1, settings, cutoffInt);
  smallComp::refineRemainderFlow(comp1Pipe, pPolyRange, trees, comp1, settings, cutoffInt);
  
  logger.logTMV("dec", comp1Pipe);
  exit(0);
  
  TaylorModelVec comp1Temp;
	vector<RangeTree *> comp1Tree;
  comp1Pipe.Picard_ctrunc_normal(comp1Temp, comp1Tree, &comp1, comp1Range, 
      step_exp_table1, paramCount, order, cutoff_threshold);
  //logger.logVRT("comp1_", comp1Tree);
  logger.logTMV("comp1", comp1Temp);
  for(int i = 0; i < 1; i++)
    comp1Temp.Picard_update_remainder(comp1Tree, &comp1, comp1Domain[0]);
  logger.logTMV("comp1", comp1Temp);
  
  
  
  
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
      
  logger.logTMV("comp2", comp2Temp);
  
  
  for(int i = 0; i < 10; i++)
    comp2Temp.Picard_update_remainder(comp2Tree, &comp2, domain[0]);
  logger.logTMV("comp2", comp2Temp);
  
  
  
  
  
  /*
  vector<Interval> rems;
  comp2Temp.Picard_only_remainder(rems, comp2Tree, &comp2, domain[0]);
  logger.logVI("rems", rems);
  */
  
  
  
  
  
	exit(0);
}



void bar() {
  logger.log("bar");
  deserializeFlows("temp.txt");
  //parseFile("temp.model");
}


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
class gps_position {
private:
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        logger.log("saving");
        ar  & degrees;
        ar  & minutes;
        
        char* abc = NULL;
        int base = 2;
        mpfr_exp_t e;
        abc = mpfr_get_str (NULL, &e, base, 0, number, MPFR_RNDN);
        
        string str(abc);
        logger.log(str);
        ar & str;
        ar & e;
        mpfr_printf ("number = %.17Rg\n", number);
  mpfr_out_str (stdout, 2, 0, number, MPFR_RNDD);
  logger.log("");
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        logger.log("loading");
        int i;
        int j;
        int base = 2;
        string str;
        mpfr_exp_t e;
        ar  & degrees;
        ar  & minutes;
        ar & str;
        ar & e;
        const char *cstr = str.c_str();
        logger.log(sbuilder() << "str: " << str);
        logger.log(sbuilder() << "cstr: " << cstr);
        logger.log(e);
        //char* abc = NULL;
        //ar & abc;
        mpfr_t back;
        mpfr_inits2(intervalNumPrecision, back, (mpfr_ptr) 0);
        mpfr_set_str(back, cstr, base, MPFR_RNDD);
        mpfr_set_exp(back, e);
        mpfr_printf ("back = %.17Rg\n", back);
  mpfr_out_str (stdout, 2, 0, back, MPFR_RNDD);
  logger.log("");
        
        
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    int degrees;
    int minutes;
    float seconds;
    mpfr_t number;
public:
    string s;
    gps_position(){};
    gps_position(int d, int m, float s) :
        degrees(d), minutes(m), seconds(s)
    {
      mpfr_inits2(intervalNumPrecision, number, (mpfr_ptr) 0);
      mpfr_set_d(number, s, MPFR_RNDD);
    }
    void print() const {
      //mpfr_printf ("mpfr = %.17Rg\n", number);
      //mpfr_out_str (stdout, 10, 0, number, MPFR_RNDD);
      logger.log(sbuilder() << "(" << degrees << "," << minutes << "," << seconds << "," << ")" << s);
    }
};

void foo2() {
  // create and open a character archive for output
    std::ofstream ofs("filename");

    // create class instance
    gps_position g(35, 59, 24.567f);
    g.s = "foo";
    Interval i1(2,3);
    logger.log(i1.toString());
    g.print();
    // save data to archive
    {
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << g;
        oa << g;
        oa << i1;
    	// archive and stream closed when destructors are called
    }

    // ... some time later restore the class instance to its orginal state
    gps_position newg, newg2;
    Interval i2;
    {
        // create and open an archive for input
        std::ifstream ifs("filename");
        boost::archive::text_iarchive ia(ifs);
        // read class state from archive
        ia >> newg;
        ia >> newg2;
        ia >> i2;
        // archive and stream closed when destructors are called
    }
    newg.print();
    newg2.print();
  
  logger.log("--------");
  
    int base = 10;
    double input = 25;

    mpfr_t number;
    mpfr_inits2(53, number, (mpfr_ptr) 0);
    mpfr_set_d(number, input, MPFR_RNDD);
    mpfr_const_pi(number, MPFR_RNDD);
    mpfr_printf ("mpfr = %.17Rg\n", number);

    char* str = NULL;
    mpfr_exp_t e;
    str = mpfr_get_str (NULL, &e, base, 0, number, MPFR_RNDN);


    
    cout << "str: " << str << endl;
    cout << "e: " << e << endl;
    cout << "stdout: ";
    mpfr_out_str(stdout, base, 0, number, MPFR_RNDD);
    cout << endl;
    

    mpfr_t back;
    mpfr_inits2(53, back, (mpfr_ptr) 0);
    
    
    char buffer[64];
    sprintf (buffer, ".%s@%ld", str, (long) e);
    cout << buffer << endl;
    mpfr_set_str (back, buffer, base, MPFR_RNDD);
    mpfr_printf ("back = %.17Rg\n", back);
    
    
    mpfr_t dif;
    mpfr_inits2(53, dif, (mpfr_ptr) 0);
    mpfr_sub(dif, number, back, MPFR_RNDU);
    mpfr_printf("dif = %.27Rg\n", dif);
    cout << "dif: ";
    mpfr_out_str(stdout, 2, 0, dif, MPFR_RNDD);
    cout << endl;
    
    
    
    //mpfr_set_str(back, str, base, MPFR_RNDD);
    //mpfr_set_exp(back, e);
    //mpfr_printf("back = %.17Rg\n", back);
    
    /*
    
    string s = "2.5000000000000000";
    const char* str2 = s.c_str();
    mpfr_t back2;
    mpfr_inits2(53, back2, (mpfr_ptr) 0);
    mpfr_set_str(back2, str2, base, MPFR_RNDD);
    mpfr_printf("back2 = %.17Rg\n", back2);
    */
    mpfr_free_str (str);
  exit(0);
  
}

void SmallCompReachability::myRun() {
  int old = logger.reset();
  logger.log("Simple Comp Run <");
  logger.inc();
  
  /*
  foo2();
  exit(0);
  /**/
  
  /*
  logger.reset();
  vector<TaylorModelVec *> plain = pDeserializeFlows("plain.txt");
  logger.log(plain.size());
  exit(0);
  
  vector<TaylorModelVec *> fcomp = pDeserializeFlows("fcomp.txt");
  
  compareFlows(fcomp, plain);
  
  //foo();
  
  //bar();
  exit(1);
  /**/
  
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
		
  logger.dec();
  logger.log("Simple Comp Run >");
  logger.restore(old);
}

void foo2(vector<MyComponent *> comps, MyComponent & all, 
      Transformer *transformer, MySettings *settings) {
  
  if(transformer->isPreconditioned == false)
    return;
  
  for(int i = 0; i < all.pipes.size(); i++) {
    cout << ".";
    TaylorModelVec & left = all.pipes[i];
    if(i == 0) {
      continue;
      all.output.push_back(left);
    }
    TaylorModelVec & right = all.pipePairs[i-1]->right;  
    TaylorModelVec composed;
    vector<Interval> rightRange;
	  right.polyRange(rightRange, settings->domain);
	  left.insert_ctrunc(composed, right, rightRange, settings->domain, settings->order, 0);
	  all.output.push_back(composed);
  }
  cout << endl;
  
}

void SmallCompSystem::my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const {
  logger.log("sc reach <");
  logger.inc();
  logger.log(sbuilder() << "# of components: " <<components.size());
  
  vector<MyComponent *> comps = createComponents(components, hfOde);
  
  
  //copy-paste from flowstar 
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
  //end of copy-paste
  
  //create output writer object
  //arguments are name and 2 indexes of variables (-1 is time)
  

  TaylorModelVec currentTMV = initialSet.tmvPre;

  //domain of the TM variable
  vector<Interval> domain = initialSet.domain;
  domain.at(0) = step_exp_table[1]; //set domain[0] to timestep
  
  settings->writer = &writer;
  settings->order = order;
  settings->step_exp_table = step_exp_table;
  settings->step_end_exp_table = step_end_exp_table;
  settings->domain = domain;
  
  vector<TaylorModelVec> pipes;
  
  if(precondition == SHRINK_WRAPPING || settings->transformer->isPreconditioned) {
    for(int i = 0; i < comps.size(); i++) {
      comps.at(i)->retainEmptyParams = true;
    }
  }
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(currentTMV, hfOde, domain, 
        settings->transformer->isPreconditioned);
  }
  MyComponent all = getSystemComponent(comps, currentTMV, hfOde, domain, 
      settings->transformer->isPreconditioned);
  

  clock_t integrClock = clock();
  double t;
  for(t = 0; t < time; t+= step) {
    //logger.log(sbuilder() << "t: " << t);
    cerr << ".";
    
    try{
      for(vector<MyComponent *>::iterator it = comps.begin(); 
          it < comps.end(); it++) {
        //logger.logTMV("init", (*it)->initSet);
        smallComp::singleStepPrepareIntegrate(**it, *settings);
        //logger.logTMV("last", (*it)->lastPipe());
      }
      settings->transformer->transform(all, comps, *settings);
      
    }catch(IntegrationException& e) {
      logger.reset();
      logger.log("IntegrationException caught");
      logger.log(e.what());
      writer.info.push_back(sbuilder() << "reason: " << e.what());
      break;
    }
    break; //only the first step, REMOVE!
  }
  writer.info.push_back(sbuilder() << "integration time: " << t);
  clock_t end = clock();
  double integrTime = double(end - integrClock) / CLOCKS_PER_SEC;
  logger.log(sbuilder() << "computation time: " << integrTime);
  writer.info.push_back(sbuilder() << "computation time: " << integrTime);
  if(settings->transformer->isWrapper) {
    //TODO
    //writer.info.push_back(sbuilder() << "shrink wraps: " << swChecker->getCount());
  } else
    writer.info.push_back(sbuilder() << "shrink wraps: 0");
  
  foo2(comps, all, settings->transformer, settings);
  writer.addComponents(comps, domain, all, settings->transformer->isPreconditioned);
  writer.writeCSV();
  writer.writeInfo();
  
  writer.finish();

  logger.dec();
  logger.log("sc reach >");
  //comps[0]->serializeFlows();
  
  
  
  serializeFlows(comps[0], "fcomp.txt");
  exit(0);
  vector<TaylorModelVec> & parsed = deserializeFlows("fcomp.txt");
  compareFlows(comps[0]->pipes, parsed);
  exit(0);
  
}
