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
    
    //logger.logTMV("pipe", pipe);
    //logger.logVI("domain", domain);
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
      
      TaylorModel added;
      
      //add the initial conditions
      tmInt.add(added, init.tms[comp[i]]);
      
      //logger.log(sbuilder() << "added: " << added.remainder.toString());
      
      Polynomial higherTerms;
      higherTerms = added.expansion - pipe.tms[comp[i]].expansion;
      
      //logger.log(sbuilder() << "tmv: " << tmv.tms[i].toString(getVNames(3)));
      //logger.log(sbuilder() << "added: " << added.toString(getVNames(3)));
        
      Interval termsBound;
      higherTerms.intEval(termsBound, domain);
      //logger.log(termsBound.toString());
      //logger.log(sbuilder() << "added rem: " << added.remainder.toString());
      Interval newRemainder = termsBound + added.remainder;
      //pipe.tms[comp[i]].remainder = newRemainder;
      remainders.at(comp[i]) = newRemainder;
      //logger.log(newRemainder.toString());
    }
    for(int i = 0; i < pipe.tms.size(); i++) {
      pipe.tms.at(i).remainder = remainders.at(i);    
    }
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
    
    // arbitrary number how many times to try to increase the remainder
    // flowstar doesn't try to increase, so they will have it as 1
    int maxTry = 40;
    bool redo = false;
    for(int j = 0; j < maxTry; j++) {
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
    for (unsigned i=0; i<comp.size(); i++) {
      pipe.tms.at(comp.at(i)) = init.tms.at(comp.at(i));
      //logger.logTM("pipe", pipe.tms.at(comp.at(i)));
    }
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
    //logger.logTMV("dec pipe", pipe);
    //contract the remainder with picard operator
    refineRemainder(comp, pipe, ode, init, domain);
    //logger.logTMV("ref pipe", pipe);
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
    
    const vector<Interval> & step_exp_table = settings.step_exp_table;
    int order = settings.order;
    const Interval & cutoff_threshold = settings.cutoff;

	  //number of taylor model parameters
    int paramCount = comp.initSet.tms[0].getParamCount();
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
	  
	  //need to find remainders more efficiently
	  p.polyRangeNormal(pPolyRange, step_exp_table);
	  
	  //initial guess for the remainder
	  vector<Interval> guess;
	  for(int i = 0; i < varCount; i++) {
	    guess.push_back(settings.estimation[0]); //TODO maybe support each variable
	  }
	  
	  //set the remainder to be initial guess
	  for(int i = 0; i < varCount; i++) {
		  p.tms[i].remainder = guess[i];
	  }
	  
	  //evaluate this one seperately to get the cutoff measures
	  TaylorModelVec tmvTemp;
	  p.Picard_ctrunc_normal(tmvTemp, trees, comp.initSet, pPolyRange, comp.odes, 
	      step_exp_table, paramCount, order, cutoff_threshold);
	  
	  //should be because of the turncated parts and uncertainties (?)
	  for(int i=0; i < varCount; i++) {
		  Polynomial polyTemp;
		  polyTemp = tmvTemp.tms[i].expansion - p.tms[i].expansion;

		  Interval intTemp;
		  polyTemp.intEvalNormal(intTemp, step_exp_table);
		  
		  cutoffInt.push_back(intTemp);
		  
      tmvTemp.tms[i].remainder += intTemp;
      p.tms[i].remainder = tmvTemp.tms[i].remainder;
	  }
	  
	  bool notSubset = false;
	  for(int i=0; i < varCount; i++) {
	    //logger.log(guess[i].toString());
	    //logger.log(tmvTemp.tms[i].remainder.toString());
		  if(tmvTemp.tms[i].remainder.subseteq(guess[i]) == false ) {
		    logger.log("not subset");
		    notSubset = true;
		    guess[i] *= 2;
		  }
	  }
	  logger.logVI("guess", guess);
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
		
		int MAX_TRIES = 40;
	  for(int j = 0; j < MAX_TRIES; j++) {
	    //logger.log(sbuilder() << "j: " << j);
	    
		  for(int i = 0; i < varCount; i++) {
		    p.tms[i].remainder = guess[i];
		  }
	    
	    
	    //logger.logVI("nrb", newRemainders);
	    //logger.logTMV("pb", p);
  		p.Picard_only_remainder(newRemainders, trees, comp.initSet, comp.odes, step_exp_table[1]);
	    //logger.logTMV("pa", p);
	    
	    logger.logVI("guess", guess);
	    logger.logVI("newre", newRemainders);
	    
	    
	    notSubset = false;
		  for(int i = 0; i < varCount; i++) {
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
		  p.Picard_only_remainder(newRemainders, trees, comp.initSet, comp.odes, step_exp_table[1]);
      
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
    TaylorModelVec p = TaylorModelVec(component.initSet);
    
    //logger.logTMV("p", p);
    
    int paramCount = p.tms[0].getParamCount();
    int varCount = p.tms.size();
    
    //find the picard polynomial
    for(int i = 1; i <= settings.order; i++)
      p.Picard_no_remainder_assign(component.initSet, component.odes, paramCount, i, settings.cutoff);
    
    
    p.cutoff(settings.cutoff);
    
    
    //logger.logTMV("p", p);
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
	  
    
	  findDecreasingRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
        
    refineRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
	  
    logger.logTMV("p", p);
    
    logger.log(component.pipes.size());
    
    component.pipes[component.pipes.size() - 1] = p;
        
    logger.restore(old);
  }
  
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
    
    if(settings.useFlow == false)
      smallComp::advance_step(component.solveIndexes, pipe, component.odes, nextInit, component.dom, settings);
    else
      smallComp::advanceFlow(component, settings);
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

void SmallCompReachability::myRun() {
  int old = logger.reset();
  logger.log("Simple Comp Run <");
  logger.inc();
  
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
  //all.log();
  
  //smallComp::foo(*comps.at(0), order, cutoff_threshold, step_exp_table, 
  //    step_end_exp_table);
  
  //single step integration
  //for(double t = 0; t < time; t+= step) {
  
  /*
  ShrinkWrapper sw = ShrinkWrapper(swChecker);
  QRTransformer qr = QRTransformer();
  NullTransformer null = NullTransformer();
  
  Transformer *transformer;
  
  transformer = new QRTransformer();
  transformer = new ShrinkWrapper(swChecker);
  //transformer = new NullTransformer();
  */
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
        //logger.logTMV("0", (*it)->pipes.at(0));
      }
      settings->transformer->transform(all, comps, *settings);
      
    }catch(IntegrationException& e) {
      logger.reset();
      logger.log("IntegrationException caught");
      logger.log(e.what());
      writer.info.push_back(sbuilder() << "reason: " << e.what());
      break;
    }
      
  }
  writer.info.push_back(sbuilder() << "integration time: " << t);
  clock_t end = clock();
  double integrTime = double(end - integrClock) / CLOCKS_PER_SEC;
  logger.log(sbuilder() << "computation time: " << integrTime);
  writer.info.push_back(sbuilder() << "computation time: " << integrTime);
  if(settings->transformer->isWrapper) {
    //TODO
    //writer.info.push_back(sbuilder() << "shrink wraps: " << swChecker->getCount());
  }
  else
    writer.info.push_back(sbuilder() << "shrink wraps: 0");
  
  foo2(comps, all, settings->transformer, settings);
  writer.addComponents(comps, domain, all, settings->transformer->isPreconditioned);
  writer.writeCSV();
  writer.writeInfo();
  
  writer.finish();

  logger.dec();
  logger.log("sc reach >");
}
