#include "SimpleComp.h"


namespace simpleComp {
  void picardIter(const vector<int> comp, TaylorModelVec & pipe, const vector<HornerForm> & ode, TaylorModelVec init, int order) {
    logger.log("picardIter <");
    logger.inc();
    TaylorModelVec integratedTMV;
    logger.log("-----------");
    
    vector<TaylorModel> newModels;
    
    for(int i = 0; i < comp.size(); i++) {
      logger.log(sbuilder() << "comp[i] = " << comp[i]);
      logger.log(ode.at(comp[i]).toString());
      
      
      TaylorModel insertedTM;
      //substitute variables with taylormodels
      //computes just the polynomial part 
      ode.at(comp[i]).insert_no_remainder_no_cutoff(insertedTM, pipe, ode.size()+1, order);
      
      TaylorModel tmInt;
      //integrate with respect to time
      insertedTM.integral_no_remainder(tmInt);
      logger.log(sbuilder() << "integratedTMV[" << comp[i] <<"]=" << tmInt.toString(getVNames(3)));
      TaylorModel currentInit = init.tms.at(comp[i]);
      
      
      TaylorModel added;
      tmInt.add(added, currentInit);
      logger.log(added.toString(getVNames(3)));
      newModels.push_back(added);
    }
    for(int i = 0; i < comp.size(); i++) {
      pipe.tms.at(comp[i]) = newModels.at(i);
    }
    logger.logTMV("pipe", pipe);
    
    logger.dec();
    logger.log("picardIter >");
  }


  void computeNewRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    logger.log("picardIterRem <");
    logger.inc();
    TaylorModelVec integratedTMV;
    //logger.logTMV("tmv: ", tmv);
    
    
    for(int i = 0; i < comp.size(); i++) {
      TaylorModel insertedTM;
      
      //logger.logTMV("tmv: ", tmv);
      vector<Interval> polyRange;
      pipe.polyRange(polyRange, domain);
      //logger.logVI("range", polyRange);
      
      Interval cutoff(-1e-55,1e-55); //don't want cutoff atm
      
      //substitute variables in TM, includes remainder
      ode.at(comp[i]).insert(insertedTM, pipe, polyRange, domain, cutoff);
      //logger.log(sbuilder() << "inserted: " << insertedTM.toString(getVNames(3)));
      
      
      TaylorModel tmInt;
      //integrate with respect to time
      insertedTM.integral(tmInt, domain.at(0));
      
      TaylorModel added;
      
      //add the initial conditions
      tmInt.add(added, init.tms[comp[i]]);
      
      Polynomial higherTerms;
      higherTerms = added.expansion - pipe.tms[comp[i]].expansion;
      
      //logger.log(sbuilder() << "tmv: " << tmv.tms[i].toString(getVNames(3)));
      //logger.log(sbuilder() << "added: " << added.toString(getVNames(3)));
      //logger.log(sbuilder() << "higherTerms: " << higherTerms.toString(getVNames(3)));
        
      Interval termsBound;
      higherTerms.intEval(termsBound, domain);
      //logger.log(termsBound.toString());
      //logger.log(sbuilder() << "added rem: " << added.remainder.toString());
      Interval newRemainder = termsBound + added.remainder;
      pipe.tms[comp[i]].remainder = newRemainder;
      //logger.log(newRemainder.toString());
    }
    logger.logTMV("end", pipe);
    logger.dec();
    logger.log("picardIterRem >");
  }

  void findDecreasingRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    logger.log("decRem <");
    logger.inc();

    vector<Interval> guess;
    for(int i=0; i<comp.size(); i++) {
      Interval guessInt(-1,1); //TODO use remainder estimation parameter
      guess.push_back(guessInt);
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
      logger.log("max increase couldn't find a remainder");
      exit(10);
    }
    logger.dec();
    logger.log("decRem >");
  }

  void refineRemainder(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, const TaylorModelVec & init, 
      const vector<Interval> & domain) {
    logger.log("refRem <");
    logger.inc();
    
    
    //store the original remainders
    vector<Interval> prevRems;
    for(int i=0; i<comp.size(); i++) {
      prevRems.push_back(pipe.tms[comp[i]].remainder);
    }
    
    
    for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
      bool redo = false;
      //apply picard operator to get new remainders
      logger.disable();
      computeNewRemainder(comp, pipe, ode, init, domain);
      logger.enable();
      
      for(int i=0; i<comp.size(); i++) {
        
        //STOP_RATIO = 0.99 (from flowstar)
        //stop if all of the remainders have almost equal widths after picard
        if(prevRems.at(i).widthRatio(pipe.tms[comp[i]].remainder) <= STOP_RATIO) {
          redo = true;
        }
        prevRems.at(i) = pipe.tms[comp[i]].remainder;
        //logger.logVI("r", prevRems);
      }
      if(redo == false) {
        break;
      }
      if(j == MAX_REFINEMENT_STEPS-1) {
        logger.reset();
        logger.log("max refinement steps");
        exit(0);
      }
    }
    //logger.logTMV("tmv", tmv);
    logger.dec();
    logger.log("refRem >");
  }

  void advance_step(const vector<int> comp, TaylorModelVec & pipe, const vector<HornerForm> & ode, const TaylorModelVec & init, const vector<Interval> & domain, int order) {
    logger.log("advancing <");
    logger.inc();
    logger.disable();
    
    for (unsigned i=0; i<comp.size(); i++) {
      pipe.tms.at(comp.at(i)) = init.tms.at(comp.at(i));
    }
    
    //apply picard operator (order times)
    for(int i = 0; i < order; i++) {
      picardIter(comp, pipe, ode, init, order);
    }
    
    //inflate the remainders until picard operator is contracting
    findDecreasingRemainder(comp, pipe, ode, init, domain);
    //contract the remainder with picard operator
    refineRemainder(comp, pipe, ode, init, domain);
    logger.enable();
    logger.logTMV("end", pipe);
    logger.dec();
    logger.log("advancing >");
  }
}

SimpleCompReachability::SimpleCompReachability()
: ContinuousReachability() {
  logger.log("simple comp reach constructor");
}



SimpleCompSystem::SimpleCompSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input)
: ContinuousSystem(ode_input, initialSet_input) {
  logger.log("simple comp system cons (ode, pipe)");
}
SimpleCompSystem::SimpleCompSystem(const ContinuousSystem & system)
: ContinuousSystem(system) {
  logger.log("simple comp system cons (sys)");
}

void SimpleCompSystem::addEmptyTM(vector<TaylorModelVec> & pipes) const {  
  TaylorModelVec temp;
  TaylorModel t;
  int dim = tmvOde.tms.size();
  for(int i = 0; i < dim; i++) {
    temp.tms.push_back(t);
  }
  pipes.push_back(temp);
}

void SimpleCompSystem::integrateComponent(const vector<int> comp, vector<TaylorModelVec> & pipes, 
    const OutputWriter writer, const vector<HornerForm> & ode, 
    const TaylorModelVec & init, const vector<Interval> & domain, 
    int order, double step, double time, vector<Interval> step_end_exp_table) const {
      
  TaylorModelVec nextInit = init;
  
  int counter = 0;
  for(double t=THRESHOLD_HIGH; t < time;) {
    logger.log("");
    logger.log(sbuilder() << "t: " << t);
    Interval stepTime = Interval(t, t + step);
    //logger.log(sbuilder() << "counter: " << counter);
    
    //add empty Taylor model in first component
    if(pipes.size() == counter) {
      addEmptyTM(pipes);
    }
    
    //logger.logTMV("nextInit", nextInit);
    TaylorModelVec & pipe = pipes.at(counter);
    //logger.logTMV("start", pipe);
    
    logger.disable();
    simpleComp::advance_step(comp, pipe, ode, nextInit, domain, order);
    logger.enable();
    
    //logger.logTMV("end", pipe);
    //evaluate TM at the end of the timestep
	  pipe.evaluate_t(nextInit, step_end_exp_table);
    //logger.logTMV("next", nextInit);
    
    //output the flowpipe for plotting
    writer.writeFlowpipe(comp, stepTime, pipe, domain);
    
    //advance the time
    t += step;
    
    if(false) {
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
    counter++;
	}
}


void SimpleCompReachability::myRun() {
  logger.log("Simple Comp Run <");
  logger.inc();
  
  
  clock_t begin, end;
	begin = clock();
	
  //copy-paste from flowstar
	compute_factorial_rec(globalMaxOrder+2);
	compute_power_4(globalMaxOrder+2);
	compute_double_factorial(2*globalMaxOrder+4);
  
  
  OutputWriter writer = OutputWriter(sbuilder() << "sc_" << outputFileName, -1, 0);
  writer.init();
  
  //analog of flowstar function
  (*pSystem).my_reach_picard(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold, writer);
  
	end = clock();
	printf("simple comp time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
  logger.dec();
  logger.log("Simple Comp Run >");
}


void SimpleCompSystem::my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const
{
  logger.log("sc reach <");
  logger.inc();
  //copy-paste from flowstar 
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	vector<PolynomialConstraint> dummy_invariant;
  //end of copy-paste
  
  //create output writer object
  //arguments are name and 2 indexes of variables (-1 is time)
  

  TaylorModelVec currentTMV = initialSet.tmvPre;

  //domain of the TM variable
  vector<Interval> domain = initialSet.domain;
  domain.at(0) = step_exp_table[1]; //set domain[0] to timestep
  
  //logger.disable();
  vector<TaylorModelVec> pipes;
  
  for(int i = 0; i < 10; i++) {
    vector<int> comp1;
    comp1.push_back(i);
    logger.log(i);
    integrateComponent(comp1, pipes, writer, hfOde, currentTMV, domain, order, step, time,  step_end_exp_table);
  }
  
  /*
  vector<int> comp1;
  comp1.push_back(0);
  integrateComponent(comp1, pipes, writer, hfOde, currentTMV, domain, order, step, time,  step_end_exp_table);
  
  vector<int> comp2;
  comp2.push_back(1);
  integrateComponent(comp2, pipes, writer, hfOde, currentTMV, domain, order, step, time, step_end_exp_table);
  
  vector<int> comp3;
  comp3.push_back(2);
  integrateComponent(comp3, pipes, writer, hfOde, currentTMV, domain, order, step, time, step_end_exp_table);
  */
  
  
  logger.log(sbuilder() << "order: " << order);
  //logger.enable();
  
  writer.finish();

  logger.dec();
  logger.log("sc reach >");
}
