#include "SimpleImpl.h"


TaylorModelVec picardIter(const vector<HornerForm> & ode, TaylorModelVec init, TaylorModelVec prev, int order) {
  mlog1("picardIter <");
  minc();
  TaylorModelVec integratedTMV;
	mlog1(sbuilder() << "last: " << prev.toString(getVNames(3)));
  mlog1("-----------");
  
  for (unsigned i=0; i<ode.size(); i++) {
		mlog1(ode.at(i).toString());
    TaylorModel insertedTM;
    
    //substitute variables with taylormodels
    //computes just the polynomial part 
    ode.at(i).insert_no_remainder_no_cutoff(insertedTM, prev, ode.size()+1, order);
    
    TaylorModel tmInt;
    //integrate with respect to time
    insertedTM.integral_no_remainder(tmInt);
    mlog1(sbuilder() << "integratedTMV[" << i <<"]=" << tmInt.toString(getVNames(3)));
    
    integratedTMV.tms.push_back(tmInt);
    
	}
  mlog1("-----------");
  
  TaylorModelVec ret;
  //add the initial conditions
  integratedTMV.add(ret, init);
  mlog("nextIter", ret);
  mdec();
  mlog1("picardIter >");
  return ret;
}


void computeNewRemainder(const vector<HornerForm> & ode, const TaylorModelVec & init, TaylorModelVec & tmv, const vector<Interval> & domain) {
  mlog1("picardIterRem <");
  minc();
  TaylorModelVec integratedTMV;
  //mlog("tmv: ", tmv);
  
  for (unsigned i=0; i<ode.size(); i++) {
    //mlog1(sbuilder() << "-------- " << i << " --------");
		//mlog1(sbuilder() << "ode: " << ode.at(i).toString());
    TaylorModel insertedTM;
    
    //mlog("tmv: ", tmv);
    vector<Interval> polyRange;
    tmv.polyRange(polyRange, domain);
    //mlog("range", polyRange);
    
    Interval cutoff(-1e-55,1e-55); //don't want cutoff atm
    
    //substitute variables in TM, includes remainder
    ode.at(i).insert(insertedTM, tmv, polyRange, domain, cutoff);
    //mlog1(sbuilder() << "inserted: " << insertedTM.toString(getVNames(3)));
    
    
    TaylorModel tmInt;
    //integrate with respect to time
    insertedTM.integral(tmInt, domain.at(0));
    
    TaylorModel added;
    
    //add the initial conditions
    tmInt.add(added, init.tms[i]);
    
		Polynomial higherTerms;
		higherTerms = added.expansion - tmv.tms[i].expansion;
    
    //mlog1(sbuilder() << "tmv: " << tmv.tms[i].toString(getVNames(3)));
    //mlog1(sbuilder() << "added: " << added.toString(getVNames(3)));
    //mlog1(sbuilder() << "higherTerms: " << higherTerms.toString(getVNames(3)));
      
    Interval termsBound;
    higherTerms.intEval(termsBound, domain);
    //mlog1(termsBound.toString());
    //mlog1(sbuilder() << "added rem: " << added.remainder.toString());
    Interval newRemainder = termsBound + added.remainder;
    tmv.tms[i].remainder = newRemainder;
    //mlog1(newRemainder.toString());
	}
  mlog("end", tmv);
  mdec();
  mlog1("picardIterRem >");
}

void findDecreasingRemainder(const vector<HornerForm> & ode, const TaylorModelVec & init, TaylorModelVec & tmv, const vector<Interval> & domain) {
  mlog1("decRem <");
  minc();

  vector<Interval> guess;
  for(int i=0; i<tmv.tms.size(); i++) {
    Interval guessInt(-1,1); //TODO use remainder estimation parameter
    guess.push_back(guessInt);
  }
  for(int i=0; i<tmv.tms.size(); i++) {
    tmv.tms[i].remainder = guess.at(i);
  }
  
  // arbitrary number how many times to try to increase the remainder
  // flowstar doesn't try to increase, so they will have it as 1
  int maxTry = 40;
  bool redo = false;
  for(int j = 0; j < maxTry; j++) {
    computeNewRemainder(ode, init, tmv, domain);
    redo = false;
    for(int i=0; i<tmv.tms.size(); i++) {
      //mlog1(tmv.tms[i].remainder.toString());
      //mlog1(guess.at(i).toString());
      
      //new remainder is not a subset of the old one - so it's bad
      if(tmv.tms[i].remainder.subseteq(guess.at(i)) == false) {
        redo = true;
        guess.at(i) *= 2; //TODO think/find a better way to increase it
      }
      //reset the remainder
      tmv.tms[i].remainder = guess.at(i);
    }
    if(redo == false) {
      break;
    }
  }
  if(redo) {
    mreset(old2);
    mlog1("max increase couldn't find a remainder");
    exit(10);
  }
  mdec();
  mlog1("decRem >");
}

void refineRemainder(const vector<HornerForm> & ode, const TaylorModelVec & init, TaylorModelVec & tmv, const vector<Interval> & domain) {
  mlog1("refRem <");
  minc();
  
  
  //store the original remainders
  vector<Interval> prevRems;
  for(int i=0; i<tmv.tms.size(); i++) {
    prevRems.push_back(tmv.tms[i].remainder);
  }
  
  
  for(int j = 0; j < MAX_REFINEMENT_STEPS; j++) {
    bool redo = false;
    
    //apply picard operator to get new remainders
    computeNewRemainder(ode, init, tmv, domain);
    
    for(int i=0; i<tmv.tms.size(); i++) {
      
      //STOP_RATIO = 0.99 (from flowstar)
      //stop if all of the remainders have almost equal widths after picard
      if(prevRems.at(i).widthRatio(tmv.tms[i].remainder) <= STOP_RATIO) {
        redo = true;
      }
      
      //update the remainder for next step
      prevRems.at(i) = tmv.tms[i].remainder;
    }
    if(redo == false) {
      break;
    }
    if(j == MAX_REFINEMENT_STEPS-1) {
      mreset(old2);
      mlog1("max refinements steps");
      exit(11);
    }
  }
  //mlog("tmv", tmv);
  mdec();
  mlog1("refRem >");
}

TaylorModelVec advance_step(const vector<HornerForm> & ode, const TaylorModelVec & init, const vector<Interval> & domain, int order) {
  mlog1("advancing <");
  minc();
  
  TaylorModelVec tmv = init;
  
  //apply picard operator (order times)
  for(int i = 0; i < order; i++) {
    tmv = picardIter(ode, init, tmv, 2);
  }
  
  //inflate the remainders until picard operator is contracting
  findDecreasingRemainder(ode, init, tmv, domain);
  
  //contract the remainder with picard operator
  refineRemainder(ode, init, tmv, domain);
  
  mdec();
  mlog1("advancing >");
  return tmv;
}

SimpleImplReachability::SimpleImplReachability()
: ContinuousReachability() {
}

void SimpleImplReachability::myRun() {
  mlog1("Simple Run <");
  minc();
  
  clock_t begin, end;
	begin = clock();
	
  //copy-paste from flowstar
	compute_factorial_rec(globalMaxOrder+2);
	compute_power_4(globalMaxOrder+2);
	compute_double_factorial(2*globalMaxOrder+4);
  
  
  OutputWriter writer = OutputWriter(sbuilder() << "si_" << outputFileName, -1, 9);
  writer.init();
  
  //analog of flowstar function
  (*pSystem).my_reach_picard(flowpipes, step, time, orders[0], precondition, estimation, bPrint, stateVarNames, cutoff_threshold, writer);
  
  end = clock();
	printf("simple impl time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
  
  mdec();
  mlog1("Simple Run >");
}


SimpleImplSystem::SimpleImplSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input)
: ContinuousSystem(ode_input, initialSet_input) {
}
SimpleImplSystem::SimpleImplSystem(const ContinuousSystem & system)
: ContinuousSystem(system) {
  mlog1("simple impl system cons (sys)");
}


void SimpleImplSystem::my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const
{
  mlog1("si reach <");
  minc();
  
  //copy-paste from flowstar 
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order);
	vector<PolynomialConstraint> dummy_invariant;
  //end of copy-paste
  
  //create output writer object
  //arguments are name and 2 indexes of variables (-1 is time)
  

  //domain of the TM variable
  TaylorModelVec currentTMV = initialSet.tmvPre;
  vector<Interval> domain = initialSet.domain;
  domain.at(0) = step_exp_table[1]; //set domain[0] to timestep
  
	for(double t=THRESHOLD_HIGH; t < time;) {
    mlog1(sbuilder() << "t: " << t);
    //flowpipe for the timestep
    TaylorModelVec step_flowpipe = advance_step(hfOde, currentTMV, domain, order);
    
    //evaluate TM at the end of the timestep
    TaylorModelVec nextModel;
	  step_flowpipe.evaluate_t(nextModel, step_end_exp_table);
    //output the flowpipe for plotting
    writer.writeFlowpipe(step_exp_table[1] + t, step_flowpipe, domain);
    
    //advance the time
    t += step;
    

    currentTMV = nextModel;
    if(false) {
			fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
			break;
		}
	}
  writer.finish();

  mdec();
  mlog1("si reach >");
}
