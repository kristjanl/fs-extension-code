#include "CompSolver.h"
#include "MyLogger.h"
#include "Utils2.h"
#include "Continuous.h"
#include "MyComponent.h"
#include "Polynomial.h"
#include "Interval.h"
#include "Transformer.h"

void Solver::setUp(MySettings *settings, IVP & ivp) {
  cout << "setting up\n";
  
  if(pSerializer == NULL) {
    //transformer name appended with .txt
    string outName = string(settings->transformer->name);
    if (settings->discardEmptyParams)
      outName.append("_dis");
    outName.append(".txt");
    mlog1(sbuilder() << "serializer name: " << outName);
    pSerializer = new TMVSerializer(outName);
  }
  
  comps = createComponents(settings, ivp.hfOde);
  printComponents(settings);


	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, settings->step, 
	    2*settings->order);
  settings->step_exp_table = step_exp_table;
  settings->step_end_exp_table = step_end_exp_table;
  
  TaylorModelVec & currentTMV = *ivp.initSet;
  //currentTMV.pushConstantToRemainder();
  
  //domain of the TM variable
  vector<Interval> & domain = ivp.domain;
  domain[0] = step_exp_table[1]; //set domain[0] to timestep
  settings->domain = domain;

  //TODO possible remove
  //vector<TaylorModelVec> pipes;
  
  //mlog1("before preparing");
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(currentTMV, ivp.hfOde, domain, 
        settings->discardEmptyParams);
  }
  //TODO refactor
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->usingPreconditioning = settings->transformer->isPreconditioned;
  }
  settings->transformer->setIntegrationMapper(comps);
  
  /*
  MyComponent all = getSystemComponent(comps, currentTMV, ivp.hfOde, domain, 
      settings->discardEmptyParams);
  */
  all = pGetSystemComponent(comps, currentTMV, ivp.hfOde, domain, 
      settings->discardEmptyParams);

  settings->domain = all->dom; //need to use all, since some params may be discarded

  //settings->log();
}

namespace compSolver{
  
  void findDecreasingRemainderFlow(TaylorModelVec & p, 
      vector<Interval> & pPolyRange, vector<RangeTree *> & trees, 
      MyComponent & comp, MySettings & settings, 
      vector<Interval> & cutoffInt) {
	  mreset(old);
    mdisable();
	  mlog1("findDecreasingRemainderFlow <");
	  minc();
	  mlog("init", comp.initSet);
    tstart1(ref_start);
    
    const vector<Interval> & step_exp_table = settings.step_exp_table;
    int order = settings.order;
    const Interval & cutoff_threshold = *settings.cutoff;

   
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
    mlog1(sbuilder() << "varCount: " << varCount);

	  //number of taylor model parameters
    int paramCount = comp.getIntergrationParamCount(); 
    mlog1(sbuilder() << "paramCount: " << paramCount);
	  //need to find remainders more efficiently
	  p.polyRangeNormal(pPolyRange, step_exp_table);
	  
	  //initial guess for the remainder
	  vector<Interval> guess;
	  for(int i = 0; i < varCount; i++) {
	    guess.push_back(settings.estimation[i]);
	  }
	  
	  //set the remainder to be initial guess
	  for(int i = 0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
		  p.tms[i].remainder = guess[i];
	  }
	  //evaluate this one seperately to get the cutoff measure
	  //and higher order terms
	  
	  //TODO make this one compositional!!!
	  
    //tstart(sc_int_noncomp);
    tend1(ref_start);
    tstart1(ref_first_picard);
	  TaylorModelVec compTemp;
    p.Picard_ctrunc_normal(compTemp, trees, &comp, pPolyRange, 
      step_exp_table, paramCount, order, cutoff_threshold);
    tend1(ref_first_picard);
    tstart1(ref_rem);
    
    mlog("comTemp", compTemp);
    mlog("tmv", compTemp);
    mlog("tm0", compTemp.tms[0]);
	  mlog("solve", comp.solveIndexes);
	  
	  //should be because of the turncated parts and uncertainties (?)
	  for(int i=0; i < varCount; i++) {
  	  //TODO make this one compositional!!! (no need to find remainders for 
  	  //solved variables
	    //bound higher order terms
		  Polynomial polyTemp;
		  polyTemp = compTemp.tms[i].expansion - p.tms[i].expansion;

		  Interval intTemp;
		  polyTemp.intEvalNormal(intTemp, step_exp_table);
		  
		  cutoffInt.push_back(intTemp);
		  
      compTemp.tms[i].remainder += intTemp;
	  }
    tend1(ref_rem);
	  tstart1(ref_subset);
    //tend(sc_int_noncomp);
	  
	  bool notSubset = false;
	  for(int i=0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
	    //mlog1(guess[i].toString());
	    //mlog1(compTemp.tms[i].remainder.toString());
		  if(compTemp.tms[i].remainder.subseteq(guess[i]) == false ) {
		    notSubset = true;
		    guess[i] *= 2;
		  }
	  }
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
    //mforce1(sbuilder() << "initial ");
		
		int counter = 0;
		//increase until you get subset remainders
	  while( notSubset ) {
	    //error when can't find decreasing
      if(counter++ == MAX_REFINEMENT_STEPS) {
        throw IntegrationException("max increase couldn't find a remainder");
      }
      //mforce1(sbuilder() << "counter to find decreasing: " << counter);
	    
	    //set remainder to guess
		  for(int i = 0; i < varCount; i++) {        
        if(comp.isSolveVar(i) == false)
          continue;
		    p.tms[i].remainder = guess[i];
		  }
		  
		  //compute new remainders
  		p.Picard_only_remainder(newRemainders, trees, &comp, step_exp_table[1]);
  		
  		//check whether any of the remainders were not subsets
	    notSubset = false;
		  for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
          
		    newRemainders[i] += cutoffInt[i];
		    if(newRemainders[i].subseteq(guess[i]) == false ) {
		      notSubset = true;
		      guess[i] *= 2;
		    }
		  }
		}
		
		//set the remainders to such that they decrease
	  for(int i = 0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
        
      //didn't need to go into while loop (newRemainders are empty)
      if(newRemainders.size() == 0) {
        //use the remainders from before loop
        p.tms[i].remainder = compTemp.tms[i].remainder;
      } else {
  	    p.tms[i].remainder = newRemainders[i];
  	  }
	  }
	  tend1(ref_subset);
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
	  //minc();
	  
    const vector<Interval> & step_exp_table = settings.step_exp_table;
    int order = settings.order;
    const Interval & cutoff_threshold = *settings.cutoff;
    
	  //number of taylor model parameters
    //int paramCount = comp.initSet.tms[0].getParamCount();
    int paramCount = comp.allTMParams.size() + 1;
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
	  
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
		
	  bool stop = false;
	  int counter = 0;
	  while( stop == false ) {
  	  if(counter++ == MAX_REFINEMENT_STEPS - 1)
	      throw IntegrationException(sbuilder() << "max refinement steps");
	    //mforce1(sbuilder() << "counter: " << counter);
      p.Picard_only_remainder(newRemainders, trees, &comp, step_exp_table[1]);
      
	    stop = true;
	    
	    for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
	      newRemainders[i] += cutoffInt[i];
	    }
	    
	    for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
	      if(p.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO) {
          //mlog1("redoing");
          stop = false;
          break;
        }
	    }
	    
	    for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
	      p.tms[i].remainder = newRemainders[i];
	    }
  	  //mforce("mr", p.getRemainders());
	  }
	  //mdec();
	  mlog1("refineRemainderFlow >");
	  mrestore(old);
	}

  void advanceFlow(MyComponent & component, MySettings & settings) {
    //tstart(sc_int_pre);
    mreset(old);
    mdisable();
    mlog1("advance flow <");
    minc();
    
    mlog("var", component.varIndexes);
    mlog("links", component.linkVars);
    mlog("init", component.initSet);
    //tstart(sc_int_all);
    
    //variable when picard approximation is stored
    TaylorModelVec & p = component.timeStepPipe;
    mlog("tsp", component.timeStepPipe);
    
    vector<int> & sIndexes = component.solveIndexes;
    
    for (unsigned i=0; i<sIndexes.size(); i++) {
      p.tms[sIndexes[i]] = component.initSet.tms[sIndexes[i]];
    }
    
    mlog("after init", p);
    
    //pSerializer->add(p, "int_start");
    
    int varCount = component.allVars.size();
    //int varCount = p.tms.size();
    
    //tend(sc_int_pre);
    
    
    tstart(sc_int_poly);
    //find the picard polynomial
    for(int i = 1; i <= settings.order; i++) {
      p.Picard_no_remainder_assign(&component, varCount + 1, i, *settings.cutoff);
    }
    mlog("poly", p);
    
    //pSerializer->add(p, "int_poly");
    
    p.cutoff(*settings.cutoff);
    tend(sc_int_poly);
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
    
    mlog("before decr", p);
    
    //pSerializer->add(p, "no_rem");
    
    tstart(sc_int_rem);
    tstart(sc_int_find_dec);
	  findDecreasingRemainderFlow(p, pPolyRange, trees, component, settings,
	      cutoffInt);
    tend(sc_int_find_dec);
    
    
    //pSerializer->add(p, "int_dec");
	  
    mlog("dec", p);
    tstart(sc_int_refine);
    refineRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
    tend(sc_int_refine);
    tend(sc_int_rem);
    
	  
    mlog("ref", p);
    //tend(sc_int_all);
    
    pSerializer->add(p, "int_end");
    
    
    mlog1(component.pipes.size());
    
    mlog("comp.tsp", component.timeStepPipe);
    mlog("tsp", p);
        
    mdec();
    mlog1("advance flow >");
    mrestore(old);    
  }

  void doSingleStep(MyComponent & component, 
    MySettings & settings) {
    //if component has already been solved return
    if(component.isSolved) {
      //mlog1("solved already");
      //mlog("component vars: ", component.compVars);
      return;
    }
    mreset(old);
    mdisable();
    mlog1("doSingleStep <");
    minc();
    for(vector<CompDependency *>::iterator it = component.dependencies.begin(); 
        it < component.dependencies.end(); it++) {
      //mlog1(sbuilder() << "link: " << (*it)->linkVar);
      MyComponent *pComp = (*it)->pComp;
      //solve all dependencies
      compSolver::doSingleStep(*pComp, settings);
    }
    mlog1("solving");
    mlog("component vars: ", component.compVars);
    
    
    //remaps previous components flowpipes
    
    component.remapTimeStepPipe();
  
    //pSerializer->add(component.timeStepPipe, "tsp");
    
    mlog1("remapping done");
    
    //component.log();
    
    advanceFlow(component, settings);
    
    mlog1("single done");
    component.isSolved = true;
    
    mdec();
    mlog1("doSingleStep >");
    mrestore(old);
  }


void createFullyCompositionalOutput(vector<MyComponent *> comps, 
      MyComponent & all, Transformer *transformer, MySettings *settings) {
  mreset(old);
  mdisable();
  mlog1("fully comp");
  
  mlog("allVars", all.allVars);
  mlog("varIndexes", all.varIndexes);
  mlog1(sbuilder() << all.dependencies.size());
  
  int pipes = all.dependencies[0]->pComp->pipePairs.size();
  TaylorModel left, right;
  
  
  for(int step = 0; step < pipes; step++) {
    mlog1(sbuilder() << "step: " << step);
    
    TaylorModelVec left, right;
    for(int var = 0; var < all.allVars.size(); var++) {
      mlog1(sbuilder() << "var: " << var);
      
      
      vector<int> lMapper;
      vector<int> rMapper;
      MyComponent *curComp = &all;
      
      
      while(true) {
        mlog1("in while");
        int pos = findPos(var, &curComp->varIndexes);
        mlog1(sbuilder() << "pos: " << pos);
        
        if(pos < curComp->varIndexes.size()) {
          mlog1("can return straight");
          
          TaylorModel singleLeft = curComp->pipePairs[step]->left.tms[pos];
          TaylorModel singleRight = curComp->pipePairs[step]->right.tms[pos];
          
          mlog("l", singleLeft);
          mlog("r", singleRight);
          
          left.tms.push_back(singleLeft.transform(lMapper));
          right.tms.push_back(singleRight.transform(rMapper));
          break;
        }
        
        
        bool foundDep = false;
        for(int j = 0; j < curComp->dependencies.size(); j++) {
          CompDependency *dep = curComp->dependencies[j];
          mlog("depAll", dep->pComp->allVars);
          
          //is variable in dependency implicit variables
          bool in = isIn(var, &dep->pComp->allVars);
          
          //if variable is not in dependency skip dependency
          if(in == false) {
            continue;
          }
          foundDep = true;
          curComp = dep->pComp;
          
          lMapper = concateMapper(dep->leftMapper, lMapper);
          rMapper = concateMapper(dep->rightMapper, rMapper);
          mlog("lm", lMapper);
          mlog("rm", rMapper);
        }
        if(foundDep == false) {
          stringstream ss;
          ss << "should never get to this point, ";
          ss << "either variable is in componenent or one of the dependencies";
          throw std::runtime_error(ss.str());
        }
      }
    }//end of var
    
    mlog("left", left);
    mlog("right", right);
    all.pipePairs.push_back(new PrecondModel(left, right));
  }//end of step
  mrestore(old);
}

//TODO move to OutputWriter
void createOutput(vector<MyComponent *> comps, MyComponent & all, 
      Transformer *transformer, MySettings *settings) {
  tstart(sc_post_composing);
  if(transformer->isPreconditioned == false)
    return;
  
  mlog1("making system flowpipes");  
  //if fully compositional  
  tstart(tr_remap3);
  
  if(transformer->getType() == TR_SINGLE_COMP) {
    //transformer preconditions single component at a time
    //need to add last flowpipe for all components
    for(int i = 0; i < comps.size(); i++) {
      comps[i]->pipePairs.push_back(
          new PrecondModel(comps[i]->timeStepPipe, comps[i]->unpairedRight));
    }
    createFullyCompositionalOutput(comps, all, transformer, settings);
  } else if(transformer->getType() == TR_ALL_COMP) {
    //tranformer maps everything to system, then precondtions
    //need to remap last integration result, add last flowpipe for system component
    all.remapTimeStepPipe();
    all.pipePairs.push_back(new PrecondModel(all.timeStepPipe, all.unpairedRight));
    
    //mforce3(old3, "all.right3", all.unpairedRight);
  } else {
    //old code might not be compatible, look into it when problems arise
    throw std::runtime_error("don't know how to make output");
  }
  tend(tr_remap3);
  
  mlog1("composing flowpipes");
  for(int i = 0; i < all.pipePairs.size(); i++) {
    cout << ".";
    //pSerializer->add(all.pipePairs[i]->left, "comp_left");
    //pSerializer->add(all.pipePairs[i]->right, "comp_right");
    TaylorModelVec composed = all.pipePairs[i]->composed(settings);
    
    //pSerializer->add(composed, "composed");
	  all.output.push_back(composed);
	  //pSerializer->add(composed, "composed");
  }
  
  cout << all.output.size();
  cout << endl;
  tend(sc_post_composing);
}


}





void Solver::solveIVP(MySettings *settings, IVP ivp) {


  exit(0);
  mforce1("solve ivp <");
  mreset(old);
  minc();
  mdisable();

  //settings->toOld()->log();
  setUp(settings, ivp);

  mforce("vars", this->all->allVars);
  mforce1(sbuilder() << "size: " << comps.size());
  
  
  clock_t integrClock = clock();
  double t;

  for(t = 0; (t + THRESHOLD_HIGH) < settings->time; t+= settings->step) {
    //mlog1(sbuilder() << "t: " << t);
    cerr << ".";
    try{
      mlog1(sbuilder() << "before transform");
      tstart(sc_transfrom);
      settings->transformer->transform(*all, comps, *settings);
      tend(sc_transfrom);
      mlog1(sbuilder() << "after transform");  
      tstart(sc_integrate);
      //solve components
      for(vector<MyComponent *>::iterator it = comps.begin(); 
          it < comps.end(); it++) {
        //mlog("init", (*it)->initSet);
        //mforce("c", (*it)->compVars);
        compSolver::doSingleStep(**it, *settings);
        //mlog("last", (*it)->lastPipe());
      }
      tend(sc_integrate);
      mlog1(sbuilder() << "after integration");
      //don't intdicate that components are solved anymore
      for(int i = 0; i < comps.size(); i++) {
        comps[i]->isSolved = false;
      }
      //make apppropriate next initial set
      //cout << "integrate: " << timeLookup["sc_integrate"] << endl;
      //cout << "transform: " << timeLookup["sc_transfrom"] << endl;
    }catch(IntegrationException& e) {
      logger.force("IntegrationException caught");
      logger.force(e.what());
      settings->writer->info.push_back(sbuilder() << "reason: " << e.what());
      break;
    }
    //break; //only the first step, REMOVE!
  }
  cerr << endl;

  settings->writer->info.push_back(sbuilder() << "int progress: " << t);
  clock_t end = clock();
  double integrTime = double(end - integrClock) / CLOCKS_PER_SEC;
  cout << "computation time: " << integrTime << endl;
  settings->writer->info.push_back(sbuilder() << "computation time: " << integrTime);
  
  //tprint("sc_transfrom");
  //tprint("tr_eval");
  //tprint("tr_remap");
  //tprint("tr_comp_pre");
  

  #ifdef no_output
    cout << "not creating my output" << endl;
  #else
    cout << "creating my output" << endl;
    compSolver::createOutput(comps, *all, settings->transformer, settings);
    //settings->transformer->addInfo(writer.info);
    addMyInfo(settings->writer->info);
    tstart(sc_post_add);
    cout << "writing output" << endl;
    //need to use all->dom, since some parameters maybe be discarded
    settings->writer->addComponents(comps, all->dom, *all, 
        settings->transformer->isPreconditioned);
    tend(sc_post_add);
    tstart(sc_post_write);
    settings->writer->writeCSV();
    settings->writer->writeInfo();
    settings->writer->finish();
    tend(sc_post_write);
  #endif
  //tprint("sc_post");
  //tprint("tr_");
  //cout << "here2" << endl;
  
  
  /*
  if(settings->transformer->isWrapper) {
    //TODO
    //writer.info.push_back(sbuilder() << "shrink wraps: " << swChecker->getCount());
  } else
    writer.info.push_back(sbuilder() << "shrink wraps: 0");
  */
  
  mdec();
  mrestore(old);

  mdec();
  mlog1("solve ivp >");
  
  if(pSerializer != NULL)
    pSerializer->serialize();

}


//TODO maybe use hfOde_centered
IVP::IVP(ContinuousSystem & system) : 
    hfOde(system.hfOde), initSet(&system.initialSet.tmvPre),
    domain(system.initialSet.domain) {
}