#include "CompSolver.h"
#include "MyLogger.h"
#include "Utils.h"
#include "Continuous.h"
#include "MyComponent.h"
#include "Polynomial.h"
#include "Transformer.h"
#include "Exceptions.h"

#ifdef DInt
  #include "DoubleInterval.h"
#else
  #include "Interval.h"
#endif


map<string,double> timeMap;
map<string,int> refMap;
map<string, map<int,int> > compStepMap;

int curStep = 0;

/*
map<string,clock_t> start;

void foo(string name) {

}*/

void foo(string name, double d) {
  timeMap[name] += d;
}
void addTime(clock_t start, string name) {
  clock_t end = clock();
  double dur = double(end - start) / CLOCKS_PER_SEC;
  foo(name, dur);
}
void addRef(int ref, string name) {
  refMap[name] += ref;
}
void addTime(clock_t start, clock_t end, string name) {
  double dur = double(end - start) / CLOCKS_PER_SEC;
  foo(name, dur);
}
void bar(string name) {
  /*
  for(map<string,double>::iterator iter = timeMap.begin(); 
      iter != timeMap.end(); ++iter) {
    //cout << "f:" << iter->first << endl;
    //cout << "s:" << iter->second << endl;
  }*/
  cout << name << ": " << timeMap[name] << endl;
}

void bar2(string name) {
  /*
  for(map<string,int>::iterator iter = refMap.begin(); 
      iter != refMap.end(); ++iter) {
    cout << "f:" << iter->first << endl;
    cout << "s:" << iter->second << endl;
  }//*/
  cout << name << ": " << refMap[name] << endl;
}

void Solver::setUp(MySettings *settings, IVP & ivp) {
  //cout << "setting up\n";
  
  if(pSerializer == NULL) {
    //transformer name appended with .txt
    string outName = string(settings->transformer->name);
    if (settings->discardEmptyParams) {
      outName.append("_dis");
    }
    outName.append(".txt");
    mlog1(sbuilder() << "serializer name: " << outName);
    pSerializer = new TMVSerializer(outName);
    //pSerializer->deactivate();
  }
  
  comps = createComponents(settings, ivp.hfOde);
  printComponents(settings);


	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, settings->step, 
	    2*settings->maxOrder);
  settings->step_exp_table = step_exp_table;
  settings->step_end_exp_table = step_end_exp_table;
  
  TaylorModelVec & currentTMV = *ivp.initSet;
  //currentTMV.pushConstantToRemainder();
  
  ivp.domain[0] = step_exp_table[1]; //set domain[0] to timestep
  settings->domain = ivp.domain;
  vector<Interval> & domain = ivp.domain;

  //mlog1("before preparing");
  int stepNum = ceil(settings->time / settings->step);
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(currentTMV, ivp.hfOde, settings);
    (*it)->pipePairs.reserve(stepNum);
  }
  
  all = pGetSystemComponent(comps, currentTMV, ivp.hfOde, settings);

  //(possibly) need to use domain from all
  //since some parameterss may be discarded
  settings->domain = all->dom; 
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
    tstart1(pic_ref_start);
    
    const vector<Interval> & step_exp_table = settings.step_exp_table;
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
    guess.reserve(varCount);
	  for(int i = 0; i < varCount; i++) {
      /*TODO decide to include or not (rem state 1/3)
      if(comp.lastRems.size() != 0) {
        guess.push_back(comp.lastRems[i] * 2);
        continue;
      }*/
      //old code (used on first iteration)
      guess.push_back(settings.estimation[i]);
	  }
	  
	  //set the remainder to be initial guess
	  for(int i = 0; i < varCount; i++) {
      if(comp.isSolveVar(i) == false)
        continue;
		  p.tms[i].remainder = guess[i];
	  }
    //pSerializer->add(p, "inspect");
	  //evaluate this one seperately to get the cutoff measure
	  //and higher order terms
	  
	  //TODO make this one compositional
	  
    //tstart(sc_int_noncomp);
    tend1(pic_ref_start);
    tstart1(pic_ref_first);
	  TaylorModelVec compTemp;
    p.Picard_ctrunc_normal(compTemp, trees, &comp, &settings, pPolyRange, 
      step_exp_table, paramCount, cutoff_threshold);
    //pSerializer->add(compTemp, "inspect");
    tend1(pic_ref_first);
    tstart1(pic_ref_rem);
    
    mlog("tmv", compTemp);
    mlog("tm0", compTemp.tms[0]);
	  mlog("solve", comp.solveIndexes);
	  
    cutoffInt.reserve(varCount);
	  //should be because of the turncated parts and uncertainties (?)
	  for(int i=0; i < varCount; i++) {
  	  //TODO make this one compositional (no need to find remainders for 
  	  //solved variables
	    //bound higher order terms
		  Polynomial polyTemp;
		  polyTemp = compTemp.tms[i].expansion - p.tms[i].expansion;

		  Interval intTemp;
		  polyTemp.intEvalNormal(intTemp, step_exp_table);
		  
		  cutoffInt.push_back(intTemp);
		  
      compTemp.tms[i].remainder += intTemp;
	  }
    //pSerializer->add(compTemp, "inspect");
    tend1(pic_ref_rem);
	  tstart1(pic_ref_subset);
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
    newRemainders.reserve(varCount);
    //mforce1(sbuilder() << "initial ");
		
		int counter = 0;
		//increase until you get subset remainders
    //mforce1(sbuilder() << "notSubset: " << notSubset);
	  while( notSubset ) {
	    //error when can't find decreasing
      if(counter++ == MAX_REFINEMENT_STEPS) {
        throw IntegrationException("max increase couldn't find a remainder (in finding contracting)");
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
    //mforce1(sbuilder() << "counter to find decreasing: " << counter);
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
	  tend1(pic_ref_subset);
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
    const Interval & cutoff_threshold = *settings.cutoff;
    
	  //number of taylor model parameters
    //int paramCount = comp.initSet.tms[0].getParamCount();
    int paramCount = comp.allTMParams.size() + 1;
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
	  
	  
	  //new remainders are stored here
		vector<Interval> newRemainders;
    newRemainders.reserve(varCount);
		
	  bool stop = false;
	  int counter = 0;
	  while( stop == false ) {
  	  if(counter++ == MAX_REFINEMENT_STEPS - 1) {
        break;
	      //throw IntegrationException(sbuilder() << "max refinement steps");
      }
	    //mforce1(sbuilder() << "counter: " << counter);
      p.Picard_only_remainder(newRemainders, trees, &comp, step_exp_table[1]);
      //pSerializer->add(p, "only_rem");
      //mforce("new", newRemainders);
      
	    stop = true;
	    
	    for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
	      newRemainders[i] += cutoffInt[i];
	    }
	    
	    for(int i = 0; i < varCount; i++) {
        if(comp.isSolveVar(i) == false)
          continue;
        //they are actually equal but mpfr libary rounds the ratio to not be 1
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
      //force n number of refinements for all components
      /*
      if(stop) {
        if(counter < 7) {
          stop = false;
          continue;
        }
        cout << counter << " ";
      }//*/
      //if(stop) {
      //  //cout << counter << endl;
      //  addRef(counter, comp.getVarName(&settings));
      //}
  	  //mforce("mr", p.getRemainders());
	  }
    compStepMap[comp.getVarName(&settings)][curStep] = counter;
    /*TODO decide to include or not (rem state 2/3)
    comp.lastRems = newRemainders;
    */
    //mforce1(sbuilder() << "counter: " << counter);
	  //mdec();
	  mlog1("refineRemainderFlow >");
	  mrestore(old);
	}

  bool print = true;
  bool print2 = false;
  void advanceFlow(MyComponent & component, MySettings & settings) {
    //tstart(sc_int_pre);
    mreset(old);
    mdisable();
    mlog1("advance flow <");
    minc();
    
    //mlog("var", component.varIndexes);
    //mlog("links", component.linkVars);
    //mlog("init", component.initSet);
    //tstart(sc_int_all);
    string compName = component.getVarName(&settings);
    /*
    if(print) {
      if(component.getVarName(&settings).find("OR_1") != string::npos) {
        print2 = true;
        mforce1("");
        mforce("init[0]", component.initSet.tms[0]);
        mforce("init[1]", component.initSet.tms[1]);
      }
    }
    */
    //mforce1(sbuilder() << "[" << component.getVarName(&settings) << "]");
    //mforce("init", component.initSet);
    
    //variable when picard approximation is stored
    //clock_t c1 = clock();
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
    //addTime(c1, compName + ".init");
    
    tstart(pic_poly);
    //find the picard polynomial
    for(int i = 1; i <= settings.maxOrder; i++) {
      p.Picard_no_remainder_assign(&component, varCount + 1, i, *settings.cutoff);
      //mforce1(sbuilder() << i);
      //mforce("poly", p);
    }
    
    mlog("poly", p);
    
    //pSerializer->add(p, "int_poly");
    
    p.cutoff(*settings.cutoff);
    tend(pic_poly);
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
	  vector<Interval> cutoffInt;
    
    mlog("before decr", p);
    
    //pSerializer->add(p, "no_rem");
    
    tstart(sc_int_rem);
    tstart(pic_decr);
	  findDecreasingRemainderFlow(p, pPolyRange, trees, component, settings,
	      cutoffInt);
    tend(pic_decr);
    
    
    //pSerializer->add(p, "int_dec");
	  
    mlog("dec", p);
    tstart(pic_ref);
    refineRemainderFlow(p, pPolyRange, trees, component, settings, cutoffInt);
    tend(pic_ref);
    tend(sc_int_rem);
    
	  
    mlog("ref", p);
    //tend(sc_int_all);
    
    pSerializer->add(p, "int_end");
    
    
    //mlog1(component.pipes.size());
    
    mlog("comp.tsp", component.timeStepPipe);
    mlog("tsp", p);
    //mforce("tsp", p.tms[1]);

    if(print) {
      //if(component.getVarName(&settings).find("OR_1") != string::npos) {
        //mforce("tsp[0]", p.tms[0]);
        //mforce("tsp[1]", p.tms[1]);
      //}
    }

    mdec();
    mlog1("advance flow >");
    mrestore(old);    
  }

  void doSingleStep(MyComponent & component, MySettings & settings) {
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
}

void Solver::solveIVP(MySettings *settings, IVP ivp) {
  //mforce1("solve ivp <");
  mreset(old);
  minc();
  mdisable();

  //settings->toOld()->log();
  setUp(settings, ivp);

  //mforce1(sbuilder() << "size: " << comps.size());
  
  
  clock_t integrClock = clock();
  double t;

  for(t = 0; (t + THRESHOLD_HIGH) < settings->time && (true || t < settings->step * 5); t+= settings->step) {
    //mforce1(sbuilder() << "t: " << t);
    if(t > 1) compSolver::print = true;
    cerr << ".";
    try{
      mlog1(sbuilder() << "before transform");
      tstart(sc_transfrom);
      settings->transformer->transform(*all, comps, *settings);
      tend(sc_transfrom);
      mlog1(sbuilder() << "after transform");  
      tstart(int_time);
      //solve components
      for(vector<MyComponent *>::iterator it = comps.begin(); 
          it < comps.end(); it++) {
        //mlog("init", (*it)->initSet);
        //mforce("c", (*it)->compVars);
        //clock_t bef = clock();
        compSolver::doSingleStep(**it, *settings);
        //clock_t aft = clock();
        //double dur = double(aft - bef) / CLOCKS_PER_SEC;
        //foo((*it)->getVarName(settings), dur);
        //mlog("last", (*it)->lastPipe());
      }
      tend(int_time);
      mlog1(sbuilder() << "after integration");
      //don't intdicate that components are solved anymore
      for(int i = 0; i < comps.size(); i++) {
        comps[i]->isSolved = false;
      }
      curStep++;
      //make apppropriate next initial set
      //cout << "integrate: " << timeLookup["int_time"] << endl;
      //cout << "transform: " << timeLookup["sc_transfrom"] << endl;
    } catch(IntegrationException& e) {
      logger.force("IntegrationException caught");
      logger.force(e.what());
      settings->writer->info.push_back(sbuilder() << "reason: " << e.what());
      break;
    }
    //break; //only the first step, REMOVE!
  }
  cerr << endl;

  settings->writer->writeRefinementInfo(comps, settings, compStepMap, curStep);
  

  clock_t end = clock();
  double integrTime = double(end - integrClock) / CLOCKS_PER_SEC;

  
  settings->writer->info.push_back(sbuilder() << "int progress: " << t);

  cout << "computation time: " << integrTime << endl;
  settings->writer->info.push_back(sbuilder() << "computation time: " << integrTime);

  post(settings);
  
  mrestore(old);
  mdec();
  mlog1("solve ivp >");
}

void Solver::post(MySettings *settings) {
  
  //tprint("sc_transfrom");
  //tprint("tr_remap");
  //tprint("tr_comp_pre");

  #ifdef output1
    settings->only1Var = true;
  #endif

  #ifdef no_output
    cout << "not creating my output" << endl;
  #else
    cout << "creating my output" << endl;
    createOutput(comps, *all, settings->transformer, settings);
    //settings->transformer->addInfo(writer.info);
    addMyInfo(settings->writer->info);
    tstart(sc_post_add);
    //cout << "writing output" << endl;
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
  
  if(pSerializer != NULL)
    pSerializer->serialize();
}


//TODO maybe use hfOde_centered
IVP::IVP(ContinuousSystem & system) : 
    hfOde(system.hfOde), initSet(&system.initialSet.tmvPre),
    domain(system.initialSet.domain) {
}
