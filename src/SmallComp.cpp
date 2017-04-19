#include "SmallComp.h"


namespace smallComp {
  struct MyException : public exception {
    const char * what () const throw () {
      return "C++ Exception";
    }
  };


  void picardIter(const vector<int> comp, TaylorModelVec & pipe, 
      const vector<HornerForm> & ode, TaylorModelVec init, int order) {
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
          paramCount, order);
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
    /*
    vector<int> pow1;
    pow1.push_back(0);
    pow1.push_back(1);
    pow1.push_back(0);
    Monomial m1(Interval(2), pow1);
    
    vector<int> pow2;
    pow2.push_back(0);
    pow2.push_back(0);
    pow2.push_back(1);
    Monomial m2(Interval(2), pow2);
        
    list<Monomial> l1;
    l1.push_back(m1);
    Polynomial p1(l1);
    list<Monomial> l2;
    l2.push_back(m2);
    Polynomial p2(l2);
    
    pipe.tms.at(0).remainder = Interval();
    pipe.tms.at(0).expansion = p1;
    pipe.tms.at(1).remainder = Interval();
    pipe.tms.at(1).expansion = p2;*/
    
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
      const vector<Interval> & domain) {
    int old = logger.reset();
    logger.disable();
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
      computeNewRemainder(comp, pipe, ode, init, domain);
      //logger.logTMV("aftnew", pipe);
      
      for(int i=0; i<comp.size(); i++) {
        
        //STOP_RATIO = 0.99 (from flowstar)
        //stop if all of the remainders have almost equal widths after picard
        if(prevRems.at(i).widthRatio(pipe.tms[comp[i]].remainder) <=
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
      const vector<Interval> & domain, int order) {
    int old = logger.reset();
    logger.disable();
    logger.log("advancing <");
    logger.inc();
    for (unsigned i=0; i<comp.size(); i++) {
      pipe.tms.at(comp.at(i)) = init.tms.at(comp.at(i));
      //logger.logTM("pipe", pipe.tms.at(comp.at(i)));
    }
    //apply picard operator (order times)
    for(int i = 0; i < order; i++) {
      picardIter(comp, pipe, ode, init, order);
      //logger.logTMV("picard pipe", pipe);
    }
    
    //inflate the remainders until picard operator is contracting
    findDecreasingRemainder(comp, pipe, ode, init, domain);
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
  
  void singleStepIntegrate(MyComponent & component, const OutputWriter writer, 
      int order, double step, double time, 
      vector<Interval> step_end_exp_table) {
    TaylorModelVec nextInit = component.initSet;
    vector<TaylorModelVec> & pipes = component.pipes;
    
    int counter = 0;
    for(double t=THRESHOLD_HIGH; t < time;) {
      //logger.log(sbuilder() << "t: " << t);
      //logger.logTMV("init", nextInit);
      Interval stepTime = Interval(t, t + step);
      //logger.log(sbuilder() << "counter: " << counter);
      //add empty Taylor model in first component
      
      //logger.logTMV("nextInit", nextInit);
      
      TaylorModelVec & pipe = pipes.at(pipes.size() - 1);
      //logger.logTMV("start", pipe);
      //logger.logTMV("nextInit", nextInit);
      
      smallComp::advance_step(component.solveIndexes, pipe, component.odes,
          nextInit, component.dom, order);
      
      //logger.logTMV("step pipe", pipe);
      
      //logger.logTMV("end", pipe);
      //evaluate TM at the end of the timestep
      pipe.evaluate_t(nextInit, step_end_exp_table);
      //logger.logTMV("next", nextInit);
      
      //output the flowpipe for plotting
      writer.writeFlowpipe(component.varIndexes, stepTime, pipe, component.dom);
      
      //advance the time
      t += step;
      
      if(false) {
        fprintf(stdout, 
            "Terminated -- The remainder estimation is not large enough.\n");
        break;
      }
      counter++;
    }
    //exit(0);
  }
  
  
  void integrateComponent2(MyComponent & component, const OutputWriter writer, 
      int order, double step, double time, vector<Interval> step_end_exp_table) {
    TaylorModelVec nextInit = component.initSet;
    vector<TaylorModelVec> & pipes = component.pipes;
    
    int counter = 0;
    for(double t=THRESHOLD_HIGH; t < time;) {
      logger.log("");
      logger.log(sbuilder() << "t: " << t);
      logger.logTMV("init2", nextInit);
      Interval stepTime = Interval(t, t + step);
      //logger.log(sbuilder() << "counter: " << counter);
      //add empty Taylor model in first component
      if(pipes.size() == counter) {
        if(component.varIndexes.size() != component.compVars.size()) {
          logger.reset();
          logger.log("not initial component, but no flowpipes");
          exit(0);
        }
        smallComp::addEmptyTM(pipes, component.varIndexes.size());
      }
      //logger.logTMV("nextInit", nextInit);
      TaylorModelVec & pipe = pipes.at(counter);
      //logger.logTMV("start", pipe);
      
      smallComp::advance_step(component.solveIndexes, pipe, component.odes, nextInit, component.dom, order);
      
      logger.logTMV("step pipe", pipe);
      
      //logger.logTMV("end", pipe);
      //evaluate TM at the end of the timestep
      pipe.evaluate_t(nextInit, step_end_exp_table);
      //logger.logTMV("next", nextInit);
      
      //output the flowpipe for plotting
      writer.writeFlowpipe(component.varIndexes, stepTime, pipe, component.dom);
      
      //advance the time
      t += step;
      
      if(false) {
        fprintf(stdout, "Terminated -- The remainder estimation is not large enough.\n");
        break;
      }
      counter++;
    }
  }
  
  TaylorModelVec shrinkWrapTMV(TaylorModelVec tmv, double factor) {
    TaylorModelVec noConst = TaylorModelVec(tmv);
    noConst.rmConstant(); 
    
    TaylorModelVec constantPart;
    noConst.mul(constantPart, -1);
    constantPart.add_assign(tmv);
    for(int i = 0; i < constantPart.tms.size(); i++) {
      constantPart.tms.at(i).remainder = Interval(0); //nullify the remainder
    }
    
    vector<TaylorModel> qTM;
    for(int i = 0; i < noConst.tms.size(); i++) {
      TaylorModel tm = noConst.tms.at(i);
      TaylorModel mult;
      tm.mul(mult, factor);
      //take only the polynomial part
      qTM.push_back(TaylorModel(mult.expansion, Interval(0)));
    }
    
    TaylorModelVec wrappedSet = TaylorModelVec(qTM);
    wrappedSet.add_assign(constantPart);
    return wrappedSet;
  }
  
  //wraps all.swInput
  void shrinkWrapSet(MyComponent & all, MyComponent * comp, double factor, 
        vector<Interval> domain) {
    int old = logger.reset();
    logger.disable();
    logger.log("sw set <");
    logger.inc();
    TaylorModelVec orig = comp->initSet;
    //comp->initSet = shrinkWrapTMV(comp->initSet, factor);
    
    logger.logTMV("input", all.swInput);
    TaylorModelVec wrapped = shrinkWrapTMV(all.swInput, factor);
    logger.log("");
    logger.logTMV("wrapped", wrapped);
    
    logger.listVi("compVars", all.compVars);
    
    TaylorModelVec newCompSet;
    
    for(int i = 0; i < comp->compVars.size(); i++) {
      //variable in component
      int varIndex = comp->compVars.at(i);
      //index of that variable in call component
      int indexInAll = find(all.compVars.begin(), all.compVars.end(), varIndex) - 
          all.compVars.begin();
      logger.log(sbuilder() << "index in all: " << indexInAll);
      logger.log(sbuilder() << "init var: " << varIndex);
      //using compVars, cause of the assumption that initial conditions 
      logger.listVi("tm paras", comp->allTMParams);
      logger.logTM("c1", wrapped.tms.at(indexInAll));
      TaylorModel tm = wrapped.tms.at(indexInAll).transform(comp->allTMParams);
      logger.logTM("c2", tm);
      logger.log(sbuilder() << "tm paramCount1: " << wrapped.tms[0].getParamCount());
      logger.log(sbuilder() << "tm paramCount2: " << tm.getParamCount());
      
      newCompSet.tms.push_back(tm);
    }
    
    
    //comp->initSet is the old set, newCompSet is the shrink wrapped one
    if(true) { // TODO add compiler flag
      if(comp->initSet.compare(newCompSet, domain) == false) {
	      //throw std::runtime_error("something decreased after shrink wrapping");
	    }
    }
    
    
    logger.logTMV("old", comp->initSet);
    comp->initSet = newCompSet;
    logger.logTMV("new", newCompSet);
    
    logger.dec();
    logger.log("sw set >");
    logger.restore(old);
  }
  
  void my_handler (const char * reason, const char * file, int line, 
      int gsl_errno) {
    throw IntegrationException(sbuilder() << reason);
  }
  
  double shrinkWrap(MyComponent & component, vector<Interval> domain, 
      vector<Interval> step_end_exp_table) {
    int old = logger.reset();
    logger.disable();
    logger.log("sw <");
    logger.inc();
    double s = -1, t = -1, d = -1, q;
    
    //TaylorModelVec last = component.pipes.at(component.pipes.size() - 1);
    //TaylorModelVec next;
    //last.evaluate_t(next, step_end_exp_table);
    
    TaylorModelVec next = component.swInput;
    
    //logger.logTMV("next", next);
    TaylorModelVec orig(next);
    logger.logTMV("orig", orig);
    
    int varSize = orig.tms.size();
    int paramCount = next.tms.at(0).getParamCount();
    
    
    //remove the constant part
    next.rmConstant(); 
    
    //logger.logTMV("noconst", next);
    
    
    //logger.log(sbuilder() << d);
    
    vector<int> vars;
    for(int i = 0; i < varSize; i++)
      vars.push_back(i);
    Matrix m(varSize);
    
    //set the elements of matrix m as the linear coefficients
    {
      for(int row = 0; row < varSize; row++) {
        Interval linearPart[varSize];
        next.tms.at(row).getLinearPart(linearPart, vars);
        for(int col = 0; col < varSize; col++) {
          m.set(linearPart[col].sup(), row, col); //TODO taking sup is not correct
        }
      }
    }
    logger.log("matrix M");
    logger.logMatrix(m);
    
    
    //find the inverse of m
    Matrix inv(varSize);
    
    gsl_error_handler_t *old_handler = gsl_set_error_handler (&my_handler);
    try{
      m.inverse(inv);
    }catch(IntegrationException& e) {
      logger.reset();
      logger.logMatrix(m);
      throw e;
    }
    gsl_set_error_handler (old_handler);
    
    logger.log("M^-1");
    logger.logMatrix(inv);
    
    
    //find the nonlinear part of m^-1 * TM
    vector<TaylorModel> linUnitTMs;
    for(int row = 0; row < varSize; row++) {
      //logger.log(sbuilder() << "row:" << row);
      
      //^m-1 * TM part
      TaylorModel result;
      
      for(int col = 0; col < varSize; col++) {
        //logger.log(sbuilder() << "(" << row << "," << col <<")");
        //logger.logd(m.get(row, col));
        TaylorModel temp;
        next.tms.at(col).mul(temp, inv.get(row,col));
        result.add_assign(temp);
      }
      logger.logTM("m^-1 * TM", result);
      
      linUnitTMs.push_back(result);
    }
    logger.logTMV("linear unit TMs", linUnitTMs);
    
    
    
    vector<TaylorModel> nlTMs;
    for(int row = 0; row < varSize; row++) {
      TaylorModel result = linUnitTMs.at(row);
      //substract the linear part (by adding -1*linear part)
      vector<Interval> coefs;
      coefs.push_back(0); //time
      for(int i = 0; i < paramCount - 1; i++) { //-1 because time is counted too
        //logger.log(i);
        if(i == row) // TODO assuming that varIndex == paramIndex
          coefs.push_back(-1);
        else
          coefs.push_back(0);
      }
      logger.logVI("coefs", coefs);
      TaylorModel tm(Polynomial(coefs), Interval(0));
      logger.logTM("-1*lin", tm);
      result.add_assign(tm);
      //logger.logTM("res2", result);
      nlTMs.push_back(result);
    }
    
    
    //nullify remainder (so you get only nonlinear part)
    for(int i = 0; i < nlTMs.size(); i++) {
      nlTMs.at(i).remainder = Interval(0);
    }
    
    TaylorModelVec nlPart(nlTMs);
    
    
    logger.logTMV("nlPart", nlPart);
    
    //logger.logTMV("C^-1 * ", tms);
    //finding d
    for(vector<TaylorModel>::iterator it = linUnitTMs.begin(); 
        it < linUnitTMs.end(); it++) {
      //logger.log(it->remainder.toString());
      //logger.log(sbuilder() << it->remainder.mag());
      if(it->remainder.mag() > d)
        d = it->remainder.mag();
    }
    logger.log(sbuilder() << "d: " << d);
    
    
    
    //finding s
    //logger.logTMV("nl", nlPart);
    //bound the nonlinear part
    vector<Interval> nlRange;
    
    nlPart.polyRange(nlRange, domain);
    
    //logger.logVI("nl", nlRange);
    for(vector<Interval>::iterator it = nlRange.begin();
        it < nlRange.end(); it++) {
      double mag = it->mag();
      s = max(s,mag);
    }
    logger.log(sbuilder() << "s: " << s);
    
    //logger.listVi("params", component.allTMParams);
    
    //find t
    for(int param = 1; param < paramCount; param++) {
      TaylorModelVec pDer;
      nlPart.derivative(pDer, param); //derivative wrt to parameter
      //logger.log(param);
      //logger.logTMV("pDer", pDer);
      
      //bound of the derivative
      vector<Interval> pDerRange;
      pDer.polyRange(pDerRange, domain);
      //logger.logVI("pd", pDerRange);
      for(vector<Interval>::iterator it2 = pDerRange.begin();
          it2 < pDerRange.end(); it2++) {
        double mag = it2->mag();
        //logger.log(sbuilder() << mag);
        t = max(t,mag);
      }
    }
    logger.log(sbuilder() << "t: " << t);
    
    int n = paramCount - 1;
    q = 1 + d*(1/((1-(n-1)*t)*(1-s))); //TODO not safe floating point rounding
    //logger.force(sbuilder() << "f(t): " << (1-(n-1)*t));
    //logger.force(sbuilder() << "f(s): " << (1-s));
    //logger.force(sbuilder() << (1/((1-s))));
    if((1 - s) <= 0) {
      throw IntegrationException(sbuilder() << 
          "Map is not shrinkable (s = " << s << ")");
    }
    if((1 - n*t) <= 0) {
      throw IntegrationException(sbuilder() << 
          "Map is not shrinkable (n = " << n << ", t = " << t << ")");
    }
    
    logger.log(sbuilder() << "1 - s: " << (1-s));
    logger.log(sbuilder() << "1 - nt: " << (1-n*t));
    logger.log(sbuilder() << "q:" << q);
    logger.dec();
    logger.log("sw >");
    logger.restore(old);
    return q;
  }
  
  void parametrizeVars(TaylorModelVec & tmv, vector<int> varsToIntroduce, 
      int paramCount) {
    int old = logger.reset();
    logger.disable();
    logger.log("parametrizing <");
    logger.inc();
    //int paramCount = tmv.tms[0].getParamCount();
    for(int i = 0; i < varsToIntroduce.size(); i++) {
      int var = varsToIntroduce[i];
      //logger.log(sbuilder() << "var: " << var);
      TaylorModel & tm = tmv.tms[var];
      
      Interval & rem = tm.remainder;
      
      Interval inf;
      rem.inf(inf);
      Interval sup;
      rem.sup(sup);
      
      Interval shift = inf + sup;
      shift.div_assign(2);
      
      //logger.log(sbuilder() << "shift: " << shift.toString());
      //logger.log(sbuilder() << "inf: " << inf.toString());
      //logger.log(sbuilder() << "sup: " << sup.toString());
      //logger.log(sbuilder() << "tm: " << tm.remainder.toString());
      
      
	    tm.expansion.add_assign(Monomial(shift, paramCount));
      Interval coef = (rem - shift).sup();
      
      //logger.log(sbuilder() << "coef: " << coef.toString());
      
      vector<int> degs;
      //logger.log(tm.getParamCount());
      for(int j = 0; j < paramCount; j++) {
        if(j == var + 1)
          degs.push_back(1);
        else
          degs.push_back(0);
      }
      
      tm.expansion.add_assign(Monomial(coef, degs));
      tm.remainder = Interval();
    }
    logger.dec();
    logger.log("parametrizing >");
    logger.restore(old);
  }
  
  void introduceParam(MyComponent & comp, vector<MyComponent *> comps, 
      vector<Interval> & domain) {
    int old = logger.reset();
    logger.disable();
    logger.log("introducing <");
    logger.inc();
    TaylorModelVec tmv = comp.swInput;
    
    
    //logger.logTMV("last", last);
    logger.logTMV("original", tmv);
    
    vector<int> varsToIntroduce;
    for(int i = 0; i < comps.size(); i++) {
      vector<int> & compVars = comps[i]->varsToBeIntroduced;
      varsToIntroduce.insert(varsToIntroduce.end(), 
          compVars.begin(), compVars.end());
      logger.listVi("vs: ", comps[i]->varsToBeIntroduced);
      comps[i]->varsToBeIntroduced.clear();
    }
    logger.listVi("introducing", varsToIntroduce);
    
    //no need to add more parameters now
    //int oParamCount = tmv.tms[0].getParamCount();
    
    //TaylorModelVec padded = tmv.addNParams(varsToIntroduce.size());
    TaylorModelVec padded = TaylorModelVec(tmv);
    int paramCount = padded.tms[0].getParamCount();
    
    
    if(varsToIntroduce.size() != 0) {
      parametrizeVars(padded, varsToIntroduce, paramCount);
    }
    
    logger.logTMV("introduced", padded);
    comp.swInput = padded;
    
    
    
    if(true) { // TODO add compiler flag
      if(tmv.compare(comp.swInput, domain) == false) {
	      throw std::runtime_error("something decreased after introducing param");
	    }
    }
    logger.dec();
    logger.log("introducing >");
    logger.restore(old);
  }
  
  double applyShrinkWrapping(MyComponent & all, vector<Interval> domain, 
      vector<Interval> step_end_exp_table, vector<MyComponent *> comps,
      OutputWriter & writer) {
    int old = logger.reset();
    logger.disable();
    logger.log("applying sw <");
    logger.inc();
    
    clock_t start = clock();
    
    all.remapLastFlowpipe();
    
    TaylorModelVec last = all.pipes.at(all.pipes.size() - 1);
    TaylorModelVec tmv;
    last.evaluate_t(tmv, step_end_exp_table);
    all.swInput = tmv;
    
    logger.logTMV("swInput", all.swInput);
    //introduces a parameter to swInput
    
    
    introduceParam(all, comps, domain);
    logger.log(all.swInput.tms[0].getParamCount());
    /*
    logger.log(all.compMappers.size());
    logger.listVi("0", all.compMappers.at(0));
    logger.listVi("1", all.compMappers.at(1));
    logger.listVi("2", all.compMappers.at(2));
    logger.listVi("3", all.compMappers.at(3));
    logger.listVi("4", all.compMappers.at(4));
    */
    logger.logTMV("swI aft", all.swInput);
    
    
    //logger.logTMV("all", all.pipes.at(all.pipes.size() - 1));
    //vector<TaylorModelVec> p2 = all.dependencies.at(0)->pComp->pipes;
    //logger.logTMV("p2", p2.at(p2.size() - 1));
    
    //calculate the shrink wrapping factor for swInput
    double swQ = shrinkWrap(all, domain, step_end_exp_table);
    
    logger.log(sbuilder() << "swQ: " << swQ);
    
    //use the computed shrink wrapping factor to modify the initial sets
    for(vector<MyComponent *>::iterator it = comps.begin(); 
        it < comps.end(); it++) {
      //sets comp.initSet from q and all.swInput
      shrinkWrapSet(all, *it, swQ, domain);
    }
    
    clock_t end = clock();
    double singleSW = double(end - start) / CLOCKS_PER_SEC;
    writer.swTime += singleSW;
    //logger.force(sbuilder() << "single: " << singleSW);
    //logger.force(sbuilder() << (end - start));
    
    logger.dec();
    logger.log("applying sw >");
    logger.restore(old);
  }
  
  
  
  void bar12() {
    
    
    vector<int> pow1;
    pow1.push_back(0);
    pow1.push_back(0);
    pow1.push_back(0);
    Monomial m1(Interval(2), pow1);
    
    vector<int> pow2;
    pow2.push_back(0);
    pow2.push_back(1);
    pow2.push_back(0);
    Monomial m2(Interval(4), pow2);
    
    vector<int> pow3;
    pow3.push_back(0);
    pow3.push_back(2);
    pow3.push_back(0);
    Monomial m3(Interval(0.5), pow3);
    
    
    
    vector<int> pow4;
    pow4.push_back(0);
    pow4.push_back(0);
    pow4.push_back(0);
    Monomial m4(Interval(1), pow4);
    
    vector<int> pow5;
    pow5.push_back(0);
    pow5.push_back(0);
    pow5.push_back(1);
    Monomial m5(Interval(3), pow5);
    
    vector<int> pow6;
    pow6.push_back(0);
    pow6.push_back(1);
    pow6.push_back(1);
    Monomial m6(Interval(1), pow6);
    
    list<Monomial> l1;
    l1.push_back(m1);
    l1.push_back(m2);
    l1.push_back(m3);
    Polynomial p1(l1);
    
    
    list<Monomial> l2;
    l2.push_back(m4);
    l2.push_back(m5);
    l2.push_back(m6);
    Polynomial p2(l2);
    
    TaylorModel t1(p1, Interval(-0.2, 0.2));
    TaylorModel t2(p2, Interval(-0.1, 0.1));
    logger.logTM("t1", t1);
    logger.logTM("t2", t2);
    
    vector<TaylorModel> tms;
    tms.push_back(t1);
    tms.push_back(t2);
    
    TaylorModelVec tmv(tms);
    
    logger.logTMV("tmv", tmv);
    
    vector<Interval> domain;
    domain.push_back(Interval(0));
    domain.push_back(Interval(-1,1));
    domain.push_back(Interval(-1,1));
    
    vector<Interval> step_exp_table, step_end_exp_table;
  	construct_step_exp_table(step_exp_table, step_end_exp_table, 0.1, 2*2); //2*order
    
    MyComponent comp;
    comp.pipes.push_back(tmv);
    comp.solveIndexes.push_back(0);
    comp.solveIndexes.push_back(1);
    comp.allTMParams.push_back(0);
    comp.allTMParams.push_back(1);
    
    
    shrinkWrap(comp, domain, step_end_exp_table);
  }
  
  void singleStepPrepareIntegrate(MyComponent & component, 
      const OutputWriter writer, int order, double step, double time, 
      vector<Interval> step_end_exp_table, TaylorModelVec tmv, 
      const vector<HornerForm> & ode, vector<Interval> domain) {
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
      smallComp::singleStepPrepareIntegrate(*pComp, writer, order, step, time, 
          step_end_exp_table, tmv, ode, domain);
    }
    //logger.log("solving");
    logger.listVi("component vars: ", component.compVars);
    
    //component.prepare(tmv, ode, domain);
    
    //remaps previous components flowpipes

    component.remapLastFlowpipe();
    
    //component.log();
    
    singleStepIntegrate(component, writer, order, step, time,
        step_end_exp_table);
    component.isSolved = true;
    
    logger.dec();
    logger.log("singleStepPrepareIntegrate >");
    logger.restore(old);
  }
  
  void integrateComponentWrapper(MyComponent & component, 
      const OutputWriter writer, int order, double step, double time, 
      vector<Interval> step_end_exp_table, TaylorModelVec tmv, const vector<HornerForm> & ode, 
      vector<Interval> domain) {
      
    //if component has already been solved return
    if(component.isSolved) {
      //logger.log("solved already");
      //logger.listVi("component vars: ", component.compVars);
      return;
    }
    for(vector<CompDependency *>::iterator it = component.dependencies.begin(); 
        it < component.dependencies.end(); it++) {
      //logger.log(sbuilder() << "link: " << (*it)->linkVar);
      MyComponent *pComp = (*it)->pComp;
      //solve all dependencies
      smallComp::integrateComponentWrapper(*pComp, writer, order, step, time, 
          step_end_exp_table, tmv, ode, domain);
    }
    //logger.log("solving");
    logger.listVi("component vars: ", component.compVars);
    
    component.prepare(tmv, ode, domain);
    component.log();
    integrateComponent2(component, writer, order, step, time, step_end_exp_table);
    component.isSolved = true;
    
    
    //bar12();
    //shrinkWrap(component, domain, step_end_exp_table);
  }
  
  
	void findDecreasingRemainder(TaylorModelVec & p, vector<Interval> & pPolyRange, 
	    vector<RangeTree *> & trees, MyComponent & comp, 
	    vector<Interval> & step_exp_table, int order, 
	    const Interval & cutoff_threshold, vector<Interval> & cutoffInt) {
	  int old = logger.reset();
	  logger.log("findDecreasingRemainder <");
	  logger.inc();
	  
	  //number of taylor model parameters
    int paramCount = comp.initSet.tms[0].getParamCount();
    //number of system variables (in the component)
    int varCount = comp.initSet.tms.size();
	  
	  //need to find remainders more efficiently
	  p.polyRangeNormal(pPolyRange, step_exp_table);
	  
	  //initial guess for the remainder
	  vector<Interval> guess;
	  for(int i = 0; i < varCount; i++)
	    guess.push_back(Interval(-1e-2,1e-2)); //TODO use remainder estimation maybe
	  
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
	    logger.log(sbuilder() << "j: " << j);
	    
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
		      logger.log("not subset");
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
	  logger.log("findDecreasingRemainder >");
	  logger.restore(old);
	}
	
	
	void refineRemainder(TaylorModelVec & p, vector<Interval> & pPolyRange, 
	    vector<RangeTree *> & trees, MyComponent & comp, 
	    vector<Interval> & step_exp_table, int order, 
	    const Interval & cutoff_threshold, vector<Interval> & cutoffInt) {
	  int old = logger.reset();
	  logger.log("refineRemainder <");
	  logger.inc();
	  
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
	    for(int i = 0; i < varCount; i++)
	      if(p.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO) {
          //logger.log("redoing");
          stop = false;
          break;
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
	  logger.log("refineRemainder >");
	  logger.restore(old);
	}
	
	void precond(TaylorModelVec & tmv, vector<Interval> & step_end_exp_table) {
	  logger.log("precond");
	  
	  logger.logTMV("tmv", tmv);
	  
	  
	  
	  //number of taylor model parameters
    int paramCount = tmv.tms[0].getParamCount();
    //number of system variables (in the component)
    int varCount = tmv.tms.size();
	  
	  logger.log(paramCount);
	  logger.log(varCount);
	  Matrix A(varCount,varCount), AT(varCount, varCount);
		preconditionQR(A, tmv, varCount, paramCount);
		
	  logger.logMatrix(A);
	  A.transpose(AT);
	  logger.logMatrix(AT);
	  
	  vector<Interval> center;
	  tmv.constant(center);
	  logger.logVI("center", center);
	  
	  TaylorModelVec c(tmv), lin, nl, rem;
	  tmv.getParts(c, lin, nl, rem);
	  
	  logger.logTMV("c", c);
	  logger.logTMV("lin", lin);
	  logger.logTMV("nl", nl);
	  logger.logTMV("rem", rem);
	  
	  nl.linearTrans_assign(AT);
	  rem.linearTrans_assign(AT);
	  lin.linearTrans_assign(AT);
	  
	  TaylorModelVec rightPart(lin);
	  rightPart.add_assign(nl);
	  rightPart.add_assign(rem);
	  
	  
	  vector<Interval> rightBound;
		rightPart.intEvalNormal(rightBound, step_end_exp_table);
	  logger.logTMV("rightPart", rightPart);
	  logger.logVI("rightBound", rightBound);
	  
	  //TODO shift this to be around 0?
	  
	  Matrix S(varCount, varCount);
  	Matrix invS(varCount, varCount);
	  for(int i = 0; i < varCount; i++) {
	    Interval intMag;
	    rightBound[i].mag(intMag);
	    double dMag = intMag.sup();
	    if(intMag.subseteq(Interval()) == false) {
			  S.set(dMag, i, i);
			  invS.set(1/dMag, i, i);
		  }
		  else {
			  S.set(0, i, i);
			  invS.set(0, i, i); //doesn't matter if 0 or 1
			}
	  }
	  
	  logger.logMatrix(S);
	  logger.logMatrix(invS);
	  
	  rightPart.linearTrans_assign(invS);
		rightPart.intEvalNormal(rightBound, step_end_exp_table);
	  logger.logTMV("rightPart", rightPart);
	  logger.logVI("rightBound", rightBound);
	  
	  
	  
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
	  
	  findDecreasingRemainder(p, pPolyRange, trees, comp, step_exp_table, order, 
	      cutoff_threshold, cutoffInt);
	  
	  refineRemainder(p, pPolyRange, trees, comp, step_exp_table, order, 
	      cutoff_threshold, cutoffInt);
	  
    logger.logTMV("p", p);
    
    //TaylorModelVec temp;
    //p.evaluate_t(temp, step_end_exp_table);
    
    parseSetting.addVar("t"); parseSetting.addVar("a"); parseSetting.addVar("b");
    TaylorModelVec temp = parseTMV(
        "my models{0.904667 + 0.0505*a + 0.005*b + [-5.09307e-5,7.86167e-5],-0.909333+0.0095*a+0.0505*b + 0.00025*a^2 + [-1.75707e-4,1.60933e-4]}");
    
    precond(temp, step_end_exp_table);
	  
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


void SmallCompSystem::my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const {
  logger.log("sc reach <");
  logger.inc();
  logger.log(sbuilder() << "# of components: " <<components.size());
  
  vector<MyComponent *> comps = createComponents(components, hfOde);
  
  
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
  
  vector<TaylorModelVec> pipes;
  
  if(precondition == SHRINK_WRAPPING) {
    for(int i = 0; i < comps.size(); i++) {
      comps.at(i)->retainEmptyParams = true;
    }
  }
  
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(currentTMV, hfOde, domain);
  }
  MyComponent all = getSystemComponent(comps, currentTMV, hfOde, domain);
  //all.log();
  
  
  MyComponent c = *comps.at(0);
  
  smallComp::foo(*comps.at(0), order, cutoff_threshold, step_exp_table, 
      step_end_exp_table);
  
  //single step integration
  //for(double t = 0; t < time; t+= step) {
  
  clock_t integrClock = clock();
  double t;
  for(t = 0; t < time; t+= step) {
    logger.log(sbuilder() << "t: " << t);
    
    
    int i = 0;
    try{
      for(vector<MyComponent *>::iterator it = comps.begin(); 
          it < comps.end(); it++) {
        //logger.logTMV("init", (*it)->initSet);
        smallComp::singleStepPrepareIntegrate(**it, writer, order, step, 
            step, step_end_exp_table, currentTMV, hfOde, domain);
        //logger.logTMV("0", (*it)->pipes.at(0));
        i++;
      }
      //prepare next initSet
      for(vector<MyComponent *>::iterator it = comps.begin(); 
          it < comps.end(); it++) {
        (*it)->isSolved = false;
        int lastIndex = (*it)->pipes.size()-1;
        
        //set the next initSet
        (*it)->pipes.at(lastIndex).evaluate_t((*it)->initSet,
            step_end_exp_table);
            
        //logger.logTMV("last", (*it)->pipes.at(lastIndex));
        //logger.listVi("vars", (*it)->varIndexes);
        //logger.logTMV("evaluated", (*it)->initSet);
      }
      
      int old2 = logger.reset();
      if(precondition == SHRINK_WRAPPING) {
        if(swChecker->checkApplicability(comps, estimation)) {
          logger.force("wrapping");
          smallComp::applyShrinkWrapping(all, domain, step_end_exp_table, comps,
              writer);
        }
      }
      logger.restore(old2);
      /*
      vector<Interval> polyRange;
      vector<Interval> range;
      comps.at(0)->initSet.polyRange(polyRange, domain);
      Interval iEv;
      comps.at(0)->initSet.tms.at(0).intEval(iEv, domain);
      
      logger.log(sbuilder() << "pRa: " << polyRange.at(0).toString());
      logger.log(sbuilder() << "iEv: " << iEv.toString());
      logger.log(sbuilder() << "rem:" << 
          comps.at(0)->initSet.tms.at(0).remainder.toString());  
      //*/
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
  writer.info.push_back(sbuilder() << "shrink wraps: " << swChecker->getCount());
  
  writer.addComponents(comps, domain);
  writer.writeCSV();
  writer.writeInfo();
  
  writer.finish();

  logger.dec();
  logger.log("sc reach >");
}
