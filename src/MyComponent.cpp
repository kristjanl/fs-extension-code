#include "MyComponent.h"

using namespace std;



MyComponent::MyComponent(vector<int> vs, vector<int> tps):varIndexes(vs),tpIndexes(tps) {
  isSolved = false;
  isPrepared = false;
  retainEmptyParams = false;
}
MyComponent::MyComponent() {
  isSolved = false;
  isPrepared = false;
  retainEmptyParams = false;
}

void MyComponent::log() {
  logger.log("component <");
  logger.inc();
  stringstream varSs;
  logger.log(varIndexes.size());
  for(int i = 0; i < varIndexes.size(); i++) {
    varSs << varIndexes.at(i) << ", ";
  }
  logger.log(sbuilder() << "varIndexes: " << varSs.str());
  logger.log(sbuilder() << "odes size: " << odes.size());
  {
    int odeIndex = 0;
    for(vector<HornerForm>::iterator it = odes.begin(); it < odes.end(); it++) {
      logger.log(sbuilder() << odeIndex++ << ": " << (*it).toString());
    }
  }
  logger.logVI("compDomain", dom);  
  logger.logTMV("compInit", initSet);
  
  for(int i = 0; i < initSet.tms.size(); i++) {
    logger.log(sbuilder() << "tm param size(" << i << "): "
        << initSet.tms.at(i).getParamCount());
  }
  logger.log(sbuilder() << "init set size: " << initSet.tms.size());
  logger.listVi("to be solved variables", varIndexes);
  logger.listVi("all variables", compVars);
  logger.listVi("solve indexes", solveIndexes);
  logger.listVi("new params", tpIndexes);
  logger.listVi("all params", allTMParams);
  logger.log(sbuilder() << "mappers size: " << compMappers.size());
  for(int i = 0; i < compMappers.size(); i++) {
    logger.listVi(
        (sbuilder() << "m" << i << 
        " (link is " << dependencies.at(i)->linkVar <<")"),
        compMappers.at(i));
  }
  logger.dec();
  logger.log("component >");
  
}

/*
creates mappers for previous components
mappers are of the form [i1,i2,...,in]
where is are either indexes to parameters in the previous components
or -2 to indicate placeholders
[-2,-2,0,-2,1] maps a1*a2 to a3*a5
*/
void MyComponent::prepareMappers() {
  int old = logger.reset();
  logger.disable();
  logger.log("preparing mappers <");
  logger.inc();

  vector< vector<int> > mappers;

  for(int i = 0; i < dependencies.size(); i++) {
    //logger.log(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    //logger.log(sbuilder() << "link: " << dependencies.at(i)->linkVar);
    //previous component parameters
    vector<int> prevCompPm = pComp->allTMParams;
    //logger.logVi("pm", prevCompPm);
    vector<int> mapper;
    //iterate all the variables in all the components
    for(int j = 0; j < allTMParams.size(); j++) {
      int var = allTMParams.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = 
          find(prevCompPm.begin(), prevCompPm.end(), var);
      
      if(it != prevCompPm.end()) {
        //find position if it was
        int pos = distance(prevCompPm.begin(), it);
        //map old position to new
        mapper.push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper.push_back(PADDING_VARIABLE);
      }
    }
    dependencies.at(i)->mapper = mapper;
    //logger.logVi("m", mapper);
    mappers.push_back(mapper);
  }
  compMappers = mappers;
  //logger.logVi("mapper0", mappers.at(0));
  //logger.logVi("mapper1", mappers.at(1));
  //return *mappers;
  //return mappers;
  logger.dec();
  logger.log("preparing mappers <");
  logger.restore(old);
}


vector< vector<int> > MyComponent::previousMappers2() {
  vector<int> allParams;
  
  //include current component parameters
  for(vector<int>::iterator it = tpIndexes.begin(); 
        it < tpIndexes.end(); it++) {
    allParams.push_back(*it);
  }
   
  logger.listVi("allParams (current)", allParams);
  //gather all indexes of all previous components
  for(int i = 0; i < dependencies.size(); i++) {
    logger.log(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    
    logger.logVi("vi", pComp->varIndexes);
    logger.logVi("pi", pComp->tpIndexes);
    allParams.insert(allParams.end(), 
        pComp->allTMParams.begin(), pComp->allTMParams.end());
  }
  logger.logVi("ai", allParams);
  //sort all indexes
  sort(allParams.begin(), allParams.end());
  //remove duplicate indexes
  allParams.erase( unique( allParams.begin(), allParams.end() ), allParams.end() );
  
  allTMParams = allParams;
  logger.logVi("allParams (with dependencies)", allParams);
  
  
  logger.logVi("ai", allParams);
  vector< vector<int> > mappers;
  
  for(int i = 0; i < dependencies.size(); i++) {
    logger.log(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    logger.log(sbuilder() << "link: " << dependencies.at(i)->linkVar);
    //previous component parameters
    vector<int> prevCompPm = pComp->allTMParams;
    logger.logVi("pm", prevCompPm);
    vector<int> mapper;
    //iterate all the variables in all the components
    for(int j = 0; j < allParams.size(); j++) {
      int var = allParams.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = find(prevCompPm.begin(), prevCompPm.end(), var);
      if(it != prevCompPm.end()) {
        //find position if it was
        int pos = distance(prevCompPm.begin(), it);
        //map old position to new
        mapper.push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper.push_back(PADDING_VARIABLE);
      }
    }
    dependencies.at(i)->mapper = mapper;
    logger.logVi("m", mapper);
    mappers.push_back(mapper);
  }
  //logger.logVi("mapper0", mappers.at(0));
  //logger.logVi("mapper1", mappers.at(1));
  
  //return *mappers;
  return mappers;
}


vector< vector<int> > MyComponent::previousMappers() {
  vector<int> allIndexes;
  
  //gather all indexes of all components
  for(int i = 0; i < previous.size(); i++) {
    allIndexes.insert(allIndexes.end(), 
        previous.at(i).varIndexes.begin(), previous.at(i).varIndexes.end());
  }
  //sort all indexes
  sort(allIndexes.begin(), allIndexes.end());
  //remove duplicate indexes
  allIndexes.erase( unique( allIndexes.begin(), allIndexes.end() ), allIndexes.end() );
  
  //logger.logVi("allIndexes", allIndexes);
  
  vector< vector<int> > *mappers = new vector< vector<int> >();
  
  for(int i = 0; i < previous.size(); i++) {
    //previous component variables
    vector<int> prevCompInd = previous.at(i).varIndexes;
    vector<int> *mapper = new vector<int>();
    //iterate all the variables in all the components
    for(int j = 0; j < allIndexes.size(); j++) {
      int var = allIndexes.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = find(prevCompInd.begin(), prevCompInd.end(), var);
      if(it != prevCompInd.end()) {
        //find position if it was
        int pos = distance(prevCompInd.begin(), it);
        //map old position to new
        mapper->push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper->push_back(PADDING_VARIABLE);
      }
    }
    mappers->push_back(*mapper);
  }
  //logger.logVi("mapper0", mappers->at(0));
  //logger.logVi("mapper1", mappers->at(1));
  
  return *mappers;
}


CompDependency::CompDependency(int link, MyComponent *pComp)
    :linkVar(link), pComp(pComp) {
  //logger.log("creating");
}


void MyComponent::addDependency(int linkVar, MyComponent *pComp) {
  dependencies.push_back(new CompDependency(linkVar, pComp));
}

void MyComponent::addVar(int var) {
  //logger.log("adding");
  varIndexes.push_back(var);
  compVars.push_back(var);
}

void addEmptyTM(vector<TaylorModelVec> & pipes, int dim) {  
  TaylorModelVec temp;
  TaylorModel t;
  for(int i = 0; i < dim; i++) {
    temp.tms.push_back(t);
  }
  pipes.push_back(temp);
}

void MyComponent::remapIVP(TaylorModelVec tmv, const vector<HornerForm> & ode, 
      vector<Interval> domain) {
  int old = logger.reset();
  logger.disable();
  logger.log("remapping IV <");
  logger.inc();

  
  logger.listVi("all params", allTMParams);
  logger.listVi("all variables", compVars);
  logger.listVi("link variables", linkVars);
  
  vector<HornerForm> compOdes;
  for(int i = 0; i < compVars.size(); i++) {
    int varIndex = compVars.at(i);
    logger.log(sbuilder() << "v: " << varIndex);
    logger.inc();
    HornerForm hf;
    if(linkVars.end() != find(linkVars.begin(), linkVars.end(), varIndex)) {
      // we don't care about the ode if the variable is a link
    } else {
      try{
        logger.listVi("compVars", compVars);
        logger.log(sbuilder() << "ode: " << ode.at(varIndex).toString());
        hf = (ode.at(varIndex)).transform(compVars);
      }
      catch (exception& e) {
        logger.reset();
        logger.log(e.what());
        logger.logVi("indexes", compVars);
        throw e;
      }
    }
    compOdes.push_back(hf);
    logger.log(hf.toString());
    logger.dec();
  }
  
  vector<Interval> compDomain;
  compDomain.push_back(domain.at(0)); //time domain
  for(int i = 0; i < allTMParams.size(); i++) {
    int paramIndex = allTMParams.at(i);
    compDomain.push_back(domain.at(paramIndex + 1)); //+1 for time
  }
  //int old2 = logger.reset();
  TaylorModelVec *compInit = new TaylorModelVec();
  for(int i = 0; i < compVars.size(); i++) {
    int varIndex = compVars.at(i);
    logger.log(sbuilder() << "init var: " << varIndex);
    //using compVars, cause of the assumption that initial conditions 
    logger.listVi("tm paras", allTMParams);
    logger.logTM("c1", tmv.tms.at(varIndex));
    TaylorModel tm = tmv.tms.at(varIndex).transform(allTMParams);
    logger.logTM("c2", tm);
    logger.log(sbuilder() << "tm paramCount: " << tm.getParamCount());
    
    (*compInit).tms.push_back(tm);
  }
  //logger.restore(old2);
  
  odes = compOdes;
  dom = compDomain;
  initSet = *compInit;
  
  logger.dec();
  logger.log("remapping IV >");
  logger.restore(old);
}
  

  
void MyComponent::prepareComponent(TaylorModelVec init, 
    const vector<HornerForm> & ode, vector<Interval> domain, 
    bool preconditioned) {
  //return if variables have already been prepared
  if(isPrepared) {
    return;
  }
  int old = logger.reset();
  logger.disable();
  logger.log("preparing component <");
  logger.inc();
  logger.listVi("vars", varIndexes);
  logger.log(sbuilder() << "dependencies: " << dependencies.size());
  //prepare variables for all the dependencies
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
        it < dependencies.end(); it++) {
    //logger.log(sbuilder() << "link: " << (*it)->linkVar);
    MyComponent *pComp = (*it)->pComp;
    pComp->prepareComponent(init, ode, domain, preconditioned);
  }

  
  prepareVariables(init, preconditioned);
  prepareMappers();
  remapIVP(init, ode, domain);
  
  isPrepared = true;
  
  logger.dec();
  logger.log("preparing component >");
  logger.restore(old);
}
  
void MyComponent::prepareVariables(TaylorModelVec init, bool preconditioned) {
  int old = logger.reset();
  logger.disable();
  logger.log("preparing variables <");
  logger.inc();

  //add all linking variables
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++) {
    compVars.push_back((*it)->linkVar);
    //TODO can remove and replace with unlinked
    linkVars.push_back((*it)->linkVar);
  }
  
  //find the indexes (relative in the component) of variables which 
  //need to be solved
  {
    vector<int> *v1 = &(varIndexes);
    vector<int> *v2 = &(compVars);
    for(vector<int>::iterator it = v1->begin(); it < v1->end(); it++) {
      //logger.log(*it);
      int pos = find(v2->begin(), v2->end(), *it) - v2->begin();
      //logger.log(sbuilder() << "pos: " << pos);
      solveIndexes.push_back(pos);
    }
  }
  
  //gather nonzero taylor model parameters from initial conditions
  for(vector<int>::iterator it = varIndexes.begin(); 
      it < varIndexes.end(); it++) {
    logger.log(sbuilder() << "var: " << *it);
    //logger.logTM("tm", init.tms.at(*it));
    vector<int> tmParams = init.tms.at(*it).getParams();
    logger.listVi("params", tmParams);
    
    
    //retain the parameter if started with point initial conditions
    if(tmParams.size() == 0 && retainEmptyParams) {
      varsToBeIntroduced.push_back(*it);
      tmParams.push_back(*it);
    }
    tpIndexes.insert(tpIndexes.end(), tmParams.begin(), tmParams.end());
  }
  sort(tpIndexes.begin(), tpIndexes.end());
  tpIndexes.erase( 
      unique(tpIndexes.begin(), tpIndexes.end()), tpIndexes.end() );
  
  
  vector<int> allParams;
      
  //include current component parameters
  for(vector<int>::iterator it = tpIndexes.begin(); 
        it < tpIndexes.end(); it++) {
    allParams.push_back(*it);
  }
  
  //logger.listVi("allParams (current)", allParams);
  //gather all indexes of all previous components
  for(int i = 0; i < dependencies.size(); i++) {
    //logger.log(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    
    //logger.logVi("vi", pComp->varIndexes);
    //logger.logVi("pi", pComp->tpIndexes);
    allParams.insert(allParams.end(), 
        pComp->allTMParams.begin(), pComp->allTMParams.end());
  }
  
  //sort all params
  sort(allParams.begin(), allParams.end());
  //remove duplicate params
  allParams.erase( 
      unique( allParams.begin(), allParams.end() ), allParams.end() );
  
  allTMParams = allParams;
  //logger.logVi("allParams (with dependencies)", allParams);
  
  logger.dec();
  logger.log("preparing variables >");
  logger.restore(old);
}  

void MyComponent::remapFlowpipes() {
  vector< vector<int> > mappers = previousMappers2();
  logger.log(sbuilder() << "mappers size: " << mappers.size());
  

  //mapping of dependent component flowpipes to current component
  int depIndex = 0;
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++, depIndex++) {
    logger.log(sbuilder() << "depIndex: " << depIndex);
    int link = (*it)->linkVar;
    logger.log(sbuilder() << "link: " << link);
    MyComponent *pComp = (*it)->pComp;
    logger.listVi("dts", pComp->tpIndexes);
    logger.listVi("mapper", (*it)->mapper);
    logger.log(pComp->pipes.size());
    logger.log("here");
    logger.logTMV("p", pComp->pipes.at(0));
    
    
    logger.log("here");
    logger.listVi("dcv", pComp->compVars);
    
    
    //position of linking variable in precomputed component
    int dLinkPos = find(pComp->compVars.begin(), pComp->compVars.end(), link) - 
        pComp->compVars.begin();
    //position of linking variable in current component
    int linkPos = find(compVars.begin(), compVars.end(), link) - 
        compVars.begin();
    if(dLinkPos >= pComp->compVars.size()) {
      logger.log(sbuilder() << "dLinkPos: " << dLinkPos);
      logger.listVi("compVars", compVars);
      throw std::invalid_argument("depend link wasn't in compVar");
    }
    if(linkPos >= compVars.size()) {
      throw std::invalid_argument("link wasn't in compVars");
    }
    logger.log(sbuilder() << "dLinkPos: " << dLinkPos);
    logger.log(sbuilder() << "linkPos: " << linkPos);
    
    
    logger.listVi("vars", varIndexes);
    logger.listVi("dvars", pComp->varIndexes); 
      
    //actual mapping of flowpipes
    int pipeIndex = 0;
    for(vector<TaylorModelVec>::iterator pipeIt = pComp->pipes.begin();
        pipeIt < pComp->pipes.end(); pipeIt++, pipeIndex++) {
      logger.log(sbuilder() << "pipeIndex: " << pipeIndex);
      logger.logTM("tm", (*pipeIt).tms.at(dLinkPos));
      logger.listVi("ma", mappers.at(depIndex));
      
      
      TaylorModel transformedLink = 
          (*pipeIt).tms.at(dLinkPos).transform(mappers.at(depIndex));
      pipes.at(pipeIndex).tms.at(linkPos) = transformedLink;
      //logger.logTMV("pipe", pipes.at(pipeIndex));
    }
    
  }
}

void MyComponent::remapLastFlowpipe() {
  int old = logger.reset();
  logger.disable();
  logger.log("remapping last <");
  logger.inc();
  int varSize = varIndexes.size() + dependencies.size();
  
  addEmptyTM(pipes, varSize);
  
  logger.log(dependencies.size());
  
  if(dependencies.size() == 0) {
    logger.restore(old);
    logger.dec();
    return;
  }
  
  //logger.disable();
  
  //mapping of dependent component flowpipes to current component
  int depIndex = 0;
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++, depIndex++) {
    logger.log(sbuilder() << "depIndex: " << depIndex);
    int link = (*it)->linkVar;
    logger.log(sbuilder() << "link: " << link);
    MyComponent *pComp = (*it)->pComp;
    logger.listVi("dts", pComp->tpIndexes);
    logger.listVi("mapper", (*it)->mapper);
    logger.log(pComp->pipes.size());
    logger.logTMV("p", pComp->pipes.at(0));
    
    
    //logger.log("here");
    //logger.listVi("dcv", pComp->compVars);
    
    
    //position of linking variable in precomputed component
    int dLinkPos = find(pComp->compVars.begin(), pComp->compVars.end(), link) - 
        pComp->compVars.begin();
    //position of linking variable in current component
    int linkPos = find(compVars.begin(), compVars.end(), link) - 
        compVars.begin();
        
        
    if(dLinkPos >= pComp->compVars.size()) {
      logger.log(sbuilder() << "dLinkPos: " << dLinkPos);
      logger.listVi("compVars", compVars);
      throw std::invalid_argument("depend link wasn't in compVar");
    }
    if(linkPos >= compVars.size()) {
      throw std::invalid_argument("link wasn't in compVars");
    }
    //logger.log(sbuilder() << "dLinkPos: " << dLinkPos);
    //logger.log(sbuilder() << "linkPos: " << linkPos);
    
    
    //logger.listVi("vars", varIndexes);
    //logger.listVi("dvars", pComp->varIndexes);
    
    
    //logger.logTMV("cpipeb", pipes.at(0));
    //logger.logTMV("dpipeb", pComp->pipes.at(0));
    
    int lastPipeIndex = pipes.size() - 1;
    int lastSourcePipeIndex = pComp->pipes.size() - 1;
    logger.log(sbuilder() << "last pipe index: " << lastPipeIndex);
    logger.logTMV("last pipe", pComp->pipes.at(lastSourcePipeIndex));
    TaylorModel transformedLink = pComp->pipes.at(lastSourcePipeIndex).tms
        .at(dLinkPos).transform(compMappers.at(depIndex));
    pipes.at(lastPipeIndex).tms.at(linkPos) = transformedLink;
    logger.logTM("remapped", transformedLink);
    
    //logger.logTMV("cpipea", pipes.at(0));
    //logger.logTMV("dpipea", pComp->pipes.at(0));
  }
  
  
  //logger.enable();
  //logger.logTMV("init", initSet);
  logger.dec();
  logger.log("remapping last >");
  logger.restore(old);
}

void MyComponent::prepare(TaylorModelVec tmv, const vector<HornerForm> & ode, 
      vector<Interval> domain) {
  logger.disable();
  prepareVariables(tmv, false);
  
  //add flowpipes if the component is not initial one
  //variables in component + link vars
  int varSize = varIndexes.size() + dependencies.size();
  
  if(dependencies.size() != 0)
    for(int i = 0; i < dependencies.at(0)->pComp->pipes.size(); i++) {
      addEmptyTM(pipes, varSize);
    }
  
  remapFlowpipes();
  logger.listVi("all tm params", allTMParams);
  
  logger.listVi("to solve variables", varIndexes);
  logger.listVi("component variables", compVars);
  logger.listVi("linking variables", linkVars);
  logger.listVi("solve indexes", solveIndexes);
  logger.log(sbuilder() << "# of pipes: " << pipes.size());
  remapIVP(tmv, ode, domain);
  logger.enable();
}

vector<MyComponent *> createComponents(vector< vector<int> > compIndexes, 
    const vector<HornerForm> & ode) {
  logger.disable();
  logger.log("creating <");
  logger.inc();
  int varSize = ode.size();
  vector<MyComponent *> components;
  MyComponent *lookup[varSize];
  
  //add variables to components
  //create a lookup from variables to components
  for(vector< vector <int> >::iterator cit = compIndexes.begin(); 
      cit < compIndexes.end(); cit++) {
    logger.listVi("componenet vars", *cit);
    MyComponent *c = new MyComponent();
    for(vector<int>::iterator it = cit->begin(); it < cit->end(); it++) {
      c->addVar(*it);
      lookup[*it] = c;
    }
    components.push_back(c);
    
  }
  
  int vars[varSize];
  for(vector<MyComponent *>::iterator cit = components.begin(); 
      cit < components.end(); cit++) {
    
    //initial vars with 0s
    memset(vars, 0, varSize*sizeof(int) );
    
    //find all the variables that the component depends on
    for(vector<int>::iterator it = (*cit)->varIndexes.begin(); 
        it < (*cit)->varIndexes.end(); it++) {
      logger.log(sbuilder() << "var: " << *it);
      logger.log(ode.at(*it).toString());
      
      //populate vars with variables that exist in the ode
      bool b = ode.at(*it).getVars(vars);
    }
    
    //nullify variables that are going to be solved in component
    for(vector<int>::iterator it = (*cit)->varIndexes.begin(); 
        it < (*cit)->varIndexes.end(); it++) {
      vars[*it] = 0;
    }
    
    //add component dependency for variables that are not in the current
    //component, but who still exist in the ode.
    for(int i = 0; i < varSize; i++) {
      if(vars[i] > 0) {
        logger.log(sbuilder() << "i: " << i);
        (*cit)->addDependency(i, lookup[i]);
      }
    }
    vector<int> v(vars, vars + sizeof vars / sizeof vars[0]);
    logger.listVi("v", v);
  }
  logger.dec();
  logger.log("creating >");
  logger.enable();
  return components;
}


MyComponent getSystemComponent(vector<MyComponent *> comps,
    TaylorModelVec init, const vector<HornerForm> & ode,
    vector<Interval> domain, bool preconditioned) {
  MyComponent allVars;
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    for(vector<int>::iterator i2 = (*it)->varIndexes.begin();
        i2 < (*it)->varIndexes.end(); i2++) {
      allVars.addDependency(*i2, *it);
    }
  }
  allVars.prepareComponent(init, ode, domain, preconditioned);
  return allVars;
}


void prepareComponents(vector<MyComponent *> & comps, TaylorModelVec init, 
    const vector<HornerForm> & ode, vector<Interval> domain, 
    bool preconditioned) {
  int old = logger.reset();
  logger.disable();
  logger.log("preparing components <");
  logger.inc();
  logger.force("shouldn't be used");
  exit(0);
  
  logger.logTMV("init", init);
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(init, ode, domain, preconditioned);
    logger.logTMV("pre Init", (*it)->initSet);
  }

  logger.dec();
  logger.log("preparing components >");
  logger.restore(old);
}

int MyComponent::nextFreeParam = -1;



PrecondModel *MyComponent::lastPre() {
  return pipePairs[pipePairs.size() - 1];
}

TaylorModelVec MyComponent::lastPipe() {
  return pipes[pipes.size() - 1];
}

bool MyComponent::isSolveVar(int var) {
  vector<int>::iterator it = find(solveIndexes.begin(), solveIndexes.end(), var);
  return it != solveIndexes.end();
}

void MyComponent::serializeFlows() {
  logger.log("writing flows");
  //exit(0);
  
  logger.reset();
  logger.log("bar");
  FILE *fpDumping = fopen("temp.txt", "w");
  logger.logTM("last", lastPipe().tms[0]);
  vector<string> names;
  names.push_back("t");
  names.push_back("a1");
  names.push_back("a2");
  
  vector<string> varNames;
  varNames.push_back("x1");
  varNames.push_back("x2");
  
  
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a1");
	parseSetting.addVar("a2");
  
  
  fprintf(fpDumping, "flowpipes{\n");
  bool first = true;
  int n = 1;
  for(vector<TaylorModelVec>::iterator pipeIt = pipes.begin();
        pipeIt < pipes.end(); pipeIt++) {
    if(first == false)
      fprintf(fpDumping, ",\n");
    first = false;
    
    //logger.logTMV("pip", *pipeIt);
    pipeIt->serialize(fpDumping, names);
    //parseTMV(sbuilder() << "my models{" << n++ << " * a1 + [-0.001,0.001]}").serialize(fpDumping, varNames, names);
    //break;
  }
  fprintf(fpDumping, "}\n");

  fclose(fpDumping);
}
void MyComponent::deserializeFlows() {
  exit(0);
  logger.log("reading flows");
  logger.reset();
  logger.log("bar");
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a1");
	parseSetting.addVar("a2");
  parseFile("temp.txt");
  vector<TaylorModelVec> & v = parseResult.pipes;
  logger.log(parseResult.pipes.size());
  for(vector<TaylorModelVec>::iterator pipeIt = v.begin();
      pipeIt < v.end(); pipeIt++) {
    logger.logTMV("p", *pipeIt);
  }
  exit(0);
}


















