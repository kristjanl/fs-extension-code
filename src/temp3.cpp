
void TaylorModelVec::Picard_ctrunc_normal(TaylorModelVec & result, 
    vector<RangeTree *> & trees, MyComponent *comp2, 
    const vector<Interval> & polyRange, const vector<Interval> & step_exp_table, 
    const int numVars, const int order, const Interval & cutoff_threshold) 
    const {
  mreset(old);
  mdisable();
  minc();
  mlog1("Picard_ctrunc_normal (2, i vec)");
  //mlog("x0", comp2->initSet);
  //mlog("ode", comp2->odes);
  //mlog("this", *this);
  //mlog("polyRange", polyRange);
  //mlog1(sbuilder() << "order: " << order);
  //mlog1(sbuilder() << "numVars: " << numVars);
  //mlog("comp", comp2->solveIndexes);
	TaylorModelVec inserted(comp2->solveIndexes.size());

	trees.clear();
	for(int i=0; i<comp2->odes.size(); ++i) {
		trees.push_back(NULL);
	}

	if(order <= 1) {
		//for(int i=0; i<ode.size(); ++i) {
	  for(int j = 0; j < comp2->solveIndexes.size(); j++) {
	    int var = comp2->solveIndexes[j];
	    //mlog1(sbuilder() << "var: " << var);
			TaylorModel tmTemp;
			//mlog1(comp2->odes[var].toString());
			comp2->odes[var].insert_ctrunc_normal(tmTemp, trees[var], *this, 
          polyRange, step_exp_table, numVars, 0, cutoff_threshold);
			inserted.tms.push_back(tmTemp);
			//mlog(sbuilder() << "tm[" << var << "]", tmTemp);
		}
	} else {
		//for(int i=0; i<ode.size(); ++i) {
	  for(int j = 0; j < comp2->solveIndexes.size(); j++) {
	    int var = comp2->solveIndexes[j];
	    mlog1(sbuilder() << "var: " << var);
			mlog1(comp2->odes[var].toString());
      mlog1(sbuilder() << numVars);
		  //order - 1, since the integration makes the degree 1 higher
			TaylorModel tmTemp;
			comp2->odes[var].insert_ctrunc_normal(tmTemp, trees[var], *this, 
          polyRange, step_exp_table, numVars, order-1, cutoff_threshold);
			inserted.tms.push_back(tmTemp);
		}
	}
	//mlog("inserted(TM)", inserted);

	TaylorModelVec integrated(inserted.tms.size());
	for(int i = 0; i < inserted.tms.size(); i++) {
		TaylorModel tmTemp;
		inserted.tms[i].integral(tmTemp, step_exp_table[1]);
		integrated.tms.push_back(tmTemp);
	}
	//mlog("integrated(TM)", integrated);
	TaylorModelVec added(comp2->solveIndexes.size());
	for(int i = 0; i < comp2->solveIndexes.size(); i++) {
	  //i is index in integrated, var is index in original
	  int var = comp2->solveIndexes[i];
		TaylorModel tmTemp;
		comp2->initSet.tms[var].add(tmTemp, integrated.tms[i]);
		added.tms.push_back(tmTemp);
	}
	
	result.clear();
	result.tms.reserve(tms.size());
	for(int i = 0; i < tms.size(); i++) {
    //mlog1(sbuilder() << "i: " << i);
    vector<int>::iterator it = find(comp2->solveIndexes.begin(), comp2->solveIndexes.end(), i);
    if(it == comp2->solveIndexes.end()) {
      result.tms.push_back(tms[i]);
    }else {
      int indexInComp = it - comp2->solveIndexes.begin();
      result.tms.push_back(added.tms[indexInComp]);
    }
	}
	//mlog("result", result);
	//mforce("result", result);
	mdec();
	mrestore(old);
}
