#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

# directory of models
modelDir = os.path.join("..", "..", "models", "compositional", "id")
flowstar = os.path.join("..", "..", "src", "flowstar")





modelTypes = ["lin", "lin_dep", "pair_dep", "sq_deg", "sq_deg_long"]
dims = [2]# + range(10, 51, 10)
dims = [2] + [10]
algos = ["flow", "comp", "nocomp"]

models = ["%s_%s_id_%s.model"%(modelType, dim, algo) for modelType in modelTypes for dim in dims for algo in algos]

for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


