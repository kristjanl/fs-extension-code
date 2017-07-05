#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

# directory of models
modelDir = os.path.join("..", "models", "bench")
flowstar = os.path.join("..", "src", "flowstar")

models = []

models = models + ['lin_3_comp.model']
models = models + ['lin_3_nocomp.model']
models = models + ['lin_20_comp.model']
models = models + ['lin_20_nocomp.model']

models = models + ['sq_deg_3_comp.model']
models = models + ['sq_deg_3_nocomp.model']
models = models + ['sq_deg_6_comp.model']
models = models + ['sq_deg_6_nocomp.model']
models = models + ['sq_deg_20_comp.model']
models = models + ['sq_deg_20_nocomp.model']


models = models + ['nl_4_comp.model']
models = models + ['nl_4_nocomp.model']
models = []
models = models + ['nl_20_comp.model']
models = models + ['nl_20_nocomp.model']

for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


