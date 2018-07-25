#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os
import argparse


# directory of scripts
scriptsDir = os.path.join("..", "..", "experiments", "scripts")
sys.path.append(scriptsDir)


import flowstar_runner


modelDir = os.path.join("..", "compositional", "id_comp_precond")
modelDir = os.path.join(".")
modelDir = os.path.join("..", "testing")
flowstar = os.path.join("..", "..", "src", "flowstar")

parser = argparse.ArgumentParser()
#parser.add_argument('modelType')
#parser.add_argument('--dim', type=int, default=2)
parser.add_argument('--only', choices=['flow', 'my'])
args = parser.parse_args()

#def runFlowstar(modelType, dim):
def runFlowstar():
  print "starting pair runner (python)"
  flowModel = "1_fl.model"
  myModel = "1_my.model"
  models = [flowModel, myModel]
  if args.only == "flow":
    models = models[:1]
  if args.only == "my":
    models = models[1:]
    
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  sys.exit(0)
  for model in models:
    print model
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()

runFlowstar()
