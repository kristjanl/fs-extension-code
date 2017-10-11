#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os
import argparse


modelDir = os.path.join("..", "compositional", "id")
flowstar = os.path.join("..", "..", "src", "flowstar")

parser = argparse.ArgumentParser()
parser.add_argument('dim', type=int)
args = parser.parse_args()

def runFlowstar(dim):
  print "starting pair runner (python)"
  modelType = "lin_dep"
  myModel = "%s_%s_id_comp.model"%(modelType, dim)
  flowModel = "%s_%s_id_flow.model"%(modelType, dim)
  
  for model in [flowModel, myModel]:
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()

runFlowstar(args.dim)
