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
parser.add_argument('modelType')
parser.add_argument('--dim', type=int, default=2)
args = parser.parse_args()

def runFlowstar(modelType, dim):
  print "starting pair runner (python)"
  myModel = "%s_%s_id_comp.model"%(modelType, dim)
  flowModel = "%s_%s_id_flow.model"%(modelType, dim)
  for model in [flowModel, myModel]:
    print model
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()

runFlowstar(args.modelType, args.dim)
