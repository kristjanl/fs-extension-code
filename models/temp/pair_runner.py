#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os
import argparse


modelDir = os.path.join("..", "compositional", "id_comp_precond")
flowstar = os.path.join("..", "..", "src", "flowstar")

parser = argparse.ArgumentParser()
parser.add_argument('modelType')
parser.add_argument('--dim', type=int, default=2)
parser.add_argument('--only', choices=['first', 'second'])
args = parser.parse_args()

def runFlowstar(modelType, dim):
  print "starting pair runner (python)"
  allModel = "%s_%s_id_comp.model"%(modelType, dim)
  singleModel = "%s_%s_cid_comp.model"%(modelType, dim)
  flowModel = "%s_%s_id_flow.model"%(modelType, dim)
  models = [flowModel, singleModel]
  if args.only == "first":
    models = models[:1]
  if args.only == "second":
    models = models[1:]
  for model in models:
    print model
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()

runFlowstar(args.modelType, args.dim)
