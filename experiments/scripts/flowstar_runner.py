#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

def runFlowstar(modelDir, flowstar, models, doLog=False):
  srcDir = os.path.dirname(flowstar)
  
  run = subprocess.Popen(["make", "-C", srcDir])
  run.wait()
  subprocess.call('date | cat >> log.txt', shell=True)
  for model in models:
    print "======%s======" %model
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()
    if doLog:
      subprocess.call('echo "%s" >> log.txt'%(model), shell=True)
    
  print "flowstar running done for %s models" %len(models)


