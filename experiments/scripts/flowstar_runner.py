#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

def runFlowstar(modelDir, flowstar, models):
  for model in models:
    print "======%s======" %model
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()
  print "flowstar running done for models: %s" %len(models)


