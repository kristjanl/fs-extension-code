import sys
import subprocess
import math
import re
import os

def foo():

  flowstar = os.path.join("..", "..", "src", "flowstar")
  srcDir = os.path.dirname(flowstar)
  
  run = subprocess.Popen(["make", "-C", srcDir])
  run.wait()
  
  
  run = subprocess.Popen([flowstar, "qr.old.txt", "qr.txt", "leftStar", "9"])
  #run = subprocess.Popen([flowstar, "f_id.old.txt", "f_id.txt", "leftStar", "9"])
  #run = subprocess.Popen([flowstar, "s_id.old.txt", "s_id.txt", "leftStar", "26"])
  run.wait()
  
  """
  for model in models:
    print "======%s======" %model
    modelFile = os.path.join(modelDir, model)
    modelStream = open(modelFile)
    run = subprocess.Popen(flowstar, stdin=modelStream)
    run.wait()
  print "flowstar running done for models: %s" %len(models)
  """
foo()