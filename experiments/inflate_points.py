#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os.path

# directory of scripts
scriptsDir = "scripts"
sys.path.append(scriptsDir)
# directory of models
modelDir = os.path.join("..", "models", "bench")

import my_functions as fs

models = []

models = models + ['moore_rot_sw_10.model']
models = models + ['moore_rot_point_sw_10.model']

models = models + ['vanderpol_sw_10.model']
models = models + ['Brusselator_sw_10.model']
models = models + ['jet_engine_sw_10.model']
models = models + ['buckling_column_sw_10.model']
models = models + ['Lotka_Volterra_sw_10.model']
models = models + ['Lorentz_sw_10.model']
models = models + ['Roessler_sw_10.model']

#slow
models = models + ['biology_I_sw_10.model']
models = models + ['biology_II_sw_10.model']

#non poly
models = models + ['lacoperon_sw_10.model']
models = models + ['coupledVanderPol_sw_10.model']

#non lin hybrid
models = models + ['nonholonomic_sw_10.model']
models = models + ['neuron_I_sw_10.model']
models = models + ['neuron_II_sw_10.model']
models = models + ['diabetic_1_sw_10.model']
models = models + ['diabetic_2_sw_10.model']

#lin hybrid
models = models + ['bouncing_ball_sw_10.model']
models = models + ['two_tanks_sw_10.model']
models = models + ['rod_reactor_sw_10.model']
models = models + ['cruise_control_sw_10.model']
models = models + ['switching_5_sw_10.model']
models = models + ['vehicle_platoon_3_sw_10.model']
models = models + ['filtered_oscillator_4_sw_10.model']
models = models + ['filtered_oscillator_8_sw_10.model']
models = models + ['filtered_oscillator_16_sw_10.model']
models = models + ['filtered_oscillator_32_sw_10.model']

def getParam(filename, param):
  with open(filename) as f:
    for line in f:
      #print line,
      m = re.search('%s (.*)' %param, line)
      if m != None:
        return m.group(1)
        

condRegex = re.compile("\[[^\]]*\]", re.IGNORECASE)
regex2 = re.compile("(.*)\n", re.IGNORECASE)
nameRegex = re.compile("(.*)\.model")

for model in models:
  if 'moore' not in model:
    continue
  modelFile = os.path.join(modelDir, model)
  est = float(getParam(modelFile, "remainder estimation"))
  output = getParam(modelFile, "output")
  content = []
  needToWrite = False
  with open(modelFile) as f:
    for line in f:
      m = re.search('([A-z0-9]*) in \[(.*),(.*)\]', line)
      if m != None and m.group(2).strip() == m.group(3).strip():
        val = float(m.group(2))
        content = content + [condRegex.sub("[%s,%s]"%(val-est,val+est), line)]
        needToWrite = True
      else:
        content = content + [line]
  
  
  
  if needToWrite:
    print model
    inflName = nameRegex.sub(r"\1_infl.model", model)
    outFile = open(os.path.join(modelDir, inflName), 'w')
    for line in content:
      if "output %s" %output in line:
        line = regex2.sub(r"\1_infl\n", line)
        #print line
      outFile.write(line)
    outFile.close()

















