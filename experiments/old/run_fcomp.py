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
models = models + ['Brusselator_fcomp.model']
models = models + ['Lorentz_fcomp.model']
models = models + ['Lotka_Volterra_fcomp.model']
models = models + ['Roessler_fcomp.model']
models = models + ['bouncing_ball_fcomp.model']
models = models + ['buckling_column_fcomp.model']
models = models + ['cruise_control_fcomp.model']
models = models + ['diabetic_1_fcomp.model']
models = models + ['diabetic_2_fcomp.model']
models = models + ['filtered_oscillator_16_fcomp.model']
models = models + ['filtered_oscillator_4_fcomp.model']
models = models + ['filtered_oscillator_8_fcomp.model']
models = models + ['jet_engine_fcomp.model']
models = models + ['moore_rot_point_fcomp.model']
models = models + ['moore_rot_fcomp.model']
models = models + ['neuron_II_fcomp.model']
models = models + ['neuron_I_fcomp.model']
models = models + ['nonholonomic_fcomp.model']
models = models + ['rod_reactor_fcomp.model']
models = models + ['switching_5_fcomp.model']
models = models + ['two_tanks_fcomp.model']
models = models + ['vanderpol_fcomp.model']
models = models + ['vehicle_platoon_3_fcomp.model']



for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


