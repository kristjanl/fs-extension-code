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



models = models + ['Brusselator_plain.model']
models = models + ['Lorentz_plain.model']
models = models + ['Lotka_Volterra_plain.model']
models = models + ['Roessler_plain.model']
models = models + ['bouncing_ball_plain.model']
models = models + ['buckling_column_plain.model']
models = models + ['cruise_control_plain.model']
models = models + ['diabetic_1_plain.model']
models = models + ['diabetic_2_plain.model']
models = models + ['filtered_oscillator_16_plain.model']
models = models + ['filtered_oscillator_4_plain.model']
models = models + ['filtered_oscillator_8_plain.model']
models = models + ['jet_engine_plain.model']
models = models + ['moore_rot_point_plain.model']
models = models + ['moore_rot_plain.model']
models = models + ['neuron_II_plain.model']
models = models + ['neuron_I_plain.model']
models = models + ['nonholonomic_plain.model']
models = models + ['rod_reactor_plain.model']
models = models + ['switching_5_plain.model']
models = models + ['two_tanks_plain.model']
models = models + ['vanderpol_plain.model']
models = models + ['vehicle_platoon_3_plain.model']

models = []
models = models + ['Brusselator_plain.model']
models = models + ['Lorentz_plain.model']
models = models + ['Lotka_Volterra_plain.model']


for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


