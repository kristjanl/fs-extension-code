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

models = models + ['Brusselator_QR.model']
models = models + ['Lorentz_QR.model']
models = models + ['Lotka_Volterra_QR.model']
models = models + ['Roessler_QR.model']
#models = models + ['biology_II_QR.model']
#models = models + ['biology_I_QR.model']
models = models + ['bouncing_ball_QR.model']
models = models + ['buckling_column_QR.model']
#models = models + ['coupledVanderPol_QR.model']
models = models + ['cruise_control_QR.model']
models = models + ['diabetic_1_QR.model']
models = models + ['diabetic_2_QR.model']
models = models + ['filtered_oscillator_16_QR.model']
#models = models + ['filtered_oscillator_32_QR.model']
models = models + ['filtered_oscillator_4_QR.model']
models = models + ['filtered_oscillator_8_QR.model']
models = models + ['jet_engine_QR.model']
#models = models + ['lacoperon_QR.model']
models = models + ['moore_rot_point_QR.model']
models = models + ['moore_rot_QR.model']
models = models + ['neuron_II_QR.model']
models = models + ['neuron_I_QR.model']
models = models + ['nonholonomic_QR.model']
models = models + ['rod_reactor_QR.model']
models = models + ['switching_5_QR.model']
models = models + ['two_tanks_QR.model']
models = models + ['vanderpol_QR.model']
models = models + ['vehicle_platoon_3_QR.model']


for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


