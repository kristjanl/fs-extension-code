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

models = models + ['Brusselator_QR_flow.model']
models = models + ['Lorentz_QR_flow.model']
models = models + ['Lotka_Volterra_QR_flow.model']
models = models + ['Roessler_QR_flow.model']
#models = models + ['biology_II_QR_flow.model']
#models = models + ['biology_I_QR_flow.model']
models = models + ['bouncing_ball_QR_flow.model']
models = models + ['buckling_column_QR_flow.model']
#models = models + ['coupledVanderPol_QR_flow.model']
models = models + ['cruise_control_QR_flow.model']
models = models + ['diabetic_1_QR_flow.model']
models = models + ['diabetic_2_QR_flow.model']
models = models + ['filtered_oscillator_16_QR_flow.model']
#models = models + ['filtered_oscillator_32_QR_flow.model']
models = models + ['filtered_oscillator_4_QR_flow.model']
models = models + ['filtered_oscillator_8_QR_flow.model']
models = models + ['jet_engine_QR_flow.model']
#models = models + ['lacoperon_QR_flow.model']
models = models + ['moore_rot_point_QR_flow.model']
models = models + ['moore_rot_QR_flow.model']
models = models + ['neuron_II_QR_flow.model']
models = models + ['neuron_I_QR_flow.model']
models = models + ['nonholonomic_QR_flow.model']
models = models + ['rod_reactor_QR_flow.model']
models = models + ['switching_5_QR_flow.model']
models = models + ['two_tanks_QR_flow.model']
models = models + ['vanderpol_QR_flow.model']
models = models + ['vehicle_platoon_3_QR_flow.model']



for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


