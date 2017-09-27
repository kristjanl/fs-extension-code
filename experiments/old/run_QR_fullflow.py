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

models = models + ['Brusselator_QR_fullflow.model']
models = models + ['Lorentz_QR_fullflow.model']
models = models + ['Lotka_Volterra_QR_fullflow.model']
models = models + ['Roessler_QR_fullflow.model']
#models = models + ['biology_II_QR_fullflow.model']
#models = models + ['biology_I_QR_fullflow.model']
models = models + ['bouncing_ball_QR_fullflow.model']
models = models + ['buckling_column_QR_fullflow.model']
#models = models + ['coupledVanderPol_QR_fullflow.model']
models = models + ['cruise_control_QR_fullflow.model']
models = models + ['diabetic_1_QR_fullflow.model']
models = models + ['diabetic_2_QR_fullflow.model']
models = models + ['filtered_oscillator_16_QR_fullflow.model']
#models = models + ['filtered_oscillator_32_QR_fullflow.model']
models = models + ['filtered_oscillator_4_QR_fullflow.model']
models = models + ['filtered_oscillator_8_QR_fullflow.model']
models = models + ['jet_engine_QR_fullflow.model']
#models = models + ['lacoperon_QR_fullflow.model']
models = models + ['moore_rot_point_QR_fullflow.model']
models = models + ['moore_rot_QR_fullflow.model']
models = models + ['neuron_II_QR_fullflow.model']
models = models + ['neuron_I_QR_fullflow.model']
models = models + ['nonholonomic_QR_fullflow.model']
models = models + ['rod_reactor_QR_fullflow.model']
models = models + ['switching_5_QR_fullflow.model']
models = models + ['two_tanks_QR_fullflow.model']
models = models + ['vanderpol_QR_fullflow.model']
models = models + ['vehicle_platoon_3_QR_fullflow.model']



for model in models:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
  
print "flowstar running done for models: %s" %models


