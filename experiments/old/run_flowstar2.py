#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

# directory of models
modelDir = os.path.join("..", "models", "bench")
flowstar = os.path.join("..", "src", "flowstar")

pairs = []

pairs = pairs + ['Brusselator_sw_10.model']
pairs = pairs + ['Lorentz_sw_10.model']
pairs = pairs + ['Lotka_Volterra_sw_10.model']
pairs = pairs + ['Roessler_sw_10.model']
#pairs = pairs + ['biology_II_sw_10.model']
#pairs = pairs + ['biology_I_sw_10.model']
pairs = pairs + ['bouncing_ball_sw_10.model']
pairs = pairs + ['buckling_column_sw_10.model']
#pairs = pairs + ['coupledVanderPol_sw_10.model']
pairs = pairs + ['cruise_control_sw_10.model']
pairs = pairs + ['diabetic_1_sw_10.model']
pairs = pairs + ['diabetic_2_sw_10.model']
pairs = pairs + ['filtered_oscillator_16_sw_10.model']
pairs = pairs + ['filtered_oscillator_32_sw_10.model']
pairs = pairs + ['filtered_oscillator_4_sw_10.model']
pairs = pairs + ['filtered_oscillator_8_sw_10.model']
pairs = pairs + ['jet_engine_sw_10.model']
pairs = pairs + ['lacoperon_sw_10.model']
pairs = pairs + ['moore_rot_point_sw_10.model']
pairs = pairs + ['moore_rot_sw_10.model']
pairs = pairs + ['neuron_II_sw_10.model']
pairs = pairs + ['neuron_I_sw_10.model']
pairs = pairs + ['nonholonomic_sw_10.model']
pairs = pairs + ['rod_reactor_sw_10.model']
pairs = pairs + ['switching_5_sw_10.model']
pairs = pairs + ['two_tanks_sw_10.model']
pairs = pairs + ['vanderpol_sw_10.model']
#pairs = pairs + ['vehicle_platoon_3_sw_10.model']


for model in pairs:
  print "======%s======" %model
  modelFile = os.path.join(modelDir, model)
  modelStream = open(modelFile)
  run = subprocess.Popen(flowstar, stdin=modelStream)
  run.wait()
print "flowstar running done for models: %s" %pairs


