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

pairs = pairs + [('moore_rot_plain.model', 'moore_rot_sw_10.model')]

#pairs = pairs + [('vanderpol_plain.model', 'vanderpol_sw_10.model')]
#pairs = pairs + [('Brusselator_plain.model', 'Brusselator_sw_10.model')]
#pairs = pairs + [('jet_engine_plain.model', 'jet_engine_sw_10.model')]
#pairs = pairs + [('buckling_column_plain.model', 'buckling_column_sw_10.model')]
#pairs = pairs + [('Lotka_Volterra_plain.model', 'Lotka_Volterra_sw_10.model')]
#pairs = pairs + [('Lorentz_plain.model', 'Lorentz_sw_10.model')]
#pairs = pairs + [('Roessler_plain.model', 'Roessler_sw_10.model')]

#slow
#pairs = pairs + [('biology_I_plain.model', 'biology_I_sw_10.model')]
#pairs = pairs + [('biology_II_plain.model', 'biology_II_sw_10.model')]

#non poly (can't run flowstar on them)
##pairs = pairs + [('lacoperon_plain.model', 'lacoperon_sw_10.model')]
##pairs = pairs + [('coupledVanderPol_plain.model', 'coupledVanderPol_sw_10.model')]

#non lin hybrid
#pairs = pairs + [('nonholonomic_plain.model', 'nonholonomic_sw_10.model')]
#pairs = pairs + [('neuron_I_plain.model', 'neuron_I_sw_10.model')]
#pairs = pairs + [('neuron_II_plain.model', 'neuron_II_sw_10.model')]
#pairs = pairs + [('diabetic_1_plain.model', 'diabetic_1_sw_10.model')]
#pairs = pairs + [('diabetic_2_plain.model', 'diabetic_2_sw_10.model')]


#lin hybrid
#pairs = pairs + [('bouncing_ball_plain.model', 'bouncing_ball_sw_10.model')]
#pairs = pairs + [('two_tanks_plain.model', 'two_tanks_sw_10.model')]
#pairs = pairs + [('rod_reactor_plain.model', 'rod_reactor_sw_10.model')]
#pairs = pairs + [('cruise_control_plain.model', 'cruise_control_sw_10.model')]
#pairs = pairs + [('switching_5_plain.model', 'switching_5_sw_10.model')]
#pairs = pairs + [('vehicle_platoon_3_plain.model', 'vehicle_platoon_3_sw_10.model')]
#pairs = pairs + [('filtered_oscillator_4_plain.model', 'filtered_oscillator_4_sw_10.model')]
#pairs = pairs + [('filtered_oscillator_8_plain.model', 'filtered_oscillator_4_sw_10.model')]
#pairs = pairs + [('filtered_oscillator_16_plain.model', 'filtered_oscillator_4_sw_10.model')]
#pairs = pairs + [('filtered_oscillator_32_plain.model', 'filtered_oscillator_4_sw_10.model')]


for (plain, sw) in pairs:
  plainFile = os.path.join(modelDir, plain)
  swFile = os.path.join(modelDir, sw)
  plainStream = open(plainFile)
  swStream = open(swFile)
  plain = subprocess.Popen(flowstar, stdin=plainStream)
  plain.wait()
  sw = subprocess.Popen(flowstar, stdin=swStream)
  sw.wait()

print "flowstar running done for models: %s" %pairs


