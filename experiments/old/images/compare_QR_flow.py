#!/usr/bin/env python 
import sys
import os

# directory of scripts
scriptsDir = "scripts"
sys.path.append(scriptsDir)

# directory of models
modelDir = os.path.join("..", "models", "bench")

import table_functions as tf

pairs = []

pairs = pairs + [('Brusselator_QR.model', 'Brusselator_QR_flow.model')]
pairs = pairs + [('Lorentz_QR.model', 'Lorentz_QR_flow.model')]
pairs = pairs + [('Lotka_Volterra_QR.model', 'Lotka_Volterra_QR_flow.model')]
pairs = pairs + [('Roessler_QR.model', 'Roessler_QR_flow.model')]
pairs = pairs + [('biology_II_QR.model', 'biology_II_QR_flow.model')]
pairs = pairs + [('biology_I_QR.model', 'biology_I_QR_flow.model')]
pairs = pairs + [('bouncing_ball_QR.model', 'bouncing_ball_QR_flow.model')]
pairs = pairs + [('buckling_column_QR.model', 'buckling_column_QR_flow.model')]
pairs = pairs + [('coupledVanderPol_QR.model', 'coupledVanderPol_QR_flow.model')]
pairs = pairs + [('cruise_control_QR.model', 'cruise_control_QR_flow.model')]
pairs = pairs + [('diabetic_1_QR.model', 'diabetic_1_QR_flow.model')]
pairs = pairs + [('diabetic_2_QR.model', 'diabetic_2_QR_flow.model')]
pairs = pairs + [('filtered_oscillator_16_QR.model', 'filtered_oscillator_16_QR_flow.model')]
#pairs = pairs + [('filtered_oscillator_32_QR.model', 'filtered_oscillator_32_QR_flow.model')]
pairs = pairs + [('filtered_oscillator_4_QR.model', 'filtered_oscillator_4_QR_flow.model')]
pairs = pairs + [('filtered_oscillator_8_QR.model', 'filtered_oscillator_8_QR_flow.model')]
pairs = pairs + [('jet_engine_QR.model', 'jet_engine_QR_flow.model')]
pairs = pairs + [('lacoperon_QR.model', 'lacoperon_QR_flow.model')]
pairs = pairs + [('moore_rot_point_QR.model', 'moore_rot_point_QR_flow.model')]
pairs = pairs + [('moore_rot_QR.model', 'moore_rot_QR_flow.model')]
pairs = pairs + [('neuron_II_QR.model', 'neuron_II_QR_flow.model')]
pairs = pairs + [('neuron_I_QR.model', 'neuron_I_QR_flow.model')]
pairs = pairs + [('nonholonomic_QR.model', 'nonholonomic_QR_flow.model')]
pairs = pairs + [('rod_reactor_QR.model', 'rod_reactor_QR_flow.model')]
pairs = pairs + [('switching_5_QR.model', 'switching_5_QR_flow.model')]
pairs = pairs + [('two_tanks_QR.model', 'two_tanks_QR_flow.model')]
pairs = pairs + [('vanderpol_QR.model', 'vanderpol_QR_flow.model')]
pairs = pairs + [('vehicle_platoon_3_QR.model', 'vehicle_platoon_3_QR_flow.model')]

tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_QR_flow")

tf.write_table("experiments_QR_flow.html", modelDir, pairs, "_QR_flow")

  
  
  










