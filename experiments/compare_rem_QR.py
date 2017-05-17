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

pairs = pairs + [('Brusselator_sw_rem.model', 'Brusselator_QR.model')]
pairs = pairs + [('Lorentz_sw_rem.model', 'Lorentz_QR.model')]
pairs = pairs + [('Lotka_Volterra_sw_rem.model', 'Lotka_Volterra_QR.model')]
pairs = pairs + [('Roessler_sw_rem.model', 'Roessler_QR.model')]
pairs = pairs + [('biology_II_sw_rem.model', 'biology_II_QR.model')]
pairs = pairs + [('biology_I_sw_rem.model', 'biology_I_QR.model')]
pairs = pairs + [('bouncing_ball_sw_rem.model', 'bouncing_ball_QR.model')]
pairs = pairs + [('buckling_column_sw_rem.model', 'buckling_column_QR.model')]
pairs = pairs + [('coupledVanderPol_sw_rem.model', 'coupledVanderPol_QR.model')]
pairs = pairs + [('cruise_control_sw_rem.model', 'cruise_control_QR.model')]
pairs = pairs + [('diabetic_1_sw_rem.model', 'diabetic_1_QR.model')]
pairs = pairs + [('diabetic_2_sw_rem.model', 'diabetic_2_QR.model')]
pairs = pairs + [('filtered_oscillator_16_sw_rem.model', 'filtered_oscillator_16_QR.model')]
#pairs = pairs + [('filtered_oscillator_32_sw_rem.model', 'filtered_oscillator_32_QR.model')]
pairs = pairs + [('filtered_oscillator_4_sw_rem.model', 'filtered_oscillator_4_QR.model')]
pairs = pairs + [('filtered_oscillator_8_sw_rem.model', 'filtered_oscillator_8_QR.model')]
pairs = pairs + [('jet_engine_sw_rem.model', 'jet_engine_QR.model')]
pairs = pairs + [('lacoperon_sw_rem.model', 'lacoperon_QR.model')]
pairs = pairs + [('moore_rot_point_sw_rem.model', 'moore_rot_point_QR.model')]
pairs = pairs + [('moore_rot_sw_rem.model', 'moore_rot_QR.model')]
pairs = pairs + [('neuron_II_sw_rem.model', 'neuron_II_QR.model')]
pairs = pairs + [('neuron_I_sw_rem.model', 'neuron_I_QR.model')]
pairs = pairs + [('nonholonomic_sw_rem.model', 'nonholonomic_QR.model')]
pairs = pairs + [('rod_reactor_sw_rem.model', 'rod_reactor_QR.model')]
pairs = pairs + [('switching_5_sw_rem.model', 'switching_5_QR.model')]
pairs = pairs + [('two_tanks_sw_rem.model', 'two_tanks_QR.model')]
pairs = pairs + [('vanderpol_sw_rem.model', 'vanderpol_QR.model')]
pairs = pairs + [('vehicle_platoon_3_sw_rem.model', 'vehicle_platoon_3_QR.model')]



tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_rem_QR")

tf.write_table("experiments_rem_QR.html", modelDir, pairs, "_rem_QR")

  
  
  










