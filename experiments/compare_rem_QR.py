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

pairs = pairs + [('Brusselator_plain.model', 'Brusselator_sw_rem.model')]
pairs = pairs + [('Lorentz_plain.model', 'Lorentz_sw_rem.model')]
pairs = pairs + [('Lotka_Volterra_plain.model', 'Lotka_Volterra_sw_rem.model')]
pairs = pairs + [('Roessler_plain.model', 'Roessler_sw_rem.model')]
pairs = pairs + [('biology_II_plain.model', 'biology_II_sw_rem.model')]
pairs = pairs + [('biology_I_plain.model', 'biology_I_sw_rem.model')]
pairs = pairs + [('bouncing_ball_plain.model', 'bouncing_ball_sw_rem.model')]
pairs = pairs + [('buckling_column_plain.model', 'buckling_column_sw_rem.model')]
pairs = pairs + [('coupledVanderPol_plain.model', 'coupledVanderPol_sw_rem.model')]
pairs = pairs + [('cruise_control_plain.model', 'cruise_control_sw_rem.model')]
pairs = pairs + [('diabetic_1_plain.model', 'diabetic_1_sw_rem.model')]
pairs = pairs + [('diabetic_2_plain.model', 'diabetic_2_sw_rem.model')]
pairs = pairs + [('filtered_oscillator_16_plain.model', 'filtered_oscillator_16_sw_rem.model')]
#pairs = pairs + [('filtered_oscillator_32_plain.model', 'filtered_oscillator_32_sw_rem.model')]
pairs = pairs + [('filtered_oscillator_4_plain.model', 'filtered_oscillator_4_sw_rem.model')]
pairs = pairs + [('filtered_oscillator_8_plain.model', 'filtered_oscillator_8_sw_rem.model')]
pairs = pairs + [('jet_engine_plain.model', 'jet_engine_sw_rem.model')]
pairs = pairs + [('lacoperon_plain.model', 'lacoperon_sw_rem.model')]
pairs = pairs + [('moore_rot_point_plain.model', 'moore_rot_point_sw_rem.model')]
pairs = pairs + [('moore_rot_plain.model', 'moore_rot_sw_rem.model')]
pairs = pairs + [('neuron_II_plain.model', 'neuron_II_sw_rem.model')]
pairs = pairs + [('neuron_I_plain.model', 'neuron_I_sw_rem.model')]
pairs = pairs + [('nonholonomic_plain.model', 'nonholonomic_sw_rem.model')]
pairs = pairs + [('rod_reactor_plain.model', 'rod_reactor_sw_rem.model')]
pairs = pairs + [('switching_5_plain.model', 'switching_5_sw_rem.model')]
pairs = pairs + [('two_tanks_plain.model', 'two_tanks_sw_rem.model')]
pairs = pairs + [('vanderpol_plain.model', 'vanderpol_sw_rem.model')]
pairs = pairs + [('vehicle_platoon_3_plain.model', 'vehicle_platoon_3_sw_rem.model')]

pairs = []
pairs = pairs + [('moore_rot_sw_rem.model', 'moore_rot_QR.model')]
pairs = pairs + [('moore_rot_point_sw_rem.model', 'moore_rot_point_QR.model')]
pairs = pairs + [('Brusselator_sw_rem.model', 'Brusselator_QR.model')]
pairs = pairs + [('Brusselator_QR.model', 'Brusselator_QR.model')]




tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_rem_QR")

tf.write_table("experiments_rem_QR.html", modelDir, pairs, "_rem_QR")

  
  
  










