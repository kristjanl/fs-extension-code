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

pairs = pairs + [('Brusselator_sw_10.model', 'Brusselator_sw_rem.model')]
pairs = pairs + [('Lorentz_sw_10.model', 'Lorentz_sw_rem.model')]
pairs = pairs + [('Lotka_Volterra_sw_10.model', 'Lotka_Volterra_sw_rem.model')]
pairs = pairs + [('Roessler_sw_10.model', 'Roessler_sw_rem.model')]
pairs = pairs + [('biology_II_sw_10.model', 'biology_II_sw_rem.model')]
pairs = pairs + [('biology_I_sw_10.model', 'biology_I_sw_rem.model')]
pairs = pairs + [('bouncing_ball_sw_10.model', 'bouncing_ball_sw_rem.model')]
pairs = pairs + [('buckling_column_sw_10.model', 'buckling_column_sw_rem.model')]
pairs = pairs + [('coupledVanderPol_sw_10.model', 'coupledVanderPol_sw_rem.model')]
pairs = pairs + [('cruise_control_sw_10.model', 'cruise_control_sw_rem.model')]
pairs = pairs + [('diabetic_1_sw_10.model', 'diabetic_1_sw_rem.model')]
pairs = pairs + [('diabetic_2_sw_10.model', 'diabetic_2_sw_rem.model')]
pairs = pairs + [('filtered_oscillator_16_sw_10.model', 'filtered_oscillator_16_sw_rem.model')]
#pairs = pairs + [('filtered_oscillator_32_sw_10.model', 'filtered_oscillator_32_sw_rem.model')]
pairs = pairs + [('filtered_oscillator_4_sw_10.model', 'filtered_oscillator_4_sw_rem.model')]
pairs = pairs + [('filtered_oscillator_8_sw_10.model', 'filtered_oscillator_8_sw_rem.model')]
pairs = pairs + [('jet_engine_sw_10.model', 'jet_engine_sw_rem.model')]
pairs = pairs + [('lacoperon_sw_10.model', 'lacoperon_sw_rem.model')]
pairs = pairs + [('moore_rot_point_sw_10.model', 'moore_rot_point_sw_rem.model')]
pairs = pairs + [('moore_rot_sw_10.model', 'moore_rot_sw_rem.model')]
pairs = pairs + [('neuron_II_sw_10.model', 'neuron_II_sw_rem.model')]
pairs = pairs + [('neuron_I_sw_10.model', 'neuron_I_sw_rem.model')]
pairs = pairs + [('nonholonomic_sw_10.model', 'nonholonomic_sw_rem.model')]
pairs = pairs + [('rod_reactor_sw_10.model', 'rod_reactor_sw_rem.model')]
pairs = pairs + [('switching_5_sw_10.model', 'switching_5_sw_rem.model')]
pairs = pairs + [('two_tanks_sw_10.model', 'two_tanks_sw_rem.model')]
pairs = pairs + [('vanderpol_sw_10.model', 'vanderpol_sw_rem.model')]
pairs = pairs + [('vehicle_platoon_3_sw_10.model', 'vehicle_platoon_3_sw_rem.model')]



tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "step_rem")

tf.write_table("experiments_step_remainder.html", modelDir, pairs, "step_rem")

  
  
  










