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

pairs = pairs + [('Brusselator_QR_flow.model', 'Brusselator_QR_fullflow.model')]
pairs = pairs + [('Lorentz_QR_flow.model', 'Lorentz_QR_fullflow.model')]
pairs = pairs + [('Lotka_Volterra_QR_flow.model', 'Lotka_Volterra_QR_fullflow.model')]
pairs = pairs + [('Roessler_QR_flow.model', 'Roessler_QR_fullflow.model')]
pairs = pairs + [('biology_II_QR_flow.model', 'biology_II_QR_fullflow.model')]
pairs = pairs + [('biology_I_QR_flow.model', 'biology_I_QR_fullflow.model')]
pairs = pairs + [('bouncing_ball_QR_flow.model', 'bouncing_ball_QR_fullflow.model')]
pairs = pairs + [('buckling_column_QR_flow.model', 'buckling_column_QR_fullflow.model')]
pairs = pairs + [('coupledVanderPol_QR_flow.model', 'coupledVanderPol_QR_fullflow.model')]
pairs = pairs + [('cruise_control_QR_flow.model', 'cruise_control_QR_fullflow.model')]
pairs = pairs + [('diabetic_1_QR_flow.model', 'diabetic_1_QR_fullflow.model')]
pairs = pairs + [('diabetic_2_QR_flow.model', 'diabetic_2_QR_fullflow.model')]
pairs = pairs + [('filtered_oscillator_16_QR_flow.model', 'filtered_oscillator_16_QR_fullflow.model')]
#pairs = pairs + [('filtered_oscillator_32_QR_flow.model', 'filtered_oscillator_32_QR_fullflow.model')]
pairs = pairs + [('filtered_oscillator_4_QR_flow.model', 'filtered_oscillator_4_QR_fullflow.model')]
pairs = pairs + [('filtered_oscillator_8_QR_flow.model', 'filtered_oscillator_8_QR_fullflow.model')]
pairs = pairs + [('jet_engine_QR_flow.model', 'jet_engine_QR_fullflow.model')]
pairs = pairs + [('lacoperon_QR_flow.model', 'lacoperon_QR_fullflow.model')]
pairs = pairs + [('moore_rot_point_QR_flow.model', 'moore_rot_point_QR_fullflow.model')]
pairs = pairs + [('moore_rot_QR_flow.model', 'moore_rot_QR_fullflow.model')]
pairs = pairs + [('neuron_II_QR_flow.model', 'neuron_II_QR_fullflow.model')]
pairs = pairs + [('neuron_I_QR_flow.model', 'neuron_I_QR_fullflow.model')]
pairs = pairs + [('nonholonomic_QR_flow.model', 'nonholonomic_QR_fullflow.model')]
pairs = pairs + [('rod_reactor_QR_flow.model', 'rod_reactor_QR_fullflow.model')]
pairs = pairs + [('switching_5_QR_flow.model', 'switching_5_QR_fullflow.model')]
pairs = pairs + [('two_tanks_QR_flow.model', 'two_tanks_QR_fullflow.model')]
pairs = pairs + [('vanderpol_QR_flow.model', 'vanderpol_QR_fullflow.model')]
pairs = pairs + [('vehicle_platoon_3_QR_flow.model', 'vehicle_platoon_3_QR_fullflow.model')]



tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_flow_fullflow")

tf.write_table("experiments_flow_fullflow.html", modelDir, pairs, "_flow_fullflow")

  
  
  










