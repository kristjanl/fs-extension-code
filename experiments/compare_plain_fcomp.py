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

pairs = pairs + [('test_comp.model', 'test_nocomp.model')]

pairs = []

pairs = pairs + [('Brusselator_plain.model', 'Brusselator_fcomp.model')]
pairs = pairs + [('Lorentz_plain.model', 'Lorentz_fcomp.model')]
pairs = pairs + [('Lotka_Volterra_plain.model', 'Lotka_Volterra_fcomp.model')]
pairs = pairs + [('Roessler_plain.model', 'Roessler_fcomp.model')]
pairs = pairs + [('biology_II_plain.model', 'biology_II_fcomp.model')]
pairs = pairs + [('biology_I_plain.model', 'biology_I_fcomp.model')]
pairs = pairs + [('bouncing_ball_plain.model', 'bouncing_ball_fcomp.model')]
pairs = pairs + [('buckling_column_plain.model', 'buckling_column_fcomp.model')]
pairs = pairs + [('coupledVanderPol_plain.model', 'coupledVanderPol_fcomp.model')]
pairs = pairs + [('cruise_control_plain.model', 'cruise_control_fcomp.model')]
pairs = pairs + [('diabetic_1_plain.model', 'diabetic_1_fcomp.model')]
pairs = pairs + [('diabetic_2_plain.model', 'diabetic_2_fcomp.model')]
pairs = pairs + [('filtered_oscillator_16_plain.model', 'filtered_oscillator_16_fcomp.model')]
#pairs = pairs + [('filtered_oscillator_32_plain.model', 'filtered_oscillator_32_fcomp.model')]
pairs = pairs + [('filtered_oscillator_4_plain.model', 'filtered_oscillator_4_fcomp.model')]
pairs = pairs + [('filtered_oscillator_8_plain.model', 'filtered_oscillator_8_fcomp.model')]
pairs = pairs + [('jet_engine_plain.model', 'jet_engine_fcomp.model')]
pairs = pairs + [('lacoperon_plain.model', 'lacoperon_fcomp.model')]
pairs = pairs + [('moore_rot_point_plain.model', 'moore_rot_point_fcomp.model')]
pairs = pairs + [('moore_rot_plain.model', 'moore_rot_fcomp.model')]
pairs = pairs + [('neuron_II_plain.model', 'neuron_II_fcomp.model')]
pairs = pairs + [('neuron_I_plain.model', 'neuron_I_fcomp.model')]
pairs = pairs + [('nonholonomic_plain.model', 'nonholonomic_fcomp.model')]
pairs = pairs + [('rod_reactor_plain.model', 'rod_reactor_fcomp.model')]
pairs = pairs + [('switching_5_plain.model', 'switching_5_fcomp.model')]
pairs = pairs + [('two_tanks_plain.model', 'two_tanks_fcomp.model')]
pairs = pairs + [('vanderpol_plain.model', 'vanderpol_fcomp.model')]
pairs = pairs + [('vehicle_platoon_3_plain.model', 'vehicle_platoon_3_fcomp.model')]

pairs = []
pairs = pairs + [('Brusselator_plain.model', 'Brusselator_fcomp.model')]
pairs = pairs + [('Lorentz_plain.model', 'Lorentz_fcomp.model')]
pairs = pairs + [('Lotka_Volterra_plain.model', 'Lotka_Volterra_fcomp.model')]


tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_plain_fcomp")

tf.write_table("experiments_plain_fcomp.html", modelDir, pairs, "_plain_fcomp")

  
  
  










