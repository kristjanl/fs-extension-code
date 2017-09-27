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
pairs = pairs + [('moore_rot_point_sw_10.model', 'moore_rot_point_sw_10_infl.model')]
pairs = pairs + [('Roessler_sw_10.model', 'Roessler_sw_10_infl.model')]
pairs = pairs + [('nonholonomic_sw_10.model', 'nonholonomic_sw_10_infl.model')]
pairs = pairs + [('bouncing_ball_sw_10.model', 'bouncing_ball_sw_10_infl.model')]
pairs = pairs + [('two_tanks_sw_10.model', 'two_tanks_sw_10_infl.model')]
pairs = pairs + [('rod_reactor_sw_10.model', 'rod_reactor_sw_10_infl.model')]
pairs = pairs + [('switching_5_sw_10.model', 'switching_5_sw_10_infl.model')]
pairs = pairs + [('filtered_oscillator_4_sw_10.model', 'filtered_oscillator_4_sw_10_infl.model')]
pairs = pairs + [('filtered_oscillator_8_sw_10.model', 'filtered_oscillator_8_sw_10_infl.model')]
pairs = pairs + [('filtered_oscillator_16_sw_10.model', 'filtered_oscillator_16_sw_10_infl.model')]


tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "step_infl")

tf.write_table("experiments_step_infl.html", modelDir, pairs, "step_infl")

  
  
  










