#!/usr/bin/env python

import sys
import os

models = []
models = models + ['Brusselator_sw_rem.model']
models = models + ['Lorentz_sw_rem.model']
models = models + ['Lotka_Volterra_sw_rem.model']
models = models + ['Roessler_sw_rem.model']
models = models + ['biology_II_sw_rem.model']
models = models + ['biology_I_sw_rem.model']
models = models + ['bouncing_ball_sw_rem.model']
models = models + ['buckling_column_sw_rem.model']
models = models + ['coupledVanderPol_sw_rem.model']
models = models + ['cruise_control_sw_rem.model']
models = models + ['diabetic_1_sw_rem.model']
models = models + ['diabetic_2_sw_rem.model']
models = models + ['filtered_oscillator_16_sw_rem.model']
models = models + ['filtered_oscillator_32_sw_rem.model']
models = models + ['filtered_oscillator_4_sw_rem.model']
models = models + ['filtered_oscillator_8_sw_rem.model']
models = models + ['jet_engine_sw_rem.model']
models = models + ['lacoperon_sw_rem.model']
models = models + ['moore_rot_point_sw_rem.model']
models = models + ['moore_rot_sw_rem.model']
models = models + ['neuron_II_sw_rem.model']
models = models + ['neuron_I_sw_rem.model']
models = models + ['nonholonomic_sw_rem.model']
models = models + ['rod_reactor_sw_rem.model']
models = models + ['switching_5_sw_rem.model']
models = models + ['two_tanks_sw_rem.model']
models = models + ['vanderpol_sw_rem.model']
models = models + ['vehicle_platoon_3_sw_rem.model']

models = models + ['moore_rot_sw_rem.model']
models = models + ['moore_rot_point_sw_rem.model']
models = models + ['Brusselator_sw_rem.model']
models = models + ['Brusselator_sw_rem.model']


def change(models):
  for m in models:
    print m
    
    inName = os.path.join('bench', m)
    
    outName = inName.replace("sw_rem", "QR")
    
    outFile = open(outName, 'w')
    with open(inName) as f:
      for line in f:
        if 'shrink wrapping' in line:
          outFile.write('    QR precondition\n')
          continue
        if 'output' in line:
          outFile.write(line.replace("sw_rem", "QR"))
          continue
        outFile.write(line)
      outFile.close()
      
change(models)
   