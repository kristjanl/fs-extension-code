#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os
import argparse


# directory of scripts
scriptsDir = os.path.join("..", "scripts")
sys.path.append(scriptsDir)

import flowstar_runner
import comparer



# directory of models
modelDir = os.path.join("..", "..", "models", "bench")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["moore_rot"]
modelTypes = ["and_fast_out_high", "and_fast_out_low"]
modelTypes = ["and_lin", "and_nl"]
modelTypes = ["bouncing_ball", "Brusselator", "buckling_column"]
modelTypes = ["test", "bouncing_ball", "Brusselator", "moore_rot", "buckling_column"]
#, "coupledVanderPol",
modelTypes = ["cruise_control", "diabetic_1", "diabetic_2", "filtered_oscillator_4", "filtered_oscillator_8", "filtered_oscillator_16", "jet_engine", "lacoperon", "Lorentz", "Lotka_Volterra", "moore_rot_point", "neuron_I", "neuron_II", "nonholonomic", "Roessler", "rod_reactor", "switching_5", "two_tanks", "vanderpol", "vehicle_platoon_3"]



modelTypes = ["moore_rot", "and_fast_out_high", "and_fast_out_low", "and_lin", "and_nl", "test", "bouncing_ball", "Brusselator", "moore_rot", "buckling_column", "cruise_control", "diabetic_1", "diabetic_2", "filtered_oscillator_4", "filtered_oscillator_8", "filtered_oscillator_16", "jet_engine", "Lorentz", "Lotka_Volterra", "moore_rot_point", "neuron_I", "neuron_II", "nonholonomic", "Roessler", "rod_reactor", "switching_5", "two_tanks", "vanderpol", "vehicle_platoon_3", "lin", "lin_dep", "pair_dep", "sq_deg_long"]

modelTypes2 = ["Brusselator", "buckling_column", "jet_engine", "Lorentz", "Lotka_Volterra", "moore_rot_point", "switching_5", "two_tanks", "vanderpol", "vehicle_platoon_3"]

modelTypes = ["Brusselator"]


algos = ["plain", "sw_10", "sw_rem"]
algos = ["plain", "sw_10", "sw_rem"]
#algos = ["qrplain"]
#algos = ["qr1", "qr2", "qr3"]


infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
#    ("int time", "Int t"), 
#    ("remap 1", "Remap1"), 
#    ("evaluate t", "Eval@t"),
#    ("precond time", "Precond t"), 
#    ("remap 2", "Remap2"),
]

models = ["%s_%s.model"%(modelType, algo) \
    for modelType in modelTypes \
    for algo in algos]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare4(modelDir, scriptsDir, modelTypes, None, algos, infoFields)




















