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
modelDir = os.path.join("..", "..", "models", "compositional", "extended", "method")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = []

modelTypes += ["and_v3"]
modelTypes += ["and_or_v2"]
modelTypes += ["Brusselator"]
modelTypes += ["buckling_column"]
modelTypes += ["jet_engine"]
modelTypes += ["Lorentz"]
modelTypes += ["Lotka_Volterra"]
modelTypes += ["moore_rot"]
modelTypes += ["Roessler"]
modelTypes += ["vanderpol"]

modelTypes += ["bouncing_ball"]
modelTypes += ["cruise_control"]
modelTypes += ["diabetic_1"]
modelTypes += ["diabetic_2"]
modelTypes += ["filtered_oscillator_4"]
modelTypes += ["filtered_oscillator_8"]
modelTypes += ["filtered_oscillator_16"]
modelTypes += ["filtered_oscillator_32"]
modelTypes += ["neuron_I"]
modelTypes += ["neuron_II"]
modelTypes += ["nonholonomic"]
modelTypes += ["rod_reactor"]
modelTypes += ["switching_5"]
modelTypes += ["two_tanks"]
modelTypes += ["vehicle_platoon_3"]

modelTypes += ["lin"]
modelTypes += ["lin_dep"]
modelTypes += ["sq_deg_long"]
modelTypes += ["pair_dep"]

#modelTypes += ["moore_rot_point"]

#modelTypes = ["rod_reactor"]

#modelTypes = ["and_v3"]

#modelTypes = ["buckling_column"]
#modelTypes += ["Brusselator"]
#modelTypes += ["Lotka_Volterra"]
#modelTypes += ["sq_deg_long"]
#modelTypes = ["filtered_oscillator_4"]


#algos = ["flowid", "flowqr"]
algos = ["id", "pa", "qr", "tqr"]
algos += ["nop", "sw1", "sw2", "sw5", "sw10"]
#algos = ["nop", "id"]

groups = [["%s_%s.model"%(t,a) for a in algos] for t in modelTypes]
    
models = [m for g in groups for m in g]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
#    ("int time", "Int t"), 
#    ("picard poly", "p poly"),
#    ("picard decreasing", "p decr"),
#    ("picard refining", "p ref"),
#    ("precond time", "Precond t"), 
]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare', 'both'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models, doLog=True)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_method", vars=[0])
elif args.action == 'both':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_method", vars=[0])




















