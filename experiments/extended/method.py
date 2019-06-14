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


modelTypes = ["Lotka_Volterra", "sq_deg_long", "lin_dep"]
modelTypes += ["and_or_v2"]
modelTypes += ["and_v3"]

#modelTypes = ["Lotka_Volterra", "sq_deg_long", "lin_dep"]

modelTypes += ["pair_dep"]
modelTypes += ["jet_engine"]
modelTypes += ["lin"]
modelTypes += ["bouncing_ball"]
modelTypes += ["Brusselator"]
modelTypes += ["buckling_column"]
modelTypes += ["moore_rot", "moore_rot_point"]

modelTypes += ["cruise_control"]
modelTypes += ["diabetic_1"]
modelTypes += ["diabetic_2"]
modelTypes += ["filtered_oscillator_4"]
modelTypes += ["filtered_oscillator_8"]
modelTypes += ["filtered_oscillator_16"]
modelTypes += ["Lorentz"]
modelTypes += ["neuron_I"]
modelTypes += ["neuron_II"]
modelTypes += ["nonholonomic"]
modelTypes += ["rod_reactor"]
modelTypes += ["Roessler"]
modelTypes += ["switching_5"]
modelTypes += ["two_tanks"]
modelTypes += ["vanderpol"]
modelTypes += ["vehicle_platoon_3"]

algos = ["flowid", "flowqr"]
algos = ["id", "pa", "qr", "tqr"]
algos += ["nop", "sw1", "sw2", "sw5", "sw10"]


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
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_method")
elif args.action == 'both':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_method")




















