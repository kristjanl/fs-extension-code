#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

import flowstar_runner as fRunner
import comparer

import argparse

# directory of scripts
scriptsDir = os.path.join("..", "scripts")
sys.path.append(scriptsDir)

# directory of models
modelDir = os.path.join("..", "..", "models", "compositional", "id")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["lin", "lin_dep", "pair_dep", "sq_deg", "sq_deg_long"]
dims = [2] + range(10, 51, 10)
algos = ["flow", "comp", "nocomp"]

models = ["%s_%s_id_%s.model"%(modelType, dim, algo) \
    for modelType in modelTypes \
    for dim in dims \
    for algo in algos]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  fRunner.runFlowstar(modelDir, flowstar, models[:3])
elif args.action == 'compare':
  print "compare"
print args

#print parser.parse_args('--foo 1'.split())
"""
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=myfoo,
                    help='sum the integers (default: find the max)')
"""
#args = parser.parse_args()

#print args.accumulate(args.integers)
#print parser.parse_args(['7', '-1', '42'])

#fRunner.runFlowstar(modelDir, flowstar, models[:3])



