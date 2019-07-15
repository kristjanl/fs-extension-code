#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os
import argparse

modelTypes = []

modelTypes += ["bouncing_ball_id"]
modelTypes += ["cruise_control_id"]
modelTypes += ["diabetic_1_id"]
modelTypes += ["diabetic_2_id"]
modelTypes += ["filtered_oscillator_4_id"]
modelTypes += ["filtered_oscillator_8_id"]
modelTypes += ["filtered_oscillator_16_id"]
modelTypes += ["filtered_oscillator_32_id"]
modelTypes += ["nonholonomic_id"]
modelTypes += ["rod_reactor_id"]
modelTypes += ["two_tanks_id"]
modelTypes += ["and_or_v2_id"]
modelTypes += ["and_v3_id"]

modelTypes = map(lambda s: s + "_fcomp.plt", modelTypes)

for m in modelTypes:
  subprocess.call(["gnuplot", m])



