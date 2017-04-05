#!/usr/bin/env python 
import sys
import subprocess
import math
import os

"""
  with open(nameLookup['basePlt']) as f3:
    for line in f3:
      if line.strip() == "$NAME$":
        outFile.write("set output './%s\n" %nameLookup['pngOut']);
"""

for m in sorted(os.listdir('../models/bench')):
  if not '_sw_10.' in m:
    continue
  outName = m.replace("sw_10.", "sw_rem.")
  print outName
  outFile = open(os.path.join('../models/bench', outName), 'w')
  with open(os.path.join('../models/bench', m)) as f:
    for line in f:
      if 'shrink wrapping' in line:
        outFile.write('  shrink wrapping rem\n')
        continue
      if 'out' in line:
        outFile.write(line.replace("sw_10", "sw_rem"))
        continue
      outFile.write(line)
    outFile.close()
    
