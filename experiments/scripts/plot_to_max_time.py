#!/usr/bin/env python 
import sys
import subprocess
import math
import os
  
import my_functions as fs
  
tend = "tend=%s" %fs.find_max_time(sys.argv[1], sys.argv[2])
#print "%s" %tend
scriptDir = sys.argv[0][:sys.argv[0].rfind('/') + 1]

#print "scriptDir: %s" %scriptDir

#../models/csvs/plot_to_max_time.py csvs/sm_vanderpol_plain.csv var1=1 var2=0 tend=0.4

#print sys.argv[1:]

args1 = [scriptDir + 'gnuplot_var.py'] + [sys.argv[1]] + sys.argv[3:] + [tend]
args2 = [scriptDir + 'gnuplot_var.py'] + [sys.argv[2]] + sys.argv[3:] + [tend]
#print "args1: %s" %args1
#print "args2: %s" %args2


no_sw_plot_p = subprocess.Popen(args1)
no_sw_plot_p.wait()

sw_plot_p = subprocess.Popen(args2)
sw_plot_p.wait()
