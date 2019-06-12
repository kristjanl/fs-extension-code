#!/usr/bin/env python 
import sys
import subprocess
import math
import os

def gnuplot_for_variable(args, nameLookup):
  var1 = int(args['var1'])
  var2 = int(args['var2'])
  
  #used global variables: inName, name, script_dir
  
  #find the directory of the csv files in case script was run from 
  #other directory
  
  outFile = open(nameLookup['pltOut'], 'w')
  
  with open(nameLookup['basePlt']) as f3:
    for line in f3:
      if line.strip() == "$NAME$":
        outFile.write("set output './%s\n" %nameLookup['pngOut']);
      elif line.strip() == "$DATA$":
        write_data(outFile, nameLookup['inFileName'], args, var1, var2)
      elif line.strip() == "$VAR_LABEL$":
        outFile.write('set ylabel "x' + str(var1) + '"\n')
      elif line.strip() == "$RANGE$":
        if 'xMin' in args:
          outFile.write('set xrange [%s:%s]\n'%(args['xMin'],args['xMax']))
          outFile.write('set yrange [%s:%s]\n'%(args['yMin'],args['yMax']))
      else:
        outFile.write(line)
  outFile.close()
  
def write_data(outFile, inFileName, args, var1, var2):  
  start = float(args['start'])
  end = float(args['end'])
  
  tend = 1000000000
  if 'tend' in args:
    tend = float(args['tend'])
    
  tstart = -1
  if 'tstart' in args:
    tstart = float(args['tstart'])
  
  varY = var2
  varYLower = 2*varY
  varYUpper = 2*varY + 1
  
  varX = var1
  varXLower = 2*varX
  varXUpper = 2*varX + 1
  
  
  with open(inFileName) as f:
    for timesteps, l in enumerate(f):
      pass
      
  startStep = int(math.floor(timesteps*start))
  endStep = int(math.floor(timesteps*end))
  
  """
  if(tend != 1000000000):
    print "printing until t=%s" %(tend)
  else:
    print "printing (%s, %s) out of %s" %(startStep, endStep, timesteps)
  """
  with open(inFileName) as inFile:
    for line in inFile:
      data = line.split(",") 
      time = float(data[0])
      
      if(time > tend):
        break
        
      
      startStep = startStep - 1
      endStep = endStep - 1
      
      
      if(time < tstart):
        print "skipping"
        continue
      
      if(startStep >= 0):
        continue
        
      #if time < 150:
      #  continue
      
      #write all corners of the rectangle
      outFile.write("%s %s\n" %(data[varYLower],data[varXLower]))
      outFile.write("%s %s\n" %(data[varYLower],data[varXUpper]))
      outFile.write("%s %s\n" %(data[varYUpper],data[varXUpper]))
      outFile.write("%s %s\n" %(data[varYUpper],data[varXLower]))
      outFile.write("%s %s\n" %(data[varYLower],data[varXLower]))
      outFile.write("\n\n")
      
      if(endStep < -1):
        break
      

def getNameMap(argv, args):
  scriptArg = argv[0]
  inArg = argv[1]

  
  scriptDir = scriptArg[:scriptArg.rfind('/') + 1]
  inDir = inArg[:inArg.rfind('/') + 1]
  inParent = inDir[:-5]
  
  pltsDir = "%splts"%inParent
  imagesDir = "%simages"%inParent
  
  if not os.path.exists(pltsDir):
    os.makedirs(pltsDir)
  if not os.path.exists(imagesDir):
    os.makedirs(imagesDir)
  
  
  inName = inArg[inArg.rfind('/')+1:-4]
  
  map = {}
  map['scriptDir'] = scriptDir
  map['basePlt'] = "%sbase.plt" %scriptDir
  map['inFileName'] = inArg
  map['inDir'] = inDir
  map['inName'] = inName
  map['pltOut'] = "%s/%s.plt" %(pltsDir, inName)
  
  if 'suffix' in args:
    map['pngOut'] = "%s/%s%s.png" %(imagesDir, inName, args['suffix'])
  else:
    map['pngOut'] = "%s/%s.png" %(imagesDir, inName)
  
  #print "reading '%s'" %map['inFileName']
  #print "writing '%s'" %map['pltOut']
  #print "writing '%s'" %map['pngOut']
  
  return map
  
  
def parseArgs(args):
# var1, var2, start, end
  map = {}
  for arg in args[2:]:
    (key,val) = arg.split('=')
    map[key] = val
  if 'var2' not in map:
    map['var2'] = '0'
  if 'start' not in map:
    map['start'] = '0.0'
  if 'end' not in map:
    map['end'] = '1.0'
  return map
  
args = parseArgs(sys.argv)

nameLookup = getNameMap(sys.argv, args)

gnuplot_for_variable(args, nameLookup)

subprocess.call(["gnuplot", nameLookup['pltOut']])
