import sys
import math
import os
import re
import subprocess

import my_functions as fs

plotAllVars = False
plotToCommon = False
#plotToCommon = True
tableCSSFile='../scripts/table.css'

def getVars(modelFile):
  csv = os.path.join("csvs", "%s.csv" %fs.getParam(modelFile, "output"))
  store = False
  s = ""
  with open(modelFile) as f:
    for line in f:
      if 'state var' in line:
        store = True
      if 'setting' in line:
        break
      if store:
        s = s + line
  return s.split(',')

def getDimension(modelFile):
  return len(getVars(modelFile))

def getVarRange2(modelFile, csvs):
  # find time that all integrations reached
  time = min( map(lambda x: fs.find_file_max_time(x), csvs) )
  
  dim = getDimension(modelFile)
  
  ranges = map(lambda c: [], csvs)
  
  for (i, csv) in enumerate(csvs):
    if not os.path.isfile(csv):
      ranges[i] = map(lambda _:None, range(dim))
      continue
    with open(csv) as f:
      for line in f:
        data = line.strip().split(',')
        if data[-1] == '':
          data = data[0:-1]
        #line with needed time
        if float(data[0]) == float(time):
          #get variable ranges
          ranges[i] = map(lambda x: float(x), data[2 + 2*dim: 2 + 3*dim])
  
  for i in range(dim):
    #find minimum value in the ranges
    minValue = ranges[0][i]
    for r in ranges[1:]:
      if minValue > r[i] and r[i] != None:
        minValue = r[i]
    
    # mark minimum value
    for (j,_) in enumerate(ranges):
      if minValue == ranges[0][i]:
        continue
      if minValue == ranges[j][i]:
        ranges[j][i] = "<b>%s</b>"%ranges[j][i]
  return ranges

def getSkippingIndexes(dim):
  if plotAllVars:
    return range(0,dim)
  
  #include first 3
  indexes = [0, 1]
  #remove if too big
  indexes = filter(lambda x: x < dim, indexes)
  
  #0, 1, 2, 5, 10, 20, 40, 80 etc
  powers = [5*2**i for i in range( int(math.ceil(math.log(dim /5.0, 2))) )]
  return indexes + powers

def filterVarRange(varRange):
  indexes = getSkippingIndexes(len(varRange))
  
  #get rid of indexes too big
  varRange = [varRange[i] for i in indexes]
  
  #discard list '[', ']'
  return str(varRange)[1:-1]

def writeData(modelFile, modelName, outFile, varRange, commonTime, nameSuffix,
    infoFields):
  outputName = fs.getParam(modelFile, "output")
  if os.path.isfile(modelFile):
    infoFile = "infos/%s.txt" %outputName
    order = fs.getParam(modelFile, "fixed orders")
    step = fs.getParam(modelFile, "fixed steps")
    time = fs.getParam(modelFile, "time")
    dim = getDimension(modelFile)
    varNames = getVars(modelFile)
  else:
    infoFile = order = time = step = dim = '-'
  
  if os.path.isfile(infoFile):
    reason = fs.getParam(infoFile, "reason:")
    intTime = fs.getParam(infoFile, "integration time:")
    compTime = fs.getParam(infoFile, "computation time:")
    swTime = fs.getParam(infoFile, "shrink wrapping time:")
    swCount = fs.getParam(infoFile, "shrink wraps:")
  else:
    reason = intTime = compTime = swTime = swCount = '-'
  
  outFile.write("  <tr align=\"center\">\n")
  outFile.write("    <td class='modelCell'><a href='%s'>%s</a></td>\n"
      %(modelFile,modelName))
  outFile.write("    <td class='dataCell'>%s</td> \n"%dim)
  outFile.write("    <td class='dataCell'>%s</td>\n"%order)
  outFile.write("    <td class='dataCell'>%s</td>\n"%step)
  outFile.write("    <td class='dataCell'>%s</td>\n"%time)
  
  for (a,_) in infoFields:
    if os.path.isfile(infoFile):
      value = fs.getParam(infoFile, "%s:"%a)
      if value == None:
        value = '-'
      #shorten reason
      if "remainder estimation" in value:
        value = "rem est"
      if "max increase" in value:
        value = "max increase"
      outFile.write("    <td class='dataCell'>%s</td>\n"%value)
    else:
      outFile.write("    <td>--</td>\n")
    
  outFile.write("    <td class='dataCell'>%s</td>\n"%filterVarRange(varRange))
  outFile.write("    <td align=\"left\">\n");
  outFile.write("    <table><tr>\n")
  
  variables = map(lambda x: x + 1, getSkippingIndexes(dim))
  for i in variables:
    imageName = "images/%s%s_%s_t_%s.png" %(outputName, nameSuffix, i, \
        commonTime)
    varName = varNames[i - 1]
    outFile.write("      <td>\n")
    if os.path.isfile(imageName):
      outFile.write("        <div>\n")
      outFile.write("          <div align='center'>%s</div>\n"%varName)
#      outFile.write("          <div align='center'>x%s</div>\n"%i)
      outFile.write("          " + 
          "<a href='%s'><img src='%s' style='width:200px;height:150px;'></a>\n"\
          %(imageName,imageName))
      outFile.write("        </div>\n")
    outFile.write("      </td>\n")
  outFile.write("    </tr></table>\n")
  outFile.write("    </td>\n")
  outFile.write("  </tr>\n")
  
  
  #print "%s \t& %s & %s" %(outputName, time, varRange[0])
  compMap = {"fcomp":"FC", "comp": "LC", "nocomp": "NC"}
#  print outputName.split("_")[-2]
  print "%s & %s & & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\" %(\
    outputName.split("_")[0],
    dim, \
    outputName.split("_")[-2], \
    compMap[outputName.split("_")[-1]], \
#    "---",\
    fs.getParam(infoFile, "computation time:"), \
    fs.getParam(infoFile, "pic_poly:"), \
    fs.getParam(infoFile, "pic_decr:"), \
    fs.getParam(infoFile, "pic_ref:"),\
    fs.getParam(infoFile, "tr_remap1:"),\
    fs.getParam(infoFile, "tr_precond:"),
    fs.getParam(infoFile, "tr_remap2:"),
    )
    
def write_table_rows(modelDir, pairs, outFile, nameSuffix, infoFields):
  for models in pairs:
    modelFiles = map(lambda m: os.path.join(modelDir, m), models)
    
    csvs = map(lambda m: \
        os.path.join("csvs", "%s.csv" %fs.getParam(m, "output") ), modelFiles)
        
    if plotToCommon:
      commonTime = min( map(lambda s: fs.find_file_max_time(s), csvs) )
    else:
      commonTime = max( map(lambda s: fs.find_file_max_time(s), csvs) )
    
    nocsvs = filter(lambda s: not os.path.isfile(s), csvs)
    
    if len(nocsvs) != 0:
      print "======== no csvs: %s" %nocsvs
    
    varRange2 = getVarRange2(modelFiles[0], csvs)
    
    for (i, model) in enumerate(models):
      #print modelFiles[i]
      #print models[i]
      #print varRange2[i]
      writeData(modelFiles[i], models[i], outFile, varRange2[i], commonTime, \
          nameSuffix, infoFields)
      
    
    outFile.write('<tr style="border-bottom:3px solid black">' + \
        '<td colspan="100%"><hr /></td></tr>\n')

def write_group_table_rows(modelDir, groups, outFile, nameSuffix, infoFields):
  for models in groups:
    modelFiles = map(lambda m: os.path.join(modelDir, m), models)
    
    csvs = map(lambda m: \
        os.path.join("csvs", "%s.csv" %fs.getParam(m, "output") ), modelFiles)
        
    if plotToCommon:
      commonTime = min( map(lambda s: fs.find_file_max_time(s), csvs) )
    else:
      commonTime = max( map(lambda s: fs.find_file_max_time(s), csvs) )
    
    nocsvs = filter(lambda s: not os.path.isfile(s), csvs)
    
    if len(nocsvs) != 0:
      print "======== no csvs: %s" %nocsvs
    
    varRange2 = getVarRange2(modelFiles[0], csvs)
    
    for (i, model) in enumerate(models):
      #print modelFiles[i]
      #print models[i]
      #print varRange2[i]
      writeData(modelFiles[i], models[i], outFile, varRange2[i], commonTime, \
          nameSuffix, infoFields)
      
    
    outFile.write('<tr style="border-bottom:3px solid black">' + \
        '<td colspan="100%"><hr /></td></tr>\n')

def write_table_start(outFile, infoFields):
  outFile.write("<html>\n")
  outFile.write("<head>\n")
  outFile.write("  <link rel='stylesheet' type='text/css' href='%s'>"%tableCSSFile)
  outFile.write("</head>\n")
  outFile.write("<body>\n")

  outFile.write("<table class='headerTable'>\n")
  outFile.write("  <tr>\n")
  outFile.write("    <th class='modelCell'>Model</th>\n")
  outFile.write("    <th class='dataCell'>Dim</th> \n")
  outFile.write("    <th class='dataCell'>Order</th>\n")
  outFile.write("    <th class='dataCell'>step</th>\n")
  outFile.write("    <th class='dataCell'>Time (goal)</th>\n")
  
  for (_,b) in infoFields:
    outFile.write("    <th class='dataCell'>%s</th>\n" %b)
  
  outFile.write("    <th class='dataCell'>W(var range)</th>\n")
  outFile.write("  </tr>\n")
  outFile.write("</table>\n")
  
  outFile.write("<table class='dataTable'>\n")


def write_table_end(outFile):
  outFile.write("</table>\n")
  outFile.write("</body>\n")
  outFile.write("</html>\n")
  outFile.close()


def write_table(tableName, modelDir, pairs, nameSuffix, infoFields):
  print "writing '%s'" %tableName
  outFile = open(tableName, 'w')
  write_table_start(outFile, infoFields)
  write_table_rows(modelDir, pairs, outFile, nameSuffix, infoFields)
  write_table_end(outFile)

def write_group_table(tableName, modelDir, groups, nameSuffix, infoFields):
  print "writing '%s'" %tableName
  outFile = open(tableName, 'w')
  write_table_start(outFile, infoFields)
  write_group_table_rows(modelDir, groups, outFile, nameSuffix, infoFields)
  write_table_end(outFile)
  
#generate a plot for a single model
def plot_variable(scriptsDir, modelFiles, var1, var2, nameSuffix):

  if var2 != 0:
    print "var2 is not time"
    sys.exit()
  csvs = map(lambda m: \
      os.path.join("csvs", "%s.csv" %fs.getParam(m, "output") ), modelFiles)

  if plotToCommon:
    time = min( map(lambda s: fs.find_file_max_time(s), csvs) )
  else:
    time = max( map(lambda s: fs.find_file_max_time(s), csvs) )
  
  tend = "tend=%s" %time
  
  #signal for no data
  if time == 100000000000000:
    return
  
  #print "time: %s" %time
  
  bounds = fs.get_range_bounds(time, csvs)
  
  #timeValues = map(lambda s: fs.get_range_up_to(time, s), csvs)
  #bounds = fs.get_bounds(timeValues)
  
  xRange = ["xMin=%s"%bounds[0][2*var2],"xMax=%s"%bounds[1][2*var2+1]]
  yRange = ["yMin=%s"%bounds[0][2*var1],"yMax=%s"%bounds[1][2*var1+1]]
  #print yRange
  
  imgSuffix = ["suffix=%s_%s_t_%s"%(nameSuffix, var1,time)]
  
  #print xRange
  #print yRange

  for csv in csvs:
    if not os.path.isfile(csv):
      continue
    args = [os.path.join(scriptsDir, 'gnuplot_var.py')] + \
        [csv] + ['var1=%s'%var1, 'var2=%s'%var2] + [tend] + xRange + yRange + imgSuffix
    p = subprocess.Popen(args)
    p.wait()  
    

#generate a plot for a single model
def group_plot_variable(scriptsDir, modelFiles, var1, var2, nameSuffix):
  if var2 != 0:
    print "var2 is not time"
    sys.exit()
  csvs = map(lambda m: \
      os.path.join("csvs", "%s.csv" %fs.getParam(m, "output") ), modelFiles)

  if plotToCommon:
    time = min( map(lambda s: fs.find_file_max_time(s), csvs) )
  else:
    time = max( map(lambda s: fs.find_file_max_time(s), csvs) )
  #for csv in csvs:
  #  print csv
  #  print fs.find_file_max_time(csv)
  
  tend = "tend=%s" %time
  
  #signal for no data
  #if time == 100000000000000:
  #  return
  
  #print "time: %s" %time
  
  bounds = fs.get_range_bounds(time, csvs)
  
  #timeValues = map(lambda s: fs.get_range_up_to(time, s), csvs)
  #bounds = fs.get_bounds(timeValues)
  
  xRange = ["xMin=%s"%bounds[0][2*var2],"xMax=%s"%bounds[1][2*var2+1]]
  yRange = ["yMin=%s"%bounds[0][2*var1],"yMax=%s"%bounds[1][2*var1+1]]
  #print xRange
  #print yRange
  
  imgSuffix = ["suffix=%s_%s_t_%s"%(nameSuffix, var1,time)]
  
  #print xRange
  #print yRange

  for csv in csvs:
    if not os.path.isfile(csv):
      continue
    args = [os.path.join(scriptsDir, 'gnuplot_var.py')] + \
        [csv] + ['var1=%s'%var1, 'var2=%s'%var2] + [tend] + xRange + yRange + imgSuffix
    p = subprocess.Popen(args)
    p.wait()      
  
#generate a plot for a single model
def plot_variable_old(scriptsDir, modelFiles, var1, var2, nameSuffix):
  csvs = map(lambda m: \
      os.path.join("csvs", "%s.csv" %fs.getParam(m, "output") ), modelFiles)
  
  time = min( map(lambda s: fs.find_file_max_time(s), csvs) )
  
  tend = "tend=%s" %time
  
  #signal for no data
  if time == 100000000000000:
    return
  
  #print "time: %s" %time
  
  timeValues = map(lambda s: fs.get_range_up_to(time, s), csvs)
  
  bounds = fs.get_bounds(timeValues)
  
  xRange = ["xMin=%s"%bounds[0][2*var2],"xMax=%s"%bounds[1][2*var2+1]]
  yRange = ["yMin=%s"%bounds[0][2*var1],"yMax=%s"%bounds[1][2*var1+1]]
  
  imgSuffix = ["suffix=%s_%s_t_%s"%(nameSuffix, var1,time)]
  
  #print xRange
  #print yRange

  for csv in csvs:
    if not os.path.isfile(csv):
      continue
    args = [os.path.join(scriptsDir, 'gnuplot_var.py')] + \
        [csv] + ['var1=%s'%var1, 'var2=%s'%var2] + [tend] + xRange + yRange + imgSuffix
    p = subprocess.Popen(args)
    p.wait()

#plotting of comparision plots

def generate_comparision_plots(scriptsDir, modelDir, pairs, nameSuffix):
  for models in pairs:
    print models
    modelFiles = map(lambda m: os.path.join(modelDir, m), models)
    
    #skip inflated model
    if len(modelFiles) == 3:
      modelFiles = modelFiles[:-1]
    dim = getDimension(modelFiles[0])
    
    for i in map(lambda x: x + 1, getSkippingIndexes(dim)):
      plot_variable(scriptsDir, modelFiles, i, 0, nameSuffix) #0 is the time


# generate all the plots in the group
def generate_group_comparision_plots(scriptsDir, modelDir, groups, nameSuffix):
  for models in groups:
    print models
    modelFiles = map(lambda m: os.path.join(modelDir, m), models)
    
    #skip inflated model
    #if len(modelFiles) == 3:
    #  modelFiles = modelFiles[:-1]
    dim = getDimension(modelFiles[0])
    
    #print dim
    
    for i in map(lambda x: x + 1, getSkippingIndexes(dim)):
      group_plot_variable(scriptsDir, modelFiles, i, 0, nameSuffix) #0 is the time

def generateHtml(scriptsDir, modelDir, modelPairs, suffix):
  oldFields = [
      ("integration time", "Time (actual)"), 
      ("computation time", "Comp time"), 
      ("shrink wrapping time", "SW time"), 
      ("shrink wraps", "SW count"),
      ("reason", "Stop reason"), 
  ]
  generateHtml(scriptsDir, modelDir, modelPairs, suffix, oldFields)

def generateHtml(scriptsDir, modelDir, modelPairs, suffix, infoFields):
  generate_comparision_plots(scriptsDir, modelDir, modelPairs, suffix)
  write_table("experiments%s.html" %suffix, modelDir, modelPairs, \
      suffix, infoFields)


def generateGroupHtml(scriptsDir, modelDir, modelGroups, suffix, infoFields):
  generate_group_comparision_plots(scriptsDir, modelDir, modelGroups, suffix)
  write_group_table("experiments%s.html" %suffix, modelDir, modelGroups, \
      suffix, infoFields)




