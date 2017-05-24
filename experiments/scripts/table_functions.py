import sys
import math
import os
import re
import subprocess

import my_functions as fs

plotToCommon = False
#plotToCommon = True

def getDimension(modelFile):
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
  return len(s.split(','))

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
        data = line.split(',')
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


def writeData(modelFile, modelName, outFile, varRange, commonTime, nameSuffix):
  outputName = fs.getParam(modelFile, "output")
  if os.path.isfile(modelFile):
    infoFile = "infos/%s.txt" %outputName
    order = fs.getParam(modelFile, "fixed orders")
    step = fs.getParam(modelFile, "fixed steps")
    time = fs.getParam(modelFile, "time")
    dim = getDimension(modelFile)
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
  
  outFile.write("  <tr>\n")
  outFile.write("    <td><a href='%s'>%s</a></td>\n"%(modelFile,modelName))
  outFile.write("    <td>%s</td> \n"%dim)
  outFile.write("    <td>%s</td>\n"%order)
  outFile.write("    <td>%s</td>\n"%step)
  outFile.write("    <td>%s</td>\n"%time)
  outFile.write("    <td>%s</td>\n"%intTime)
  outFile.write("    <td>%s</td>\n"%compTime)
  outFile.write("    <td>%s</td>\n"%swTime)
  outFile.write("    <td>%s</td>\n"%swCount)
  outFile.write("    <td>%s</td>\n"%reason)
  outFile.write("    <td>%s</td>\n"%str(varRange)[1:-1])
  outFile.write("    <td>\n");
  outFile.write("    <table><tr>\n")
  for i in range(1, dim+1):
    imageName = "images/%s%s_%s_t_%s.png" %(outputName, nameSuffix, i, commonTime)
    outFile.write("      <td>\n")
    if os.path.isfile(imageName):
      outFile.write("        <div>\n")
      outFile.write("          <div align='center'>%s</div>\n"%varRange[i-1])
      outFile.write("          <a href='%s'><img src='%s' style='width:200px;height:150px;'></a>\n" %(imageName,imageName))
      outFile.write("        </div>\n")
    outFile.write("      </td>\n")
  outFile.write("    </tr></table>\n")
  outFile.write("    </td>\n")
  outFile.write("  </tr>\n")

def write_table_rows(modelDir, pairs, outFile, nameSuffix):
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
          nameSuffix)
      
    
    outFile.write('<tr style="border-bottom:3px solid black">' + \
        '<td colspan="100%"><hr /></td></tr>\n')


def write_table_start(outFile):
  outFile.write("<html>\n")
  outFile.write("<body>\n")

  outFile.write("<table >\n")
  outFile.write("  <tr>\n")
  outFile.write("    <th>Model</th>\n")
  outFile.write("    <th>Dim</th> \n")
  outFile.write("    <th>Order</th>\n")
  outFile.write("    <th>step</th>\n")
  outFile.write("    <th>Time (goal)</th>\n")
  outFile.write("    <th>Time (actual)</th>\n")
  outFile.write("    <th>Computation time</th>\n")
  outFile.write("    <th>Shrink wrapping time</th>\n")
  outFile.write("    <th>Shrink wrapping count</th>\n")
  outFile.write("    <th>Stop reason</th>\n")
  outFile.write("    <th>Size of the variables range</th>\n")
  outFile.write("    <th>Plots</th>\n")
  outFile.write("  </tr>\n")


def write_table_end(outFile):
  outFile.write("</table>\n")
  outFile.write("</body>\n")
  outFile.write("</html>\n")
  outFile.close()


def write_table(tableName, modelDir, pairs, nameSuffix):
  outFile = open(tableName, 'w')
  write_table_start(outFile)
  write_table_rows(modelDir, pairs, outFile, nameSuffix)
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
    modelFiles = map(lambda m: os.path.join(modelDir, m), models)
    
    #skip inflated model
    if len(modelFiles) == 3:
      modelFiles = modelFiles[:-1]
    print modelFiles
    dim = getDimension(modelFiles[0])
    for i in range(1, dim+1):
      plot_variable(scriptsDir, modelFiles, i, 0, nameSuffix) #0 is the time











