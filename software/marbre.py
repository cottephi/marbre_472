import glob
import sys
import os
import getopt
import shutil
import matplotlib.pyplot as plt
import collections
import numpy as np
import math
from sklearn import linear_model
from sortedcontainers import SortedDict
from numbers import Number

def usage():
  prfloat("Usage is:")
  prfloat("-h : print this message")
  prfloat("-i path/to/file : input data file")
  prfloat("-c : optional config file, where one can specify x or z cuts")
  
  
def load_cuts(my_cut_file):
  if not os.path.isfile(my_cut_file):
    prfloat("ERROR: cut file " + my_cut_file + " not found.")
    exit(0)
  cutfile = open(my_cut_file,"r")
  cuts = cutfile.readlines()
  cutfile.close()
  cutx = []
  cutz = []
  for cut in cuts:
    cut = cut.split("\n")[0]
    if len(cut.split(" ")) != 3:
      prfloat("load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut).")
      exit(0)
    if cut.split(" ")[0] == "x":
      cutx.append([cut.split(" ")[1],cut.split(" ")[2]])
    elif cut.split(" ")[0] == "z":
      cutz.append([cut.split(" ")[1],cut.split(" ")[2]])
    else:
      prfloat("2 load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut).")
      exit(0)
  return [cutx,cutz]
      
      
def apply_cuts(rawdata, rawdatacuts, cuts):
  cutdata = []
  cutdata.append(SortedDict())
  new_set = False
  ignore = False
  for x in rawdata:    #loop over data
    for cut in cuts[0]:    #for each x-z, loop over the x-cuts to see if the x-z pair must be ignored 
      if cut[0] == cut[1] and float(x) == float(cut[0]):    #if the cut is a single point and not a range and x or y matches a single cut, ignore the value but does not create a new data set
        rawdatacuts[x] = 0
        ignore = True
        break
      elif cut[0] == "-inf" and float(x) < float(cut[1]):
        new_set = True
        ignore = True
        rawdatacuts[x] = 0
        break
      elif float(cut[1]) == "inf" and float(x) > float(cut[0]):
        new_set = True
        ignore = True
        rawdatacuts[x] = 0
        break
      elif float(x) >= float(cut[0]) and float(x) <= float(cut[1]):
        new_set = True
        ignore = True
        rawdatacuts[x] = 0
        break
      else:
        ignore = False
        
    if ignore:
      continue
      
    for cut in cuts[1]:    #for each x-z, loop over the z-cuts to see if the x-z pair must be ignored 
      if cut[0] == cut[1] and float(rawdata[x]) == float(cut[0]):    #if the cut is a single point and not a range and x or y matches a single cut, ignore the value but does not create a new data set
        ignore = True
        rawdatacuts[x] = 0
        break
      elif cut[0] == "-inf" and float(rawdata[x]) < float(cut[1]):
        new_set = True
        ignore = True
        rawdatacuts[x] = 0
        break
      elif cut[1] == "inf" and float(rawdata[x]) > float(cut[0]):
        new_set = True
        ignore = True
        rawdatacuts[x] = 0
        break
      elif float(rawdata[x]) >= float(cut[0]) and float(rawdata[x]) <= float(cut[1]):
        new_set = True
        ignore = True
        rawdatacuts[x] = 0
        break
      else:
        ignore = False
          
    if ignore:
      continue
      
    else:
      if new_set and len(cutdata[0]) != 0:    #if previous data matched a cut ranged, create a new data range 
        cutdata.append(SortedDict())
      cutdata[-1][x] = rawdata[x]
      new_set = False
      
  return cutdata
  
def main(argv):

  #Booleans to ID the arguments
  opt_cut_file = False
  cut_file_arg = ""
  opt_data_file = False
  data_file_arg = ""
  
  #Read the arguments
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:c:", [])
  except getopt.GetoptError as err:
    prfloat(str(err))
    usage()
    exit(0)
  for opt, arg in opts:
    if opt == "-h":
      usage()
      exit(1)
    elif opt == "-c":
      opt_cut_file = True
      cut_file_arg = arg
    elif opt == "-i":
      opt_data_file = True
      data_file_arg = arg
    else:
      prfloat("Unrecognized arguments")
      usage()
  
  if not opt_data_file:
    prfloat("ERROR: need a data file.")
  if not os.path.isfile(data_file_arg):
    prfloat("ERROR: data file " + data_file_arg + " file not found.")
    exit(0)
  if opt_cut_file:
    if not os.path.isfile(cut_file_arg):
      prfloat("ERROR: cut file " + data_file_arg + " not found.")
      exit(0)

  data_file = open(data_file_arg, "r")
  lines = data_file.readlines()[1:]
  data_file.close()
  raw_data = SortedDict()
  for line in lines:
    raw_data[float(line.split("\n")[0].split(";")[4])]=float(line.split("\n")[0].split(";")[0])
  
  #for x in sorted(raw_data.iterkeys()):
    #prfloat(x + " " + raw_data[x])
  cut_data = []
  raw_data_cuts = raw_data
  if opt_cut_file:
    cut_data = apply_cuts(raw_data, raw_data_cuts, load_cuts(cut_file_arg))
    
  marble_x = np.array(0)
  marble_z = np.array(0)
  marble_x, marble_z = zip(*sorted(cut_data[0].items()))
  marble_data = SortedDict()
  for x in cut_data[0]:
    if not abs(np.mean(marble_z)-cut_data[0][x]) > 4*math.sqrt(np.var(marble_z)):
      marble_data[x] = cut_data[0][x]
  marble_x, marble_z = zip(*sorted(cut_data[-1].items()))
  for x in cut_data[-1]:
    if not abs(np.mean(marble_z)-cut_data[-1][x]) > 4*math.sqrt(np.var(marble_z)):
      marble_data[x] = cut_data[-1][x]
  
  marble_x, marble_z = zip(*sorted(marble_data.items()))
  marble_x = np.array(marble_x)
  marble_z = np.array(marble_z)
  marble_regr = linear_model.LinearRegression()
  marble_regr.fit(marble_x[:,np.newaxis], marble_z)
  x_test = np.linspace(np.min(marble_x), np.max(marble_x), 100)
  plt.figure(0)
  plt.title("Marble Height vs length")
  plt.xlabel('x(micrometer)', fontsize=18)
  plt.ylabel('z(micrometer)', fontsize=16)
  plt.scatter(marble_x, marble_z, color='black')
  plt.plot(x_test, marble_regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=3)
  plt.savefig("marble.pdf")
  
  holes_data = [SortedDict()]
  for hole in cut_data[1:-1]:
    for x in hole:
      holes_data[-1][x] = hole[x]-marble_regr.predict(hole[x])[0]
    holes_data.append(SortedDict())
  holes_data = holes_data[:-1]
  plt.figure(1)
  plt.title("Holes Height vs length")
  plt.xlabel('x(micrometer)', fontsize=18)
  plt.ylabel('z(micrometer)', fontsize=16)
  for hole in holes_data:
    print(len(hole))
    x, z = zip(*sorted(hole.items()))
    plt.scatter(x, z, color='black')
  plt.show()
  plt.savefig("holes.pdf")
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
