#Fait par Philippe Cotte 
#Décembre 2017.
#CEA Saclay / Irfu / DPhP


import sys
import os
import getopt
from sortedcontainers import SortedDict

from my_analysis import my_analysis

def usage():
  print("Usage is:")
  print("-h : print this message")
  print("-i path/to/file : input data file")
  print("-c : optional config file, where one can specify x or z cuts")
  
  
def load_cuts(my_cut_file):
  if not os.path.isfile(my_cut_file):
    print("ERROR: cut file " + my_cut_file + " not found.")
    exit(0)
  cutfile = open(my_cut_file,"r")
  cuts = cutfile.readlines()
  cutfile.close()
  cutx = []
  cutz = []
  for cut in cuts:
    cut = cut.split("\n")[0]
    if len(cut.split(" ")) != 3:
      print("load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut).")
      exit(0)
    if cut.split(" ")[0] == "x":
      cutx.append([cut.split(" ")[1],cut.split(" ")[2]])
    elif cut.split(" ")[0] == "z":
      cutz.append([cut.split(" ")[1],cut.split(" ")[2]])
    else:
      print("2 load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut).")
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
    print(str(err))
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
      print("Unrecognized arguments")
      usage()
  
  if not opt_data_file:
    print("ERROR: need a data file.")
  if not os.path.isfile(data_file_arg):
    print("ERROR: data file " + data_file_arg + " file not found.")
    exit(0)
  if opt_cut_file:
    if not os.path.isfile(cut_file_arg):
      print("ERROR: cut file " + data_file_arg + " not found.")
      exit(0)

  data_file = open(data_file_arg, "r")
  lines = data_file.readlines()[1:]
  data_file.close()
  raw_data = SortedDict()
  for line in lines:
    raw_data[float(line.split("\n")[0].split(";")[4])]=float(line.split("\n")[0].split(";")[0])
  
  cut_data = []
  if opt_cut_file:
    raw_data_cuts = raw_data
    cut_data = apply_cuts(raw_data, raw_data_cuts, load_cuts(cut_file_arg))
  else:
    cut_data.append(raw_data)
  
  my_analysis(cut_data)
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
