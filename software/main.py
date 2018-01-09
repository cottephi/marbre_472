#Fait par Philippe Cotte 
#Decembre 2017.
#CEA Saclay / Irfu / DPhP

from glob import glob
import sys
import os
import getopt
from marble_fit import marble
from sortedcontainers import SortedDict
import matplotlib.pyplot as plt

from my_analysis import my_analysis
from my_analysis import plot_holes
import pandas

def usage():
  print("Usage is:")
  print("-h : print this message")
  print("-i path/to/file : input data file")
  print("-c : optional config file, where one can specify x or z cuts")
  
  
def load_cuts(my_cut_file):
  if not os.path.isfile(my_cut_file):
    print("ERROR: cut file " + my_cut_file + " not found.")
    exit(1)
  cutfile = open(my_cut_file,"r")
  cuts = cutfile.readlines()
  cutfile.close()
  cutx = []
  cutz = []
  l_marble_file = []
  l_calle_file = []
  for cut in cuts:
    cut = cut.split("\n")[0]
    if cut[0] == "#":
      continue
    if len(cut.split(" ")) != 3 and cut.split(" ")[0] != "MARBLE" and cut.split(" ")[0] != "CALLE":
      print("load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut).")
      exit(1)
    if cut.split(" ")[0] == "x":
      cutx.append([cut.split(" ")[1],cut.split(" ")[2]])
    elif cut.split(" ")[0] == "z":
      cutz.append([cut.split(" ")[1],cut.split(" ")[2]])
    elif cut.split(" ")[0] == "MARBLE":
      if glob(cut.split(" ")[1]):
        l_marble_file = sorted(glob(cut.split(" ")[1]))
        print("   Found marble files :" , l_marble_file)
      else:
        print("ERROR : could not find Marble reference file")
        exit(1)
    elif cut.split(" ")[0] == "CALLE":
      if glob(cut.split(" ")[1]):
        marble_ref = cut.split(" ")[1].replace("calle*","marbre_1000Hz.csv")
        if os.path.isfile(marble_ref):
          l_calle_file.append(marble_ref)
        l_calle_file = l_calle_file + sorted(glob(cut.split(" ")[1]))
        print("   Found calle files :", l_calle_file)
      else:
        print("ERROR : could not find calle file")
        exit(1)
    else:
      print("Load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut), or by MARBLE or CALLE.")
      exit(1)
  return [[cutx,cutz], l_marble_file, l_calle_file]
      
      
def sort_data(lines, cut_file_arg, ID = 0):
  cutdata = [[[],[]]]
  rawdata = [[],[]]
  ignore = True
  do_cuts = False
  cutfile_content = []
  cuts = []
  l_marble_file = ""
  l_calle_file = ""
  if cut_file_arg != "":
    print("  Loading cutfile " + cut_file_arg + "...")
    cutfile_content = load_cuts(cut_file_arg)
    do_cuts = True
  ################for line in lines:#####################
    cuts = cutfile_content[0]
    l_marble_file = cutfile_content[1]
    l_calle_file = cutfile_content[2]
  for line in lines:
    x = float(line.split("\n")[0].split(";")[4])
    z = float(line.split("\n")[0].split(";")[0])
    rawdata[0].append(x)
    rawdata[1].append(z)
    addto = 0
    ################if do_cuts#####################
    if do_cuts:
      for i in range(0,len(cuts[0])):    #for each x-z, loop over the x-cuts to see if the x-z pair must be ignored 
        if (cuts[0][i][0] == "-inf" and x < float(cuts[0][i][1])) \
        or (x > float(cuts[0][i][0]) and cuts[0][i][1] == "inf") \
        or (x > float(cuts[0][i][0]) and cuts[0][i][1] == "+inf") \
        or (x > float(cuts[0][i][0]) and x < float(cuts[0][i][1])):
          addto = i
          ignore = False
      if ignore:
        continue
        
      for i in range(0,len(cuts[1])):    #for each x-z, loop over the z-cuts to see if the x-z pair must be ignored
        if ( cuts[1][i][0] == "-inf" and z < float(cuts[1][i][1]) ) \
        or ( cuts[1][i][0] == cuts[1][i][1] and z == float(cuts[1][i][0]) ) \
        or ( z > float(cuts[1][i][0]) and (cuts[1][i][1] == "inf" or cuts[1][i][1] == "+inf") )\
        or ( z > float(cuts[1][i][0]) and z < float(cuts[1][i][1]) ):    
          ignore = True
      
      if not ignore:
        while len(cutdata) <= addto:
          cutdata.append([[],[]])
        cutdata[addto][0].append(x)
        cutdata[addto][1].append(z)
    ################if do_cuts#####################
  ################for line in lines:#####################
  if not do_cuts:
    cutdata.append(rawdata)
  
  return [cutdata, rawdata, l_marble_file, l_calle_file]
  
def main(argv):

  #Booleans to ID the arguments p
  opt_cut_file = False
  cut_file_arg = ""
  opt_data = False
  data_arg = ""
  opt_marble_fit_file = False
  marble_fit_file_arg = ""
  opt_separate_in_holes = False
  
  #Doit etre Vrai si des mesures du marbres ont ete prises avant ET apres la mesure, ou si que le marbre a ete mesure.
  #Set to True if you want to fit the marble
  do_marble_fit = True
  #Set to true if your first and last groups of measurements are marble that are not cut out of the data by the cut you chose
  l_marble_file = []
  l_calle_file = []
  
  #Read the arguments
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:c:m:s", [])
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
      opt_data = True
      data_arg = arg
    elif opt == "-m":
      opt_marble_fit_file = True
      marble_fit_file_arg = arg
    elif opt == "-s":
      opt_separate_in_holes = True
    else:
      print("Unrecognized arguments")
      usage()
  
  if not opt_data:
    print("ERROR: need a data input.")
  if not os.path.isfile(data_arg) and not os.path.isdir(data_arg):
    print("ERROR: " + data_arg + " not found.")
    exit(0)
  if opt_cut_file:
    if not os.path.isfile(cut_file_arg):
      print("ERROR: cut file " + cut_file_arg + " not found.")
      exit(0)
  if opt_marble_fit_file:
    if not os.path.isfile(marble_fit_file_arg):
      print("ERROR: marble fit file " + marble_fit_file_arg + " not found.")
      exit(0)
      

  l_datafiles = []
  if os.path.isfile(data_arg):
    l_datafiles.append(data_arg)
  elif os.path.isdir(data_arg):
    l_datafiles = sorted(glob(data_arg + '/*merged*'))
  row = len(l_datafiles)
  l_l_cutdata = []
  l_l_raw_data = []
  l_l_glued_cutdata = [[[],[]]]
  col = 0
  
  for i in range(0,len(l_datafiles)):
    print("Opening datafile ",l_datafiles[i],"...")
    datafiles = open(l_datafiles[i], "r")
    lines = datafiles.readlines()
    datafiles.close()
    if "Distance" in lines[0]:
      lines = lines[1:]
    print(" Sorting data...")
    tmp_l_l_cut_data, tmp_l_raw_data, l_marble_file, l_calle_file = sort_data(lines, cut_file_arg)
    print(" ...done")
    l_l_raw_data.append(tmp_l_raw_data)
    if opt_marble_fit_file:
      print(" Found specified marble file " + marble_fit_file_arg)
      l_marble_file = [marble_fit_file_arg]
    if do_marble_fit and (opt_cut_file or opt_marble_fit_file):
      print(" Correcting row " + str(i+1) + " with marble...")
      if i == 0:
        tmp_l_l_cut_data = marble(tmp_l_l_cut_data, l_marble_file, i+1, l_calle_file)
      else:
        tmp_l_l_cut_data = marble(tmp_l_l_cut_data, l_marble_file, i+1, [])
      print(" ...done")
      print("")
    l_l_cutdata = l_l_cutdata + tmp_l_l_cut_data
    if i == 0:
      col = len(l_l_cutdata)
      
      
  print("Data has ",row," rows and ",col," columns")
  print("Plotting raw data...")
  plot_holes(l_l_raw_data, row, col, "raw_")
  print("...done")
  if opt_separate_in_holes:
    print("Analysing data...")
    my_analysis(l_l_cutdata, row, col)
    print("done")
  else:
    for l_data in l_l_cutdata:
      l_l_glued_cutdata[0][0] = l_l_cutdata[0][0] + l_data[0]
      l_l_glued_cutdata[0][1] = l_l_cutdata[0][1] + l_data[1]
    print("Analysing non-separated data...")
    my_analysis(l_l_glued_cutdata,1,1)
    print("...done")
  
  #plt.show()
  #plt.clf()
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
