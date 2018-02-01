#Fait par Philippe Cotte 
#Decembre 2017.
#CEA Saclay / Irfu / DPhP

from glob import glob
import sys
import os
import getopt
from marble_fit import marble
from marble_fit import plot_other_marble_file
from sortedcontainers import SortedDict
import matplotlib.pyplot as plt

from my_analysis import my_analysis
from my_analysis import plot_holes
import pandas

def usage():
  print("Usage is:")
  print("-h : print this message")
  print("-i path/to/file : input data file")
  print("-c : cut file, where one can specify x or z cuts")
  print("-a : cali files (marbles and, possibly, cales)")
  print("-m : to specify a single marble file to be used as reference")
  print("-s : to apply the separation from the cut file")
  print("-t : only plot raw data and cali data. Usefull to quickly check the aspect of the measurements")
  
  
def load_cali(data_arg, opt_marble_fit_file = False, marble_fit_file_arg = ""):
  l_marble_file = []
  l_cale_file = []
  l_other_marble_file = []
    
  if glob(data_arg + "/cale/*"):
    l_cale_file = glob(data_arg + "/cale/*")
    print("   Found cale files :" , l_cale_file)
  elif glob(os.path.dirname(data_arg) + "/cale/*"):
    l_cale_file = glob(os.path.dirname(data_arg) + "/cale/*")
    print("   Found cale files :" , l_cale_file)
  elif glob(os.path.dirname(data_arg) + "/../cale/*"):
    l_cale_file = glob(os.path.dirname(data_arg) + "/../cale/*")
    print("   Found cale files :" , l_cale_file)
    
  if glob(data_arg + "/other_files/*"):
    l_other_marble_file = glob(data_arg + "/other_files/*")
    print("   Found other files :" , l_other_marble_file)
  elif glob(os.path.dirname(data_arg) + "/other_files/*"):
    l_other_marble_file = glob(os.path.dirname(data_arg) + "/other_files/*")
    print("   Found other files :" , l_other_marble_file)
  elif glob(os.path.dirname(data_arg) + "/../other_files/*"):
    l_other_marble_file = glob(os.path.dirname(data_arg) + "/../other_files/*")
    print("   Found other files :" , l_other_marble_file)
    
  if opt_marble_fit_file:
    print(" Found specified marble file " + marble_fit_file_arg)
    l_marble_file = [marble_fit_file_arg]
  else:
    if glob(data_arg + "/cali/*merged*"):
      l_marble_file = sorted(glob(data_arg + "/cali/*merged*"))
      print("   Found marble files :" , l_marble_file)
    elif glob(os.path.dirname(data_arg) + "/cali/*merged*"):
      l_marble_file = sorted(glob(os.path.dirname(data_arg) + "/cali/*merged*"))
      print("   Found marble files :" , l_marble_file)
    elif glob(os.path.dirname(data_arg) + "/../cali/*merged*"):
      l_marble_file = sorted(glob(os.path.dirname(data_arg) + "/../cali/*merged*"))
      print("   Found marble files :" , l_marble_file)
    else:
      print("**WARNING: No califiles found. Ignoring marble calibration.")
    
  return[l_marble_file, l_cale_file, l_other_marble_file]
  
def load_cuts(my_cut_file):
  cutfile = open(my_cut_file,"r")
  cuts = cutfile.readlines()
  cutfile.close()
  cutx = []
  cutz = []
  for cut in cuts:
    cut = cut.split("\n")[0]
    if cut[0] == "#":
      continue
    if len(cut.split(" ")) != 3 and len(cut.split(" ")) != 11:
      print("load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut).")
      exit(1)
    if cut.split(" ")[0] == "x":
      cutx.append([])
      for i in range(1,len(cut.split(" "))):
        cutx[-1].append(cut.split(" ")[i])
    elif cut.split(" ")[0] == "z":
      cutz.append([cut.split(" ")[1],cut.split(" ")[2]])
    else:
      print("Load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut) for each measurement row.")
      exit(1)
  return [cutx, cutz]
      
      
def sort_data(lines, data_arg, row = 0):
  cutdata = [[[],[]]]
  rawdata = [[],[]]
  ignore = True
  do_cuts = False
  cutfile_content = []
  cutx = []
  cutz = []
  cut_file_arg = ""
  if os.path.isfile(data_arg + "/cuts.txt"):
    cut_file_arg = data_arg + "/cuts.txt"
    print("   Found cuts file :" , cut_file_arg)
  elif os.path.isfile(os.path.dirname(data_arg) + "/cuts.txt"):
    cut_file_arg = os.path.dirname(data_arg) + "/cuts.txt"
    print("   Found cuts files :" , cut_file_arg)
  elif os.path.isfile(os.path.dirname(data_arg) + "/../cuts.txt"):
    cut_file_arg = os.path.dirname(data_arg) + "/../cuts.txt"
    print("   Found cuts files :" , cut_file_arg)
  else:
    print("**WARNING: No cuts file found. Ignoring cuts.")
  
  if cut_file_arg != "":
    print("  Loading cutfile " + cut_file_arg + "...")
    cutfile_content = load_cuts(cut_file_arg)
    do_cuts = True
  ################for line in lines:#####################
    cutx = cutfile_content[0]
    cutz = cutfile_content[1]
  for line in lines:
    x = float(line.split("\n")[0].split(";")[4])
    z = float(line.split("\n")[0].split(";")[0])
    rawdata[0].append(x)
    rawdata[1].append(z)
    addto = 0
    ################if do_cuts#####################
    if do_cuts:
      for i in range(0,len(cutx)):    #for each x-z, loop over the x-cuts to see if the x-z pair must be ignored 
        if (cutx[i][0 + 2*row] == "-inf" and x < float(cutx[i][1 + 2*row])) \
        or (x > float(cutx[i][0 + 2*row]) and cutx[i][1 + 2*row] == "inf") \
        or (x > float(cutx[i][0 + 2*row]) and cutx[i][1 + 2*row] == "+inf") \
        or (x > float(cutx[i][0 + 2*row]) and x < float(cutx[i][1 + 2*row])):
          addto = i
          ignore = False
      if ignore:
        continue
        
      for i in range(0,len(cutz)):    #for each x-z, loop over the z-cuts to see if the x-z pair must be ignored
        if ( cutz[i][0] == "-inf" and z < float(cutz[i][1]) ) \
        or ( cutz[i][0] == cutz[i][1] and z == float(cutz[i][0]) ) \
        or ( z > float(cutz[i][0]) and (cutz[i][1] == "inf" or cutz[i][1] == "+inf") )\
        or ( z > float(cutz[i][0]) and z < float(cutz[i][1]) ):    
          ignore = True
      
      if not ignore:
        ignore = True
        while len(cutdata) <= addto:
          cutdata.append([[],[]])
        cutdata[addto][0].append(x)
        cutdata[addto][1].append(z)
    ################if do_cuts#####################
  ################for line in lines:#####################
  if not do_cuts:
    cutdata.append(rawdata)
  cutdata = [ xz for xz in cutdata if xz != [[],[]] ]
  return [cutdata, rawdata]
  
def main(argv):

  #Booleans to ID the arguments p
  opt_cali_file = False
  opt_cut_file = False
  opt_data = False
  data_arg = ""
  opt_marble_fit_file = False
  marble_fit_file_arg = ""
  opt_separate_in_holes = False
  opt_is_test = False
  
  l_marble_file = []
  l_cale_file = []
  l_other_marble_file = []
  sigmaCopperLaser = 0
  
  #Read the arguments
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:cam:st", [])
  except getopt.GetoptError as err:
    print(str(err))
    usage()
    exit(0)
  for opt, arg in opts:
    if opt == "-h":
      usage()
      exit(0)
    elif opt == "-a":
      opt_cali_file = True
    elif opt == "-c":
      opt_cut_file = True
    elif opt == "-i":
      opt_data = True
      data_arg = arg
    elif opt == "-m":
      opt_marble_fit_file = True
      marble_fit_file_arg = arg
    elif opt == "-s":
      opt_separate_in_holes = True
    elif opt == "-t":
      opt_is_test = True
    else:
      print("Unrecognized arguments")
      usage()
      
  if not opt_data:
    print("ERROR: need a data input.")
  if not os.path.isfile(data_arg) and not os.path.isdir(data_arg):
    print("ERROR: " + data_arg + " not found.")
    exit(1)
  if opt_marble_fit_file:
    if not os.path.isfile(marble_fit_file_arg):
      print("ERROR: marble fit file " + marble_fit_file_arg + " not found.")
      exit(1)
      

  
  l_datafiles = []
  while data_arg[-1] == '/':
    data_arg = data_arg[:-1]
  outdirectory = data_arg
  if os.path.isfile(data_arg):
    l_datafiles.append(data_arg)
    outdirectory = outdirectory.replace(".csv","") + "_plots"
  elif os.path.isdir(data_arg):
    if not  os.path.isdir(data_arg + "/data"):
      print("ERROR: need data directory INSIDE the LEM directory",data_arg)
      exit(1)
    outdirectory = outdirectory + "/plots"
    l_datafiles = sorted(glob(data_arg + '/data/*merged*'))
  if not os.path.isdir(outdirectory):
    os.mkdir(outdirectory)
  outdirectory = outdirectory + "/"
  
  if opt_cali_file:
    l_marble_file, l_cale_file, l_other_marble_file = load_cali(data_arg, opt_marble_fit_file, marble_fit_file_arg)
  row = len(l_datafiles)
  l_l_cutdata = []
  l_l_raw_data = []
  l_l_glued_cutdata = [[[],[]]]
  col = 0
  sigmarble = 0
  
  for i in range(0,len(l_datafiles)):
    print("Opening datafile ",l_datafiles[i],"...")
    datafiles = open(l_datafiles[i], "r")
    title = os.path.basename(l_datafiles[i]).replace("merged_","").replace(".csv","")
    lines = [ line for line in datafiles.readlines() if not "Distance" in line ]
    datafiles.close()
    print(" Sorting data...")
    tmp_l_l_cut_data, tmp_l_raw_data = sort_data(lines, data_arg, i)
    print(" ...done")
    l_l_raw_data.append(tmp_l_raw_data)
    if l_marble_file != []:
      print(" Correcting row " + str(i+1) + " with marble...")
      if i == 0:
        tmp_l_l_cut_data, sigmarble = marble(tmp_l_l_cut_data, l_marble_file, i+1, l_cale_file, outdirectory)
      else:
        tmp_l_l_cut_data, sigmarble = marble(tmp_l_l_cut_data, l_marble_file, i+1, [], outdirectory)
      print(" ...done")
      print("")
    l_l_cutdata = l_l_cutdata + tmp_l_l_cut_data
    if i == 0:
      col = len(l_l_cutdata)
      
      
  if l_other_marble_file != []:
    sigmaCopperLaser = plot_other_marble_file(l_other_marble_file, outdirectory)
    
  print("Data has ",row," rows and ",col," columns")
  print("Plotting raw data...")
  plot_holes(l_l_raw_data, row, col, "raw_", outdirectory)#, ["raw " + ti for ti in title])
  if opt_is_test:
    plt.show()
    #plt.clf()
    #plt.close()
    exit(0)
  print("...done")
  if opt_separate_in_holes:
    print("Analysing data...")
    my_analysis(l_l_cutdata, row, col, sigmarble, sigmaCopperLaser, outdirectory)
    print("done")
  else:
    for l_data in l_l_cutdata:
      l_l_glued_cutdata[0][0] = l_l_cutdata[0][0] + l_data[0]
      l_l_glued_cutdata[0][1] = l_l_cutdata[0][1] + l_data[1]
    print("Analysing non-separated data...")
    my_analysis(l_l_glued_cutdata,1,1, sigmarble, sigmaCopperLaser, outdirectory)
    print("...done")
  
  #plt.show()
  #plt.clf()
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
