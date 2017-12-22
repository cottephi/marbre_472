#Fait par Philippe Cotte 
#Décembre 2017.
#CEA Saclay / Irfu / DPhP

from glob import glob
import sys
import os
import getopt
from marble_fit import marble
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
      
      
def apply_cuts(rawdata, cuts):
  cutdata = []
  cutdata.append(SortedDict())
  new_set = False
  ignore = False
  for x in rawdata:    #loop over data
    for cut in cuts[0]:    #for each x-z, loop over the x-cuts to see if the x-z pair must be ignored 
      if cut[0] == cut[1] and float(x) == float(cut[0]):    #if the cut is a single point and not a range and x or y matches a single cut, ignore the value but does not create a new data set
        ignore = True
        break
      elif cut[0] == "-inf" and float(x) < float(cut[1]):
        new_set = True
        ignore = True
        break
      elif float(cut[1]) == "inf" and float(x) > float(cut[0]):
        new_set = True
        ignore = True
        break
      elif float(x) >= float(cut[0]) and float(x) <= float(cut[1]):
        new_set = True
        ignore = True
        break
      else:
        ignore = False
        
    if ignore:
      continue
      
    for cut in cuts[1]:    #for each x-z, loop over the z-cuts to see if the x-z pair must be ignored 
      if cut[0] == cut[1] and float(rawdata[x]) == float(cut[0]):    #if the cut is a single point and not a range and x or y matches a single cut, ignore the value but does not create a new data set
        ignore = True
        break
      elif cut[0] == "-inf" and float(rawdata[x]) < float(cut[1]):
        new_set = True
        ignore = True
        break
      elif cut[1] == "inf" and float(rawdata[x]) > float(cut[0]):
        new_set = True
        ignore = True
        break
      elif float(rawdata[x]) >= float(cut[0]) and float(rawdata[x]) <= float(cut[1]):
        new_set = True
        ignore = True
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
  opt_data = False
  data_arg = ""
  
  #Doit etre Vrai si des mesures du marbres ont ete prises avant ET apres la mesure, ou si que le marbre a ete mesure.
  do_marble_fit = True
  
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
      opt_data = True
      data_arg = arg
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
      print("ERROR: cut file " + data_arg + " not found.")
      exit(0)

  l_data = []
  if os.path.isfile(data_arg):
    l_data.append(data_arg)
  elif os.path.isdir(data_arg):
    l_data = sorted(glob(data_arg + '/*merged*'))

  l_sd_data = []
  for i in range(0,len(l_data)):
    tmp_l_sd_data = []
    data = open(l_data[i], "r")
    lines = data.readlines()
    data.close()
    if "Distance" in lines[0]:
      lines = lines[1:]
    sd_raw_data = SortedDict() #SortedDict va automatiquement trier les donnees par ordre croissant en x. Donc peut importe dans quelle ordre elles ont ete prise, tout sera comme il faut.
    for line in lines:
      sd_raw_data[float(line.split("\n")[0].split(";")[4])]=float(line.split("\n")[0].split(";")[0])

    if opt_cut_file:
      tmp_l_sd_data = apply_cuts(sd_raw_data, load_cuts(cut_file_arg))
    else:
      tmp_l_sd_data.append(sd_raw_data)
    if do_marble_fit:
      tmp_l_sd_data = marble(tmp_l_sd_data, i+1)
    l_sd_data = l_sd_data + tmp_l_sd_data
  
  my_analysis(l_sd_data, 5)
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
