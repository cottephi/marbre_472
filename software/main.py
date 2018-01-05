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
      else:
        print("ERROR : could not find Marble reference file")
        exit(1)
    elif cut.split(" ")[0] == "CALLE":
      if glob(cut.split(" ")[1]):
        marble_ref = cut.split(" ")[1].replace("calle*","marbre_1000Hz.csv")
        if os.path.isfile(marble_ref):
          l_calle_file.append(marble_ref)
        l_calle_file = l_calle_file + sorted(glob(cut.split(" ")[1]))
      else:
        print("ERROR : could not find calle file")
        exit(1)
    elif cut[0] == "#":
      continue
    else:
      print("Load_cuts ERROR : bad cut format " + cut + ". Must start by x or z, followed by two values (range of the cut), or by MARBLE or CALLE.")
      exit(1)
  return [[cutx,cutz], l_marble_file, l_calle_file]
      
      
def apply_cuts(rawdata, cuts, ID = 0):
  cutdata = []
  cutdata.append(SortedDict())
  new_set = False
  ignore = False
  l_marble_file = cuts[1]
  l_calle_file = cuts[2]
  cuts = cuts[0]
  i = 0
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
      
  return [cutdata,l_marble_file, l_calle_file]
  
def main(argv):

  #Booleans to ID the arguments p
  opt_cut_file = False
  cut_file_arg = ""
  opt_data = False
  data_arg = ""
  opt_marble_fit_file = False
  marble_fit_file_arg = ""
  
  #Doit etre Vrai si des mesures du marbres ont ete prises avant ET apres la mesure, ou si que le marbre a ete mesure.
  #Set to True if you want to fit the marble
  do_marble_fit = True
  #Set to true if your first and last groups of measurements are marble that are not cut out of the data by the cut you chose
  edges_are_marble = True
  l_marble_file = []
  l_calle_file = []
  
  #Read the arguments
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:c:m:", [])
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
      

  l_data = []
  if os.path.isfile(data_arg):
    l_data.append(data_arg)
  elif os.path.isdir(data_arg):
    l_data = sorted(glob(data_arg + '/*merged*'))
  row = len(l_data)
  l_sd_data = []
  l_sd_raw_data = []
  col = 0
  
  for i in range(0,len(l_data)):
    tmp_l_sd_data = []
    data = open(l_data[i], "r")
    lines = data.readlines()
    data.close()
    if "Distance" in lines[0]:
      lines = lines[1:]
    sd_raw_data = SortedDict() #SortedDict va automatiquement trier les donnees par ordre croissant en x. Donc peut importe dans quelle ordre elles ont ete prise, tout sera comme il faut.
    for line in lines:
      conti = False
      x = float(line.split("\n")[0].split(";")[4])
      x_tmp = x
      z = float(line.split("\n")[0].split(";")[0])
      while x_tmp in sd_raw_data: #evite d'avoir deux fois la meme valeur de x, ce qui ne fonctionne pas pour les dictionnaires
        x_tmp = x_tmp+0.01
        if x_tmp > x + 1:
          conti = True
          break
      if conti:
        continue
      sd_raw_data[x_tmp] = z
    l_sd_raw_data.append(sd_raw_data)
    if opt_cut_file:
      tmp_l_sd_data,l_marble_file, l_calle_file = apply_cuts(sd_raw_data, load_cuts(cut_file_arg), i)
    else:
      tmp_l_sd_data.append(sd_raw_data)
    if opt_marble_fit_file:
      l_marble_file = [marble_fit_file_arg]
    if do_marble_fit and opt_cut_file:
      if i == 0:
        tmp_l_sd_data = marble(tmp_l_sd_data, l_marble_file, i+1, l_calle_file, edges_are_marble)
      else:
        tmp_l_sd_data = marble(tmp_l_sd_data, l_marble_file, i+1, [], edges_are_marble)
    l_sd_data = l_sd_data + tmp_l_sd_data
    if i == 0:
      col = len(l_sd_data)
  plot_holes(l_sd_raw_data, row, col, "raw_")
  my_analysis(l_sd_data, row, col)
  
  #plt.show()
  #plt.clf()
  
if __name__ == '__main__':
  if len(sys.argv) == 1:
    usage()
    exit(0)
  main(sys.argv[1:])
