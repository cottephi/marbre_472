import sys
import os
import pandas
import matplotlib.pyplot as plt
from toolbox import *

if len(sys.argv) != 3:
  print("Error: need 2 input files")
  exit(1)
if not os.path.isfile(sys.argv[1]):
  print("Error: file ",sys.argv[1]," not found")
  exit(1)
if not os.path.isfile(sys.argv[2]):
  print("Error: file ",sys.argv[1]," not found")
  exit(1)
title = os.path.basename(os.path.dirname(os.path.dirname(sys.argv[2]))) + "-" + os.path.basename(os.path.dirname(os.path.dirname(sys.argv[1]))) + "_" +  os.path.basename(sys.argv[1])
outdirectory = os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[1])))
infile1 = open(sys.argv[1], "r")
infile2 = open(sys.argv[2], "r")
z1 = infile1.readline().split(";")
z2 = infile2.readline().split(";")
if len(z1) != len(z2):
  print("Error: both files do not have the same number of holes.")
  exit(1)
infile1.close()
infile2.close()

df_dz = pandas.DataFrame({'dz':[ float(l2)-float(l1) for l1,l2 in zip(z1,z2) ]})
binsize = 1. #microns
myRange = [min(df_dz['dz'])*0.99,max(df_dz['dz'])*1.01]
nbins=int((myRange[1]-myRange[0])/binsize)
fig = plt.figure(0,figsize=(9, 9))
sb = fig.add_subplot(111)
sb.set_title(title)
sb.set_xlabel('Thickness difference(micrometer)')
sb.set_ylabel('count / 5 microns')
i_count, binned_dz , _ = sb.hist(df_dz['dz'], bins=nbins, range = myRange)
#binned_dz = binned_dz[:-1]
#print("Fitting thick diff...")
#fit_par = [len(df_dz),0,5]
#gaussians_param, chisquare, result = singlegaussfit(binned_dz, i_count, fit_par)
#sb.plot(result[0], result[1], 'r-', linewidth = .5)
#fitbox = FormFitBox(gaussians_param, df_dz['dz'], [chisquare])
#sb.text(plt.xlim()[0], 0.5*i_count.max(), fitbox, horizontalalignment='left', color = 'r', fontsize = 20)
statbox = FormStatBox(df_dz['dz'])
sb.text(plt.xlim()[0], 0.5*i_count.max(), statbox, horizontalalignment='left', color = 'r', fontsize = 20)
fig.savefig(outdirectory + "/" + title + '.pdf')
print("thickness diff histo saved in " + outdirectory + "/" + title + '.pdf')
