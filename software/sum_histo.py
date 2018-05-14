import sys
import os
import pandas
import matplotlib.pyplot as plt
from decimal import Decimal
from toolbox import *
import math
from glob import glob

LEMs = sorted(glob("../data/eltos/CFR-35/*"))
#["CFR-35/A-059", "CFR-35/A-053", "CFR-35/A-054", "CFR-35/A-055", "CFR-35/A-057", "CFR-35/A-062", "CFR-35/A-065", "CFR-35/A-067", "CFR-35/A-068", "CFR-35/A-069", "CFR-35/A-070", "CFR-35/A-072", "CFR-35/A-066", "CFR-35/A-076", "CFR-35/A-075"]
analysis = ["LEM", "Copper"]
figID = 0

for ana in analysis:
  df_data = [ pandas.DataFrame({'z':list(map(float,open(lem + "/plots/"+ana+".txt","r").readline().split("\n")[0].split(";")))}) if os.path.isfile(lem + "/plots/"+ana+".txt") else print(lem + "/plots/"+ana+".txt not found") for lem in LEMs ]

  std_devs = []
  all_data = pandas.DataFrame()
  for i in range(0,len(df_data)):
    df_data[i] = df_data[i].ix[ df_data[i]['z'] > 1e-3 ]
    all_data = all_data.append(df_data[i], ignore_index=True)
    std_devs.append(df_data[i]['z'].describe()['std'])
    lemi = LEMs[i].split("/")[-1]
    column = "\\mathrm{" + lemi + ": }" + str(int(df_data[i]['z'].describe()['mean'])) + " \pm " + str(int(std_devs[-1]/math.sqrt(len(df_data[i])))) + "\;\mu m \\mathrm{; RMS: }" + str(int(std_devs[-1])) + "\;\mu m"
    df_data[i].rename(columns={'z':column}, inplace = True)

  #data_ranges_mins, data_ranges_maxs = zip(*[ [min(df_z[df_z.columns[0]]), max(df_z[df_z.columns[0]])] for df_z in df_data ])
  std_devs, df_data, LEMs = map(list,zip(*sorted(zip(std_devs, df_data, LEMs))))
  binsize = 5. #microns
  #myRange = [min(data_ranges_mins)*0.99,max(data_ranges_maxs)*1.01]
  myRange = [950,1250]
  if ana == "Copper":
    myRange = [20,100]
  nbins=int((myRange[1]-myRange[0])/binsize)
  fig = plt.figure(figID,figsize=(9, 9))
  sb = fig.add_subplot(111)
  sb.set_xlabel('Thickness ($\mu m$)')
  sb.set_ylabel('count / ' + str(int(binsize)) + ' $\mu m$')
  cm = plt.cm.get_cmap('plasma')

  for i in range(0,len(df_data)):
    #sb.hist(df_data[i][df_data[i].columns[0]], bins=nbins, range = myRange, facecolor = "none", alpha = 1, ls = "solid", edgecolor = cm(i/len(df_data)), linewidth = 3)
    sb.hist(df_data[i][df_data[i].columns[0]], bins=nbins, range = myRange, facecolor = cm(i/len(df_data)), alpha = 1, ls = "solid", edgecolor = 'black', linewidth = 1, label=r'$'+df_data[i].columns[0]+'$')
  sb.legend(loc='best')
  fig.savefig("./"+ana+"_sum_histo.pdf")
  figID = figID + 1
  
  statbox = FormStatBox(all_data['z'])
  fig2 = plt.figure(figID, figsize=(9, 9))
  sb2 = fig2.add_subplot(111)
  sb2.set_xlabel('Thickness ($\mu m$)')
  sb2.set_ylabel('count / ' + str(int(binsize)) + ' $\mu m$')
  i_count, binned_z , _  = sb2.hist(all_data[all_data.columns[0]], bins=nbins, range = myRange, ls = "solid", edgecolor = 'black', linewidth = 1, label=r'$'+all_data.columns[0]+'$')
  sb2.text(myRange[0], 0.3*i_count.max(), statbox, horizontalalignment='left', color = 'r')#, fontsize = .7*fontsize)
  figID = figID + 1
  fig2.savefig("./"+ana+"_sum_all_histo.pdf")
  
