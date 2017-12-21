import numpy as np
from gaussfit import *
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas
import scipy
from marble_fit import marble
import matplotlib.gridspec as gridspec

def my_analysis(l_sd_data, do_marble_fit = True):
  if do_marble_fit:
    l_sd_data = marble(l_sd_data)
  int_nrow = 0
  int_ncol = 0
  if int(math.sqrt(len(l_sd_data))) == math.sqrt(len(l_sd_data)):
    int_nrow = math.sqrt(len(l_sd_data))
    int_ncol = math.sqrt(len(l_sd_data))
  elif int(len(l_sd_data)/int(math.sqrt(len(l_sd_data)))) == len(l_sd_data)/int(math.sqrt(len(l_sd_data))):
    int_nrow = int(math.sqrt(len(l_sd_data)))
    int_ncol = len(l_sd_data)/int_nrow
  else:
    int_nrow = int(math.sqrt(len(l_sd_data)))
    int_ncol = int(math.sqrt(len(l_sd_data))) + 1
  l_sbID = [str(int_nrow) + str(int_ncol) + str(i+int_ncol) for i in range(1,len(l_sd_data)+1)]
  sb = [[],[]]
  fig = plt.figure(1,figsize=(24, 24))
  outer = gridspec.GridSpec(int_nrow, int_ncol, wspace=0.2, hspace=0.2)
  sd_all_data = l_sd_data[0].copy()
    
  
  for i in range(1,len(l_sd_data)):
    z_sup = [z for z in list(l_sd_data[i].values()) if z > (np.array(list(l_sd_data[i].values())).max()-np.array(list(l_sd_data[i].values())).min())/2]
    z_inf = [z for z in list(l_sd_data[i].values()) if z <= (np.array(list(l_sd_data[i].values())).max()-np.array(list(l_sd_data[i].values())).min())/2]
    inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[i-1], wspace=0.2, hspace=0.1)
    df_z_sup = pandas.DataFrame({'z':z_sup})
    df_z_inf = pandas.DataFrame({'z':z_inf})
    
    sb[0].append(plt.Subplot(fig, inner[0]))
    #sb.append(fig.add_subplot(l_sbID[i-1]))
    sb[0][i-1].set_title('Marble', fontsize=14)
    sb[0][i-1].set_xlabel('z(micrometer)', fontsize=14)
    sb[0][i-1].set_ylabel('count', fontsize=12)
    i_count_sup, _ , _ = sb[0][i-1].hist(df_z_inf['z'], bins='auto')
    #gaussfit(df_z['z'], sb[0][i-1])
    #plt.text(df_z['z'].min()+(df_z['z'].max()-df_z['z'].min())/4, i_count_sup.max()-(i_count_sup.max()-i_count_sup.min())/3, df_z['z'].describe(),horizontalalignment='center')
    fig.add_subplot(sb[0][i-1])
    
    sb[1].append(plt.Subplot(fig, inner[1]))
    #sb.append(fig.add_subplot(l_sbID[i-1]))
    sb[1][i-1].set_title('LEM', fontsize=14)
    sb[1][i-1].set_xlabel('z(micrometer)', fontsize=14)
    i_count_sup, _ , _ = sb[1][i-1].hist(df_z_sup['z'], bins='auto')
    #gaussfit(df_z['z'], sb[1][i-1])
    #plt.text(df_z['z'].min()+(df_z['z'].max()-df_z['z'].min())/4, i_count_sup.max()-(i_count_sup.max()-i_count_sup.min())/3, df_z['z'].describe(),horizontalalignment='center')
    fig.add_subplot(sb[1][i-1])
    sd_all_data.update(l_sd_data[i])
  plt.savefig("holes_histo.pdf")
  
  fig = plt.figure(2,figsize=(12, 8))
  tpl_x, tpl_z = zip(*(sd_all_data.items()))
  plt.scatter(tpl_x, tpl_z, color='black')
  plt.title("Holes Height vs length")
  plt.xlabel('x(micrometer)', fontsize=14)
  plt.ylabel('z(micrometer)', fontsize=12)
  plt.savefig("holes.pdf")
  #plt.show()
