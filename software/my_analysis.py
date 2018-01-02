import numpy as np
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas
import scipy
import matplotlib.gridspec as gridspec
from toolbox import *

def my_analysis(l_sd_data, row = 1):
    plot_data(l_sd_data, row)
  
def plot_data(l_sd_data, row = 1):
  int_nrow, int_ncol = GetNcolNrow(l_sd_data)
  sb = [[],[]]
  fig = plt.figure(1,figsize=(36, 36))
  plt.subplots_adjust(left=0.01, bottom=0.02, right=0.99, top=0.99, wspace=0.2, hspace=0.2)
  outer = gridspec.GridSpec(int_nrow, int_ncol, wspace=0.2, hspace=0.2)
  l_sd_all_data = []
  k = 0
  plot_range_sup = [0,0]
  plot_range_inf = [0,0]
  for i in range(0,len(l_sd_data)):
    if i % row == 0:#j'utilise le range de la hauteur du premier trou comme reference, pour avoir le meme range pour tous les trous
      l_sd_all_data.append(l_sd_data[0].copy())
      k = k + 1
    sd_z_sup = SortedDict((x,l_sd_data[i][x]) for x in l_sd_data[i] if l_sd_data[i][x] > (np.array(list(l_sd_data[i].values())).max()-np.array(list(l_sd_data[i].values())).min())/2)
    sd_z_inf = SortedDict((x,l_sd_data[i][x]) for x in l_sd_data[i] if l_sd_data[i][x] <= (np.array(list(l_sd_data[i].values())).max()-np.array(list(l_sd_data[i].values())).min())/2)
    df_z_sup = None
    df_z_inf = None
    if i == 0:
      sd_z_sup = cut(sd_z_sup, [4])
      sd_z_inf = cut(sd_z_inf, [4])
      df_z_sup = pandas.DataFrame({'z':list(sd_z_sup.values())})
      df_z_inf = pandas.DataFrame({'z':list(sd_z_inf.values())})
      plot_range_sup[0]=float(df_z_sup.mean())-5*math.sqrt(float(df_z_sup.var()))
      plot_range_sup[1]=float(df_z_sup.mean())+5*math.sqrt(float(df_z_sup.var()))
      plot_range_inf[0]=float(df_z_inf.mean())-3*math.sqrt(float(df_z_inf.var()))
      plot_range_inf[1]=float(df_z_inf.mean())+3*math.sqrt(float(df_z_inf.var()))
    else:
      sd_z_sup = cut(sd_z_sup, plot_range_sup)
      sd_z_inf = cut(sd_z_inf, plot_range_inf)
      df_z_sup = pandas.DataFrame({'z':list(sd_z_sup.values())})
      df_z_inf = pandas.DataFrame({'z':list(sd_z_inf.values())})
      
    inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[i], wspace=0.2, hspace=0.1)
    sb[0].append(plt.Subplot(fig, inner[0]))
    sb[0][i].set_title('hole ' + str(i+1) + ' Marble', fontsize=14)
    sb[0][i].set_xlabel('z(micrometer)', fontsize=14)
    sb[0][i].set_ylabel('count', fontsize=12)
    i_count_inf, _ , _ = sb[0][i].hist(df_z_inf['z'], bins='auto', range = [plot_range_inf[0],plot_range_inf[1]])#, color = 'C1')
    #gaussfit(df_z_inf['z'], sb[0][i])
    statbox = FormStatBox(df_z_inf['z'])
    sb[0][i].text(plot_range_inf[0], 0.9*i_count_inf.max(), statbox,horizontalalignment='left')
    fig.add_subplot(sb[0][i])
    
    sb[1].append(plt.Subplot(fig, inner[1]))
    sb[1][i].set_title('hole  ' + str(i+1) + ' LEM', fontsize=14)
    sb[1][i].set_xlabel('z(micrometer)', fontsize=14)
    i_count_sup, binned_z , _ = sb[1][i].hist(df_z_sup['z'], bins='auto', range = [plot_range_sup[0],plot_range_sup[1]])#, color = 'C0')
    binned_z = binned_z[:-1]
    fit_par = [len(l_sd_data[i])/2,1050,math.sqrt(float(df_z_sup.var()))/2,len(l_sd_data[i])/2,1100,math.sqrt(float(df_z_sup.var()))/2]
    doublegaussfit(binned_z, i_count_sup, fit_par, sb[1][i])
    #gaussfit(df_z_sup['z'], sb[1][i])
    statbox = FormStatBox(df_z_sup['z'])
    sb[1][i].text(plot_range_sup[0], 0.9*i_count_sup.max(), statbox,horizontalalignment='left')
    fig.add_subplot(sb[1][i])
    l_sd_all_data[k-1].update(l_sd_data[i])
  plt.savefig("holes_histo.pdf")
  plt.clf()
  
  sb = []
  fig = plt.figure(2,figsize=(12, 18))
  grid = gridspec.GridSpec(row, 1, wspace=0.2, hspace=0.5)
  for i in range(0,len(l_sd_all_data)):
    sb.append(fig.add_subplot(int(str(row)+'1'+str(i+1))))
    tpl_x, tpl_z = zip(*(l_sd_all_data[i].items()))
    plt.scatter(tpl_x, tpl_z, color='black')
    sb[i].set_title("Holes of row " + str(i+1))
    if i == row-1:
      sb[i].set_xlabel('x(micrometer)', fontsize=14)
    sb[i].set_ylabel('z(micrometer)', fontsize=12)
  plt.savefig("holes.pdf")
  #plt.show()
  plt.clf()

def GetNcolNrow(l_sd_data):
  if int(math.sqrt(len(l_sd_data))) == math.sqrt(len(l_sd_data)):
    int_nrow = int(math.sqrt(len(l_sd_data)))
    int_ncol = int(math.sqrt(len(l_sd_data)))
  elif int(len(l_sd_data)/int(math.sqrt(len(l_sd_data)))) == len(l_sd_data)/int(math.sqrt(len(l_sd_data))):
    int_nrow = int(math.sqrt(len(l_sd_data)))
    int_ncol = int(len(l_sd_data)/int_nrow)
  else:
    int_nrow = int(math.sqrt(len(l_sd_data)))
    int_ncol = int(math.sqrt(len(l_sd_data))) + 1
  return int_nrow, int_ncol
