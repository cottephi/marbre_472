import numpy as np
import matplotlib.cm as cm
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas
import scipy
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tk
import itertools as it
from toolbox import *

def my_analysis(l_sd_data, row = 1):
    plot_data(l_sd_data, row)
  
def plot_data(l_sd_data, row = 1):
  int_nrow, int_ncol = GetNcolNrow(l_sd_data)
  sb_plot_data = [[],[]]
  fig = plt.figure(1,figsize=(36, 36))
  plt.subplots_adjust(left=0.01, bottom=0.02, right=0.99, top=0.99, wspace=0.2, hspace=0.2)
  outer = gridspec.GridSpec(int_nrow, int_ncol, wspace=0.2, hspace=0.2)
  l_sd_all_data = []
  k = 0
  x_y_thick_sigmathick_rim_sigmarim = [[],[],[],[],[],[]]
  for i in range(0,len(l_sd_data)):
    plot_color = 'C0'
    plot_range_sup = [1000,1300]
    plot_range_inf = [-50,50]
    underflow = len([z for z in l_sd_data[i].values() if z < plot_range_sup[0]])
    overflow = len([z for z in l_sd_data[i].values() if z > plot_range_sup[1]])
    if overflow > len(l_sd_data[i].values())/4:
      plot_range_sup[1]=1300
      plot_color = 'g'
    if underflow > len(l_sd_data[i].values())/4:
      plot_range_sup[0]=800
      plot_color = 'g'
    binsize=5 #microns
    nbin_sup=int((plot_range_sup[1]-plot_range_sup[0])/binsize)
    nbin_inf=int((plot_range_inf[1]-plot_range_inf[0])/binsize)
    if i % row == 0:#j'utilise le range de la hauteur du premier trou comme reference, pour avoir le meme range pour tous les trous
      l_sd_all_data.append(l_sd_data[0].copy())
      k = k + 1
    #sd_z_sup = SortedDict((x,l_sd_data[i][x]) for x in l_sd_data[i]) if l_sd_data[i][x] > (np.array(list(l_sd_data[i].values())).max()-np.array(list(l_sd_data[i].values())).min())/2)
    sd_z_sup = SortedDict((x,l_sd_data[i][x]) for x in l_sd_data[i] if l_sd_data[i][x] > plot_range_sup[0] and l_sd_data[i][x] < plot_range_sup[1])
    sd_z_inf = SortedDict((x,l_sd_data[i][x]) for x in l_sd_data[i] if l_sd_data[i][x] > plot_range_inf[0] and l_sd_data[i][x] < plot_range_inf[1])
    #sd_z_inf = SortedDict((x,l_sd_data[i][x]) for x in l_sd_data[i])# if l_sd_data[i][x] <= (np.array(list(l_sd_data[i].values())).max()-np.array(list(l_sd_data[i].values())).min())/2)
    if len(sd_z_inf) == 0:
      skip_inf = True
    df_z_sup = None
    df_z_inf = None
    sd_z_sup = cut(sd_z_sup, [5])
    sd_z_inf = cut(sd_z_inf, [5])
    df_z_sup = pandas.DataFrame({'z':list(sd_z_sup.values())})
    df_z_inf = pandas.DataFrame({'z':list(sd_z_inf.values())})
    
    inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[i], wspace=0.2, hspace=0.1)
    sb_plot_data[0].append(plt.Subplot(fig, inner[0]))
    sb_plot_data[0][i].set_title('hole ' + str(i+1) + ' Marble', fontsize=14)
    sb_plot_data[0][i].set_xlabel('z(micrometer)', fontsize=14)
    sb_plot_data[0][i].set_ylabel('count', fontsize=12)
    i_count_inf, _ , _ = sb_plot_data[0][i].hist(df_z_inf['z'], bins=nbin_inf, range = [plot_range_inf[0],plot_range_inf[1]], color = 'C1')
    #gaussfit(df_z_inf['z'], sb_plot_data[0][i])
    statbox = FormStatBox(df_z_inf['z'])
    sb_plot_data[0][i].text(plot_range_inf[0], 0.9*i_count_inf.max(), statbox,horizontalalignment='left')
    fig.add_subplot(sb_plot_data[0][i])
    
    sb_plot_data[1].append(plt.Subplot(fig, inner[1]))
    sb_plot_data[1][i].set_title('hole  ' + str(i+1) + ' LEM', fontsize=14)
    sb_plot_data[1][i].set_xlabel('z(micrometer)', fontsize=14)
    i_count_sup, binned_z , _ = sb_plot_data[1][i].hist(df_z_sup['z'], bins=nbin_sup, range = [plot_range_sup[0],plot_range_sup[1]], color = plot_color)
    sb_plot_data[1][i].set_ylim([0,1.1*i_count_sup.max()])
    binned_z = binned_z[:-1]
    Max = find_local_max(i_count_sup.copy(),np.array(binned_z.copy()), binsize)
    if Max[0] == Max[1]:
      Max[0] = Max[1]-50
    possible_maximums = list(it.permutations(Max,2))
    attempt = 0
    rsquares = []
    gaussians_params = []
    fitted_func = []
    while attempt < len(possible_maximums):
      rsquare = 0
      result = []
      gaussians_param = []
      fit_par = [len(l_sd_data[i])/2,possible_maximums[attempt][0],math.sqrt(float(df_z_sup.var()))/2,len(l_sd_data[i])/2,possible_maximums[attempt][1],math.sqrt(float(df_z_sup.var()))/2]
      fit_par_range = [1000*fit_par[0],-1000*fit_par[0],plot_range_sup[0],plot_range_sup[1],0,40,1000*fit_par[3],-1000*fit_par[3],plot_range_sup[0],plot_range_sup[1],0,40]
      gaussians_param, rsquare, result = doublegaussfit(binned_z, i_count_sup, fit_par, sb_plot_data[1][i],fit_par_range)
      rsquares.append(rsquare)
      fitted_func.append(result)
      gaussians_params.append(gaussians_param)
      attempt = attempt + 1
    best_fit = [ir for ir,r in enumerate(rsquares) if r == max(rsquares)][0]
    sb_plot_data[1][i].plot(fitted_func[best_fit][0], fitted_func[best_fit][1], 'r-')
    gaussians_params = gaussians_params[best_fit]
    z1 = min(gaussians_params[1], gaussians_params[4])
    if z1 == gaussians_params[1]:
      sigma1 = gaussians_params[2]
      sigma2 = gaussians_params[5]
    else:
      sigma1 = gaussians_params[5]
      sigma2 = gaussians_params[2]
    z2 = max(gaussians_params[1], gaussians_params[4])
    x_y_thick_sigmathick_rim_sigmarim[0].append(1+i%row)
    x_y_thick_sigmathick_rim_sigmarim[1].append(5-int((i-i%row)/row))
    x_y_thick_sigmathick_rim_sigmarim[2].append(2*z1-z2)
    x_y_thick_sigmathick_rim_sigmarim[3].append(math.sqrt(4*sigma1**2+sigma2**2))
    x_y_thick_sigmathick_rim_sigmarim[4].append(z2-z1)
    x_y_thick_sigmathick_rim_sigmarim[5].append(math.sqrt(sigma1**2+sigma2**2))
    statbox = FormStatBox(df_z_sup['z'])
    fitbox = FormFitBox(gaussians_params, df_z_sup['z'])
    #sb_plot_data[1][i].text(plot_range_sup[0], 0.9*i_count_sup.max(), statbox,horizontalalignment='left')
    sb_plot_data[1][i].text(Max[0]*0.83, 0.75*i_count_sup.max(), fitbox,horizontalalignment='left', color = 'r')
    fig.add_subplot(sb_plot_data[1][i])
    l_sd_all_data[k-1].update(l_sd_data[i])
  fig.savefig("holes_histo.pdf")
  plot_holes(l_sd_all_data, row)
  plot_thicknesses_map(x_y_thick_sigmathick_rim_sigmarim)

def plot_holes(l_sd_all_data, row, title = ""):
  sb_plot_holes = []
  ID = 2
  if title == "":
    ID = 3
  fig = plt.figure(ID,figsize=(6, 9))
  grid = gridspec.GridSpec(row, 1, wspace=0.2, hspace=0.5)
  for i in range(0,len(l_sd_all_data)):
    sb_plot_holes.append(fig.add_subplot(int(str(row) + '1' + str(i+1))))
    tpl_x, tpl_z = zip(*(l_sd_all_data[i].items()))
    sb_plot_holes[i].scatter(tpl_x, tpl_z, s=1, color='black')
    sb_plot_holes[i].set_title("Holes of row " + str(i+1))
    if i == row-1:
      sb_plot_holes[i].set_xlabel('x(micrometer)', fontsize=14)
    sb_plot_holes[i].set_ylabel('z(micrometer)', fontsize=12)
  fig.savefig(title + "holes.pdf")
  
def plot_thicknesses_map(x_y_thick_sigmathick_rim_sigmarim):
  x,y,thick,sigmathick,rim,sigmarim = np.array(x_y_thick_sigmathick_rim_sigmarim)
  plot_2D_map(x,y,thick,4,"map of LEM thickness. Each pixel is a measurement hole","2D_LEM_thickness_distri.pdf")
  x,y,thick,sigmathick,rim,sigmarim = np.array(x_y_thick_sigmathick_rim_sigmarim)
  plot_2D_map(x,y,rim,5,"map of rim thickness. Each pixel is a measurement hole","2D_rim_thickness_distri.pdf")
  
def plot_2D_map(x,y,z,figID,title, savename):
  fig = plt.figure(figID,figsize=(12, 12))
  sb_plot_2D_map = fig.add_subplot(111)
  marker_size = good_marker_size(x,y,fig,sb_plot_2D_map)
  p = sb_plot_2D_map.scatter(x,y,c=z, s=marker_size**2, marker='.', cmap=cm.plasma, linewidth=0)
  #sb_plot_2D_map.axis([min(x)-1., max(x)+1., min(y)-1., max(y)+1.])
  sb_plot_2D_map.set_yticks([])
  sb_plot_2D_map.set_xticks([])
  sb_plot_2D_map.set_title(title)
  for i,j in zip(x,y):
    ID = int(i+(5-j)*5)
    mytext = str(ID) + "\n" + str(int(z[int(i-1+(5-j)*5)]))
    sb_plot_2D_map.annotate(mytext,xy=(i,j), color='g')
  fig.colorbar(p)
  fig.savefig(savename)

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
