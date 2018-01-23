import numpy as np
import matplotlib.cm as cm
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib
import pandas
import scipy
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tk
import itertools as it
from toolbox import *
from decimal import Decimal

def my_analysis(l_l_cut_data, row = 1, col = 1, sigmarble = 0, sigmaCopperLaser = 0, outdirectory = "./"):
  int_nrow, int_ncol = row, col #GetNcolNrow(l_l_cut_data)
  sb_plot_data = [[],[]]
  fig = plt.figure(1,figsize=(18*col, 18*row))
  plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
  outer = gridspec.GridSpec(int_nrow, int_ncol, wspace=0.2, hspace=0.3)
  l_l_all_data = []
  k = 0
  thick_sigmathick_rim_sigmarim = [[],[],[],[],[],[],[],[],[]]
  fontsize = 2.5
  if row == 1:
    fontsize = 40
  for i in range(0,len(l_l_cut_data)):
    print("   Analysing hole ",i+1,"...")
    #skip_inf = True
    plot_color = 'b'
    plot_range_sup = [900,1300]
    #plot_range_inf = [-50,50]
    underflow = len([z for z in l_l_cut_data[i][1] if z < plot_range_sup[0]])
    overflow = len([z for z in l_l_cut_data[i][1] if z > plot_range_sup[1]])
    if overflow > len(l_l_cut_data[i])/4:
      plot_range_sup[1]=1300
      plot_color = 'g'
    if underflow > len(l_l_cut_data[i])/4:
      plot_range_sup[0]=600
      plot_color = 'g'
    binsize = 5 #microns
    nbin_sup=int((plot_range_sup[1]-plot_range_sup[0])/binsize)
    #nbin_inf=int((plot_range_inf[1]-plot_range_inf[0])/binsize)
    if i % col == 0:
      l_l_all_data.append(l_l_cut_data[i].copy())
      k = k + 1
    #lzinf = [ z for z in l_l_cut_data[i][1] if z > plot_range_inf[0] and z < plot_range_inf[1] ]
    lzsup = [ z for z in l_l_cut_data[i][1] if z > plot_range_sup[0] and z < plot_range_sup[1] ]
    #df_z_inf = pandas.DataFrame({'z':lzinf})
    df_z_sup = pandas.DataFrame({'z':lzsup})
    #for j in range(0,len(l_l_cut_data[i][0])):
      #x = l_l_cut_data[i][0][j]
      #z = l_l_cut_data[i][1][j]
      #if z > plot_range_sup[0] and z < plot_range_sup[1]:
        #df_z_sup.loc[len(df_z_sup)] = z
      #if z > plot_range_inf[0] and z < plot_range_inf[1]:
        #df_z_inf.loc[len(df_z_inf)] = z
    #df_z_inf = cut(df_z_inf, [5])
    df_z_sup = cut(df_z_sup, [5])
    #if not skip_inf:
      #inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[i], wspace=0.2, hspace=0.1)
      #sb_plot_data[0].append(plt.Subplot(fig, inner[0]))
      #sb_plot_data[0][i].set_title('hole ' + str(i+1) + ' Marble', fontsize=14)
      #sb_plot_data[0][i].set_xlabel('z(micrometer)', fontsize=14)
      #i_count_inf, _ , _ = sb_plot_data[0][i].hist(df_z_inf['z'], bins=nbin_inf, range = [plot_range_inf[0],plot_range_inf[1]], color = 'b')
      #statbox = FormStatBox(df_z_inf['z'])
      #sb_plot_data[0][i].text(plot_range_inf[0], 0.9*i_count_inf.max(), statbox,horizontalalignment='left')
      #fig.add_subplot(sb_plot_data[0][i])
      #sb_plot_data[1].append(plt.Subplot(fig, inner[1]))
    #else:
    inner = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[i], wspace=0.2, hspace=0.1)
    sb_plot_data[1].append(plt.Subplot(fig, outer[i]))
    i_count_sup, binned_z , _ = sb_plot_data[1][i].hist(df_z_sup['z'], bins=nbin_sup, range = [plot_range_sup[0],plot_range_sup[1]], color = plot_color)
    y1, y2 = sb_plot_data[1][i].get_window_extent().get_points()[:, 1]
    yscale = (y2-y1)/(1.1*i_count_sup.max())
    #fontsize = 5
    sb_plot_data[1][i].set_title('hole  ' + str(i+1) + ' LEM', fontsize = 1.2*fontsize)
    sb_plot_data[1][i].set_xlabel('z(micrometer)', fontsize = fontsize, labelpad = .1)
    sb_plot_data[1][i].set_ylabel('count', fontsize = fontsize, labelpad = 1)
    sb_plot_data[1][i].tick_params(axis='both', which='major', labelsize = fontsize, length = .8*fontsize, width = .1*fontsize, pad = 0*fontsize)
    sb_plot_data[1][i].set_ylim([0,1.1*i_count_sup.max()])
    sb_plot_data[1][i].title.set_position([.5,.9])
    [i.set_linewidth(0.1) for i in sb_plot_data[1][i].spines.values()]
    binned_z = binned_z[:-1]
    Max = find_local_max(i_count_sup.copy(),np.array(binned_z.copy()), binsize)
    if len(Max) == 0:
      Max = [1140,1070]
    if len(Max) == 1:
      Max.append(Max[0]-60)
    possible_maximums = [ z for z in Max if (z < Max[0] and Max[0]-z < 90 and Max[0]-z > 40) or z == Max[0] ]
    if len(possible_maximums) == 0:
      possible_maximums.append([1140,1070])
    
    print("   fitting data...")
    
    tmp_thick_sigmathick_rim_sigmarim, gaussians_params, fitted_func1, fitted_func2, chisquare = mydoublefit(df_z_sup, binned_z, i_count_sup, possible_maximums, plot_range_sup, sigmarble)
    
    for j in range(0,len(thick_sigmathick_rim_sigmarim)):
      thick_sigmathick_rim_sigmarim[j].append(tmp_thick_sigmathick_rim_sigmarim[j])
    if len(fitted_func1) != 0:
      sb_plot_data[1][i].plot(fitted_func1[0], fitted_func1[1], 'r-', linewidth = .5)
    if len(fitted_func2) != 0:
      sb_plot_data[1][i].plot(fitted_func2[0], fitted_func2[1], 'r-', linewidth = .5)
    statbox = FormStatBox(df_z_sup['z'])
    fitbox = FormFitBox(gaussians_params, df_z_sup['z'], chisquare)
    print(fitbox)
    sb_plot_data[1][i].text(plot_range_sup[0], 0.3*i_count_sup.max(), fitbox, horizontalalignment='left', color = 'r', fontsize = .7*fontsize)
    fig.add_subplot(sb_plot_data[1][i])
    for xz in range(0,len(l_l_all_data[k-1])):
      l_l_all_data[k-1][xz] = l_l_all_data[k-1][xz] + l_l_cut_data[i][xz]
    print("")
  fig.savefig(outdirectory + "/holes_histo.pdf")
  print("   Holes histograms saved in " + outdirectory + "/holes_histo.pdf")
  plot_holes(l_l_all_data, row, col, "", outdirectory)
  plot_thicknesses_map(thick_sigmathick_rim_sigmarim, row, col, sigmaCopperLaser, outdirectory)

def plot_holes(l_l_all_data, row, col, name = "", outdirectory = "./"):
  ID = 2
  if name == "":
    ID = 3
  fig = plt.figure(ID,figsize=(5*col, 2*row))
  grid = gridspec.GridSpec(row, 1, wspace=0.2, hspace=0.5)
  sb_plot_holes = [ fig.add_subplot(int(str(row) + '1' + str(i+1))) for i in range(0,len(l_l_all_data)) ]
  for i in range(0,len(l_l_all_data)):
    tpl_x = np.array(l_l_all_data[i][0])
    tpl_z = np.array(l_l_all_data[i][1])
    sb_plot_holes[i].scatter(tpl_x, tpl_z, s=1, color='black')
    sb_plot_holes[i].set_title("row "+str(i+1))
    if i == row-1:
      sb_plot_holes[i].set_xlabel('x(micrometer)', fontsize=14)
    sb_plot_holes[i].set_ylabel('z(micrometer)', fontsize=12)
  fig.savefig(outdirectory + "/" + name + "holes.pdf")
  print("   Histo saved in " + outdirectory + "/" + name + "holes.pdf")
  
def plot_thicknesses_map(thick_sigmathick_rim_sigmarim, row, col, sigmaCopperLaser = 0, outdirectory = "./"):
  thick, sigmathick, FR4, sigmaFR4, Cu, NLem, NFR4, N, sigmarble = np.array(thick_sigmathick_rim_sigmarim)
  sigmaMEANthick = [ math.sqrt(sigthick**2/nLem**2 + sigmar**2/n**2) for sigthick,nLem,sigmar,n in zip(sigmathick,NLem,sigmarble,N) ]
  sigmaMEANFR4 = [ math.sqrt(sigfr4**2/nFR4**2 + sigmar**2/n**2) for sigfr4,nFR4,sigmar,n in zip(sigmaFR4,NFR4,sigmarble,N) ]
  sigmaMEANCU = [ math.sqrt(sigfr4**2 + sigthick) for sigfr4,sigthick in zip(sigmaMEANFR4,sigmaMEANthick) ]
  sigmaMEANFR4 = [ math.sqrt(4*sigfr4**2 + sigthick**2) for sigfr4,sigthick in zip(sigmaMEANFR4,sigmaMEANthick) ]
  print("Result : ")
  for i in range(0,len(thick)):
    print(" Hole ",i,": ",thick[i],"+/-",sigmaMEANthick[i]," microns thick. Thickness standard deviation: ",math.sqrt(sigmathick[i]**2-sigmaCopperLaser**2))
    print("    Cu: ",Cu[i],"+/-",sigmaMEANCU[i]," microns, FR4: ", FR4[i],"+/-",sigmaMEANFR4[i]," microns")
  plot_2D_map(thick,sigmaMEANthick,4,"map of LEM thickness. Each pixel is a measurement hole","2D_LEM_thickness_distri.pdf", row, col, 1050, 1250, outdirectory)
  plot_2D_map(FR4,sigmaMEANFR4,5,"map of FR4 thickness. Each pixel is a measurement hole","2D_FR4_thickness_distri.pdf", row, col, 800, 1300, outdirectory)
  plot_2D_map(Cu,sigmaMEANCU,6,"map of rim thickness. Each pixel is a measurement hole","2D_rim_thickness_distri.pdf", row, col, 30,90, outdirectory)
  
def plot_2D_map(z,sigma,figID,title, savename, row, col, v0 = 800, v1 = 1300, outdirectory = "./"):
  fig = plt.figure(figID,figsize=(3*col, 3*row))
  sb_plot_2D_map = fig.add_subplot(111)
  
  x = [i%(col)+1 for i in range(0,len(z)) if z[i] != 0]
  y = [row-(int(i/col)) for i in range(0,len(z)) if z[i] != 0]
  z = z[z != 0]
  
  marker_size = good_marker_size(x,y,fig,sb_plot_2D_map)
  p = sb_plot_2D_map.scatter(x,y,c=z, s=marker_size[0]*marker_size[1], marker='.', cmap=cm.plasma, linewidth=0, vmin=v0, vmax=v1)
  #sb_plot_2D_map.axis([min(x)-1., max(x)+1., min(y)-1., max(y)+1.])
  sb_plot_2D_map.set_yticks([])
  sb_plot_2D_map.set_xticks([])
  sb_plot_2D_map.set_title(title)
  for i in range(0,len(z)):
    mytext = str(i+1) + "\n" + str(round(Decimal(z[i]),2)) + "+/-" + str(round(Decimal(sigma[i]),4))
    sb_plot_2D_map.annotate(mytext,xy=(x[i]-0.25,y[i]), color='white', path_effects=[PathEffects.withStroke(linewidth=2, foreground="black")])
  fig.colorbar(p)
  fig.savefig(outdirectory + "/" + savename)
  print("   2D map saved in " + outdirectory + "/" + savename)
  
def mydoublefit(df_z, binned_z, i_count, maxima, plot_range, sigmarble):
  thick_sigmathick_rim_sigmarim = []
  attempt = 1
  fit_par1 = [len(df_z)/2,maxima[0],6]
  fit_par_range1 = [0,1000*fit_par1[0],plot_range[0],plot_range[1],0,30]
  range_z1 = np.array([z for z in binned_z if z < maxima[0]+30 and z > maxima[0]-30])
  range_count1 = np.array([i for [z,i] in zip(binned_z, i_count) if z < maxima[0]+30 and z > maxima[0]-30])
  print("    Fitting whole LEM...")
  gaussians_param1, rsquare1, result1 = singlegaussfit(range_z1, range_count1, fit_par1, fit_par_range1)
  chisquare = [rsquare1]
  if rsquare1 < 0.8:
    print("   Fit failed")
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    return thick_sigmathick_rim_sigmarim, [], [], [], []
  z1 = gaussians_param1[1]
  sigma1 = gaussians_param1[2]
  rsquares = []
  fitted_func = []
  gaussians_params = []
  gaussians_param2 = []
  ###########while attempt < len(possible_maximums):###############
  while attempt < len(maxima):
    print("    attempting to fit FR4...")
    fit_par2 = [len(df_z)/2,maxima[attempt],12]
    fit_par_range2 = [0,1000*fit_par2[0],plot_range[0],plot_range[1],0,30]
    range_z2 = np.array([z for z in binned_z if z < maxima[attempt]+50 and z > maxima[attempt]-50])
    range_count2 = np.array([i for [z,i] in zip(binned_z, i_count) if z < maxima[attempt]+50 and z > maxima[attempt]-50])
    gaussians_param2, rsquare2, result2 = singlegaussfit(range_z2, range_count2, fit_par2, fit_par_range2)
    if rsquare2 > 0 and abs(gaussians_param1[1]-gaussians_param2[1]) < 90 and abs(gaussians_param1[1]-gaussians_param2[1]) > 40:# and gaussians_param1[0]/gaussians_param2[0] < 9:
      rsquares.append(rsquare2)
      fitted_func.append(result2)
      gaussians_params.append(gaussians_param2)
    elif rsquare2 <= 0:
      print("    discarding result (R**2 < 0)")
    elif abs(gaussians_param1[1]-gaussians_param2[1]) > 90 and abs(gaussians_param1[1]-gaussians_param2[1]) < 40:
      print("    discarding result (means too close or too far away)")
    else:
      print("    discarding result (FR4 fit too small)")
    attempt = attempt + 1
  ###########while attempt < len(possible_maximums):###############
  if len(rsquares) == 0:
    print("    ...FR4 fit failed")
    NLem = df_z['z'].describe()['count']
    thick_sigmathick_rim_sigmarim.append(z1)
    thick_sigmathick_rim_sigmarim.append(sigma1)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(NLem)
    thick_sigmathick_rim_sigmarim.append(0)
    thick_sigmathick_rim_sigmarim.append(NLem)
    thick_sigmathick_rim_sigmarim.append(sigmarble)
    return thick_sigmathick_rim_sigmarim, gaussians_param1 + [0,0,0], result1, [], chisquare
      
  else:
    best_fit = [ir for ir,r in enumerate(rsquares) if r == min(rsquares)][0]
    gaussians_param2 = gaussians_params[best_fit]
    NLem = df_z['z'].describe()['count']*gaussians_param1[0]/(gaussians_param1[0]+gaussians_param2[0])
    NFR4 = df_z['z'].describe()['count']*gaussians_param2[0]/(gaussians_param1[0]+gaussians_param2[0])
    chisquare.append(rsquare2)
    print("    best fit : ",gaussians_param2[1]," with Chi**2=",min(rsquares))
    z2 = gaussians_param2[1]
    sigma2 = gaussians_param2[2]
    thick_sigmathick_rim_sigmarim.append(z1)
    thick_sigmathick_rim_sigmarim.append(sigma1)
    thick_sigmathick_rim_sigmarim.append(2*z2-z1)
    thick_sigmathick_rim_sigmarim.append(sigma2)
    thick_sigmathick_rim_sigmarim.append(z1-z2)
    thick_sigmathick_rim_sigmarim.append(NLem)
    thick_sigmathick_rim_sigmarim.append(NFR4)
    thick_sigmathick_rim_sigmarim.append(df_z['z'].describe()['count'])
    thick_sigmathick_rim_sigmarim.append(sigmarble)
    return thick_sigmathick_rim_sigmarim, list(gaussians_param1) + list(gaussians_param2), result1, fitted_func[best_fit], chisquare

