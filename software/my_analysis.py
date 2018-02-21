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
import csv

def my_analysis(l_l_cut_data, row = 1, col = 1, sigmarble = 0, sigmaCopperLaser = 0, outdirectory = "./"):
  fig = plt.figure(1,figsize=(3*col, 3*row))
  plt.subplots_adjust(left=-0.1*col/4 + 0.165, bottom=-0.05*row/4 + 0.105, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
  outer = gridspec.GridSpec(row, col, wspace=0.2, hspace=0.3)
  l_l_all_data = []
  k = 0
  thick_sigmathick_rim_sigmarim = [[],[],[],[],[],[],[],[],[]]
  fontsize = 10
  #if row == 1:
    #fontsize = 40
  for i in range(0,len(l_l_cut_data)):
    print("   Analysing hole ",i+1,"...")
    plot_color = 'b'
    plot_range_sup = [900,1300]
    binsize = 5 #microns
    if i % col == 0:
      l_l_all_data.append(l_l_cut_data[i].copy())
      k = k + 1
    lzsup = cut([ z for z in l_l_cut_data[i][1] if z > 500 ],[5])
    in_range = len([z for z in lzsup if z > plot_range_sup[0] and z < plot_range_sup[1] ])/len(lzsup)
    underflow = len([z for z in lzsup if z < plot_range_sup[0]])
    overflow = len([z for z in lzsup if z > plot_range_sup[1]])
    while in_range < .9:
      if underflow/len(lzsup) > .15:
        plot_range_sup[0] = plot_range_sup[0] - 100
        print("   Decreasing min range to",plot_range_sup[0])
      if overflow/len(lzsup) > .15:
        plot_range_sup[1] = plot_range_sup[1] + 100
        print("   Increasing max range to",plot_range_sup[1])
      underflow = len([z for z in lzsup if z < plot_range_sup[0]])
      overflow = len([z for z in lzsup if z > plot_range_sup[1]])
      in_range = len([z for z in lzsup if z > plot_range_sup[0] and z < plot_range_sup[1] ])
      plot_color = 'g'

    nbin_sup=int((plot_range_sup[1]-plot_range_sup[0])/binsize)
    lzsup = [ z for z in lzsup if z > plot_range_sup[0] and z < plot_range_sup[1] ]
    df_z_sup = pandas.DataFrame({'z':lzsup})
   
    sb_plot_data = plt.Subplot(fig, outer[i])
    i_count_sup, binned_z , _ = sb_plot_data.hist(df_z_sup['z'], bins=nbin_sup, range = [plot_range_sup[0],plot_range_sup[1]], color = plot_color)
    y1, y2 = sb_plot_data.get_window_extent().get_points()[:, 1]
    yscale = (y2-y1)/(1.1*i_count_sup.max())
    sb_plot_data.set_title('hole  ' + str(i+1) + ' LEM', fontsize = 1.2*fontsize)
    sb_plot_data.set_xlabel('z($\mu m$)', fontsize = .8*fontsize, labelpad = .1)
    sb_plot_data.set_ylabel('count / '+str(binsize)+'$\; \mu m$', fontsize = .8*fontsize, labelpad = 1)
    sb_plot_data.tick_params(axis='both', which='major', labelsize = .7*fontsize, length = .8*fontsize, width = .1*fontsize, pad = 0*fontsize)
    sb_plot_data.set_ylim([0,1.1*i_count_sup.max()])
    sb_plot_data.title.set_position([.5,.9])
    [i.set_linewidth(0.1) for i in sb_plot_data.spines.values()]
    binned_z = binned_z[:-1]
    Max = find_local_max(i_count_sup.copy(),np.array(binned_z.copy()), binsize)
    if len(Max) == 0:
      Max = [1140,1070]
    if len(Max) == 1:
      Max.append(Max[0]-60)
    possible_maximums = [ z for z in Max if (z < Max[0] and Max[0]-z < 90 and Max[0]-z > 30) or z == Max[0] ]
    if len(possible_maximums) == 0:
      possible_maximums.append([1140,1070])
    
    print("   fitting data...")
    
    tmp_thick_sigmathick_rim_sigmarim, gaussians_params, fitted_func1, fitted_func2, chisquare = mydoublefit(df_z_sup, binned_z, i_count_sup, possible_maximums, plot_range_sup, sigmarble)
    
    for j in range(0,len(thick_sigmathick_rim_sigmarim)):
      thick_sigmathick_rim_sigmarim[j].append(tmp_thick_sigmathick_rim_sigmarim[j])
    if len(fitted_func1) != 0:
      sb_plot_data.plot(fitted_func1[0], fitted_func1[1], 'r-', linewidth = .5)
    if len(fitted_func2) != 0:
      sb_plot_data.plot(fitted_func2[0], fitted_func2[1], 'r-', linewidth = .5)
    fitbox = FormFitBox(gaussians_params, df_z_sup['z'], chisquare, binsize)
    sb_plot_data.text(plot_range_sup[0], 0.3*i_count_sup.max(), fitbox, horizontalalignment='left', color = 'r', fontsize = .7*fontsize)
    fig.add_subplot(sb_plot_data)
    for xz in range(0,len(l_l_all_data[k-1])):
      l_l_all_data[k-1][xz] = l_l_all_data[k-1][xz] + l_l_cut_data[i][xz]
    print("")
  fig.savefig(outdirectory + "holes_histo.pdf")
  print("   Holes histograms saved in " + outdirectory + "holes_histo.pdf")
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
      sb_plot_holes[i].set_xlabel('x($\mu m$)', fontsize=14)
    sb_plot_holes[i].set_ylabel('z($\mu m$)', fontsize=12)
  fig.savefig(outdirectory + name + "holes.pdf")
  print("   Histo saved in " + outdirectory + name + "holes.pdf")
  
def plot_thicknesses_map(thick_sigmathick_rim_sigmarim, row, col, sigmaCopperLaser = 0, outdirectory = "./"):
  thick, sigmathick, FR4, sigmaFR4, Cu, NLem, NFR4, N, sigmarble = np.array(thick_sigmathick_rim_sigmarim)
  
  plot_2D_map(thick,[],4,"map of LEM thickness. Each pixel is a measurement hole","2D_LEM_thickness_distri.pdf", row, col, 1050, 1250, outdirectory)
  #plot_2D_map(FR4,sigmaMEANFR4,5,"map of FR4 thickness. Each pixel is a measurement hole","2D_FR4_thickness_distri.pdf", row, col, 800, 1300, outdirectory)
  plot_2D_map(Cu,[],7,"map of copper thickness. Each pixel is a measurement hole","2D_rim_thickness_distri.pdf", row, col, 30,90, outdirectory)
  mean_thick = np.mean(thick)
  if mean_thick == 0:
    print("Mean LEM thickness is 0")
    return
  thick_relat = np.array([ thi - mean_thick for thi in thick ])
  #sigma_mean = math.sqrt(sum([sig**2 for sig in sigmaMEANthick]))/len(thick)
  #sigma_relat = np.array([ 100*math.sqrt(sig**2+thi**2*sigma_mean**2/mean_thick**2)/mean_thick for sig,thi in zip(sigmaMEANthick, thick) ])
  plot_2D_map(thick_relat,[],6,"map of relative LEM thickness compared to mean thickness. Each pixel is a measurement hole","2D_relat_LEM_thickness_distri.pdf", row, col, -50,50, outdirectory)
  myRangeLEM = [1050,1200]
  myRangeCu = [30,90]
  if col*row != 25:
    myRangeLEM = []
    myRangeCu = []
  if col*row != 1:
    plot_thickness_histo(thick,'LEM thickness distribution of the holes', outdirectory, 8, myRangeLEM)
    plot_thickness_histo(Cu,'Copper thickness distribution of the holes', outdirectory, 9, myRangeCu)
  
def plot_thickness_histo(thicknesses, title, outdirectory, figID, myRange):
  with open(outdirectory + title.split(' ')[0] + '.txt', 'w') as outfile:
    spamwriter = csv.writer(outfile, delimiter=';')
    spamwriter.writerow(thicknesses)
  binsize = 5 #microns
  if myRange == []:
    myRange = [min(thicknesses)*0.99,max(thicknesses)*1.01]
  nbins=int((myRange[1]-myRange[0])/binsize)
  df_thicknesses = pandas.DataFrame({'thick':thicknesses})
  statbox = FormStatBox(df_thicknesses['thick'])
  fig = plt.figure(figID,figsize=(9, 9))
  sb = fig.add_subplot(111)
  sb.set_title(title)
  sb.set_xlabel('Thickness($\mu m$)')
  sb.set_ylabel('count / '+str(binsize)+'$\; \mu m$')
  i_count, binned_thick , _ = sb.hist(df_thicknesses['thick'], bins=nbins, range = myRange)
  sb.text(sb.get_xlim()[0], 0.5*i_count.max(), statbox, horizontalalignment='left', color = 'r', fontsize = 20)
  fig.savefig(outdirectory  + title.split(' ')[0] + '.pdf')
  print("   thickness histo saved in " + outdirectory + title.split(' ')[0] + '.pdf')
  
  
def plot_2D_map(z,sigma,figID,title, savename, row, col, v0 = 800, v1 = 1300, outdirectory = "./"):
  fig = plt.figure(figID,figsize=(3*col, 3*row))
  sb_plot_2D_map = fig.add_subplot(111)
  
  x = [i%(col)+1 for i in range(0,len(z)) if z[i] != 0]
  y = [row-(int(i/col)) for i in range(0,len(z)) if z[i] != 0]
  z = z[z != 0]
  
  marker_size = good_marker_size(x,y,fig,sb_plot_2D_map)
  p = sb_plot_2D_map.scatter(x,y,c=z, s=marker_size[0]*marker_size[1], marker='.', cmap=cm.plasma, linewidth=0, vmin=v0, vmax=v1)
  sb_plot_2D_map.set_yticks([])
  sb_plot_2D_map.set_xticks([])
  sb_plot_2D_map.set_title(title)
  for i in range(0,len(z)):
    if len(sigma) == len(z):
      mytext = str(i+1) + "\n" + str(round(Decimal(z[i]),2)) + "+/-" + str(round(Decimal(sigma[i]),4))
    else:
      mytext = str(i+1) + "\n" + str(round(Decimal(z[i]),2))
    sb_plot_2D_map.annotate(mytext,xy=(x[i]-0.25,y[i]), color='white', path_effects=[PathEffects.withStroke(linewidth=2, foreground="black")], fontsize = 20)
  fig.colorbar(p)
  fig.savefig(outdirectory + savename)
  print("   2D map saved in " + outdirectory + savename)
  
def mydoublefit(df_z, binned_z, i_count, maxima, plot_range, sigmarble):
  thick_sigmathick_rim_sigmarim = []
  attempt = 1
  fit_par1 = [len(df_z)/2,maxima[0],6]
  fit_par_range1 = [0,1000*fit_par1[0],plot_range[0],plot_range[1],0,30]
  range_z1 = np.array([z for z in binned_z if z < maxima[0]+40 and z > maxima[0]-40])
  range_count1 = np.array([i for [z,i] in zip(binned_z, i_count) if z < maxima[0]+40 and z > maxima[0]-40])
  print("    Fitting whole LEM...")
  gaussians_param1, chisquare1, result1 = singlegaussfit(range_z1, range_count1, fit_par1, fit_par_range1)
  chisquare = [chisquare1]
  #if chisquare1 > 50:
    #print("   Fit failed")
    #thick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #hick_sigmathick_rim_sigmarim.append(0)
    #thick_sigmathick_rim_sigmarim.append(0)
    #return thick_sigmathick_rim_sigmarim, [], [], [], []
  z1 = gaussians_param1[1]
  sigma1 = gaussians_param1[2]
  chisquares = []
  fitted_func = []
  gaussians_params = []
  gaussians_param2 = []
  ###########while attempt < len(possible_maximums):###############
  while attempt < len(maxima):
    print("    attempting to fit FR4...")
    fit_par2 = [len(df_z)/2,maxima[attempt],12]
    fit_par_range2 = [0,1000*fit_par2[0],plot_range[0],plot_range[1],0,15]
    range_z2 = np.array([z for z in binned_z if z < maxima[attempt]+30 and z > maxima[attempt]-30])
    range_count2 = np.array([i for [z,i] in zip(binned_z, i_count) if z < maxima[attempt]+30 and z > maxima[attempt]-30])
    gaussians_param2, chisquare2, result2 = singlegaussfit(range_z2, range_count2, fit_par2, fit_par_range2)
    if chisquare2 > 0  and chisquare2 < 50 and abs(gaussians_param1[1]-gaussians_param2[1]) < 90 and abs(gaussians_param1[1]-gaussians_param2[1]) > 40:# and gaussians_param1[0]/gaussians_param2[0] < 9:
      chisquares.append(chisquare2)
      fitted_func.append(result2)
      gaussians_params.append(gaussians_param2)
    elif chisquare2 <= 0:
      print("    discarding result (R**2 < 0)")
    elif abs(gaussians_param1[1]-gaussians_param2[1]) > 90 and abs(gaussians_param1[1]-gaussians_param2[1]) < 40:
      print("    discarding result (means too close or too far away)")
    else:
      print("    discarding result (FR4 fit too small)")
    attempt = attempt + 1
  ###########while attempt < len(possible_maximums):###############
  if len(chisquares) == 0:
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
    best_fit = [ir for ir,r in enumerate(chisquares) if r == min(chisquares)][0]
    gaussians_param2 = gaussians_params[best_fit]
    NLem = df_z['z'].describe()['count']*gaussians_param1[0]/(gaussians_param1[0]+gaussians_param2[0])
    NFR4 = df_z['z'].describe()['count']*gaussians_param2[0]/(gaussians_param1[0]+gaussians_param2[0])
    chisquare.append(chisquare2)
    print("    best fit : ",gaussians_param2[1]," with Chi**2=",min(chisquares))
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

