import numpy as np
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas
import scipy
import os
import matplotlib.gridspec as gridspec
from toolbox import *
import re

def marble(l_l_cut_data, file_marble_file = [], ID = 1, file_cale_file = [], outdirectory = "./"):
  l_marble_data = [[],[]]
  if ID > len(file_marble_file):
    print("ERROR: should have as many marble measurements than holes measurements")
  print("  Opening marble file " + file_marble_file[ID-1] + "...")
  #the rows of the real data and the rows of the marble data should correspond
  marble_file = open(file_marble_file[ID-1], "r")
  content_marble_file = [ line for line in marble_file.readlines() if not "Distance" in line ]
  marble_file.close()
  #beginning of measurment usually is crap
  #startcut = -900000
  #remove header line
  #read the file
  nx = 0
  do_fit = True
  for line in content_marble_file:
    x = float(line.split("\n")[0].split(";")[4])
   # if x < startcut:
    #  continue
    if not x in l_marble_data[0]:
      nx = nx + 1
    l_marble_data[0].append(x)
    l_marble_data[1].append(float(line.split("\n")[0].split(";")[0]))
  if nx < 10:
    do_fit = False
  #remove the marble data too far from the mean
  if do_fit:
    l_marble_data = cut(l_marble_data,[2])
  l_l_cut_data, lr_fit, sigmarble = marble_fit(l_marble_data, l_l_cut_data, ID, outdirectory, do_fit)
  if len(file_cale_file) != 0:
    print("  Plotting cale data...")
    plot_cale(file_cale_file, outdirectory)
  return l_l_cut_data, sigmarble
    
def plot_cale(file_cale_file, outdirectory = "./"):
  i = 0
  cale = 0
  sb_plot_cale = []
  marble_ref = []
  l_cale_data = []
  means = []
  title = []
  binsize = 5
  is_marble_ref = False
  marble_file_position = [ i for i,fname in enumerate(file_cale_file) if "marbre" in os.path.basename(fname) ]
  if len(marble_file_position) > 1:
    print("   ERROR: found several marble files amoung cale files. I do not know which one to use for reference.")
    exit(1)
  elif len(marble_file_position) == 1:
    marble_file_position = marble_file_position[0]
    is_marble_ref = True
    print("   Found marble reference file for cale")
    cale_file = open(file_cale_file[marble_file_position], "r")
    content_cale_file = [ line for line in cale_file.readlines() if not "Distance" in line ]
    cale_file.close()
    marble_ref = np.array([float(line.split("\n")[0].split(";")[0]) for line in content_cale_file])
  for i in range(0,len(file_cale_file)):
    print("   Opening cale file " + file_cale_file[i] + "...")
    cale_file = open(file_cale_file[i], "r")
    content_cale_file = [ line for line in cale_file.readlines() if not "Distance" in line ]
    cale_file.close()
    if is_marble_ref and i == marble_file_position:
      plot_title = "marble_ref"
    else:
      plot_title = os.path.basename(file_cale_file[i])
    if is_marble_ref:
      l_cale_data = [float(float(line.split("\n")[0].split(";")[0])-marble_ref.mean()) for line in content_cale_file if float(line.split("\n")[0].split(";")[0]) > 0]
    else:
      l_cale_data = [ float(line.split("\n")[0].split(";")[0]) for line in content_cale_file if float(line.split("\n")[0].split(";")[0]) > 0 ]
    if max(l_cale_data)-min(l_cale_data) > 500:
      l_cale_data = [data for data in l_cale_data if data > (max(l_cale_data)-min(l_cale_data))/2 + min(l_cale_data)]
      l_cale_data = cut(l_cale_data,[4])
    nbin = int((max(l_cale_data) - min(l_cale_data))/binsize)
    if nbin < 50:
      nbin = 'auto'
    df_z = pandas.DataFrame({'z':l_cale_data})
    fig = plt.figure(100+i,figsize=(12, 12))
    #sb_plot_cale.append(fig.add_subplot(211))
    #sb_plot_cale[-1].set_xlabel("time (ms)")
    #sb_plot_cale[-1].set_ylabel("z (micrometer)")
    #sb_plot_cale[-1].set_ylim([0.9*min(l_cale_data),1.1*max(l_cale_data)])
    #sb_plot_cale[-1].set_title("z vs time")
    #sb_plot_cale[-1].plot([ [i,x] for i,x in enumerate(l_cale_data)], 'g-')
    sb_plot_cale.append(fig.add_subplot(111))
    i_count, binned_z, _ = sb_plot_cale[-1].hist(df_z['z'], bins = nbin)
    binned_z = binned_z[:-1]
    maxz = [l_cale_data[i] for i,c in enumerate(i_count) if c == i_count.max()][0]
    fit_par = [len(df_z), maxz, math.sqrt(float(df_z.var()))]
    gaussians_param, rsquare, result = singlegaussfit(binned_z, i_count, fit_par)
    if len([param for param in gaussians_param if param > 1e20]) != 0:
      print("   fit failed, parameter too high :'(")
      gaussians_param = []
    sb_plot_cale[-1].plot(result[0], result[1], 'r-')
    fitbox = FormFitBox(gaussians_param, df_z['z'], [rsquare], binsize)
    sb_plot_cale[-1].text(df_z['z'].min(), 0.5*i_count.max(), fitbox, horizontalalignment='left', fontsize = 36)
    sb_plot_cale[-1].set_xlabel("Height ($\mu m$)")
    sb_plot_cale[-1].set_ylabel('Count / %.2g $\; \mu m$' % binsize)
    sb_plot_cale[-1].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    sb_plot_cale[-1].set_title("Height of " + plot_title + " reference")
    fig.savefig(outdirectory + "cale_" + str(plot_title.replace('.csv','.pdf')))
    print("   ...plot saved in " + outdirectory + "cale_" + str(plot_title) + ".pdf")
    plt.close('all')
    if i == marble_file_position:
      cale = 0
    else:
      cale = re.findall(r'\d+',plot_title)
      if len(cale) == 2:
        cale = float(cale[0] + "." + cale[1])
      elif len(cale) == 1:
        cale = float(cale[0])
      else:
        print(   "ERROR: wrong cale filename format")
        exit(1)
    title.append(cale)
    means.append(np.mean(df_z['z']))
  title = np.array(title)
  means = np.array(means)
  print("   Computing laser slope...")
  lr_fit = linear_model.LinearRegression()
  lr_fit.fit(title[:,np.newaxis], means)
  print("   cale cali plot saved in " + outdirectory + "cale_cale_cali.pdf")
  npl_x_test = np.linspace(np.min(title), np.max(title), 100)
  plt.plot(npl_x_test, lr_fit.predict(npl_x_test[:,np.newaxis]), color='blue', linewidth=3)
  line_equation = 'z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2))
  plt.text(3*(title.max()-title.min())/4,1.005*means.max(),line_equation,horizontalalignment='center', color='r', fontsize = 36)
  plt.scatter(title,means)
  plt.savefig(outdirectory + "cale_cali.pdf")
  plt.close()

def marble_fit(l_marble_data, l_l_cut_data = None, ID = 1, outdirectory = "./", do_fit = True):
  l_l_corrected_data = []
  npa_marble_x = np.array(l_marble_data[0])
  npa_marble_z = np.array(l_marble_data[1])
  if do_fit:
    lr_fit = linear_model.LinearRegression()
    print("   Computing marble slope...")
    lr_fit.fit(npa_marble_x[:,np.newaxis], npa_marble_z)
    print('  ...z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2)))
    print("   substracting marble from data...")
    if not l_l_cut_data is None:
      for l_cut_data in l_l_cut_data:
        if l_cut_data == [[],[]]:
          l_l_corrected_data.append([[],[]])
          continue
        lx,lz = zip(*[ [x,(z - lr_fit.predict(x))[0]] for x,z in zip(l_cut_data[0],l_cut_data[1]) ])
        l_l_corrected_data.append([list(lx),list(lz)])
    print("   plotting marble data...")
    sigmarble = plot_fit(npa_marble_x, npa_marble_z, lr_fit, ID, outdirectory)
    return [l_l_corrected_data, lr_fit, sigmarble]
  else:
    if not l_l_cut_data is None:
      print("   substracting marble from data...")
      for l_cut_data in l_l_cut_data:
        lx,lz = zip(*[ [x,(z - npa_marble_z.mean())] for x,z in zip(l_cut_data[0],l_cut_data[1]) ])
        l_l_corrected_data.append([list(lx),list(lz)])
    return [l_l_corrected_data, None, 0]


def plot_fit(npa_x, npa_z, lr_fit, ID = 1, outdirectory = "./"):
  npl_x_test = np.linspace(np.min(npa_x), np.max(npa_x), 100)
  #df_z = pandas.DataFrame({'z':npa_z})
  npa_z_corr = [(z - lr_fit.predict(x))[0] for x,z in zip(npa_x, npa_z)]
  df_z_corr = pandas.DataFrame({'z':npa_z_corr})
  fig = plt.figure(10+ID,figsize=(12, 12))
  #plt.tight_layout()
  sb1 = fig.add_subplot(211)
  sb1.set_title("Marble Height vs length (z vs x, z histo before and after slope corr.)")
  sb1.set_xlabel('x($\mu m$)', fontsize=14)
  sb1.set_ylabel('z($\mu m$)', fontsize=12)
  sb1.scatter(npa_x, npa_z, color='black')
  sb1.axis([min(npa_x)-1., max(npa_x)+1., 0.99*min(npa_z), max(npa_z)*1.01])
  sb1.plot(npl_x_test, lr_fit.predict(npl_x_test[:,np.newaxis]), color='blue', linewidth=3)
  line_equation = 'z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2))
  plt.text(3*(npa_x.max()-npa_x.min())/4,1.005*npa_z.max(),line_equation,horizontalalignment='center', color='r', fontsize = 36)
  #sb2 = fig.add_subplot(312)
  #sb2.set_xlabel('z(micrometer)', fontsize=14)
  #sb2.set_ylabel('count', fontsize=12)
  #i_count, _, _ = sb2.hist(df_z['z'], bins='auto')
  #statbox = FormStatBox(df_z['z'])
  #sb2.text(df_z['z'].min(), 0.9*i_count.max(), statbox, horizontalalignment='left')
  sb3 = fig.add_subplot(212)
  sb3.set_xlabel('z($\mu m$)', fontsize=14)
  i_count, binned_z, _ = sb3.hist(df_z_corr['z'], bins='auto')
  binsize = binned_z[1]-binned_z[0]
  sb3.set_ylabel('Count / %.2g $\; \mu m$' % binsize, fontsize=12)
  binned_z = binned_z[:-1]
  fit_par = [len(df_z_corr), 0, math.sqrt(float(df_z_corr.var()))]
  gaussians_param, rsquare, result = singlegaussfit(binned_z, i_count, fit_par)
  sb3.plot(result[0], result[1], 'r-')
  fitbox = FormFitBox(gaussians_param, df_z_corr['z'], [rsquare], binsize)
  sb3.text(df_z_corr['z'].min(), 0.5*i_count.max(), fitbox, horizontalalignment='left', fontsize = 24, color = 'r')
  fig.savefig(outdirectory + "marble_" + str(ID) + ".pdf",bbox_inches = "tight")
  print("   ...plot saved in marble_" + outdirectory + "marble_" + str(ID) + ".pdf")
  return gaussians_param[2]
  plt.close('all')
  
def plot_other_marble_file(other_marble_file, outdirectory = "./"):
  sb_plot_marble = []
  plt.subplots_adjust(left=0.01, bottom=0.1, right=0.99, top=0.99, wspace=0.2, hspace=0.2)
  binsize = 5
  sigma = 0
  for i in range(0,len(other_marble_file)):
    fig = plt.figure(1000+i,figsize=(5, 8))
    File = other_marble_file[i]
    l_marble_data = [[],[]]
    print("  Opening other marble file " + File + "...")
  #the rows of the real data and the rows of the marble data should correspond
    marble_file = open(File, "r")
    content_marble_file = [ line for line in marble_file.readlines() if not "Distance" in line ]
    marble_file.close()
    name = os.path.basename(File).replace(".csv",".pdf")
    #beginning of measurment usually is crap
    #startcut = -900000
    #remove header line
    if "Distance" in content_marble_file[0]:
      content_marble_file = content_marble_file[1:]
    #read the file
    for line in content_marble_file:
      #if float(line.split("\n")[0].split(";")[4])<startcut:
        #continue
      l_marble_data[0].append(float(line.split("\n")[0].split(";")[4]))
      l_marble_data[1].append(float(line.split("\n")[0].split(";")[0]))
    print("   plotting marble data...")
    l_marble_data = cut(l_marble_data,[4])

    nbin = int((max(l_marble_data[1]) - min(l_marble_data[1]))/binsize)
    if nbin < 50:
      nbin = 'auto'
    plot_title = os.path.basename(File)
    df_z = pandas.DataFrame({'z':l_marble_data[1]})
    #sb_plot_marble.append(fig.add_subplot(211))
    #sb_plot_marble[-1].set_xlabel("time (ms)")
    #sb_plot_marble[-1].set_ylabel("z (micrometer)")
    #sb_plot_marble[-1].set_ylim([0.9*min(l_marble_data[1]),1.1*max(l_marble_data[1])])
    #sb_plot_marble[-1].set_title("z vs time")
    #sb_plot_marble[-1].plot([ [i,x] for i,x in enumerate(l_marble_data[1])], 'g-')
    sb_plot_marble.append(fig.add_subplot(111))
    i_count, binned_z, _ = sb_plot_marble[-1].hist(df_z['z'], bins='auto')
    binned_z = binned_z[:-1]
    maxz = [l_marble_data[1][i] for i,c in enumerate(i_count) if c == i_count.max()][0]
    fit_par = [len(df_z), maxz, math.sqrt(float(df_z.var()))]
    gaussians_param, rsquare, result = singlegaussfit(binned_z, i_count, fit_par)
    sb_plot_marble[-1].plot(result[0], result[1], 'r-')
    fitbox = FormFitBox(gaussians_param, df_z['z'], [rsquare], binsize)
    if "copper" or "cuivre" in name:
      sigma = gaussians_param[2]
    sb_plot_marble[-1].text(df_z['z'].min(), 0.5*i_count.max(), fitbox, horizontalalignment='left', fontsize = 36, color = 'r')
    sb_plot_marble[-1].set_xlabel("Height ($\mu m$)")
    sb_plot_marble[-1].set_ylabel('Count / %.2g $\; \mu m$' % binsize)
    sb_plot_marble[-1].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    sb_plot_marble[-1].set_title(plot_title)
    fig.savefig(outdirectory + name)
    print("   ...plot saved in " + outdirectory + name)
  plt.close()
  return sigma
    
