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

def marble(l_l_cut_data, file_marble_file = [], ID = 1, file_calle_file = []):
  l_marble_data = [[],[]]
  if ID > len(file_marble_file):
    print("ERROR: should have as many marble measurements than holes measurements")
  print("  Opening marble file " + file_marble_file[ID-1] + "...")
  #the rows of the real data and the rows of the marble data should correspond
  marble_file = open(file_marble_file[ID-1], "r")
  content_marble_file = marble_file.readlines()
  marble_file.close()
  #beginning of measurment usually is crap
  startcut = -900000
  #remove header line
  if "Distance" in content_marble_file[0]:
    content_marble_file = content_marble_file[1:]
  #read the file
  for line in content_marble_file:
    if float(line.split("\n")[0].split(";")[4])<startcut:
      continue
    l_marble_data[0].append(float(line.split("\n")[0].split(";")[4]))
    l_marble_data[1].append(float(line.split("\n")[0].split(";")[0]))
  #remove the marble data too far from the mean

  l_marble_data = cut(l_marble_data,[4])
  
  l_l_cut_data, lr_fit = marble_fit(l_marble_data, l_l_cut_data, ID)
  if len(file_calle_file) != 0:
    print("  Plotting calle data...")
    plot_calle(file_calle_file)
  return l_l_cut_data
    
def plot_calle(file_calle_file):
  i = 0
  sb_plot_calle = []
  marble_ref = []
  l_calle_data = []
  binsize = 5
  is_marble_ref = False
  if "marbre" in file_calle_file[0]:
    is_marble_ref = True
    calle_file = open(file_calle_file[0], "r")
    content_calle_file = calle_file.readlines()
    calle_file.close()
    if "Distance" in content_calle_file[0]:
      content_calle_file = content_calle_file[1:]
    marble_ref = np.array([float(line.split("\n")[0].split(";")[0]) for line in content_calle_file])
  for i in range(0,len(file_calle_file)):
    if is_marble_ref and i == 0:
      plot_title = "marble_ref"
      l_calle_data = marble_ref
    else:
      plot_title = os.path.basename(file_calle_file[i])
      print("   Opening calle file " + file_calle_file[i] + "...")
      calle_file = open(file_calle_file[i], "r")
      content_calle_file = calle_file.readlines()
      calle_file.close()
      if "Distance" in content_calle_file[0]:
        content_calle_file = content_calle_file[1:]
      if is_marble_ref:
        l_calle_data = [float(float(line.split("\n")[0].split(";")[0])-marble_ref.mean()) for line in content_calle_file if float(line.split("\n")[0].split(";")[0]) > 0]
      else:
        l_calle_data = [float(line.split("\n")[0].split(";")[0]) for line in content_calle_file if float(line.split("\n")[0].split(";")[0]) > 0]
    if max(l_calle_data)-min(l_calle_data) > 500:
      l_calle_data = [data for data in l_calle_data if data > (max(l_calle_data)-min(l_calle_data))/2 + min(l_calle_data)]
      l_calle_data = cut(l_calle_data,[4])
      print(max(l_calle_data))
    nbin = int((max(l_calle_data) - min(l_calle_data))/binsize)
    if nbin < 50:
      nbin = 'auto'
    df_z = pandas.DataFrame({'z':l_calle_data})
    fig = plt.figure(100+i,figsize=(12, 12))
    #sb_plot_calle.append(fig.add_subplot(211))
    #sb_plot_calle[-1].set_xlabel("time (ms)")
    #sb_plot_calle[-1].set_ylabel("z (micrometer)")
    #sb_plot_calle[-1].set_ylim([0.9*min(l_calle_data),1.1*max(l_calle_data)])
    #sb_plot_calle[-1].set_title("z vs time")
    #sb_plot_calle[-1].plot([ [i,x] for i,x in enumerate(l_calle_data)], 'g-')
    sb_plot_calle.append(fig.add_subplot(111))
    i_count, binned_z, _ = sb_plot_calle[-1].hist(df_z['z'], bins = nbin)
    binned_z = binned_z[:-1]
    maxz = [l_calle_data[i] for i,c in enumerate(i_count) if c == i_count.max()][0]
    fit_par = [len(df_z), maxz, math.sqrt(float(df_z.var()))]
    gaussians_param, rsquare, result = singlegaussfit(binned_z, i_count, fit_par)
    if len([param for param in gaussians_param if param > 1e20]) != 0:
      print("   fit failed, parameter too high :'(")
      gaussians_param = []
    sb_plot_calle[-1].plot(result[0], result[1], 'r-')
    fitbox = FormFitBox(gaussians_param, df_z['z'])
    sb_plot_calle[-1].text(df_z['z'].min(), 0.9*i_count.max(), fitbox, horizontalalignment='left')
    sb_plot_calle[-1].set_xlabel("Height (micrometer)")
    sb_plot_calle[-1].set_ylabel("Count")
    sb_plot_calle[-1].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    sb_plot_calle[-1].set_title("Height of " + plot_title + " reference")
    fig.savefig("calle_" + str(plot_title) + ".pdf")
    plt.close()
    print("   ...plot saved in calle_" + str(plot_title) + ".pdf")


def marble_fit(l_marble_data, l_l_cut_data = None, ID = 1):
  npa_marble_x = np.array(l_marble_data[0])
  npa_marble_z = np.array(l_marble_data[1])
  lr_fit = linear_model.LinearRegression()
  print("   Computing marble slope...")
  lr_fit.fit(npa_marble_x[:,np.newaxis], npa_marble_z)
  print('  ...z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2)))
  l_l_corrected_data = []
  print("   Applying fit to data...")
  if not l_l_cut_data is None:
    for l_cut_data in l_l_cut_data:
      l_l_corrected_data.append([[],[]])
      for x,z in zip(l_cut_data[0],l_cut_data[1]):
        l_l_corrected_data[-1][0].append(x)
        l_l_corrected_data[-1][1].append((z - lr_fit.predict(x))[0])
  print("   plotting marble data...")
  plot_fit(npa_marble_x, npa_marble_z, lr_fit, ID)
  return [l_l_corrected_data, lr_fit]


def plot_fit(npa_x, npa_z, lr_fit, ID = 1):
  npl_x_test = np.linspace(np.min(npa_x), np.max(npa_x), 100)
  #df_z = pandas.DataFrame({'z':npa_z})
  npa_z_corr = [(z - lr_fit.predict(x))[0] for x,z in zip(npa_x, npa_z)]
  df_z_corr = pandas.DataFrame({'z':npa_z_corr})
  fig = plt.figure(10+ID,figsize=(12, 12))
  #plt.tight_layout()
  sb1 = fig.add_subplot(211)
  sb1.set_title("Marble Height vs length (z vs x, z histo before and after slope corr.)")
  sb1.set_xlabel('x(micrometer)', fontsize=14)
  sb1.set_ylabel('z(micrometer)', fontsize=12)
  sb1.scatter(npa_x, npa_z, color='black')
  sb1.axis([min(npa_x)-1., max(npa_x)+1., 0.99*min(npa_z), max(npa_z)*1.01])
  sb1.plot(npl_x_test, lr_fit.predict(npl_x_test[:,np.newaxis]), color='blue', linewidth=3)
  line_equation = 'z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2))
  plt.text(3*(npa_x.max()-npa_x.min())/4,1.005*npa_z.max(),line_equation,horizontalalignment='center', color='r')
  #sb2 = fig.add_subplot(312)
  #sb2.set_xlabel('z(micrometer)', fontsize=14)
  #sb2.set_ylabel('count', fontsize=12)
  #i_count, _, _ = sb2.hist(df_z['z'], bins='auto')
  #statbox = FormStatBox(df_z['z'])
  #sb2.text(df_z['z'].min(), 0.9*i_count.max(), statbox, horizontalalignment='left')
  sb3 = fig.add_subplot(212)
  sb3.set_xlabel('z(micrometer)', fontsize=14)
  sb3.set_ylabel('count', fontsize=12)
  i_count, binned_z, _ = sb3.hist(df_z_corr['z'], bins='auto')
  binned_z = binned_z[:-1]
  fit_par = [len(df_z_corr), 0, math.sqrt(float(df_z_corr.var()))]
  gaussians_param, rsquare, result = singlegaussfit(binned_z, i_count, fit_par)
  sb3.plot(result[0], result[1], 'r-')
  fitbox = FormFitBox(gaussians_param, df_z_corr['z'])
  sb3.text(df_z_corr['z'].min(), 0.9*i_count.max(), fitbox, horizontalalignment='left')
  fig.savefig("marble_" + str(ID) + ".pdf",bbox_inches = "tight")
  plt.close()
  print("   ...plot saved in marble_" + str(ID) + ".pdf")
  
def plot_other_marble_file(other_marble_file):
  sb_plot_marble = []
  plt.subplots_adjust(left=0.01, bottom=0.1, right=0.99, top=0.99, wspace=0.2, hspace=0.2)
  binsize = 5
  for i in range(0,len(other_marble_file)):
    fig = plt.figure(1000+i,figsize=(5, 8))
    File = other_marble_file[i]
    l_marble_data = [[],[]]
    print("  Opening other marble file " + File + "...")
  #the rows of the real data and the rows of the marble data should correspond
    marble_file = open(File, "r")
    content_marble_file = marble_file.readlines()
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
    fitbox = FormFitBox(gaussians_param, df_z['z'])
    sb_plot_marble[-1].text(df_z['z'].min(), 0.9*i_count.max(), fitbox, horizontalalignment='left')
    sb_plot_marble[-1].set_xlabel("Height (micrometer)")
    sb_plot_marble[-1].set_ylabel("Count")
    sb_plot_marble[-1].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    sb_plot_marble[-1].set_title(plot_title)
    fig.savefig(name)
  plt.close()
  print("   ...plot saved in other_marble.pdf")
    
