import numpy as np
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas
import scipy
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
    #note that sorted dict can not have twice the same key (here, the x coordinate). A workaround has been put in place for the real data (see main.py), but not here.
    l_marble_data[0].append(float(line.split("\n")[0].split(";")[4]))
    l_marble_data[1].append(float(line.split("\n")[0].split(";")[0]))
  #remove the marble data too far from the mean

  l_marble_data = cut(l_marble_data,[4])
  
  l_l_cut_data, lr_fit = marble_fit(l_marble_data, l_l_cut_data, ID)
  if len(file_calle_file) != 0:
    print("  Plotting calle data...")
    plot_calle(file_calle_file, lr_fit)
  return l_l_cut_data
    
def plot_calle(file_calle_file,lr_fit):
  i = 0
  sb_plot_calle = []
  marble_ref = []
  l_calle_data = []
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
      plot_title = file_calle_file[i].replace("1000Hz","")
      plot_title = remove_letters(plot_title)
      plot_title = plot_title.replace("..","")
      plot_title = plot_title.replace("2017","")
      plot_title = plot_title.replace("2018","")
      plot_title = plot_title.replace("/","")
      plot_title = plot_title.replace("_","")
      plot_title = plot_title.replace(".csv","")
      print("   Opening calle file " + file_calle_file[i] + "...")
      calle_file = open(file_calle_file[i], "r")
      content_calle_file = calle_file.readlines()
      calle_file.close()
      if "Distance" in content_calle_file[0]:
        content_calle_file = content_calle_file[1:]
      if is_marble_ref:
        l_calle_data = [float(float(line.split("\n")[0].split(";")[0])-marble_ref.mean()) for line in content_calle_file]
      else:
        l_calle_data = [float(float(line.split("\n")[0].split(";")[0])-lr_fit.predict(float(line.split("\n")[0].split(";")[0]))) for line in content_calle_file]
      #-lr_fit.predict(float(line.split("\n")[0].split(";")[0]))
    df_z = pandas.DataFrame({'z':l_calle_data})
    fig = plt.figure(100+i,figsize=(12, 12))
    sb_plot_calle.append(fig.add_subplot(111))
    i_count, _, _ = sb_plot_calle[-1].hist(df_z['z'], bins='auto')
    statbox = FormStatBox(df_z['z'])
    sb_plot_calle[-1].text(df_z['z'].min(), 0.9*i_count.max(), statbox,horizontalalignment='left')
    sb_plot_calle[-1].set_xlabel("Height (micrometer)")
    sb_plot_calle[-1].set_ylabel("Count")
    sb_plot_calle[-1].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    sb_plot_calle[-1].set_title("Height of " + plot_title + " reference")
    fig.savefig("calle_" + str(plot_title) + ".pdf")
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
  plot_fit(npa_marble_x, npa_marble_z,l_marble_data, lr_fit, ID)
  return [l_l_corrected_data, lr_fit]


def plot_fit(npa_x, npa_z, l_xz, lr_fit, ID = 1):
  npl_x_test = np.linspace(np.min(npa_x), np.max(npa_x), 100)
  df_z = pandas.DataFrame({'z':l_xz[1]})
  fig = plt.figure(10+ID,figsize=(12, 12))
  #plt.tight_layout()
  sb1 = fig.add_subplot(211)
  sb1.set_title("Marble Height vs length")
  sb1.set_xlabel('x(micrometer)', fontsize=14)
  sb1.set_ylabel('z(micrometer)', fontsize=12)
  sb1.scatter(npa_x, npa_z, color='black')
  sb1.axis([min(npa_x)-1., max(npa_x)+1., 0.99*min(npa_z), max(npa_z)*1.01])
  sb1.plot(npl_x_test, lr_fit.predict(npl_x_test[:,np.newaxis]), color='blue', linewidth=3)
  line_equation = 'z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2))
  plt.text(3*(npa_x.max()-npa_x.min())/4,1.005*npa_z.max(),line_equation,horizontalalignment='center', color='r')
  sb2 = fig.add_subplot(212)
  sb2.set_xlabel('z(micrometer)', fontsize=14)
  sb2.set_ylabel('count', fontsize=12)
  i_count, _, _ = sb2.hist(df_z['z'], bins='auto')
  statbox = FormStatBox(df_z['z'])
  sb2.text(df_z['z'].min(), 0.9*i_count.max(), statbox,horizontalalignment='left')
  fig.savefig("marble_" + str(ID) + ".pdf",bbox_inches = "tight")
  print("   ...plot saved in marble_" + str(ID) + ".pdf")
  

