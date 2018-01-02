import numpy as np
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas
import scipy
from toolbox import *

def marble(l_sd_data, file_marble_file = [], ID = 1, file_calle_file = []):
  l_sd_marble_data = [SortedDict()]
  if len(file_marble_file) != 0:
    if ID > len(file_marble_file):
      print("ERROR: should have as many marble measurements than holes measurements")
    marble_file = open(file_marble_file[ID-1], "r")
    content_marble_file = marble_file.readlines()
    marble_file.close()
    startcut = -900000
    if "Distance" in content_marble_file[0]:
      content_marble_file = content_marble_file[1:]
    for line in content_marble_file:
      if float(line.split("\n")[0].split(";")[4])<startcut:
        continue
      l_sd_marble_data[0][float(line.split("\n")[0].split(";")[4])]=float(line.split("\n")[0].split(";")[0])
    #l_sd_marble_data[0]=cut(l_sd_marble_data[0],[4])
  else:
    l_sd_marble_data[0]=cut(l_sd_data[0],[4])
  if len(l_sd_data) >1 and file_marble_file == "":
    l_sd_marble_data.append(SortedDict())
    l_sd_marble_data[1] = cut(l_sd_data[-1],[4])
        
  if len(l_sd_data) >2:
    l_sd_data, lr_fit = marble_fit(l_sd_marble_data, l_sd_data, ID)
    if len(file_calle_file) != 0:
      plot_calle(file_calle_file, lr_fit)
    return l_sd_data
  else:
    _,lr_fit = marble_fit(l_sd_marble_data, ID)
    if len(file_calle_file) != 0:
      plot_calle(file_calle_file, lr_fit)
    
def plot_calle(file_calle_file,lr_fit):
  i = 0
  for File in file_calle_file:
    plot_title = File.replace("../data/calle/calle","")
    plot_title = plot_title.replace("_","")
    plot_title = plot_title.replace("1000Hz","")
    plot_title = plot_title.replace(".csv","")
    calle_file = open(File, "r")
    content_calle_file = calle_file.readlines()
    calle_file.close()
    if "Distance" in content_calle_file[0]:
      content_calle_file = content_calle_file[1:]
    sd_calle_data = [float(float(line.split("\n")[0].split(";")[0])-lr_fit.predict(float(line.split("\n")[0].split(";")[0]))) for line in content_calle_file]
    #-lr_fit.predict(float(line.split("\n")[0].split(";")[0]))
    df_z = pandas.DataFrame({'z':sd_calle_data})
    fig = plt.figure(5,figsize=(12, 12))
    sb = fig.add_subplot(111)
    i_count, _, _ = sb.hist(df_z['z'], bins='auto')
    statbox = FormStatBox(df_z['z'])
    sb.text(df_z['z'].min(), 0.9*i_count.max(), statbox,horizontalalignment='left')
    sb.set_xlabel("Height micrimeter)")
    sb.set_ylabel("Count")
    sb.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    sb.set_title("Height of " + plot_title + " reference")
    plt.savefig("calle_" + str(plot_title) + ".pdf")
    plt.clf()
    i = i + 1


def marble_fit(l_sd_marble_data, l_sd_data = None, ID = 1):
  sd_marble = l_sd_marble_data[0].copy()
  if len(l_sd_marble_data) == 2:
    sd_marble.update(l_sd_marble_data[1])
  l_marble_x, l_marble_z = zip(*(sd_marble.items()))
  npa_marble_x = np.array(l_marble_x)
  npa_marble_z = np.array(l_marble_z)
  lr_fit = linear_model.LinearRegression()
  lr_fit.fit(npa_marble_x[:,np.newaxis], npa_marble_z)
  if not l_sd_data is None:
    l_sd_corrected_data = [SortedDict()]
    for sd_data in l_sd_data[1:-1]:
      for float_x in sd_data:
        l_sd_corrected_data[-1][float_x] = sd_data[float_x]-lr_fit.predict(sd_data[float_x])[0]
      l_sd_corrected_data.append(SortedDict())
    l_sd_corrected_data = l_sd_corrected_data[:-1]
    l_sd_data = l_sd_corrected_data
  
  plot_fit([npa_marble_x, npa_marble_z,l_sd_marble_data], lr_fit, ID)
  return [l_sd_data, lr_fit]


def plot_fit(l_marble_xzdata, lr_fit, ID = 1):
  npa_x = l_marble_xzdata[0]
  npa_z = l_marble_xzdata[1]
  npl_x_test = np.linspace(np.min(npa_x), np.max(npa_x), 100)
  df_z = [pandas.DataFrame({'z':list(l_marble_xzdata[2][0].values())})]
  fig = plt.figure(0,figsize=(12, 12))
  #plt.tight_layout()
  sb1 = fig.add_subplot(211)
  sb1.set_title("Marble Height vs length")
  sb1.set_xlabel('x(micrometer)', fontsize=14)
  sb1.set_ylabel('z(micrometer)', fontsize=12)
  sb1.scatter(npa_x, npa_z, color='black')
  sb1.plot(npl_x_test, lr_fit.predict(npl_x_test[:,np.newaxis]), color='blue', linewidth=3)
  line_equation = 'z(x)=' + str('%.2E' % Decimal(lr_fit.coef_[0])) + 'x+' + str(round(Decimal(lr_fit.intercept_),2))
  plt.text((npa_x.max()-npa_x.min())/2,npa_z.max()-(npa_z.max()-npa_z.min())/4,line_equation,horizontalalignment='center')
  if len(l_marble_xzdata[2]) == 2:
    df_z.append(pandas.DataFrame({'z':list(l_marble_xzdata[2][1].values())}))
    sb2 = fig.add_subplot(223)
    sb2.set_xlabel('z(micrometer)', fontsize=14)
    sb2.set_ylabel('count', fontsize=12)
    i_count, _ , _ = sb2.hist(df_z[0]['z'], bins='auto')
    #gaussfit(df_z[0]['z'], sb2)
    statbox = FormStatBox(df_z[0]['z'])
    sb2.text(df_z[0]['z'].min(), 0.9*i_count.max(), statbox,horizontalalignment='left')
    sb3 = fig.add_subplot(224)
    sb3.set_xlabel('z(micrometer)', fontsize=14)
    sb3.set_ylabel('count', fontsize=12)
    i_count, _, _ = sb3.hist(df_z[1]['z'], bins='auto')
    #gaussfit(df_z[1]['z'], sb3)
    statbox = FormStatBox(df_z[1]['z'])
    sb3.text(df_z[1]['z'].min(), 0.9*i_count.max(), statbox,horizontalalignment='left')
  elif len(l_marble_xzdata[2]) == 1:
    sb2 = fig.add_subplot(212)
    sb2.set_xlabel('z(micrometer)', fontsize=14)
    sb2.set_ylabel('count', fontsize=12)
    i_count, _, _ = sb2.hist(df_z[0]['z'], bins='auto')
    statbox = FormStatBox(df_z[0]['z'])
    sb2.text(df_z[0]['z'].min(), 0.9*i_count.max(), statbox,horizontalalignment='left')
  else:
    print("marble_fit::plot_fit ERROR : wrong size for marble z data")
    exit(1)
  plt.savefig("marble_" + str(ID) + ".pdf",bbox_inches = "tight")
  plt.clf()
  #plt.show()
