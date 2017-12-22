from sortedcontainers import SortedDict
import numpy as np
import scipy
import math
import pandas
from decimal import Decimal
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def cut(sd_z, myrange):
  sd_cut_z = SortedDict()
  _,l_z=zip(*(sd_z.items()))
  if len(myrange) == 1:
    for x in sd_z:
      if not abs(np.mean(l_z)-sd_z[x]) > myrange[0]*math.sqrt(np.var(l_z)):
        sd_cut_z[x] = sd_z[x]
  elif len(myrange) == 2:
    for x in sd_z:
      if sd_z[x] < myrange[1] and sd_z[x] > myrange[0]:
        sd_cut_z[x] = sd_z[x]
  else:
    sd_cut_z = sd_z
  return sd_cut_z
  
def FormStatBox(df_z_selected):
  return 'count : ' + str(int(df_z_selected.describe()['count'])) + '\nmean : ' + str(round(Decimal(df_z_selected.describe()['mean']),2)) + '\nstd : ' + str(round(Decimal(df_z_selected.describe()['std']),2))
  
def gaussfit(df_z, sb):
  l_fit_range = []
  float_var_dataframe = df_z.var()
  float_mean_dataframe = df_z.mean()
  for float_z in df_z:
    if float_z > float_mean_dataframe-3*math.sqrt(float_var_dataframe) and float_z < float_mean_dataframe+3*math.sqrt(float_var_dataframe):
      l_fit_range.append(float_z)
  npa_fit_range = np.array(l_fit_range)
  npl_x_plot = np.linspace(float_mean_dataframe-3*math.sqrt(float_var_dataframe), float_mean_dataframe+3*math.sqrt(float_var_dataframe), 1000)
  float_m, float_s = scipy.stats.norm.fit(npa_fit_range)
  npa_result = df_z.count()*scipy.stats.norm.pdf(npl_x_plot, float_m, float_s)
  sb.plot(npl_x_plot, npa_result)


def doublegaussfit(x,proba,par):
  gauss_fit = lambda p, x: p[0]*(1/np.sqrt(2*math.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2))+p[3]*(1/np.sqrt(2*math.pi*(p[5]**2)))*np.exp(-(x-p[4])**2/(2*p[5]**2)) #1d Gaussian func
  e_gauss_fit = lambda p, x, y: (gauss_fit(p,x) -y) #1d Gaussian fit
  out = leastsq(e_gauss_fit, par[:], args=(x, proba), maxfev=100000, full_output=1)
  xxx = np.arange(min(x),max(x),x[1]-x[0])
  ccc = gauss_fit(par,xxx) # this will only work if the units are pixel and not wavelength
  fig = plt.figure(0,figsize=(9, 9)) #make a plot
  ax1 = fig.add_subplot(111)
  ax1.plot(x,proba,'gs') #spectrum
  ax1.plot(xxx,ccc,'b-') #fitted spectrum
  plt.show()
