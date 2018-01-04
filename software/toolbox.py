from sortedcontainers import SortedDict
import numpy as np
import scipy
import math
import pandas
from decimal import Decimal
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import re

def cut(sd_z, myrange):
  if len(sd_z) == 0:
    return sd_z
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
  return 'count : ' + str(int(df_z_selected.describe()['count'])) + '\nmean : ' + str(round(Decimal(df_z_selected.describe()['mean']),0)) + '\nstd : ' + str(round(Decimal(df_z_selected.describe()['std']),2))
  
def FormFitBox(param, df_z_selected):
  return 'count : ' + str(int(df_z_selected.describe()['count'])) + '\nFit result:\n First gaussian:\n  mean=' + str(round(Decimal(param[1]),2)) + '\n  sigma=' + str(abs(round(Decimal(param[2]),2))) + '\n Second gaussian:\n  mean=' + str(round(Decimal(param[4]),2)) + '\n  sigma=' + str(abs(round(Decimal(param[5]),2))) 
  
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

def gauss_fit(p,x):
  return p[0]*(1/np.sqrt(2*math.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2))+p[3]*(1/np.sqrt(2*math.pi*(p[5]**2)))*np.exp(-(x-p[4])**2/(2*p[5]**2))

def e_gauss_fit(p, x, y):
  range_p = p[6:]
  if len(range_p) != 2*(len(p)-len(range_p)):
    print("ERROR: should have twice as many ranges than parameter in doublegaussfit")
    exit(1)
  #if p[0] > range_p[0] and p[0] < range_p[1] and p[1] > range_p[2] and p[1] < range_p[3] and p[2] > range_p[4] and p[2] < range_p[5] and p[3] > range_p[6] and p[3] < range_p[7] and p[4] > range_p[8] and p[4] < range_p[9] and p[5] > range_p[10] and p[5] < range_p[11]:
  if p[2] > range_p[4] and p[2] < range_p[5] and p[5] > range_p[10] and p[5] < range_p[11]:
    return gauss_fit(p,x) -y
  else:
    ret = [1e6 for i in [None]*(len(p))]
    return ret
  
def doublegaussfit(x,proba,par, sb, range_p = []):
  out = []
  if range_p == []:
    p,cov,infodict,mesg,ier = leastsq(e_gauss_fit, par[:], args=(x, proba), maxfev=100000, full_output=1)
  else:
    p,cov,infodict,mesg,ier = leastsq(e_gauss_fit, par[:]+range_p[:], args=(x, proba), maxfev=100000, full_output=1)
  xxx = np.arange(min(x),max(x),x[1]-x[0])
  ccc = gauss_fit(p,xxx) # this will only work if the units are pixel and not wavelength
  ss_err = (infodict['fvec']**2).sum()
  ss_tot = ((proba-proba.mean())**2).sum()
  rsquare = 1-(ss_err/ss_tot)
  return p,rsquare, [xxx,ccc]
  
def find_local_max(counts, values, binsize):
  posmax = [i for i,x in enumerate(counts) if x == counts.max()]
  Max = [values[posmax[0]]]
  counts2 = counts.copy()
  counts2[posmax[0]] = 0
  attempt = 0
  
  while len(Max) <= 3 and attempt < len(counts):
    attempt = attempt + 1
    rec = True
    posmax.append([i for i,x in enumerate(counts2) if x == counts2.max()][0])
    if posmax[-1] == 0 or posmax[-1] == len(counts2)-1:
      counts2[posmax[-1]] = 0
      continue
    for i in range(0,len(posmax)-1):
      if (posmax[-1] <= posmax[i]+int(5/binsize) and posmax[-1] >= posmax[i]-int(15/binsize)) or counts[posmax[-1]+1] == 0 or counts[posmax[-1]-1] == 0:
        counts2[posmax[-1]] = 0
        rec = False
        break
    if rec:
      Max.append(values[posmax[-1]])
      counts2[posmax[-1]] = 0
  return Max
  
def remove_letters(r):
    return re.sub(r'[a-zA-Z]', r'', r)
  
def good_marker_size(x,y,fig,sb):

  # initialize a plot to determine the distance between the data points in pixel:    
  s = 0.0
  points = sb.scatter(x,y,s=s,marker='s')
  #sb.axis([min(x)-1., max(x)+1., min(y)-1., max(y)+1.])

  # retrieve the pixel information:
  xy_pixels = sb.transData.transform(np.vstack([x,y]).T)
  xpix, ypix = xy_pixels.T

  # In matplotlib, 0,0 is the lower left corner, whereas it's usually the upper 
  # right for most image software, so we'll flip the y-coords
  width, height = fig.canvas.get_width_height()
  ypix = height - ypix

  # this assumes that your data-points are equally spaced
  s1 = xpix[1]-xpix[0]
  return s1
