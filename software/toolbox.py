from sortedcontainers import SortedDict
import numpy as np
import scipy
import math
import pandas
from decimal import Decimal
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import re

def cut(xz, myrange):
  if type(xz) is pandas.DataFrame:
    if len(xz) == 0:
      return xz
    if len(myrange) == 1:
      return pandas.DataFrame({'z':[ z for z in xz['z'] if not abs(xz['z'].mean()-z) > myrange[0]*math.sqrt(xz['z'].var()) ]})
    elif len(myrange) == 2:
      return pandas.DataFrame({'z':[ z for z in xz['z'] if z > myrange[0] and z < myrange[1] ]})
    else:
     return xz
  
  elif type(xz) is list:
    if len(xz) == 0:
      return xz
    if len(myrange) == 1:
      if len(xz) == 2:
        np_xz = np.array(xz[1])
        lx,lz = zip(*[ [x,z] for x,z in zip(xz[0],xz[1]) if not abs(np_xz.mean()-z) > myrange[0]*math.sqrt(np_xz.var()) ])
        return [list(lx),list(lz)]
      else:
        np_z = np.array(xz)
        return [ z for z in xz if not abs(np_z.mean()-z) > myrange[0]*math.sqrt(np_z.var()) ]
    elif len(myrange) == 2:
      if len(xz) == 2:
        lx,lz = zip(*[ [x,z] for x,z in zip(xz[0],xz[1]) if z > myrange[0] and z < myrange[1] ])
        return [list(lx),list(lz)]
      else:
        np_z = np.array(xz)
        return [ z for z in xz if z > myrange[0] and z < myrange[1] ]
    else:
      return xz
  
def FormStatBox(df_z_selected):
  return 'count : ' + str(int(df_z_selected.describe()['count'])) + '\nmean : ' + str(round(Decimal(df_z_selected.describe()['mean']),0))+'$\; \mu m$' + '\nstd : ' + str(round(Decimal(df_z_selected.describe()['std']),0))+'$\; \mu m$'
  
def FormFitBox(param, df_z_selected, chisquare, binsize = 0):
  if len(param) == 0:
    return 'fit failed'
  if len(param) == 3:
    if binsize != 0:
      param[0] = param[0]/binsize;
    return 'count : ' + str(int(df_z_selected.describe()['count'])) + '\nFit result:\n Int=' + str(round(Decimal(param[0]),0)) + '\n mean=' + str(round(Decimal(param[1]),0))+'$\; \mu m$' + '\n sigma=' + str(abs(round(Decimal(param[2]),0)))+'$\; \mu m$' + '\n  chi2=' + str(abs(round(Decimal(chisquare[0]),2)))
  elif len(param) == 6:
    if binsize != 0:
      param[0] = param[0]/binsize;
      param[3] = param[3]/binsize;
    return 'count : ' + str(int(df_z_selected.describe()['count'])) + '\nFit result:\n First gaussian:\n  Int=' + str(round(Decimal(param[0]),0)) + '\n  mean=' + str(round(Decimal(param[1]),0))+'$\; \mu m$' + '\n  sigma=' + str(abs(round(Decimal(param[2]),0)))+'$\; \mu m$' + '\n  chi2=' + str(abs(round(Decimal(chisquare[0]),2))) + '\n Second gaussian:\n  Int=' + str(round(Decimal(param[3]),0)) + '\n  mean=' + str(round(Decimal(param[4]),0))+'$\; \mu m$' + '\n  sigma=' + str(abs(round(Decimal(param[5]),0)))+'$\; \mu m$' + '\n  chi2=' + str(abs(round(Decimal(chisquare[1]),2)))
  else:
    print("ERROR: unknown number of parameter in FormFitBox")
    exit(1)


def singlegaussfit(x,proba,par, range_p = []):
  out = []
  chisquare = 100
  if len(par)+len(range_p) > len(x):
    print("   You gave me less x than parameters, can not fit in those conditions!")
    return None, chisquare, None
  print("    Fit attempt around possible mean ", par[1], "...")
  if range_p == []:
    p,cov,infodict,mesg,ier = leastsq(e_single_gauss_fit, par[:], args=(x, proba), maxfev=100000, full_output=1)
  else:
    p,cov,infodict,mesg,ier = leastsq(e_single_gauss_fit, par[:]+range_p[:], args=(x, proba), maxfev=100000, full_output=1)
  xxx = np.arange(min(x), max(x), x[1]-x[0])
  ccc = single_gauss_fit(p,xxx) 
  if proba.sum() != 0:
    chisquare = np.array([ (c-pro)*(c-pro) for c,pro in zip(ccc,proba) ]).sum()/proba.sum()
  xxx = np.arange(p[1]-(max(x)-min(x))/2, p[1]+(max(x)-min(x))/2, x[1]-x[0])
  ccc = single_gauss_fit(p,xxx) # this will only work if the units are pixel and not wavelengths
  
  print("     Mean found : ", p[1], " chi**2=", chisquare)
  return p[:len(p)-len(range_p)], chisquare, [xxx,ccc]

def e_single_gauss_fit(p, x, y):
  range_p = p[3:]
  binsize = x[2]-x[1]
  if len(range_p) != 2*(len(p)-len(range_p)) and len(range_p) != 0:
    print("ERROR: should have twice as many ranges than parameters in singlegaussfit")
    exit(1)
  #if p[0] > range_p[0] and p[0] < range_p[1] and p[1] > range_p[2] and p[1] < range_p[3] and p[2] > range_p[4] and p[2] < range_p[5] and p[3] > range_p[6] and p[3] < range_p[7] and p[4] > range_p[8] and p[4] < range_p[9] and p[5] > range_p[10] and p[5] < range_p[11]:
  if len(range_p) != 0:
    if p[2] > range_p[4] and p[2] < range_p[5]:
      return single_gauss_fit(p,x+binsize/2) -y
    else:
      ret = [1e6 for i in [None]*(len(p))]
      return ret
  else:
    return single_gauss_fit(p,x) -y

def single_gauss_fit(p, x):
  return p[0]*(1/np.sqrt(2*math.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2))
  
def doublegaussfit(x,proba,par, range_p = []):
  out = []
  print("    Fit attempt around possible means ",par[1]," and ",par[4],"...")
  if range_p == []:
    p,cov,infodict,mesg,ier = leastsq(e_double_gauss_fit, par[:], args=(x, proba), maxfev=100000, full_output=1)
  else:
    p,cov,infodict,mesg,ier = leastsq(e_double_gauss_fit, par[:]+range_p[:], args=(x, proba), maxfev=100000, full_output=1)
  xxx = np.arange(min(x),max(x),x[1]-x[0])
  ccc = double_gauss_fit(p,xxx) # this will only work if the units are pixel and not wavelength
  ss_err = (infodict['fvec']**2).sum()
  ss_tot = ((proba-proba.mean())**2).sum()
  chisquare = 1-(ss_err/ss_tot)
  print("     Means found : ",p[1]," and ",p[4],", R**2=",chisquare)
  return p[:len(p)-len(range_p)],chisquare, [xxx,ccc]

def e_double_gauss_fit(p, x, y):
  range_p = p[6:]
  binsize = x[2]-x[1]
  if len(range_p) != 2*(len(p)-len(range_p)):
    print("ERROR: should have twice as many ranges than parameter in doublegaussfit")
    exit(1)
  #if p[0] > range_p[0] and p[0] < range_p[1] and p[1] > range_p[2] and p[1] < range_p[3] and p[2] > range_p[4] and p[2] < range_p[5] and p[3] > range_p[6] and p[3] < range_p[7] and p[4] > range_p[8] and p[4] < range_p[9] and p[5] > range_p[10] and p[5] < range_p[11]:
  if p[2] > range_p[4] and p[2] < range_p[5] and p[5] > range_p[10] and p[5] < range_p[11]:
    return double_gauss_fit(p,x+binsize/2) -y
  else:
    ret = [1e6 for i in [None]*(len(p))]
    return ret

def double_gauss_fit(p, x):
  return p[0]*(1/np.sqrt(2*math.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2))+p[3]*(1/np.sqrt(2*math.pi*(p[5]**2)))*np.exp(-(x-p[4])**2/(2*p[5]**2))
  
def find_local_max(counts, values, binsize):
  posmax = [i for i,x in enumerate(counts) if x == counts.max()]
  Max = [values[posmax[0]]]
  counts2 = counts.copy()
  counts2 = np.array([ z if i < posmax[0]-int(25/binsize) or i > posmax[0]+int(25/binsize) else 0 for i,z in enumerate(counts2)])
  attempt = 0
  
  while len(Max) < 4 and attempt < len(counts):
    attempt = attempt + 1
    rec = True
    tmp_posmax = [i for i,x in enumerate(counts2) if x == counts2.max()][0]
    if tmp_posmax == 0 or tmp_posmax == len(counts2)-1:
      counts2[tmp_posmax] = 0
      continue
    for i in range(0,len(posmax)):
      if tmp_posmax >= posmax[i]-int(30/binsize) or counts[tmp_posmax+1] == 0 or counts[tmp_posmax-1] == 0:
        counts2[tmp_posmax] = 0
        rec = False
        break
    if rec:
      Max.append(values[tmp_posmax])
      posmax.append(tmp_posmax)
      counts2 = np.array([ z if i < tmp_posmax-int(30/binsize) or i > tmp_posmax+int(30/binsize) else 0 for i,z in enumerate(counts2)])
  print("   Local max found : ",Max)
  return Max
  

def gainfit(x, proba, par, range_p = []):
  out = []
  chisquare = 100
  if len(par)+len(range_p) > len(x):
    print("   You gave me less x than parameters, can not fit in those conditions!")
    return None, chisquare, None
  print("    Fiting gain...")
  if range_p == []:
    p,cov,infodict,mesg,ier = leastsq(e_gain_fit, par[:], args=(x, proba), maxfev=100000, full_output=1)
  else:
    p,cov,infodict,mesg,ier = leastsq(e_gain_fit, par[:]+range_p[:], args=(x, proba), maxfev=100000, full_output=1)
  xxx = np.arange(0, 4000,10)
  ccc = gain_fit(p,xxx) 
  if proba.sum() != 0:
    chisquare = np.array([ (c-pro)*(c-pro) for c,pro in zip(ccc,proba) ]).sum()/proba.sum()
  
  print("Parameters found: ", p, ", chi**2=", chisquare)
  return p[:len(p)-len(range_p)], chisquare, [xxx,ccc]

def e_gain_fit(p, x, y):
  range_p = p[3:]
  binsize = x[2]-x[1]
  if len(range_p) != 2*(len(p)-len(range_p)) and len(range_p) != 0:
    print("ERROR: should have twice as many ranges than parameters in singlegaussfit")
    exit(1)
  #if p[0] > range_p[0] and p[0] < range_p[1] and p[1] > range_p[2] and p[1] < range_p[3] and p[2] > range_p[4] and p[2] < range_p[5] and p[3] > range_p[6] and p[3] < range_p[7] and p[4] > range_p[8] and p[4] < range_p[9] and p[5] > range_p[10] and p[5] < range_p[11]:
  if len(range_p) != 0:
    if p[2] > range_p[4] and p[2] < range_p[5]:
      return gain_fit(p,x) -y
    else:
      ret = [1e6 for i in [None]*(len(p))]
      return ret
  else:
    return gain_fit(p,x) -y

def gain_fit(p, x):
  return p[0]*np.exp(p[1]*np.exp(-p[2]/(x/0.001))*0.001)
  
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
  s1 = xpix
  s2 = ypix
  if len(xpix) > 1:
    s1 = xpix[1]-xpix[0]
  if len(ypix) > 1:
    y2 = ypix[0]
    test = 1
    while y2 == ypix[0] and test < len(ypix):
      y2 = ypix[test]
      test = test + 1
    if y2 == ypix[0]:
      s2 = y2
    else:
      s2 = y2-ypix[0]
  return [s1,s2]
  
def GetNcolNrow(l_df_data):
  if int(math.sqrt(len(l_df_data))) == math.sqrt(len(l_df_data)):
    int_nrow = int(math.sqrt(len(l_df_data)))
    int_ncol = int(math.sqrt(len(l_df_data)))
  elif int(len(l_df_data)/int(math.sqrt(len(l_df_data)))) == len(l_df_data)/int(math.sqrt(len(l_df_data))):
    int_nrow = int(math.sqrt(len(l_df_data)))
    int_ncol = int(len(l_df_data)/int_nrow)
  else:
    int_nrow = int(math.sqrt(len(l_df_data)))
    int_ncol = int(math.sqrt(len(l_df_data))) + 1
  return int_nrow, int_ncol
