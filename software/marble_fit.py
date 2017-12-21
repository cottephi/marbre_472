import numpy as np
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas
import scipy

def marble_fit(my_data):
  marble_x, marble_z = zip(*sorted(my_data[0].items()))
  marble_z_forplot = []
  marble_z_forplot.append(marble_z)
  marble_data = SortedDict()
  for x in my_data[0]:
    if not abs(np.mean(marble_z)-my_data[0][x]) > 4*math.sqrt(np.var(marble_z)):
      marble_data[x] = my_data[0][x]
  if len(my_data) >1:
    marble_x, marble_z = zip(*sorted(my_data[-1].items()))
    marble_z_forplot.append(marble_z)
    for x in my_data[-1]:
      if not abs(np.mean(marble_z)-my_data[-1][x]) > 4*math.sqrt(np.var(marble_z)):
        marble_data[x] = my_data[-1][x]
  
  marble_x, marble_z = zip(*sorted(marble_data.items()))
  marble_x = np.array(marble_x)
  marble_z = np.array(marble_z)
  my_marble_regr = linear_model.LinearRegression()
  my_marble_regr.fit(marble_x[:,np.newaxis], marble_z)
  corr_data = apply_marble_fit(my_marble_regr, my_data)
  
  plot_fit(marble_x, marble_z, my_marble_regr, marble_z_forplot)
  
  return corr_data
  
def apply_marble_fit(fit, data):
  corrected_data = [SortedDict()]
  for data in data[1:-1]:
    for x in data:
      corrected_data[-1][x] = data[x]-fit.predict(data[x])[0]
    corrected_data.append(SortedDict())
  corrected_data = corrected_data[:-1]
  return corrected_data
  
def plot_fit(x, z, regr, z_forplot):
  x_test = np.linspace(np.min(x), np.max(x), 100)
  z_df = []
  z_df.append(pandas.DataFrame({'z':z_forplot[0]}))
  fig = plt.figure(0,figsize=(12, 12))
  fig.tight_layout()
  ax1 = fig.add_subplot(211)
  ax1.set_title("Marble Height vs length")
  ax1.set_xlabel('x(micrometer)', fontsize=14)
  ax1.set_ylabel('z(micrometer)', fontsize=12)
  ax1.scatter(x, z, color='black')
  ax1.plot(x_test, regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=3)
  plt.text((x.max()-x.min())/2,z.max()-(z.max()-z.min())/4,"z(x)=" + str(regr.coef_) + "x+" + str(regr.intercept_),horizontalalignment='center')
  if len(z_forplot) == 2:
    z_df.append(pandas.DataFrame({'z':z_forplot[1]}))
    ax2 = fig.add_subplot(223)
    ax2.set_xlabel('z(micrometer)', fontsize=14)
    ax2.set_ylabel('count', fontsize=12)
    binned_z, binned_x, _ = ax2.hist(z_df[0]['z'], bins='auto')
    ax2.hist(z_df[0]['z'], bins='auto')
    plt.text(z_df[0]['z'].min()+(z_df[0]['z'].max()-z_df[0]['z'].min())/4, binned_z.max()-(binned_z.max()-binned_z.min())/3, z_df[0]['z'].describe(),horizontalalignment='center')
    ax3 = fig.add_subplot(224)
    ax3.set_xlabel('z(micrometer)', fontsize=14)
    ax3.set_ylabel('count', fontsize=12)
    binned_z, binned_x, _ = ax3.hist(z_df[1]['z'], bins='auto')
    x_plot = np.linspace(z_df[1]['z'].min(), z_df[1]['z'].max(), 1000)
    m = scipy.stats.poisson.fit(z_df[0]['z'])
    result = len(z_forplot[1])*scipy.stats.poisson.pdf(x_plot, m)
    ax3.plot(x_plot, result)
    plt.text(z_df[1]['z'].min()+(z_df[1]['z'].max()-z_df[1]['z'].min())/4, binned_z.max()-(binned_z.max()-binned_z.min())/3, z_df[1]['z'].describe(),horizontalalignment='center')
  elif len(z_forplot) == 1:
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel('z(micrometer)', fontsize=14)
    ax2.set_ylabel('count', fontsize=12)
    binned_z, binned_x, _ = ax2.hist(z_df[0]['z'], bins='auto')
    ax2.hist(z_df[0]['z'], bins='auto')
    plt.text(z_df[0]['z'].min()+(z_df[0]['z'].max()-z_df[0]['z'].min())/4, binned_z.max()-(binned_z.max()-binned_z.min())/3, z_df[0]['z'].describe(),horizontalalignment='center')
  else:
    print("marble_fit::plot_fit ERROR : wrong size for marble z data")
    exit(1)
  plt.savefig("marble.pdf",bbox_inches = "tight")
  plt.show()
  
def poisson(k, lamb):
  return (lamb**k/factorial(k)) * np.exp(-lamb)

