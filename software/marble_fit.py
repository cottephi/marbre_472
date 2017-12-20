import numpy as np
from sortedcontainers import SortedDict
import math
from sklearn import linear_model
import matplotlib.pyplot as plt


def marble_fit(my_data):
  marble_x = np.array(0)
  marble_z = np.array(0)
  marble_x, marble_z = zip(*sorted(my_data[0].items()))
  marble_data = SortedDict()
  for x in my_data[0]:
    if not abs(np.mean(marble_z)-my_data[0][x]) > 4*math.sqrt(np.var(marble_z)):
      marble_data[x] = my_data[0][x]
  marble_x, marble_z = zip(*sorted(my_data[-1].items()))
  if len(my_data) >1:
    for x in my_data[-1]:
      if not abs(np.mean(marble_z)-my_data[-1][x]) > 4*math.sqrt(np.var(marble_z)):
        marble_data[x] = my_data[-1][x]
  
  marble_x, marble_z = zip(*sorted(marble_data.items()))
  marble_x = np.array(marble_x)
  marble_z = np.array(marble_z)
  my_marble_regr = linear_model.LinearRegression()
  my_marble_regr.fit(marble_x[:,np.newaxis], marble_z)
  x_test = np.linspace(np.min(marble_x), np.max(marble_x), 100)
  plt.figure(0)
  plt.title("Marble Height vs length")
  plt.xlabel('x(micrometer)', fontsize=18)
  plt.ylabel('z(micrometer)', fontsize=16)
  plt.scatter(marble_x, marble_z, color='black')
  plt.plot(x_test, my_marble_regr.predict(x_test[:,np.newaxis]), color='blue', linewidth=3)
  plt.savefig("marble.pdf")
  corr_data = apply_marble_fit(my_marble_regr, my_data)
  return corr_data
  
def apply_marble_fit(fit, data):
  corrected_data = [SortedDict()]
  for data in data[1:-1]:
    for x in data:
      corrected_data[-1][x] = data[x]-fit.predict(data[x])[0]
    corrected_data.append(SortedDict())
  corrected_data = corrected_data[:-1]
  return corrected_data
