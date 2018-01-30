
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np

x = np.linspace(0,100, 10)
row = 5
col = 5
ft = 10
fig = plt.figure(1,figsize=(3*row, 3*col))
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
outer = gridspec.GridSpec(row, col, wspace=0.2, hspace=0.3)
for i in range(0,col*row ):
  sb = plt.Subplot(fig, outer[i])
  sb.plot(x,x)
  sb.text(0.6, 0.2, "text", horizontalalignment='left', color = 'r', fontsize = ft)
  fig.add_subplot(sb)
plt.show()

