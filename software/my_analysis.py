from marble_fit import marble
from sortedcontainers import SortedDict
import matplotlib.pyplot as plt

def my_analysis(l_sd_data, do_marble_fit = True):
  if do_marble_fit:
    marble(l_sd_data)
    print(len(l_sd_data))
  plt.figure(2)
  plt.title("Holes Height vs length")
  plt.xlabel('x(micrometer)', fontsize=18)
  plt.ylabel('z(micrometer)', fontsize=16)
  for hole in l_sd_data:
    x, z = zip(*sorted(hole.items()))
    plt.scatter(x, z, color='black')
  plt.savefig("holes.pdf")
  #plt.show()
