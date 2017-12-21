from marble_fit import marble_fit
from sortedcontainers import SortedDict
import matplotlib.pyplot as plt

def my_analysis(my_data_to_analyse, do_marble_fit = True):
  holes_data = None
  if do_marble_fit:
    holes_data = marble_fit(my_data_to_analyse)
  else:
    holes_data = my_data_to_analyse
  plt.figure(2)
  plt.title("Holes Height vs length")
  plt.xlabel('x(micrometer)', fontsize=18)
  plt.ylabel('z(micrometer)', fontsize=16)
  for hole in holes_data:
    x, z = zip(*sorted(hole.items()))
    plt.scatter(x, z, color='black')
  plt.savefig("holes.pdf")
  #plt.show()
