import numpy as np
import scipy
import matplotlib as plt
import math

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
