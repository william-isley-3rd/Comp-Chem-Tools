#! /usr/bin/env python
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

# This part reads in the excel file, and parses it into dataframse
path_in = 'C:/Users/isle132/Documents/NGDE/HRUS Data/From Forrest Test Data'
exp_data_file = 'FTIR Plate Run.xlsx'

xl_exp_data = pd.ExcelFile(path_in+'/'+exp_data_file)
df_xl_in = xl_exp_data.parse("Nucleation Data", index_col=0)
df_xl_in.columns = df_xl_in.columns.str.replace(' minutes', '')

time = df_xl_in.columns
freq_sel = (df_xl_in.index.get_level_values(0) <= 1601) & (df_xl_in.index.get_level_values(0) >= 1575)
test_signal = df_xl_in.loc[freq_sel].T
normalized_min = 0.0001
normalized_data = (test_signal - test_signal.min()) / (test_signal.max() - test_signal.min()+ normalized_min) \
                  + normalized_min
avrami_transform = normalized_data.ix[2:].apply(lambda x: np.log(np.log(1/(1-x))))
avrami_transform.index = list(map(lambda x: np.log(x), avrami_transform.index.astype(float)))

#x_time = avrami_transform.index.values.reshape(len(avrami_transform.index), 1)
x_time = avrami_transform.index.values
# Do Linear Regression Individually for each Frequency
avrami_coef = []
avrami_inter = []
avrami_r = []
for column in avrami_transform:
    #y = avrami_transform[column].values.reshape(len(avrami_transform.index), 1)
    y = avrami_transform[column].values
    slope, intercept, r_value, p_value, std_err = linregress(x_time, y)
    avrami_coef.append(slope)
    avrami_inter.append(intercept)
    avrami_r.append(r_value)

# print out the Regression Information
print('average is ', np.average(avrami_coef), ' x + ', np.average(avrami_inter),
      '\nSTDev is ', np.std(avrami_coef), ' x + ', np.std(avrami_inter),
      '\nMax is ', np.max(avrami_coef), ' x + ', np.max(avrami_inter),
      '\nMin is ', np.min(avrami_coef), ' x + ', np.min(avrami_inter),
      '\nR**2 is ', np.average(r_value))

reg_plot_y = x_time * np.average(avrami_coef) + np.average(avrami_inter)

# Plot Difference Spectrum in time
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
df_xl_in.plot(ax=ax1, legend=False, xlim=(1000,1800), ylim=(0,0.12))
ax1.set_xlabel('Frequency / cm$^{-1}$')
ax1.set_ylabel('Difference Absorbance')
# ax1.set_xlim(ax1.get_xlim()[::-1])
ax1.xaxis.label.set_fontsize(20)
ax1.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
# Plot Difference Spectrum for 1575 - 1600 cm-1
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
test_signal.plot(ax=ax2, legend=False, ylim=(0, 0.12))
ax2.set_xlabel('Time / min')
ax2.set_ylabel('Difference Absorbance')
ax2.xaxis.label.set_fontsize(20)
ax2.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
# Plot Extent Transformation
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
normalized_data.plot(ax=ax3, legend=False, ylim=(0, 1))
ax3.set_xlabel('Time / min')
ax3.xaxis.label.set_fontsize(20)
ax3.yaxis.label.set_fontsize(20)
ax3.set_ylabel('Extent Transformation (%)')
plt.tight_layout(0.05)
# Plot Avrami Kinetics
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
avrami_transform.plot(ax=ax4, legend=False, alpha=0.3)
plt.plot(x_time, reg_plot_y, color='black')
ax4.set_xlabel('LN[Time] / LN[min]')
ax4.set_ylabel('LN[LN[1/(1/X)]]')
ax4.xaxis.label.set_fontsize(20)
ax4.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
plt.show()


