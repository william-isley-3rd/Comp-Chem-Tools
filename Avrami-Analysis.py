#! /usr/bin/env python
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.fftpack import rfft, irfft, fftfreq


def normalize(df):
    result = df.copy()
    eps = 0.0000001
    for feature_name in df.columns:
        max_value = df[feature_name].max()
        min_value = df[feature_name].min()
        result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value + eps) + eps
    return result


def extent_trans_X(x, n, k):
    return 1 - np.exp(-1 * np.exp(k) * pow(x, n))


def running_mean(x, n=5):
    return np.convolve(x, np.ones((n,))/n)[(n-1):]


# This part reads in the excel file, and parses it into dataframse
# path_in = 'C:/Users/isle132/Documents/NGDE/Data/From Forrest Test Data'
# exp_data_file = 'FTIR Plate Run.xlsx'
# exp_data_file = 'Copy of FTIR Nucleation Run.xlsx'
path_in = 'C:/Users/isle132/Documents/NGDE/Data/'
exp_data_file = 'FTIR Nucleation Run v2.xlsx'

xl_exp_data = pd.ExcelFile(path_in+'/'+exp_data_file)
df_xl_in = xl_exp_data.parse("High Conc Run Py", index_col=0)
# df_xl_in = xl_exp_data.parse("Medium Conc Run Py", index_col=0)
# df_xl_in = xl_exp_data.parse("Low Conc Run Py", index_col=0)
# df_xl_in = xl_exp_data.parse("Nucleation Data", index_col=0)
# df_xl_in.columns = df_xl_in.columns.str.replace(' minutes', '')


injection_time = 5.0
# high concentration growth, time after injection
growth_onset = 4.0

# medium concentration growth, time after injection
# growth_onset = 10

time = df_xl_in.columns
all_t_signal = df_xl_in.copy()
all_t_smooth_signal = df_xl_in.copy()

# FFT FILTER
sample_freq = all_t_signal.index[1] - all_t_signal.index[0]  # units of frequency cm**-1
fft_freqs = fftfreq(len(all_t_signal), d=sample_freq)        # bins into cycles / (cm**-1)
for time_column in all_t_signal:
    ft_signal_data = rfft(all_t_signal[time_column])
    ft_signal_data[fft_freqs < -0.10] = 0                    # Cut out all freqs larger than 0.05 cycles / (cm**-1)
    ft_signal_smooth = irfft(ft_signal_data)
    all_t_smooth_signal[time_column] = ft_signal_smooth

# Subset the frequencies near the peak
freq_sel = (all_t_smooth_signal.index.get_level_values(0) <= 1600) & \
           (all_t_smooth_signal.index.get_level_values(0) >= 1590)
test_signal = all_t_smooth_signal.loc[freq_sel].T
test_signal = test_signal[test_signal.index >= injection_time]  # injection time at 5 minutes
test_signal.index = test_signal.index - injection_time  # injection time at 5 minutes


normalized_data_unproc = normalize(test_signal)
# SAVITZKY GOLAY FILTER
# apply a Savitzky Golay filter to smooth the data, using a 9 point window with order 3 polynomial
normalized_data = normalized_data_unproc.apply(savgol_filter, args=(9, 3,))

# take values for time greater than 1, to prevent singularity in log
avrami_transform = normalized_data.ix[1:].apply(lambda x: np.log(np.log(1/(1-x))))
avrami_transform.index = list(map(lambda x: np.log(x), avrami_transform.index.astype(float)))

impingement_factor = 1
austin_rickett_transform = normalized_data.ix[1:].apply(lambda x: np.log(pow(1-x, -1 * impingement_factor)-1) +
                                                        np.log(impingement_factor))
austin_rickett_transform.index = list(map(lambda x: np.log(x), austin_rickett_transform.index.astype(float)))

x_ln_time = avrami_transform.index.values
avrami_slopes = pd.DataFrame(index=x_ln_time[1:])  # make the empty array

# set up the time selectors
nucl_ln_time_sel = (avrami_transform.index.get_level_values(0) <= np.log(growth_onset-0.8)) & \
                    (avrami_transform.index.get_level_values(0) >= 0.2)
growth_ln_time_sel = (avrami_transform.index.get_level_values(0) >= np.log(growth_onset)) & \
                     (avrami_transform.index.get_level_values(0) <= np.log(40.0))
x_ln_time_nucl_sub = avrami_transform.loc[nucl_ln_time_sel].index.values
x_ln_time_growth_sub = avrami_transform.loc[growth_ln_time_sel].index.values

# Do Linear Regression Individually for each Frequency
avrami_coef_w = []
avrami_inter_w = []
avrami_r_w = []
avrami_coef_n = []
avrami_inter_n = []
avrami_r_n = []
avrami_coef_n_s = []
avrami_inter_n_s = []
avrami_r_n_s = []
avrami_coef_g = []
avrami_inter_g = []
avrami_r_g = []
for column in avrami_transform:
    y_w = avrami_transform[column].values
    y_n = avrami_transform.loc[nucl_ln_time_sel, column].values
    y_g = avrami_transform.loc[growth_ln_time_sel, column].values
    slope_w, intercept_w, r_value_w, p_value_w, std_err_w = linregress(x_ln_time, y_w)
    avrami_coef_w.append(slope_w)
    avrami_inter_w.append(intercept_w)
    avrami_r_w.append(r_value_w)
    slope_n, intercept_n, r_value_n, p_value_n, std_err_n = linregress(x_ln_time_nucl_sub, y_n)
    avrami_coef_n.append(slope_n)
    avrami_inter_n.append(intercept_n)
    avrami_r_n.append(r_value_n)
    slope_g, intercept_g, r_value_g, p_value_g, std_err_g = linregress(x_ln_time_growth_sub, y_g)
    avrami_coef_g.append(slope_g)
    avrami_inter_g.append(intercept_g)
    avrami_r_g.append(r_value_g)
    avrami_slopes[column] = np.diff(y_w)/np.diff(x_ln_time)

# Get the slopes from the avrami plots - time exponent
avrami_slopes.fillna(method='backfill')
# apply a Savitzky Golay filter to smooth the data, using a 7 point window with order 3 polynomial
# avrami_slopes_smooth = avrami_slopes.apply(savgol_filter, args=(5, 2,), mode='mirror')
avrami_slopes_smooth = avrami_slopes.copy()
normalized_data['mean'] = normalized_data.mean(axis=1)
avrami_slopes['mean'] = avrami_slopes.mean(axis=1)
avrami_slopes_smooth['mean'] = avrami_slopes_smooth.mean(axis=1)

# Get the n vs X data
transformed_fraction = pd.DataFrame(index=normalized_data.ix[1:, 'mean'].values[1:])
transformed_fraction['mean'] = avrami_slopes_smooth['mean'].values
transformed_fraction = transformed_fraction.sort_index()

# print out the Regression Information
print('For the entire Data set')
print('average is ', np.average(avrami_coef_w), ' x + ', np.average(avrami_inter_w),
      '\nSTDev is ', np.std(avrami_coef_w), ' x + ', np.std(avrami_inter_w),
      '\nMax is ', np.max(avrami_coef_w), ' x + ', np.max(avrami_inter_w),
      '\nMin is ', np.min(avrami_coef_w), ' x + ', np.min(avrami_inter_w),
      '\nR**2 is ', np.average(avrami_r_w))

print('For the nucleation subset of the full Data set')
print('average is ', np.average(avrami_coef_n), ' x + ', np.average(avrami_inter_n),
      '\nSTDev is ', np.std(avrami_coef_n), ' x + ', np.std(avrami_inter_n),
      '\nMax is ', np.max(avrami_coef_n), ' x + ', np.max(avrami_inter_n),
      '\nMin is ', np.min(avrami_coef_n), ' x + ', np.min(avrami_inter_n),
      '\nR**2 is ', np.average(avrami_r_n))

print('For the growth subset of the full Data set')
print('average is ', np.average(avrami_coef_g), ' x + ', np.average(avrami_inter_g),
      '\nSTDev is ', np.std(avrami_coef_g), ' x + ', np.std(avrami_inter_g),
      '\nMax is ', np.max(avrami_coef_g), ' x + ', np.max(avrami_inter_g),
      '\nMin is ', np.min(avrami_coef_g), ' x + ', np.min(avrami_inter_g),
      '\nR**2 is ', np.average(avrami_r_g))

# Prepare data for plotting
reg_plot_y_w = x_ln_time * np.average(avrami_coef_w) + np.average(avrami_inter_w)
reg_plot_y_n = x_ln_time * np.average(avrami_coef_n) + np.average(avrami_inter_n)
reg_plot_y_g = x_ln_time * np.average(avrami_coef_g) + np.average(avrami_inter_g)

x_time_extent = normalized_data.index.values
X_n = extent_trans_X(x_time_extent, np.average(avrami_coef_n), np.average(avrami_inter_n))
X_w = extent_trans_X(x_time_extent, np.average(avrami_coef_w), np.average(avrami_inter_w))

abs_max = np.max(df_xl_in.loc[freq_sel].values)
# Plot Difference Spectrum in time
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
# df_xl_in.plot(ax=ax1, legend=False, xlim=(1350, 1750), ylim=(0, abs_max))
all_t_smooth_signal.plot(ax=ax1, legend=False, xlim=(1350, 1750), ylim=(0, abs_max))
ax1.set_xlabel('Frequency / cm$^{-1}$')
ax1.set_ylabel('Difference Absorbance')
# ax1.set_xlim(ax1.get_xlim()[::-1])
ax1.xaxis.label.set_fontsize(20)
ax1.yaxis.label.set_fontsize(20)
ax1.tick_params(top='on', direction='in')
plt.tight_layout(0.05)
# Plot Difference Spectrum for 1575 - 1600 cm-1
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
test_signal.plot(ax=ax2, legend=False, xlim=(0, growth_onset+15), ylim=(0, abs_max))
ax2.set_xlabel('Time / min')
ax2.set_ylabel('Difference Absorbance')
ax2.xaxis.label.set_fontsize(20)
ax2.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
# Plot Extent Transformation
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
normalized_data.plot(ax=ax3, legend=False, xlim=(0, growth_onset+15), ylim=(0, 1), alpha=0.3)
plt.plot(x_time_extent, X_n, color='blue')
# plt.plot(x_time_extent, X_w, color='black')
ax3.set_xlabel('Time / min')
ax3.xaxis.label.set_fontsize(20)
ax3.yaxis.label.set_fontsize(20)
ax3.set_ylabel('Extent Transformation (X)')
plt.tight_layout(0.05)
# Plot Avrami Kinetics
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
avrami_transform.plot(ax=ax4, legend=False, ylim=(-4, 2), alpha=0.3)
# plt.plot(x_ln_time, reg_plot_y_w, color='black')
plt.plot(x_ln_time, reg_plot_y_n, color='blue')
plt.plot(x_ln_time, reg_plot_y_g, color='red')
ax4.set_xlabel('LN[Time]')
ax4.set_ylabel('LN[-LN[1-X]]')
ax4.xaxis.label.set_fontsize(20)
ax4.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
# Plot Avrami slope
fig5 = plt.figure(5)
ax5 = fig5.add_subplot(111)
avrami_slopes_smooth.plot(ax=ax5, legend=False, alpha=0.3)
avrami_slopes_smooth.plot(y='mean', color='black', ax=ax5, legend=False, alpha=0.5)
ax5.set_xlabel('LN[Time]')
ax5.set_ylabel('n (Avrami Exponent)')
ax5.set_ylim(0, 5)
ax5.xaxis.label.set_fontsize(20)
ax5.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
# Plot Fraction Transformed vs n
fig6 = plt.figure(6)
ax6 = fig6.add_subplot(111)
transformed_fraction.plot(ax=ax6, legend=False)
ax6.set_xlabel('Fraction Transformed')
ax6.set_ylabel('n Avrami Exponent')
ax6.set_ylim(0, 5)
ax6.xaxis.label.set_fontsize(20)
ax6.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
fig7 = plt.figure(7)
ax7 = fig7.add_subplot(111)
austin_rickett_transform.plot(ax=ax7, legend=False)
ax7.set_xlabel('LN[Time]')
ax7.set_ylabel(r'LN$[\frac{1}{(1-X)}-1]$')
# ax7.set_ylabel(r'LN$[\frac{1}{(1-X)^2}-1]$')
ax7.xaxis.label.set_fontsize(20)
ax7.yaxis.label.set_fontsize(20)
plt.tight_layout(0.05)
plt.show()


