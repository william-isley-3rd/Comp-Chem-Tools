from scipy import integrate
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import lognorm
from scipy.signal import argrelextrema
import glob, os
import pandas as pd

#############################################################################################
#                                                                                           #
# VERSION :: 0.2.1  July  2017                                                              #
#                                                                                           #
# AUTHORS :: William Christian Isley III                                                    #
#                                                                                           #
#                                                                                           #
# DISCLAIMER :: This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by the Free Software    #
# Foundation, either version 3 of the License, or any later version                         #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY  #
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE  #
# See the GNU General Public License for more details.                                      #
#                                                                                           #
#############################################################################################
#                                       NOTE ON UNITS                                       #
#                         SI Units are used throughout the code.                            #
#############################################################################################
#                                    ------------                                           #
# 				                     NOMENCLATURE                                           #
# 				                     ------------                                           #
# __________________________________________________________________________________________#
# Variables            	        Meaning (Units)                                             #
# __________________________________________________________________________________________#
# angle                         The angle at which the DLS detects the scattering event     #
#                               Is set by the instrument, select from 173, 90, 45 degrees   #
# wavelength                    Wavelength of the incident laser light (in angstroms)       #
# T                             Temperature (K)                                             #
# eta                           viscosity of solution (Pa s)                                #
# n_ref                         refractive index (unitless)                                 #
# rho                           density of solution ()                                      #
# g_scaler                      rescales the predicted autocorrelations by this factor      #
#                               this is used to match experimentally measured g(r)'s        #
#                               which are typically less than 1                             #
# concentration                 sets the units on the initial concentration.                #
# distribution_list             sets the concentration, mean radius, and sigma for a        #
#                               lognormal distribution. Should be a list of 3 element lists #
#                               To make more than 1 distribution, add additional elements   #
#                               to the list.                                                #
#############################################################################################

# constants definition
pi = np.pi
boltz = 1.3806 * 10**(-23)     # Boltzmann Constant J/(K mol)
N_avo = 6.022 * 10**23         # Avogadro's Number

# Options from experimental apparatus
angle_back = 173*pi/180        # backscatter angle
angle_90 = pi/2                # 90 degree scatter angle
angle_45 = pi/4                # 45 degree scatter angle
wave_532 = 532*10              # wavelength at 532 nm
wave_488 = 488*10              # wavelength at 488 nm
# Options for water solvent
eta_water = 8.90 * 10 ** (-4)  # viscosity of water at 25 degrees C (units Pa s)
n_ref_water = 1.33             # refractive index of water
# Options for ethanol solvent
eta_EtOH = 1.083 * 10 ** (-3)  # viscosity of ethanol at 25 degrees C (units Pa s)
n_ref_EtOH = 1.361             # refractive index of ethanol

# set the experimental parameters
T = 298.1                  # 25 degree C
rho = 2.165 * 1000         # density
wavelength = wave_532      # 532 nm in angstroms
angle = angle_90           # in backscatter mode
eta = eta_EtOH             # viscosity of ethanol at 25 degrees C
n_refraction = n_ref_EtOH  # using ethanol
beta = 1.0                 # correction factor dependent on the geometry, and alignment of laser

g_scaler = 0.95
concentration = 1.0
N_0 = concentration
distribution_list = [[N_0 / np.log(1.25), 1500/2, 1.25]]
# distribution_list = [[N_0 / np.log(1.2), 2500, 1.2],
#                      [1*10**14*N_0 / np.log(1.5), 5, 1.5]]

# q is computed from the experimental properties
# refraction index
# wavelength - choose between wave_532 and wave_488
# scattering angle - choose between angle_45, angle_90, angle_180
q = (4 * pi * n_refraction / wavelength) * np.sin(angle / 2)            # in 1 / Angstroms
c_gamma = 2 * boltz * T * pow(q * (10**10), 2) / (3 * pi * eta) * (10**10)  # in Angstroms / s


def lognorm_f(x, triple_list):
    [pf, r_mean, sigma] = triple_list
    a = 1 / (2 * pow(np.log(sigma), 2))
    mu = np.log(r_mean)
    return pf / x * np.exp(-pow(np.log(x) - mu, 2) * a)


def nucleation_exp_decay_f(x, pf_exp, theta, rho_exp):
    exp_theta = np.exp(theta)
    exp_pf = - theta * pow(4 * pi * rho_exp / 3, 2/3)
    return pf_exp * exp_theta * np.exp(exp_pf * pow(x, 2))


def lognorm_sum(x, triple_list):
    return sum(map(lambda triple: lognorm_f(x, triple), triple_list))


def compute_denom_lognorm(x, triple_list):
    return pow(np.sin(q*x) - q*x*np.cos(q*x), 2) * lognorm_sum(x, triple_list)


def compute_denom_lognorm_val(x, triple):
    return pow(np.sin(q * x) - q * x * np.cos(q * x), 2) * lognorm_f(x, triple)


def compute_num_lognorm(x, time, triple_list):
    return pow(np.sin(q * x) - q * x * np.cos(q * x), 2) * np.exp(- c_gamma * time / x) *\
           lognorm_sum(x, triple_list)


def compute_num_lognorm_val(x, time, triple):
    return pow(np.sin(q * x) - q * x * np.cos(q * x), 2) * np.exp(- c_gamma * time / x) *\
           lognorm_f(x, triple)

time_list = np.append(np.linspace(10**-7, 10**-6, num=50), [np.linspace(10**-6, 10**-5, num=50),
                                                            np.linspace(10**-5, 10**-4, num=50),
                                                            np.linspace(10**-4, 10**-3, num=50),
                                                            np.linspace(10**-3, 10**-2, num=50),
                                                            np.linspace(10**-2, 10**-1, num=50),
                                                            np.linspace(10**-1, 1, num=50)])

# in order to computed the integrals, we need to separate them into different subsets of space
# this helps improve the numerical convergence of the integration.  No analytic solution to the integrals
# for DLS signals is available, and must be done numerically.
# If issues are had with numerical convergence, play with the x limits for the integrations.
g_t_list = []        # This is a g_1
g_t_error_list = []  # this is a g_1
g_t_previous = 0.0   # this is a g_1
for t in time_list:
    seg_denom = []
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 0, 1))
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 1, 50))
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 50, 200))
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 200, 500))
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 500, 1000))
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 1000, 2000))
    seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm(x, distribution_list), 2000, np.inf))
    seg_denom_sum = np.sum(seg_denom, axis=0)

    seg_num = []
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 0, 1))
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 1, 50))
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 50, 200))
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 200, 500))
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 500, 1000))
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 1000, 2000))
    seg_num.append(integrate.quad(lambda x: compute_num_lognorm(x, t, distribution_list), 2000, np.inf))
    seg_num_sum = np.sum(seg_num, axis=0)
    g_t = seg_num_sum[0] / seg_denom_sum[0]
    g_t_error = seg_num_sum[1] / seg_denom_sum[1]
    if g_t > 1.0:
        g_t = g_t_previous
    else:
        g_t_previous = g_t
    g_t_list.append(g_t)
    g_t_error_list.append(g_t_error)

# g_2 = 1+beta*g_1(q,t)**2 , assume beta=1
g_2_list = [1 + beta * pow(x, 2) for x in g_t_list]

g_t_scale_list = [g_scaler * x for x in g_t_list]

time_list_us = [x * 10**6 for x in time_list]
# print(list(zip(time_list, g_t_scale_list)))

# This section of the code reads in the experimental data, and plots it to compare results.
# Import the data from the .correlation.csv file
# Plot the delay time vs correlation
path_in = 'C:/Users\isle132/Desktop/For Python/140 nm Silica Data/173 deg/Cell Center'

os.chdir(path_in)  # go to the path specified
cor_files = glob.glob('*.correlation.csv')

cor_dfs = [pd.read_csv(cor_fp, encoding="ISO-8859-1").assign(Data=os.path.splitext(os.path.splitext(cor_fp)[0])[0])
           for cor_fp in cor_files]
cor_df = pd.concat(cor_dfs, ignore_index=True)

cor_groups = cor_df.groupby('Data')
fig1, ax1 = plt.subplots()
ax1.margins(0.05)
ax1.set_xscale('log')
ax1.set_xlim(cor_df['Delay Time       (µs)'].min(), cor_df['Delay Time       (µs)'].max())
for name, group in cor_groups:
    ax1.plot(group['Delay Time       (µs)'], group['Correlation'], label=name, linewidth=2)
# ax1.legend(loc='upper right')

# plt.figure(1)
# plt.subplot(211)
ax1.plot(time_list_us, g_t_scale_list, linestyle='None', marker='o',  alpha=0.3)
# plt.axis([10**(-7), 1, 0, 1])
plt.title('Predicted Autocorrelation using Lognormal Distributions')
plt.xlabel('Time (us)')
plt.ylabel('$g^{(1)}$(t)')
plt.xscale('log')
plt.tight_layout(0.2)
# plt.grid(True, which='both')

# This section plots the Mass squared * scattering intensity, vs r
pf2 = pow(4*pi*rho, 2) / (q**6)
r_space = np.arange(1, 10000, step=1)
msquared_pop = []
pop = []
for r in r_space:
    pop.append(lognorm_sum(r, distribution_list))
    msquared_pop.append(pf2 * compute_denom_lognorm(r, distribution_list))

cdf_np_mat = np.zeros((len(distribution_list), r_space.size))
norm_factor = 0
for dist_index, dist_curve in enumerate(distribution_list):
    cdf_np_mat[dist_index] = dist_curve[0] * lognorm.cdf(r_space, s=np.log(dist_curve[2]), scale=dist_curve[1], loc=0)
    norm_factor += dist_curve[0]

# determine the r quantiles for the cdfs 
cdf = cdf_np_mat.sum(axis=0) / norm_factor
cdf_10 = (np.abs(cdf-0.1)).argmin()
cdf_25 = (np.abs(cdf-0.25)).argmin()
cdf_50 = (np.abs(cdf-0.5)).argmin()
cdf_75 = (np.abs(cdf-0.75)).argmin()
cdf_90 = (np.abs(cdf-0.9)).argmin()

# set the locations of the labels for the cdf graph
pop_label = []
msquared_label = []
label = [[r_space[cdf_10]/10, cdf[cdf_10], '10% CDF at {0} nm'.format(r_space[cdf_10]/10)],
         [r_space[cdf_25]/10, cdf[cdf_25], '25% CDF at {0} nm'.format(r_space[cdf_25]/10)],
         [r_space[cdf_50]/10, cdf[cdf_50], '50% CDF at {0} nm'.format(r_space[cdf_50]/10)],
         [r_space[cdf_75]/10, cdf[cdf_75], '75% CDF at {0} nm'.format(r_space[cdf_75]/10)],
         [r_space[cdf_90]/10, cdf[cdf_90], '90% CDF at {0} nm'.format(r_space[cdf_90]/10)]]

pop_modes_list = argrelextrema(np.array(pop), np.greater)[0]
msquared_pop_modes_list = argrelextrema(np.array(msquared_pop), np.greater)[0]
# label the peaks of the population graph
for mode in pop_modes_list:
    label.append([r_space[mode]/10, cdf[mode],
                  'Pop Mode at {0} nm'.format(r_space[mode]/10)])
    pop_label.append([r_space[mode]/10, pop[mode],
                      'Pop Mode at {0} nm'.format(r_space[mode]/10)])
# label the peaks of the msquared pop graph
for mode in msquared_pop_modes_list:
    if msquared_pop[mode] < msquared_pop[msquared_pop_modes_list[0]] / 2:
        continue
    msquared_label.append([r_space[mode] / 10, msquared_pop[mode],
                           'M**2 * Pop Mode at {0} nm'.format(r_space[mode] / 10)])
    label.append([r_space[mode] / 10, cdf[mode],
                  'M**2 * Pop Mode at {0} nm'.format(r_space[mode] / 10)])

plt.figure(2)
plt.plot(r_space/10, pop, linestyle='None', marker='o', alpha=0.3)
plt.xscale('log')
plt.title('Population Distribution Radii Distributions')
plt.xlabel('Particle Radius (r) (nm)')
plt.ylabel('N(r)')
plt.ylim(0, 1.2*np.max(pop))
for val in pop_label:
    plt.annotate(
        val[2],
        xy=(val[0], val[1]), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
    )
plt.tight_layout(0.2)

plt.figure(3)
plt.plot(r_space/10, msquared_pop, linestyle='None', marker='o', alpha=0.3)
plt.xscale('log')
plt.title('Mass-squared weighted Radii Distributions')
plt.xlabel('Particle Radius (r) (nm)')
plt.ylabel('$M^2P(r)N(r)$')
plt.ylim(0, 1.2*np.max(msquared_pop))
for val in msquared_label:
    plt.annotate(
        val[2],
        xy=(val[0], val[1]), xytext=(-20, -20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
    )
plt.tight_layout(0.2)

plt.figure(4)
plt.plot(r_space/10, cdf, marker='o', alpha=0.3)
plt.xscale('log')
for val in label:
    plt.annotate(
        val[2],
        xy=(val[0], val[1]), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
    )
plt.title('Cumulative Distribution Function')
plt.xlabel('Particle Radius (r) (nm)')
plt.ylabel('Cumulative Distribution Function')
plt.ylim(0, 1.1)
plt.tight_layout(0.2)
plt.show()

