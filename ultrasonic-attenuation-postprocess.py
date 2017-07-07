#! /usr/bin/env python
import pandas as pd
import numpy as np
import cmath
from scipy.stats import lognorm
import matplotlib.pyplot as plt

__author__ = 'William Isley III'

#############################################################################################
#                                                                                           #
# VERSION :: 0.1.2  July  2017                                                              #
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
# This software will compute the High Resolution Ultra-sonic Spectra (HRUS) for a given     #
# particle size distribution. This is currently specified as a sum of log-normal            #
# distributions that are specifed by the user in the code below. This is then compared to   #
# experimentally measured HRUS data.  The code will randomly sample about the specified     #
# distribution and print out the set of best fit particle size distibutions.                #
#############################################################################################
#                                    ------------                                           #
# 				                     NOMENCLATURE                                           # 
# 				                     ------------                                           # 
# __________________________________________________________________________________________#
# Variables            	            Meaning (Units)                                         #
# __________________________________________________________________________________________#
# alpha_exp[i]         	Attenuation Coefficient
# c_exp[i]              Ultrasonic Velocity
# alpha_tot[i]          Attenuation Coefficient for total Predicted
# alpha_I[i]            Intrinsic Absorption Attenuation Coeff
# alpha_S[i]            Scattering Losses (Monopole / Dipole) Attenuation Coeff
# alpha_T[i]            Thermal Absorption Attenuation Coeff
# alpha_VI[i]           Viscoinertial Absorption Attenuation Coeff
# phi                   Disperse Phase Volume Fraction
# f[i]                  Frequency of ultrasonic Probe
# r_j[j]                Radius of particle j
# rho_susp              Density of Suspended Phase  (kg/m3)
# rho_solv              Density of Solvent Phase    (kg/m3)
# c_solv                Ultrasonic Velocity of solvent phase (m/s)
# c_susp                Ultrasonic Velocity of suspended phase (m/s)
# eta_solv              Shear Viscosity of solvent Phase (Pa*s)
# Cp_susp               Heat Capacity of Suspended Phase (J/[kg*K])
# Cp_solv               Heat Capacity of Solvent Phase   (J/[kg*K])
# tau_susp              Thermal Conductivity of suspended Phase (J/[m*sec*K])
# tau_solv              Thermal Conductivity of Solvent Phase   (J/[m*sec*K])
# beta_susp             Thermal Expansion Factor of Suspended Phase (1/K) Coefficient of Volume Expansion
# beta_solv             Thermal Expansion Factor of Solvent Phase   (1/K)
# alpha_susp            Attenuation Coefficient for Suspended Phase
# alpha_solv            Attenuation Coefficient for Solvent Phase
# k_solv                Complex Propogation for Solvent Phase (1/m)
# k_susp                Complex Propogation for Suspended Phase (1/m)
# delta_v_susp          Viscous skin depth (m)
# delta_t_susp          Thermal skin depth (m)
# prop_const            Complex Propogation Constant
# A_0, A_1              Scattering Coefficients

pi = cmath.pi
np.random.seed(42)
ppf_list = np.arange(0.01, 0.99, step=0.01)


# Function definitions up here, explanation of how code works in main section below
def define_parameters():
    """ No Input
        Defines all global input parameters. This function must be called in order to work
        But only after the DataFrame containing the input Parameters has been read. """
    global rho_susp, rho_solv, c_susp, c_solv, eta_solv, Cp_susp, Cp_solv, tau_susp, tau_solv, \
        beta_susp, beta_solv, alpha_susp, alpha_solv, alpha_freq_power_susp, alpha_freq_power_solv, \
        kappa_solv, kappa_susp, T
    rho_susp = df_const['Suspended']['Density']
    rho_solv = df_const['Solvent']['Density']
    c_susp = df_const['Suspended']['Ultrasonic Velocity']
    c_solv = df_const['Solvent']['Ultrasonic Velocity']
    eta_solv = df_const['Solvent']['Viscosity']
    Cp_susp = df_const['Suspended']['Specific Heat Capacity']
    Cp_solv = df_const['Solvent']['Specific Heat Capacity']
    tau_susp = df_const['Suspended']['Thermal Conductivity']
    tau_solv = df_const['Solvent']['Thermal Conductivity']
    beta_susp = df_const['Suspended']['Thermal Expansion Coefficient']
    beta_solv = df_const['Solvent']['Thermal Expansion Coefficient']
    alpha_susp = df_const['Suspended']['Attenuation Coefficient']
    alpha_solv = df_const['Solvent']['Attenuation Coefficient']
    alpha_freq_power_susp = df_const['Suspended']['Attenuation Power Dependence']
    alpha_freq_power_solv = df_const['Solvent']['Attenuation Power Dependence']
    kappa_susp = 1 / (c_susp ** 2 * rho_solv)
    kappa_solv = 1 / (c_solv ** 2 * rho_solv)


# Computes the Ultrasonic Attenuation as a function of particle radius, frequency of probe and population of particle
def compute_alpha_of_r_freq(r_j, freq, num_j):
    """ r_j: float, freq: float, num_j: int
        radius of particle j, frequency probed, number / volume of particle j"""
    phi_j = 4 / 3 * pi * r_j**3 * float(num_j)

    k_solv = 2 * pi * freq / c_solv + alpha_solv * pow(freq, alpha_freq_power_solv) * 1j

    # For Computing Intrinsic Absorption Attentuation
    alpha_I = phi_j * alpha_susp * pow(freq, alpha_freq_power_susp) + \
              (1 - phi_j) * alpha_solv * pow(freq, alpha_freq_power_solv)

    # For Computing Scatting Loss Attenuation
    temp_alpha_S = 0.333333 * pow((kappa_susp-kappa_solv)/kappa_solv, 2) \
                    + pow((rho_susp - rho_solv)/(2*rho_susp+rho_solv), 2)
    alpha_S = 0.5 * phi_j * pow(k_solv, 4) * pow(r_j, 3) * temp_alpha_S

    # For Computing Thermal Absorption Attenuation
    temp_alpha_T = pow(1 - (beta_susp * rho_solv * Cp_solv)/(beta_solv * rho_susp * Cp_susp), 2)
    delta_t_susp = pow(2 * tau_susp / (rho_susp * Cp_susp * 2 * pi * freq), 0.5)
    delta_t_solv = pow(2 * tau_solv / (rho_solv * Cp_solv * 2 * pi * freq), 0.5)
    z_susp = (1 + 1j)*r_j / delta_t_susp
    z_solv = (1 + 1j)*r_j / delta_t_solv
    H_inv = 1 / (1 - 1j * z_solv) - (tau_solv / tau_susp) * cmath.tan(z_susp) / (cmath.tan(z_susp) - z_susp)
    H_actual = 1 / H_inv
    gamma = Cp_susp / Cp_solv
    alpha_T = 3 * phi_j * k_solv * H_actual * (gamma - 1) / (2 * pow(z_solv, 2)) * temp_alpha_T

    # For computing Viscoinertial Attenuation
    delta_v_solv = pow(2 * eta_solv / (rho_solv * 2 * pi * freq), 0.5)
    temp_s = 9 * delta_v_solv / (4 * r_j) * (1 + delta_v_solv / r_j)
    T_v = 0.5 + 9 * delta_v_solv / (4 * r_j)
    alpha_VI_denom = pow(rho_susp + T_v * rho_solv, 2) + pow(temp_s, 2) * pow(rho_solv, 2)
    alpha_VI = 0.5 * phi_j * k_solv * temp_s * (rho_susp - rho_solv) / alpha_VI_denom

    # Sum total alphas
    alpha_total = alpha_I + alpha_S + alpha_T + alpha_VI
    return alpha_I.real, alpha_S.real, alpha_T.real, alpha_VI.real, alpha_total.real,


def compute_scattering_coefficient(r_j, freq, T=293.1):
    """ r_j: float, freq: float
        Computes A0 and A1 scattering coefficients. This function takes a mean radius and
        ultrasonic frequency as input parameters. """

    # Compute skin depths and z factors
    k_solv = 2 * pi * freq / c_solv + alpha_solv * pow(freq, alpha_freq_power_solv) * 1j
    k_susp = 2 * pi * freq / c_susp + alpha_susp * pow(freq, alpha_freq_power_susp) * 1j
    delta_t_susp = pow(2 * tau_susp / (rho_susp * Cp_susp * 2 * pi * freq), 0.5)
    z_susp = (1 + 1j)*r_j / delta_t_susp
    delta_t_solv = pow(2 * tau_solv / (rho_solv * Cp_solv * 2 * pi * freq), 0.5)
    z_solv = (1 + 1j)*r_j / delta_t_solv
    delta_v_solv = pow(2 * eta_solv / (rho_solv * 2 * pi * freq), 0.5)

    # Unitless Factors
    H_inv = 1 / (1 - 1j * z_solv) - (tau_solv / tau_susp) * cmath.tan(z_susp) / (cmath.tan(z_susp) - z_susp)
    H_actual = 1 / H_inv
    temp_s = 9 * delta_v_solv / (4 * r_j) * (1 + delta_v_solv / r_j)
    T_v = 0.5 + 9 * delta_v_solv / (4 * r_j)

    # Compute A0 scattering Factor
    term1_pf = 0.33333333 * 1j * pow(k_solv, 3) * pow(r_j, 3)
    term1_1 = rho_solv * pow(k_susp, 2) / (rho_susp * pow(k_solv, 2))
    # Sonas
    zeta = T * pow(beta_solv, 2) * pow(c_solv, 2) / Cp_solv
    term2_pf = pow(k_solv, 3) * pow(r_j, 3) * zeta * H_actual / pow(z_solv, 2)
    term2_2 = beta_susp * rho_solv * Cp_solv / (beta_solv * rho_susp * Cp_susp)
    # McClements
    # term2_pf = pow(k_solv, 2) * r_j * c_solv * T_v * rho_solv * tau_solv * H_actual
    # term2_1 = beta_solv / (rho_solv * Cp_solv)
    # term2_2 = beta_susp / (rho_susp * Cp_susp)
    # A_0 = term1_pf * (term1_1 - 1) + term2_pf * pow(term2_1 - term2_2, 2)
    A_0 = term1_pf * (term1_1 - 1) + term2_pf * pow(1 - term2_2, 2)

    # Compute A1 scattering factor
    term3_nom = (rho_susp - rho_solv) * (1 + T_v + 1j * temp_s)
    term3_denom = (rho_susp + T_v * rho_solv + 1j * temp_s * rho_solv)
    A_1 = 0.33333333 * term1_pf * term3_nom / term3_denom

    return A_0, A_1


def compute_propagation_ratio(r_pop_list, r_factor, freq_sample):
    """ pop_r_dist: [(float, float)], freq: float
        Computes K, and returns c and alpha scattering coefficients. This function takes a list of tuples
        where the values are radius of cluster size k and population of cluster size k as first input, and the
        ultrasonic frequency as a second input parameter. """

    sum_N_f0 = 0.0
    sum_N2_f02_fpi2 = 0.0
    k_solv = 2 * pi * freq_sample / c_solv + alpha_solv * pow(freq_sample, alpha_freq_power_solv) * 1j
    prefactor = 1 / (1j * k_solv)
    for row in r_pop_list:
        (r_k, N_j) = (row[0] * r_factor, row[1])
        A_0, A_1 = compute_scattering_coefficient(r_k, freq_sample)
        f_0_j = prefactor * (A_0 + 3 * A_1)
        f_pi_j = prefactor * (A_0 - 3 * A_1)
        sum_N_f0 += N_j * f_0_j
        sum_N2_f02_fpi2 += (N_j ** 2) * (f_0_j ** 2 - f_pi_j ** 2)

    sum_1 = (4 * pi / (k_solv ** 2)) * sum_N_f0
    sum_2 = (4 * pi ** 2 / (k_solv ** 4)) * sum_N2_f02_fpi2
    ratio_prop_to_baseline = 1 + sum_1 + sum_2
    return ratio_prop_to_baseline


def compute_propagation_constant(prop_ratio, freq_sample, c_base, alphaf2_base, freq_base):
    k_baseline = 2 * pi * freq_base / c_base + alphaf2_base * 1j
    prop_const = pow(k_baseline ** 2 * prop_ratio, 0.5)
    speed_theory = 2 * pi * freq_sample / prop_const.real
    atten_theory = prop_const.imag
    return speed_theory, atten_theory


def compute_lognorm_dist(concentration, mean_r, sigma_r):
    list_r = lognorm.ppf(ppf_list, s=np.log(sigma_r), scale=mean_r, loc=0)
    list_numdensity = []
    for i in list_r:
        list_numdensity.append(lognorm.pdf(i, np.log(sigma_r), loc=0, scale=mean_r))

    list_numdensity /= np.sum(list_numdensity)
    list_numdensity *= concentration
    list_r_num = np.column_stack([list_r, list_numdensity])
    return list_r_num


def make_random_gaussian_list(mean, mean_width, sigma, sigma_width, samples=1000):
    random_list = np.column_stack((np.random.normal(mean, scale=mean_width, size=samples),
                                   np.random.normal(sigma, scale=sigma_width, size=samples)))
    # remove list values with sigma's <= 1, these will cause NaN errors in log normal
    random_list = random_list[
        np.logical_not(np.logical_or(random_list[:, 1] <= 1.0, random_list[:, 0] <= 1.0))
    ]
    return random_list


# This part reads in the excel file, and parses it into dataframse
path_in = 'C:/Users/isle132/Documents/NGDE/HRUS Data/From Forrest Test Data'
exp_data_file = 'gold_50nm_test_input.xls'
# exp_data_file = 'silica_140nm_test_input.xlsx'
path_out = 'C:/Users/isle132/Documents/NGDE/HRUS Data/From Forrest Test Data/'
outfile = 'gold_50nm_test_out3.xls'
# outfile = 'silica_140nm_test_out.xls'

# This part parses the input excel file, defines the parameters set
# and reads in the measurements
xl_exp_data = pd.ExcelFile(path_in+'/'+exp_data_file)
df_const = xl_exp_data.parse("constants", index_col=0)
df_const.columns = [col.strip() for col in df_const.columns]
define_parameters()
df_exp_us_data = xl_exp_data.parse("measurements")
df_exp_freq_data = xl_exp_data.parse("frequencies", header=0)
df_all_exp_data = pd.merge(df_exp_us_data, df_exp_freq_data, on='Peak Number')
df_all_exp_data.columns = [col.strip() for col in df_all_exp_data.columns]
df_all_exp_data['Peak Number'].apply(np.int64)
df_all_exp_data['Baseline Cell'].apply(np.int64)

# Here you set the initial test distribution to sample
# The code currently samples a gaussian distribution around the
# probe mean for both mean_radius, and sigma_radius
# the width determines how wide the sample is
# A broad scan will serve well as a first pass, then try to hone in on hotspots
# that have a low SSD value
test_mean = 25
test_mean_width = 10
test_sigma = 1.0
test_sigma_width = 1
col_for_df = ['r_mean', 'sigma']
seed_values = make_random_gaussian_list(test_mean, test_mean_width, test_sigma, test_sigma_width, samples=10000)
df_rmean_sigma_sample = pd.DataFrame(seed_values, columns=col_for_df)
test_concentration = 50 * 10 ** -3 * 10 ** -3 * 6.022 * 10 ** 23
r_min = 1
r_max = 4000
r_sample = int(4000)
r_to_meters = 10 ** (-10)

for index, measurement in df_all_exp_data.iterrows():
    baseline_flag = measurement['Baseline Cell']
    if baseline_flag == 1:
        freq_base = measurement['Cell 1 (kHz)']*1000
        freq_test = measurement['Cell 2 (kHz)']*1000
        c_base = measurement['U1 [m/s]']
        c_test = measurement['U2 [m/s]']
        alphaf2_base = measurement['N1 [1/m]']
        alphaf2_test = measurement['N2 [1/m]']
    elif baseline_flag == 2:
        freq_test = measurement['Cell 1 (kHz)']*1000
        freq_base = measurement['Cell 2 (kHz)']*1000
        c_test = measurement['U1 [m/s]']
        c_base = measurement['U2 [m/s]']
        alphaf2_test = measurement['N1 [1/m]']
        alphaf2_base = measurement['N2 [1/m]']
    time = measurement['TIME [s]']
    peak_flag = int(measurement['Peak Number'])
    c_temp_list = []
    alpha_temp_list = []
    ssd_temp_list = []
    time_list = []
    ratio_list = []
    for index2, sample in df_rmean_sigma_sample.iterrows():
        if 'ratio p'+str(peak_flag) not in df_rmean_sigma_sample.columns:
            r_num_list_all = compute_lognorm_dist(test_concentration, sample['r_mean'], sample['sigma'])
            prop_ratio_theory = compute_propagation_ratio(r_num_list_all, r_to_meters, freq_test)
            ratio_list.append(str(prop_ratio_theory))
        else:
            prop_ratio_theory = complex(sample['ratio p'+str(peak_flag)])
        # compute speed and attenuation
        c_theory, alpha_theory = compute_propagation_constant(prop_ratio_theory, freq_test,
                                                              c_base, alphaf2_base, freq_base)
        ssd = pow((c_theory - c_test), 2) + pow((alpha_theory - alphaf2_test), 2)
        time_list.append(time)
        c_temp_list.append(c_theory)
        alpha_temp_list.append(alpha_theory)
        ssd_temp_list.append(ssd)

    if 'ratio p'+str(peak_flag) not in df_rmean_sigma_sample.columns:
        df_rmean_sigma_sample['ratio p'+str(peak_flag)] = pd.Series(ratio_list)
    df_rmean_sigma_sample['time' + str(index)] = pd.Series(time_list)
    df_rmean_sigma_sample['c ' + str(index)] = pd.Series(c_temp_list)
    df_rmean_sigma_sample['alpha ' + str(index)] = pd.Series(alpha_temp_list)
    df_rmean_sigma_sample['SSD ' + str(index)] = pd.Series(ssd_temp_list)

ssd_cols = [col for col in df_rmean_sigma_sample if 'SSD' in col]
df_rmean_sigma_sample['SSD AVG'] = df_rmean_sigma_sample[ssd_cols].mean(axis=1)
reorder_cols = ['r_mean', 'sigma', 'SSD AVG'] + [col for col in df_rmean_sigma_sample if col not in ('r_mean', 'sigma', 'SSD AVG')]
df_rmean_sigma_sample = df_rmean_sigma_sample[reorder_cols]
df_rmean_sigma_sample.sort_values('SSD AVG').to_csv(path_out+outfile, sep='\t')

(x, y, z) = (df_rmean_sigma_sample['r_mean'], df_rmean_sigma_sample['sigma'], df_rmean_sigma_sample['SSD AVG'])
plt.hexbin(x, y, C=z, gridsize=200, cmap='gnuplot')
plt.ylabel('sigma (nm)')
plt.xlabel('r_mean (nm)')
cb = plt.colorbar()
cb.set_label('SSD AVG')
plt.show()

# r_num_list_all = compute_lognorm_dist(test_concentration, mean_radius, std_dev)
# plt.plot(r_num_list_all[:, 0], r_num_list_all[:, 1])
# plt.axis([0, 500, 0, 0.1])
# plt.show()
