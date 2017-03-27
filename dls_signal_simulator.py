from scipy import integrate
import numpy as np
from matplotlib import pyplot as plt

boltz = 1.3806 * 10**(-23)
N_avo = 6.022 * 10**23
T = 300
eta = 0.890 * 10**(-3)
n_refraction = 1.33
wavelength = 532*10  # in angstroms

pi = np.pi
q = (2 * pi * n_refraction / wavelength) * np.sin(pi/4)
c_gamma = boltz * T * pow(q * (10**10), 2) / (6 * pi * eta) * (10**10)  # in angstroms / s


def lognorm_f(x, triple_list):
    [pf, r_mean, sigma] = triple_list
    a = 1 / (2 * pow(np.log(sigma), 2))
    mu = np.log(r_mean)
    return pf / x * np.exp(-pow(np.log(x) - mu, 2) * a)


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
g_t_list = []
g_t_error_list = []

concentration = 1.0
N_0 = concentration
# distribution_list = [[N_0 / np.log(1.1), 5, 1.2]]
distribution_list = [[N_0 / np.log(1.3), 2500, 1.3],
                     [2*10**13*N_0 / np.log(1.5), 5, 1.5]]

g_t_previous = 0.0
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

plt.figure(1)
# plt.subplot(211)
plt.plot(time_list, g_t_list, linestyle='None', marker='o',  alpha=0.3)
plt.axis([10**(-7), 1, 0, 1])
plt.title('Predicted Autocorrelation using Lognormal Distributions')
plt.xlabel('Time (s)')
plt.ylabel('$g^(1)$(t)')
plt.xscale('log')
plt.tight_layout(0.2)

# g_t_list_bypart = []
# g_t_error_list_bypart = []
# for t in time_list:
#     t_dist_num = []
#     t_dist_denom = []
#     for distribution in distribution_list:
#         seg_denom = []
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 0, 1))
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 1, 50))
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 50, 200))
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 200, 500))
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 500, 1000))
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 1000, 2000))
#         seg_denom.append(integrate.quad(lambda x: compute_denom_lognorm_val(x, distribution), 2000, np.inf))
#         seg_denom_sum = np.sum(seg_denom, axis=0)
#         t_dist_denom.append(seg_denom_sum)
#
#         seg_num = []
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 0, 1))
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 1, 50))
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 50, 200))
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 200, 500))
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 500, 1000))
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 1000, 2000))
#         seg_num.append(integrate.quad(lambda x: compute_num_lognorm_val(x, t, distribution), 2000, np.inf))
#         seg_num_sum = np.sum(seg_num, axis=0)
#         t_dist_num.append(seg_num_sum)
#     t_dist_denom_sum = np.sum(t_dist_denom, axis=0)
#     t_dist_num_sum = np.sum(t_dist_num, axis=0)
#     g_t = t_dist_num_sum[0] / t_dist_denom_sum[0]
#     g_t_error = t_dist_num_sum[1] / t_dist_denom_sum[1]
#     if t_dist_num_sum[0] > t_dist_denom_sum[0]:
#         g_t = g_t_previous
#     else:
#         g_t_previous = g_t
#     g_t_list_bypart.append(g_t)
#     g_t_error_list_bypart.append(g_t_error)
#
# plt.subplot(212)
# plt.plot(time_list, g_t_list_bypart, linestyle='None', marker='o',  alpha=0.3)
# plt.axis([10**(-7), 1, 0, 1])
# plt.title('Predicted Autocorrelation using Lognormal Distributions by parts')
# plt.xlabel('Time (s)')
# plt.ylabel('$g^(1)$(t)')
# plt.xscale('log')
# plt.tight_layout(0.2)

rho = 2.165 * 1000
pf2 = pow(4*pi*rho, 2) / (q**6)
r_space = np.arange(1, 10000, step=1)
signal = []
pop = []
for r in r_space:
    pop.append(lognorm_sum(r, distribution_list))
    signal.append(pf2 * compute_denom_lognorm(r, distribution_list))

plt.figure(2)
plt.plot(r_space/10, signal, linestyle='None', marker='o', alpha=0.3)
plt.xscale('log')
plt.title('Mass-squared weighted Radii Distributions')
plt.xlabel('Particle Radius (r) (nm)')
plt.ylabel('$M^2P(r)N(r)$')
plt.tight_layout(0.2)
plt.show()

