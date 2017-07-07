#! /usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'William Isley III'

############################################################################################
#                                                                                          #
# VERSION :: 0.1.1  March 2017                                                             #
#                                                                                          #
# AUTHORS :: William Christian Isley III                                                   #
#                                                                                          #
#                                                                                          #
# DISCLAIMER :: This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by the Free Software   #
# Foundation, either version 3 of the License, or any later version                        #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY #
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE #
# See the GNU General Public License for more details.                                     #
#                                                                                          #
############################################################################################

# constants
MW = 0.748                 # kg/mol
rho = 2241                 # kg/m3
Na = 6.0225 * pow(10, 23)  # avogadro's number
MAX = 28                   # number of bins
pi = np.pi

# setting up the volumes of the bins
vol_1 = MW/(rho*Na)
vol_p = np.empty([MAX])
diameter_p = np.empty([MAX])
mass_p = np.empty([MAX])
q = pow(10.0, (8.0/(MAX-2)))  # Calculation of geometric spacing factor which depends on the number of nodes.

for i in range(1, MAX):
    vol_p[i] = vol_1*pow(q, i)
    diameter_p[i] = pow((6.0*vol_p[i]/pi), 0.333333)
    mass_p[i] = rho * vol_p[i]

vol_p[0] = vol_1
diameter_p[0] = pow((6.0*vol_1/pi), 0.333333)
mass_p[0] = rho * vol_1

print(df.columns.values)
const_index_list = ['N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10' 'N11' 'N12' 'N13' 'N14'
                    'N15' 'N16' 'N17' 'N18' 'N19' 'N20' 'N21' 'N22' 'N23' 'N24' 'N25' 'N26'
                    'N27' 'N28']
df_particle_const = pd.DataFrame(index=const_index_list)
df_particle_const['Volume'] = vol_p
df_particle_const['Diameter'] = diameter_p
list_to_broadcast = np.empty([1])
list_to_broadcast = np.append(np.append(np.empty([1]), df_particle_const.loc['N2':'N28', 'Volume'].values, axis=0),
                              np.empty([1]), axis=0)


file_path_name = 'C:/Users/isle132/Documents/NGDE/Data/post-process/nucl-coag-surf-bins.dat'
df = pd.read_csv(file_path_name, header=0, index_col=0, sep='\s+')
# df['Vavg'] = df.mul(list_to_broadcast, axis=1).sum(axis=1) / df['Ntot']
df['Ntot'] = df.iloc[:, 2:].sum(axis=1)
df['Diameter_avg'] = pow((6.0*(df.mul(list_to_broadcast, axis=1).sum(axis=1) / df['Ntot'])/pi), 0.333333) * 10**9
print(df.head())

x_time = df.index.values
y_diametervtime = df['Diameter_avg'].values

plt.plot(x_time, y_diametervtime,  alpha=0.8)
plt.xlabel('Time (s)')
plt.ylabel('Avg Particle Diameter (nm)')
plt.title('Particle Diameter Growth over Time')
plt.show()

#xs = df_transpose['Diameter'].loc['N2':'N28'].values
#ys = df_transpose.iloc[1:28, :100].values
#print(xs, ys)
#plt.plot(xs, ys,  alpha=0.8)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('Particle Diameter')
#plt.ylabel('Particle Population')
#plt.title('Particle Growth over Time')
#plt.show()
