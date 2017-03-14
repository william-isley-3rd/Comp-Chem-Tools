#! /usr/bin/env python
import argparse
import re
import math

__author__ = 'William Isley III'

############################################################################################
#                                                                                          #
# VERSION :: 1.4.2  October 2015                                                           #
#                                                                                          #
# AUTHORS :: William Christian Isley III                                                   #
#                                                                                          #
#                                                                                          #
# SOFTWARE :: This script was created using Python 2.7                                     #
#             Depends on argparse, re, and math packages.
#             This script is available at github.com/william-isley-3rd/Comp-Chem-Tools     #
#                                                                                          #
# INPUT :: This script has been tested to function with Gaussian 09 frequency output       #
#          files, as of G09.d01\                                                           #
#                                                                                          #
# VARIABLE SPECIFICATIONS: Default Temperature is set to 298.15 K                          #
# The Quasiharmonic frequency corrections are performed for frequencies below 50 cm^-1     #
# Both the entropic and Zero Point Energy (ZPE) terms are corrected                        #
#                                                                                          #
# Update v1.4.2 Change Log: Updated print statements, Improved code readability, added     #
# additional commentary text to script                                                     #
#                                                                                          #
# Update v1.4.1 Change Log: Updated the printing formatting, no functionality change.      #
#                                                                                          #
# Update v1.4 Change Log: Bug Fix: Fixed the temperature change parameter to function      #
# for all frequencies in the output file, not just the ones being corrected                #
# Additionally included the previously missing thermodynamic changes to enthalpy and       #
# entropy for translations, and rotation. Added missing R(T2 - T1) to enthalpy and         #
# missing S_old*(T2- T1) term to free energy.                                              #
#                                                                                          #
# Update v1.3 Change Log: The option to perform Quasiharmonic Correction on Transition     #
# state files was added, use flag --isTS . The option for a constant factor of frequency   #
# scaling has been added, see --fs option                                                  #
#                                                                                          #
# Update v1.2: Change log: update the parsing of frequencies so undefined frequencies are  #
# not appended to the array                                                                #
#                                                                                          #
# Update v1.1: Change log: updated the reporting of variables to correctly reflect gibbs,  #
# enthalpy and entropy                                                                     #
#                                                                                          #
# DISCLAIMER :: This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by the Free Software   #
# Foundation, either version 3 of the License, or any later version                        #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY #
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE #
# See the GNU General Public License for more details.                                     #
#                                                                                          #
############################################################################################

# CONSTANTS
RCONST = 1.9858555 / 1000
KBOLTZ = 0.69503476
HARTREE_TO_KCAL = 627.5095

# USER SPECIFIED PARAMETERS
# Electronic Spin
SPIN_2 = 0.0
SPIN_2_IDEAL = 0.0

# Thermal Variables
TEMP_K = 298.15
TEMP_K_DEFAULT = 298.15

# the new replacement frequency
# default set to 50 wave numbers
# can be changed by parameter
TARGET_FREQ = 50.00000

# Frequency Scaling factor for all frequencies (can be changed by a flag)
FREQ_SCALING = 1.000

# READ IN PARAMETERS FROM OUTPUT
# Thermodynamic Quantities needed for making corrections
FREQ_ALL = []  # contains all frequencies found in outfile, wavenumber units
THERMAL_VIB_F_CORR = 0.0  # tracks changes in vib E_v from low freq, kcal/mol
THERMAL_VIB_T_CORR = 0.0  # tracks changes in vib E_v from new T, kcal/mol
ZPE_CORR = 0.0  # tracks changes in ZPE, a.u.
ENTHALPY_T_CORR = 0.0  # tracks changes in enthalpy from T, kcal/mol
ENTROPY_F_CORR = 0.0  # tracks changes in ettropy from low freq corr, kcal/mol
ENTROPY_T_CORR = 0.0  # tracks changes in entropy from T changes, kcal/mol

# variables read in from gaussian output file
# Units of energy are in a.u.
ZPE_READ = None
ZPE_N_ELECTRONIC = None
ENTHALPY_CORR_READ = None
THERMAL_ENTHALPY = None
GIBBS_READ = None
THERMAL_FREE_ENERGY = None
# Flags from G09 Outfile
TS_ON = False
IMAG_FLAG = False

# variables for final corrections
# values are defined below
# GIBBS_CORR = 0.0  # tracks changes in gibbs corr, kcal/mol
# ENTHALPY_TOTAL = 0.0  # final enthalpy correction, a.u.
# ENTHALPY_ALL_CORR = 0.0  # all corrections to enthalpy, kcal/mol
# ENTHALPY_T_CORR_ALL = 0.0  # all correction to enthalpy from T, kcal/mol
# ENTROPY_OLD = 0.0  # tracks old entropy for changes in T, needed for final gibbs corr.
# GIBBS_TOTAL = 0.0  # final gibbs energy, a.u.

# Read Arguments Passed into Script
HELP_TEXT = """\nThis script corrects for low-frequencies.
It is intended for use to correct frequencies which are unphyiscally low
as a result of improper treatment with the Quasiharmonic Oscillator
Approximation. Use the options below to parse your output."""
parser = argparse.ArgumentParser(description=HELP_TEXT)
parser.add_argument('--temp', "-t",
                    help='Correct Freq to specified Temperature (K), '
                         '298.15 K is the Default',
                    nargs='?',
                    type=float,
                    default=298.15)
parser.add_argument('--isTS',
                    help='Flag freq file as TS calculation, ignore first frequency',
                    action="store_true",
                    default=False)
parser.add_argument("--corrFreq", "-f",
                    help='Change Target Corrected Freq (cm**-1) to new value, '
                         'Default target is 50 cm**-1 (recommended)',
                    type=float,
                    default=50.00000)
parser.add_argument("--freqScale", "-s",
                    help='Change Frequency scaling factor. Scales all frequencies by s. Default 1.0000 (no scaling)'
                         ' scaled frequency corrections are reported with changing temperature',
                    type=float,
                    default=1.0000)
parser.add_argument("freq_file",
                    help='G09 File Name with Frequencies')
parser.parse_args()
args = parser.parse_args()
if args.temp != '298.15':
    TEMP_K = args.temp
if args.isTS:
    TS_ON = True
if args.corrFreq != '50.00000':
    TARGET_FREQ = args.corrFreq
if args.freqScale != '1.0000':
    FREQ_SCALING = args.freqScale

# Open up G09 OUTFILE for reading
FILE = open(args.freq_file, mode='r')
# Regular Expressions to find in G09 OUTFILE
for line in FILE:
    if re.search('Frequencies\s+--\s+', line):
        line_freq = re.split('\s+', line)
        for x in range(3, len(line_freq)):
            if line_freq[x]:
                FREQ_ALL.append(float(line_freq[int(x)]))
        # if line_freq[3]:
        #     FREQ_ALL.append(float(line_freq[3]))
        # if line_freq[4]:
        #     FREQ_ALL.append(float(line_freq[4]))
        # if line_freq[5]:
        #     FREQ_ALL.append(float(line_freq[5]))
    elif re.search('Zero-point correction=\s*([-\w\.]*)', line):
        ZPE_READ = float(re.search('Zero-point correction=\s+([-\w\.]*)', line).group(1))
    elif re.search('zero-point Energies=\s*([-\w\.]*)', line):
        ZPE_N_ELECTRONIC = float(re.search('zero-point Energies=\s+([-\w\.]*)', line).group(1))
    elif re.search('Thermal correction to Enthalpy', line):
        ENTHALPY_CORR_READ = float(re.search('Thermal correction to Enthalpy=\s+([-\w\.]*)', line).group(1))
    elif re.search('Thermal correction to Gibbs Free Energy=\s+([-\w\.]*)', line):
        GIBBS_READ = float(re.search('Thermal correction to Gibbs Free Energy=\s+([-\w\.]*)', line).group(1))
    elif re.search('thermal Enthalpies=\s+([-\w\.]*)', line):
        THERMAL_ENTHALPY = float(re.search('thermal Enthalpies=\s+([-\w\.]*)', line).group(1))
    elif re.search('thermal Free Energies=\s+([-\w\.]*)', line):
        THERMAL_FREE_ENERGY = float(re.search('thermal Free Energies=\s+([-\w\.]*)', line).group(1))
    elif re.search(' <Sx>= [\w\.]* <Sy>= [\w\.]* <Sz>= [\w\.]* <S\*\*2>= ([\w\.]*) S= [\w\.]', line):
        SPIN_2 = float(re.search(' <Sx>= [\w\.]* <Sy>= [\w\.]* <Sz>= [\w\.]* <S\*\*2>= ([\d\.]*) S= [\w\.]', line)
                       .group(1))
    elif re.search('Charge\s*=\s*[\-\w]*\s*Multiplicity\s*=\s*(\d*)', line):
        SPIN = float(re.search('Charge\s*=\s*[\-\w]*\s*Multiplicity\s*=\s*(\d*)', line).group(1))
        SPIN_2_IDEAL = (SPIN - 1) / 2 * ((SPIN - 1) / 2 + 1)
    elif re.search('Temperature\s*([\d\.]*)\s*Kelvin', line):
        TEMPORARY_TEMP = float(re.search('Temperature\s*([\d\.]*)\s*Kelvin', line).group(1))
        if TEMPORARY_TEMP.__eq__(TEMP_K_DEFAULT):
            continue
        else:
            TEMP_K_DEFAULT = TEMPORARY_TEMP
            print("\nTemperature Default Changed to {} \n\n", format(TEMP_K_DEFAULT))
FILE.close()
print("Done Reading {} \n\n".format(args.freq_file))

if not THERMAL_FREE_ENERGY:
    print("Frequencies Not Found")
    exit()

ELECTRONIC_E = ZPE_N_ELECTRONIC - ZPE_READ

print("Assuming NON-LINEAR MOLECULE for Rotational Temperature Corrections")
print("Assuming RIGID ROTOR Approximation\n")
print("Before Low Frequency Corrections")
print("The electronic energy (Hartrees):          {}", format(ELECTRONIC_E))
print("The ZPE and Electronic E (Hartrees):       {}", format(ZPE_N_ELECTRONIC))
print("The Thermal Enthalpy (Hartrees):           {}", format(THERMAL_ENTHALPY))
print("The Thermal Free Energy (Hartrees):        {}", format(THERMAL_FREE_ENERGY))
print("The S**2 value is:                         {}", format(SPIN_2))
print("with the ideal:                            {}", format(SPIN_2_IDEAL))
print("\n")
print("The Zero Point Corrections (Hartrees):     {}", format(ZPE_READ))
print("Thermal correction to Enthalpy(Hartrees):  {}", format(ENTHALPY_CORR_READ))
print("Thermal correction to Gibbs(Hartrees):     {}", format(GIBBS_READ))
print("\n")

if TS_ON:
    # if a transition state calculation is not then, remove first frequency from array
    print("Frequency not corrected for TS optimization:")
    print("Frequency (cm**-1): {}", format(FREQ_ALL[0]))
    FREQ_ALL.pop(0)

# begin frequency replacements for frequencies less than TARGET_FREQ
print("Frequencies in need of correction: ")
# see gaussian white pages for equations of entropy and vibrational computations
for freq in FREQ_ALL:
    if 0.0 < freq < TARGET_FREQ:
        print("Frequency (cm**-1): {}", format(freq))

        # compute frequency factors used for both entropy and enthalpy corrections
        # F variables are for changing frequencies
        # T variables are for scaled temperatures
        PHI_OVER_T_OLD_FREQ = freq * (1 / KBOLTZ) * (1 / TEMP_K_DEFAULT)
        PHI_OVER_T_NEW_FREQ = TARGET_FREQ * (1 / KBOLTZ) * (1 / TEMP_K_DEFAULT)
        PHI_OVER_T_OLD_TEMP = PHI_OVER_T_NEW_FREQ
        PHI_OVER_T_NEW_TEMP = FREQ_SCALING * TARGET_FREQ * (1 / KBOLTZ) * (1 / TEMP_K)

        # compute corrections to entropy for low frequencies (most difference should be from this term)
        S_Vib_OLD_F = RCONST * (PHI_OVER_T_OLD_FREQ / (math.exp(PHI_OVER_T_OLD_FREQ) - 1) -
                                math.log(1 - math.exp(-PHI_OVER_T_OLD_FREQ)))
        S_Vib_NEW_F = RCONST * (PHI_OVER_T_NEW_FREQ / (math.exp(PHI_OVER_T_NEW_FREQ) - 1) -
                                math.log(1 - math.exp(-PHI_OVER_T_NEW_FREQ)))
        S_FREQ_CORR = S_Vib_NEW_F - S_Vib_OLD_F
        ENTROPY_F_CORR += S_FREQ_CORR

        # compute corrections to entropy for low frequencies at different T (most difference should be from this term)
        S_Vib_OLD_TEMP = S_Vib_NEW_F
        S_Vib_NEW_TEMP = RCONST * (PHI_OVER_T_NEW_TEMP / (math.exp(PHI_OVER_T_NEW_TEMP) - 1) -
                                   math.log(1 - math.exp(-PHI_OVER_T_NEW_TEMP)))
        S_TEMP_CORR = S_Vib_NEW_TEMP - S_Vib_OLD_TEMP
        ENTROPY_T_CORR += S_TEMP_CORR

        # compute correction to ZPE for having low frequencies
        E_Vib_OLD_FREQ = RCONST * (freq * (1 / KBOLTZ) * (0.5 + 1 / (math.exp(PHI_OVER_T_OLD_FREQ) - 1)))
        E_Vib_NEW_FREQ = RCONST * (TARGET_FREQ * (1 / KBOLTZ) * (0.5 + 1 / (math.exp(PHI_OVER_T_NEW_FREQ) - 1)))
        E_Vib_FREQ_CORR = E_Vib_NEW_FREQ - E_Vib_OLD_FREQ
        ZPE_CORR += RCONST * (FREQ_SCALING * TARGET_FREQ * (1 / KBOLTZ) * 0.5) - (RCONST * (freq * (1 / KBOLTZ) * 0.5))
        THERMAL_VIB_F_CORR += E_Vib_FREQ_CORR

        # compute correction to ZPE for having new Temperature
        E_VIB_OLD_TEMP = E_Vib_NEW_FREQ
        E_VIB_NEW_TEMP = RCONST * (FREQ_SCALING * TARGET_FREQ * (1 / KBOLTZ) *
                                   (0.5 + 1 / (math.exp(PHI_OVER_T_NEW_TEMP) - 1)))
        E_Vib_TEMP_CORR = E_VIB_NEW_TEMP - E_VIB_OLD_TEMP
        THERMAL_VIB_T_CORR += E_Vib_TEMP_CORR

    # Loops over imaginary Frequencies 
    # Replaces with TARGET_FREQ unless isTS is specified
    elif freq < 0.0:
        # imaginary frequency corrections added here
        print("Frequency (cm**-1): {}  WARNING Imaginary Frequecy Found! Recommended further optimization. WARNING",
              format(freq))

        # compute values needed for both entropy and enthalpy corrections
        PHI_OVER_T_NEW_FREQ = TARGET_FREQ * (1 / KBOLTZ) * (1 / TEMP_K_DEFAULT)
        PHI_OVER_T_NEW_TEMP = FREQ_SCALING * TARGET_FREQ * (1 / KBOLTZ) * (1 / TEMP_K)

        # compute correction to entropy for imaginary freq 
        S_Vib_OLD_F = 0.0
        S_Vib_NEW_F = RCONST * (PHI_OVER_T_NEW_FREQ / (math.exp(PHI_OVER_T_NEW_FREQ) - 1) -
                                math.log(1 - math.exp(-PHI_OVER_T_NEW_FREQ)))
        S_FREQ_CORR = S_Vib_NEW_F - S_Vib_OLD_F
        ENTROPY_F_CORR += S_FREQ_CORR

        # compute corrections to entropy for low frequencies at different T (most difference should be from this term)
        S_Vib_NEW_TEMP = RCONST * (PHI_OVER_T_NEW_TEMP / (math.exp(PHI_OVER_T_NEW_TEMP) - 1) -
                                   math.log(1 - math.exp(-PHI_OVER_T_NEW_TEMP)))
        S_TEMP_CORR = S_Vib_NEW_TEMP - S_Vib_NEW_F
        ENTROPY_T_CORR += S_TEMP_CORR

        # compute correction to ZPE for having low frequencies
        E_Vib_OLD_FREQ = 0.0
        E_Vib_NEW_FREQ = RCONST * (TARGET_FREQ * (1 / KBOLTZ) * (0.5 + 1 / (math.exp(PHI_OVER_T_NEW_FREQ) - 1)))
        E_Vib_FREQ_CORR = E_Vib_NEW_FREQ - E_Vib_OLD_FREQ
        ZPE_CORR += RCONST * (FREQ_SCALING * TARGET_FREQ * (1 / KBOLTZ) * 0.5) - (RCONST * (freq * (1 / KBOLTZ) * 0.5))
        THERMAL_VIB_F_CORR += E_Vib_FREQ_CORR

        # compute correction to ZPE for having new Temperature
        E_VIB_OLD_TEMP = E_Vib_NEW_FREQ
        E_VIB_NEW_TEMP = RCONST * (FREQ_SCALING * TARGET_FREQ * (1 / KBOLTZ) *
                                   (0.5 + 1 / (math.exp(PHI_OVER_T_NEW_TEMP) - 1)))
        E_Vib_TEMP_CORR = E_VIB_NEW_TEMP - E_VIB_OLD_TEMP
        THERMAL_VIB_T_CORR += E_Vib_TEMP_CORR

        IMAG_FLAG = True

    # For all Other Frequencies 
    elif freq >= TARGET_FREQ and ((TEMP_K != TEMP_K_DEFAULT) or (FREQ_SCALING != 1.0)):
        # compute needed for both entropy and enthalpy corrections
        PHI_OVER_T_OLD_TEMP = freq * (1 / KBOLTZ) * (1 / TEMP_K_DEFAULT)
        PHI_OVER_T_NEW_TEMP = FREQ_SCALING * freq * (1 / KBOLTZ) * (1 / TEMP_K)

        # compute corrections to entropy for low frequencies at different T (most difference should be from this term)
        S_Vib_OLD_TEMP = RCONST * (PHI_OVER_T_OLD_TEMP / (math.exp(PHI_OVER_T_OLD_TEMP) - 1) -
                                   math.log(1 - math.exp(-PHI_OVER_T_OLD_TEMP)))
        S_Vib_NEW_TEMP = RCONST * (PHI_OVER_T_NEW_TEMP / (math.exp(PHI_OVER_T_NEW_TEMP) - 1) -
                                   math.log(1 - math.exp(-PHI_OVER_T_NEW_TEMP)))
        S_TEMP_CORR = S_Vib_NEW_TEMP - S_Vib_OLD_TEMP
        ENTROPY_T_CORR += S_TEMP_CORR
        # compute correction to ZPE for having new Temperature
        E_VIB_OLD_TEMP = RCONST * (freq * (1 / KBOLTZ) * (0.5 + 1 / (math.exp(PHI_OVER_T_OLD_TEMP) - 1)))
        E_VIB_NEW_TEMP = RCONST * (FREQ_SCALING * freq * (1 / KBOLTZ) * (0.5 + 1 / (math.exp(PHI_OVER_T_NEW_TEMP) - 1)))
        E_Vib_TEMP_CORR = E_VIB_NEW_TEMP - E_VIB_OLD_TEMP
        ZPE_CORR += (RCONST * (FREQ_SCALING * freq * (1 / KBOLTZ) * 0.5) - (RCONST * (freq * (1 / KBOLTZ) * 0.5)))
        THERMAL_VIB_T_CORR += E_Vib_TEMP_CORR

####################################################################
# FINAL THERMAL COMPUTATIONS
# compute the enthalpy, entropy, and free energies final corrections
####################################################################

# Compute corrections from change in T flag
if TEMP_K != TEMP_K_DEFAULT:
    # if new temperature, we have to correct enthalpies and entropies of translations 
    # S_trans = R/2 * (5ln TEMP_K + 3lnM) - const
    translation_entropy_Tcorr = (RCONST * 5 / 2) * (math.log(TEMP_K) - math.log(TEMP_K_DEFAULT))
    # H_trans = 3/2 R TEMP_K + pV
    translation_enthalpy_Tcorr = 3 / 2 * RCONST * (TEMP_K - TEMP_K_DEFAULT)
    # add corrections to global variable
    ENTROPY_T_CORR += translation_entropy_Tcorr
    ENTHALPY_T_CORR += translation_enthalpy_Tcorr
    # if new temperature, we have to correct entropies of translations
    # only part that scales with temperature is 3/2RT
    rotation_entropy_Tcorr = (RCONST * 3 / 2) * (math.log(TEMP_K) - math.log(TEMP_K_DEFAULT))
    # H_rot = 3/2 R TEMP_K 
    rotation_enthalpy_Tcorr = 3 / 2 * RCONST * (TEMP_K - TEMP_K_DEFAULT)
    # add corrections to global variable
    ENTROPY_T_CORR += rotation_entropy_Tcorr
    ENTHALPY_T_CORR += rotation_enthalpy_Tcorr

    # Add RdTEMP_K corr for term (Internal Energy -> Enthalpy)
    ENTHALPY_T_CORR += RCONST * (TEMP_K - TEMP_K_DEFAULT)

#################################################
# ENTHALPY FINAL VALUES
# all thermal corrections (kcal/mol)
ENTHALPY_ALL_CORR = (THERMAL_VIB_F_CORR + THERMAL_VIB_T_CORR + ENTHALPY_T_CORR)
# All enthalpy corrections from changing T, kcal/mol
ENTHALPY_T_CORR_ALL = THERMAL_VIB_T_CORR + ENTHALPY_T_CORR
# Final Enthalpy total (a.u.)
ENTHALPY_TOTAL = THERMAL_ENTHALPY + (ENTHALPY_ALL_CORR / HARTREE_TO_KCAL)
# recompute the total internal thermal "correction"
THERMAL_CONTR_CORR = ENTHALPY_TOTAL - ELECTRONIC_E
####################################################
# ENTROPY
# Compute old Entropy in case of TEMP_K flag requested
ENTROPY_OLD = HARTREE_TO_KCAL * (-1 * (THERMAL_FREE_ENERGY - THERMAL_ENTHALPY) / TEMP_K_DEFAULT)
###################################################
# FREE ENERGY
# final gibbs correction (kcal/mol)
GIBBS_CORR = ENTHALPY_ALL_CORR - (TEMP_K * (ENTROPY_F_CORR + ENTROPY_T_CORR) + ENTROPY_OLD * (TEMP_K - TEMP_K_DEFAULT))
# recompute the total Gibbs Free Energy
GIBBS_TOTAL = THERMAL_FREE_ENERGY + (GIBBS_CORR / HARTREE_TO_KCAL)
# Final Gibbs correction (a.u.)
GIBBS_CONTR_CORR = GIBBS_TOTAL - ELECTRONIC_E
#####################################################
# ZPE
# recompute ZPE with low freq corrections
ZPE_TOTAL = ZPE_READ + (ZPE_CORR / HARTREE_TO_KCAL)
################################################# 
######################################################################
# END COMPUTATIONS
######################################################################
# Format numbers for output (fixed 6 decimal places, 12 spaces)
print("All Frequencies Analyzed\n")
print("THERMAL CORRECTION AT {} K\n", format(TEMP_K))
print("The magnitude of changes from low frequencies is the following:")
print("The low Freq ZPE Correction (kcal/mol):                {:18.6f}", format(ZPE_CORR))
print("The Low Freq E_v Correction (kcal/mol):                {:18.6f}", format(THERMAL_VIB_F_CORR))
print("The Low Freq Entropic Correction (kcal/mol*K):         {:18.6f}", format(ENTROPY_F_CORR))
print("\n")
print("Corrections to thermodynamics to reach {} K from {} K are:", format(TEMP_K, TEMP_K_DEFAULT))
print("The new Temperature E_v Correction (kcal/mol):           {:18.6f}", format(THERMAL_VIB_T_CORR))
print("The new Temperature H Correction (kcal/mol):             {:18.6f}", format(ENTHALPY_T_CORR_ALL))
print("The new Temperature Entropic Correction (kcal/mol*K):    {:18.6f}", format(ENTROPY_T_CORR))
print("\n")
print("Total Requested Thermal and Low Freq Corrections\n")
print("The Total Enthalpy Correction (kcal/mol):                 {:18.6f}", format(ENTHALPY_ALL_CORR))
print("The Total Free Energy Correction (kcal/mol):              {:18.6f}", format(GIBBS_CORR))
print("\n")
print("The Thermal and Low Freq Corrected Energies:\n")
print("The electronic energy (Hartrees):                         {:18.6f}", format(ELECTRONIC_E))
print("The Corrected Thermal Enthalpy (Hartrees):                {:18.6f}", format(ENTHALPY_TOTAL))
print("The Corrected Thermal Free Energy (Hartrees):             {:18.6f}", format(GIBBS_TOTAL))
print("\n")
print("The Corrected Zero Point Vibrational Energy (Hartrees):   {:18.6f}", format(ZPE_TOTAL))
print("The Corrected Thermal Enthalpy Contribution (Hartrees):   {:18.6f}", format(THERMAL_CONTR_CORR))
print("The Corrected Thermal Free Energy Contribution (Hartrees):{:18.6f}", format(GIBBS_CONTR_CORR))

if IMAG_FLAG:
    print("WARNING IMAGINARY FREQUENCIES WERE FOUND!!")

print("\n")
