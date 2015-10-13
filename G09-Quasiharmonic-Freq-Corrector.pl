#!/usr/bin/perl


############################################################################################
#                                                                                          #
# VERSION :: 1.4.2  October 2015                                                           #
#                                                                                          #
# AUTHORS :: William Christian Isley III                                                   #
#                                                                                          #
# DISCLAIMER :: This script is provided "as-is", without any warranties or support.        #
# The authors assume no liability or responsiblity for the use of this software.           #
#                                                                                          #
# SOFTWARE :: This script was created using Perl v 5.20                                    #
# This script has dependency on the pacakges, Getopt::Long                                 #
# This script is available at github.com/william-isley-3rd/Comp-Chem-Tools                 #
#                                                                                          #
# INPUT :: This script has been tested to function with Gaussian 09 frequency output       #
# files, as of G09.d01                                                                     #
#                                                                                          #
# VARIABLE SPECIFICATIONS: Default Temperature is set to 298.15 K                          #
# The Quasiharmonic frequency corrections are performed for frequencies below 50 cm^-1     #
# Both the entropic and Zero Point Eenergy (ZPE) terms are corrected                       #
#                                                                                          #
# Update v1.4.2 Change Log: Updated print statments, Improved code readability, added      #
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
# Update v1.2: Change log: update the parsing of frequencies so undefined frequences are   #
# not appended to the array                                                                #
#                                                                                          #
# Update v1.1: Change log: updated the reporting of variables to correctly reflect gibbs,  #
# enthalpy and entropy                                                                     #
#                                                                                          #
############################################################################################


use strict;
use warnings;
use Getopt::Long;

# constants
my $Rconst = 1.9858775/1000;
my $kb = 0.69503476;
my $hartreetokcal = 627.5095;
my $imagflag = 0;
my $help = "\nThis script corrects for low-frequencies
The default cutoff wavenumber is 50 cm^-1
The default temperature is 298.15 K
The default frequency scaling factor is set to 1.0

Options available for the script (include these before your the filename):
--temp=X or --t=X        ; Change the temperature used to X
--corrFreq=X or --cf=X   ; Change the cutoff frequency used to X
--freqScale=X or --fs=X  ; Scale all frequencies used in the calcualtion by X
--isTS                   ; Frequencies are for a transitions state, so the the first negative frequency will be ignored

Examples:
Change cutoff wavenumber to 100 cm^-1       ; perl freq_replacement_g09.pl --cf=100 filename.out
Change cutoff to 75 and temperature to 300  ; perl freq_replacement_g09.pl --cf=75 --t=300 filename.out
Run the script for a transition state       ; perl freq_replacement_g09.pl --isTS filename.out
Change frequency scaling to 0.978           ; perl freq_replacement_g09.pl --fs=0.978 filename.out

Any number of options can be combined as shown in the 2nd example";
my $printHelp=0;

# electronic spin 
my $S_2;
my $S_2_ideal;

# Thermal Variables (or constants)
my $T = -1.0;
my $Tdefault = 298.15;

# Frequency scaling factor for all frequencies (changed by flag)
my $freqScalingFactor = 1.000;

# variables for making the corrections
my @frequencies;                # contains all frequencies found in outfile, wavenumber units
my $Thermal_vib_Fcorr=0.0;      # tracsk changes in vib E_v from low freq, kcal/mol
my $Thermal_vib_Tcorr=0.0;      # tracks changes in vib E_v from new T, kcal/mol
my $ZPE_corr=0.0;               # tracks changes in ZPE, a.u.
my $enthalpyTcorr = 0.0;        # tracks changes in enthalpy from T, kcal/mol
my $entropyFcorr = 0.0;         # tracks changes in ettropy from low freq corr, kcal/mol
my $entropyTcorr = 0.0;         # tracks changes in entropy from T changes, kcal/mol
my $gibbscorr = 0.0;            # tracks changes in gibbs corr, kcal/mol

# the new replacement frequency 
# default set to 50 wave numbers
# can be changed by paramater
my $corrFreq = 50.00000;

# variables read in from gaussian output file
# Units of energy are in a.u.
my $electronicEnergy;
my $ZPEread;
my $ZPEelectronic;
my $enthalpyread;
my $thermalEnthalpy;
my $gibbsread;
my $thermalFreeEnergy;
# Flags from G09 Outfile
my $thermalflag=0;
my $TS_on=0;

# variables for final corrections
my $enthalpy_total;     # final enthalpy correction, a.u.
my $enthalpy_all_corr;  # all corrections to enthalpy, kcal/mol
my $enthalpy_T_corr;    # all correction to enthalpy from T, kcal/mol
my $oldentropy;         # tracks old entropy for changes in T, needed for final gibbs corr.
my $gibbs_total;        # final gibbs energy, a.u.


GetOptions ('help|h' => \$printHelp,
            'temp|t=f' => \$T,
            'corrFreq|cf=f' => \$corrFreq,
            'freqScale|fs=f' => \$freqScalingFactor,
            'isTS' => \$TS_on)
or die("Error in command line arguments");

if ($printHelp){
    print "$help\n";
    exit 0;
} 

# input file
my $file=$ARGV[0]; #this reads the output file

open(INFILE, "<$file") or die "Could not Read ".$file."\n";
print "Reading ".$file."\n";

while (my $line = <INFILE>){
   if ($line =~ /Frequencies/) { #search for all frequencies 
      my @linesplit = split(/\s+/,$line);
      if (defined $linesplit[3] && $linesplit[3] ne '') {
           push(@frequencies, $linesplit[3]);
      }
      if (defined $linesplit[4] && $linesplit[4] ne '') {
            push(@frequencies, $linesplit[4]);
      }
      if (defined $linesplit[5] && $linesplit[5] ne '') {
           push(@frequencies, $linesplit[5]);
      }
      
   }
   elsif($line =~/Zero-point correction\=\s+([-\w\.]*)/) { # searches for ZPE 
     $ZPEread = $1;
   }
   elsif($line =~ /zero-point Energies\=\s+([-\w\.]*)/) { # searches for ZPE + electronic energy
     $ZPEelectronic = $1;
   }
   elsif ($line =~ /Thermal correction to Enthalpy\=\s+([-\w\.]*)/) {
     $enthalpyread = $1;
   }
   elsif ($line =~ /Thermal correction to Gibbs Free Energy\=\s+([-\w\.]*)/) {
     $gibbsread = $1;
   }
   elsif($line =~ /thermal Free Energies\=\s+([-\w\.]*)/) { # searches for thermal free energy
     $thermalFreeEnergy = $1; 
     $thermalflag = 1;
   }
   elsif($line =~ /thermal Enthalpies\=\s+([-\w\.]*)/) { # searches for thermal enthalpy
      $thermalEnthalpy = $1;
   }
   elsif($line =~ / \<Sx\>\= [\w\.]* \<Sy\>\= [\w\.]* \<Sz\>\= [\w\.]* \<S\*\*2\>\= ([\w\.]*) S\= [\w\.]/ ) {
      $S_2 = $1;
   }
   elsif($line =~ /Charge\s*\=\s*[\-\w]*\s*Multiplicity\s*\=\s*(\w*)/ ){
      $S_2_ideal = ($1-1)/2*(($1-1)/2+1);
   }
   elsif($line =~ /Temperature\s*([\w\.]*)\s*Kelvin/) {
      my $temporaryTemperature = $1; 
      if ($temporaryTemperature == $Tdefault) {next;}
      else {
         $Tdefault = $1;
         print "\nTemperature Default Changed to ".$Tdefault."\n\n";
      }
   } 
}
close(INFILE);
print "Done reading ".$file."\n\n";

if ($thermalflag == 0) { die "Frequencies not found\n";}

$electronicEnergy = $ZPEelectronic - $ZPEread; # obtain electronic energy

if ($T==-1.0) { $T = $Tdefault;} # If no temperature was specified by User, set $T to temperature of freq file (or 298.15 K)

{ no warnings 'uninitialized';
print "Assuming NON-LINEAR MOLECULE for Rotational Temperature Corrections\n\n";

print "Assuming RIGID ROTOR Approximation\n\n";

print "Before Low Frequency Corrections\n";
print "The electronic energy (Hartrees):          ".$electronicEnergy."\n";
print "The ZPE and Electronic E (Hartrees):       ".$ZPEelectronic."\n";
print "The Thermal Enthalpy (Hartrees):           ".$thermalEnthalpy."\n";
print "The Thermal Free Energy (Hartrees):        ".$thermalFreeEnergy."\n";
print "The S**2 value is:                         ".$S_2."\n".
      "with the ideal:                            ".$S_2_ideal."\n";
print "\n";

print "The Zero Point Corrections (Hartrees):     ".$ZPEread."\n";
print "Thermal correction to Enthalpy(Hartrees):  ".$enthalpyread."\n";
print "Thermal correction to Gibbs(Hartrees):     ".$gibbsread."\n";
print "\n";

} # end no warnings initialized

if ($TS_on eq 1) { # if a transition state calculation is not then, remove first frequency from array
    print "Frequency not corrected for TS optimization:\n";
    print "Frequency (cm**-1): ".$frequencies[0]."\n";
    shift(@frequencies);
}
#begin frequency replacements for frequencies less than $corrFreq
print "Frequencies in need of correction: \n";

{ no warnings 'uninitialized';  # ignore unitialized variables, so it doesn't throw warnings in the output
for my $i (0 .. $#frequencies){   # see gaussian white pages for equations of entropy and vibrational computations

    # Loops over real frequencies less than $corrFreq
    if (($frequencies[$i] < $corrFreq) and ($frequencies[$i] > 0.0) ) {
       print "Frequency (cm**-1): ".$frequencies[$i]."\n";

       # compute frequency factors used for both entropy and enthalpy corrections
       # F variables are for changing frequencies
       # T variables are for scaled temperatures
       my $phi_div_T_oldf = $frequencies[$i]*(1/$kb)*(1/$Tdefault);
       my $phi_div_T_newf = $corrFreq*(1/$kb)*(1/$Tdefault);
       my $phi_div_T_oldT = $phi_div_T_newf;
       my $phi_div_T_newT = ($freqScalingFactor)* $corrFreq*(1/$kb)*(1/$T);


       # compute corrections to entropy for low frequencies (most difference should be from this term)
       my $S_v_oldf = $Rconst*( $phi_div_T_oldf / (exp($phi_div_T_oldf) - 1) - log(1 - exp( -$phi_div_T_oldf) ) );
       my $S_v_newf = $Rconst*( $phi_div_T_newf / (exp($phi_div_T_newf) - 1) - log(1 - exp( -$phi_div_T_newf) ) );
       my $S_Fcorr = $S_v_newf - $S_v_oldf;
       $entropyFcorr += $S_Fcorr;

       # compute corrections to entropy for low frequencies at different T (most difference should be from this term)
       my $S_v_oldT = $S_v_newf;
       my $S_v_newT = $Rconst*( $phi_div_T_newT / (exp($phi_div_T_newT) - 1) - log(1 - exp( -$phi_div_T_newT) ) );
       my $S_Tcorr = $S_v_newT - $S_v_oldT;
       $entropyTcorr += $S_Tcorr;


       # compute correction to ZPE for having low frequencies
       my $E_v_oldf = $Rconst * ( $frequencies[$i]*( 1/$kb) * ( 0.5 + 1/( exp($phi_div_T_oldf) - 1) ) );
       my $E_v_newf = $Rconst * ( $corrFreq * (1/$kb) * ( 0.5 + 1/( exp($phi_div_T_newf) - 1) ) );
       my $E_v_Fcorr = $E_v_newf - $E_v_oldf;
       $ZPE_corr += ( $Rconst * ( ($freqScalingFactor)* $corrFreq * (1/$kb) * (0.5) ) - ($Rconst * ( $frequencies[$i] * (1/$kb) * (0.5) ) ) ); 
       $Thermal_vib_Fcorr += $E_v_Fcorr;

       # compute correction to ZPE for having new Temperature
       my $E_v_oldT = $E_v_newf;
       my $E_v_newT = $Rconst * ( ($freqScalingFactor)* $corrFreq * (1/$kb) * ( 0.5 + 1/( exp($phi_div_T_newT) - 1) ) );
       my $E_v_Tcorr = $E_v_newT - $E_v_oldT;
       $Thermal_vib_Tcorr += $E_v_Tcorr;


    }
    # Loops over imaginary Frequencies 
    # Replaces with corrFreq unless isTS is specified
    elsif (($frequencies[$i] < $corrFreq) and ($frequencies[$i] < 0.0)  ){ #imaginary frequency corrections added here
       print "Frequency (cm**-1): ".$frequencies[$i]."  WARNING Imaginary Frequecy Found! Recommended further optimization. WARNING\n";

       # compute values needed for both entropy and enthalpy corrections
       my $phi_div_T_newf = $corrFreq*(1/$kb)*(1/$Tdefault);
       my $phi_div_T_newT = ($freqScalingFactor)*$corrFreq*(1/$kb)*(1/$T);

       # compute correction to entropy for imaginary freq 
       my $S_v_oldf = 0.0;
       my $S_v_newf = $Rconst*( $phi_div_T_newf / (exp($phi_div_T_newf) - 1) - log(1 - exp( -$phi_div_T_newf) ) );
       my $S_Fcorr = $S_v_newf - $S_v_oldf;
       $entropyFcorr += $S_Fcorr;

       # compute corrections to entropy for low frequencies at different T (most difference should be from this term)
       my $S_v_newT = $Rconst*( $phi_div_T_newT / (exp($phi_div_T_newT) - 1) - log(1 - exp( -$phi_div_T_newT) ) );
       my $S_Tcorr = $S_v_newT - $S_v_newf;
       $entropyTcorr += $S_Tcorr;


       # compute correction to ZPE for having low frequencies
       my $E_v_oldf = 0.0;
       my $E_v_newf = $Rconst * ( $corrFreq * (1/$kb) * ( 0.5 + 1/( exp($phi_div_T_newf) - 1) ) );
       my $E_v_Fcorr = $E_v_newf - $E_v_oldf;
       $ZPE_corr += ( $Rconst * ( ($freqScalingFactor)*$corrFreq * (1/$kb) * (0.5) ) - ($Rconst * ( $frequencies[$i] * (1/$kb) * (0.5) ) ) );
       $Thermal_vib_Fcorr += $E_v_Fcorr;

       # compute correction to ZPE for having new Temperature
       my $E_v_oldT = $E_v_newf;
       my $E_v_newT = $Rconst * ( ($freqScalingFactor)*$corrFreq * (1/$kb) * ( 0.5 + 1/( exp($phi_div_T_newT) - 1) ) );
       my $E_v_Tcorr = $E_v_newT - $E_v_oldT;
       $Thermal_vib_Tcorr += $E_v_Tcorr;

       $imagflag = 1;

    }
    elsif (($frequencies[$i] >= $corrFreq ) and ( ($T != $Tdefault ) or ($freqScalingFactor != 1.0)) ){ 
       # compute needed for both entropy and enthalpy corrections
       my $phi_div_T_oldT = $frequencies[$i]*(1/$kb)*(1/$Tdefault);
       my $phi_div_T_newT = ($freqScalingFactor)*$frequencies[$i]*(1/$kb)*(1/$T);

       # compute corrections to entropy for low frequencies at different T (most difference should be from this term)
       my $S_v_oldT = $Rconst*( $phi_div_T_oldT / (exp($phi_div_T_oldT) - 1) - log(1 - exp( -$phi_div_T_oldT) ) );
       my $S_v_newT = $Rconst*( $phi_div_T_newT / (exp($phi_div_T_newT) - 1) - log(1 - exp( -$phi_div_T_newT) ) );
       my $S_Tcorr = $S_v_newT - $S_v_oldT;
       $entropyTcorr += $S_Tcorr;
       # compute correction to ZPE for having new Temperature
       my $E_v_oldT = $Rconst * ( $frequencies[$i]*( 1/$kb) * ( 0.5 + 1/( exp($phi_div_T_oldT) - 1) ) );
       my $E_v_newT = $Rconst * ( ($freqScalingFactor)*$frequencies[$i] * (1/$kb) * ( 0.5 + 1/( exp($phi_div_T_newT) - 1) ) );
       my $E_v_Tcorr = $E_v_newT - $E_v_oldT;
       $ZPE_corr += ( $Rconst * ( ($freqScalingFactor)*$frequencies[$i] * (1/$kb) * (0.5) ) - ($Rconst * ( $frequencies[$i] * (1/$kb) * (0.5) ) ) );
       $Thermal_vib_Tcorr += $E_v_Tcorr;

    }

}
####################################################################
# FINAL THERMAL COMPUTATIONS
# compute the enthalpy, entropy, and free energies final corrections
####################################################################

# Compute corrections from change in T flag
if ( $T != $Tdefault) {

   # if new temperature, we have to correct enthalpies and entropies of translations 
   my $translation_entropy_Tcorr = ($Rconst * 5 / 2) * ( log($T)  -  log($Tdefault) ); # S_trans = R/2 * (5ln T + 3lnM) - const
   my $translation_enthalpy_Tcorr = 3/2 * $Rconst * ($T - $Tdefault);  # H_trans = 3/2 R T + pV
   # add corrections to global variable
   $entropyTcorr += $translation_entropy_Tcorr;
   $enthalpyTcorr += $translation_enthalpy_Tcorr;
   # if new temperature, we have to correct entropies of translations
   my $rotation_entropy_Tcorr = ($Rconst * 3 / 2) * ( log($T) - log($Tdefault)); # only part that scales with temperature is 3/2RT
   my $rotation_enthalpy_Tcorr = 3/2 * $Rconst * ($T - $Tdefault);  # H_rot = 3/2 R T 
   # add corrections to global variable
   $entropyTcorr += $rotation_entropy_Tcorr;
   $enthalpyTcorr += $rotation_enthalpy_Tcorr;
   
   # Add RdT corr for term (Internal Energy -> Enthalpy )
   $enthalpyTcorr += $Rconst * ($T - $Tdefault);
}



#################################################
# ENTHALPY FINAL VALUES
# all thermal corrections (kcal/mol)
$enthalpy_all_corr = (($Thermal_vib_Fcorr + $Thermal_vib_Tcorr + $enthalpyTcorr) );
# All enthalpy corrections from changing T, kcal/mol
$enthalpy_T_corr = ($Thermal_vib_Tcorr + $enthalpyTcorr);
# Final Enthalpy total (a.u.)
$enthalpy_total = $thermalEnthalpy + ($enthalpy_all_corr / $hartreetokcal);
#recompute the total internal thermal "correction"
my $corrThermalCont = $enthalpy_total - $electronicEnergy;
####################################################
# ENTROPY
# Compute old Entropy in case of T flag requested
$oldentropy = $hartreetokcal * (-1 * ($thermalFreeEnergy - $thermalEnthalpy) / $Tdefault);
###################################################
# FREE ENERGY
# final gibbs correction (kcal/mol)
$gibbscorr = ( ( $enthalpy_all_corr ) - ( ($entropyFcorr + $entropyTcorr)*($T) + $oldentropy*($T - $Tdefault) ));
# recompute the total Gibbs Free Energy
$gibbs_total = $thermalFreeEnergy + ($gibbscorr / $hartreetokcal );
# Final Gibbs correction (a.u.)
my $corrGibbsCont = $gibbs_total - $electronicEnergy;
#####################################################
# ZPE
# recompute ZPE with low freq corrections
my $ZPE = $ZPEread + ($ZPE_corr / $hartreetokcal);
$ZPE = sprintf("%.6f",$ZPE);
################################################# 
# Format numbers for output (fixed 6 and 5 decimal places)
$ZPE_corr = sprintf("%.5f",$ZPE_corr);
$Thermal_vib_Fcorr = sprintf("%.5f",$Thermal_vib_Fcorr);
$Thermal_vib_Tcorr = sprintf("%.5f",$Thermal_vib_Tcorr);

$entropyFcorr = sprintf("%.5f",$entropyFcorr);
$entropyTcorr = sprintf("%.5f",$entropyTcorr);

$enthalpyTcorr = sprintf("%.5f",$enthalpyTcorr);
$enthalpy_all_corr = sprintf("%.5f",$enthalpy_all_corr);
$enthalpy_T_corr = sprintf("%.5f",$enthalpy_T_corr);
$enthalpy_total = sprintf("%.6f",$enthalpy_total);
$corrThermalCont = sprintf("%.6f",$corrThermalCont);

$gibbscorr = sprintf("%.5f",$gibbscorr);
$gibbs_total = sprintf("%.6f",$gibbs_total);
$corrGibbsCont = sprintf("%.6f",$corrGibbsCont);

######################################################################
# END COMPUTATIONS
######################################################################

print "All Frequencies Analyzed\n\n";
print "THERMAL CORRECTION AT ".$T." K\n\n";
print "The magnitude of changes from low frequencies is the following:\n";
print "The low Freq ZPE Correction (kcal/mol):                  ".$ZPE_corr."\n";
print "The Low Freq E_v Correction (kcal/mol):                  ".$Thermal_vib_Fcorr."\n";
print "The Low Freq Entropic Correction (kcal/mol*K):           ".$entropyFcorr."\n";
print "\n";

print "Corrections to thermodynamics to reach ".$T." K from ".$Tdefault." K are:\n";
print "The new Temperature E_v Correction (kcal/mol):           ".$Thermal_vib_Tcorr."\n";
print "The new Temperature H Correction (kcal/mol):             ".$enthalpy_T_corr."\n";
print "The new Temperature Entropic Correction (kcal/mol*K):    ".$entropyTcorr."\n";
print "\n";

print "Total Requested Thermal Corrections\n";
print "The Total Enthalpy Correction (kcal/mol):              ".$enthalpy_all_corr."\n";
print "The Total Free Energy Correction (kcal/mol):          ".$gibbscorr."\n";
print "\n";

print "The Thermally Corrected Results:\n";
print "The electronic energy (Hartrees):                            ".$electronicEnergy."\n";
print "The Corrected Thermal Enthalpy (Hartrees):                   ".$enthalpy_total."\n";
print "The Corrected Thermal Free Energy (Hartrees):                ".$gibbs_total."\n\n";

print "The Corrected Zero Point Vibrational Energy (Hartrees):      ".$ZPE."\n";
print "The Corrected Thermal Enthalpy Contribution (Hartrees):      ".$corrThermalCont."\n";
print "The Corrected Thermal Free Energy Contribution (Hartrees):   ".$corrGibbsCont."\n";

if ($imagflag == 1) { print "WARNING IMAGINARY FREQUENCIES WERE FOUND!!  \n "; }

print "\n \n"
} # end of section where it ignores undeclared variables

