#!/usr/bin/perl

############################################################################################
#                                                                                          #
# VERSION :: 2.0.1                                                                         #
#                                                                                          #
# AUTHORS :: William Christian Isley III                                                   #
#                                                                                          #
#                                                                                          #
# SOFTWARE :: This script was created using Perl v 5.20                                    #
# This script has dependency on the pacakges, Getopt::Long,                                #
# This script is available via github.com/william-isley-3rd/Comp-Chem-Tools                #
#                                                                                          #
# INPUT :: This script has been tested to function with regularly formated ADF, ORCA, and  #
# Gaussian 09 Output files. Gaussian files require "Standard Coordinates" to function.     #
#                                                                                          #       
# IMPORTANT :: This the ADF lowest energy indexer and specific step indexer are exclusive. #
#                                                                                          #
# DISCLAIMER :: This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by the Free Software   #
# Foundation, either version 3 of the License, or any later version                        #
# This program is distributed in the hope that it will be useful, but WIHOUT ANY WARRANTY  #
# without even the implied warranty of MERCANTABILITY or FITNESS FOR A PARTICULAR PURPOSE  #
# See the GNU General Public License for more details.                                     #
#                                                                                          #
############################################################################################


use strict;
use warnings;
use Getopt::Long;

# global flags for input
my $targetGeoStep = -1;
my $g09File = 0;
my $adfFile = 0;
my $orcaFile= 0;
my $adfLowEnergyFile = 0;
my $printHelp = 0;

# get the options for input
GetOptions ('help|h' => \$printHelp,
            'targetStep|step=i' => \$targetGeoStep,
            'G09file|g09=s' => \$g09File,
            'Orcafile|orca=s' => \$orcaFile,
            'ADFfile|adf=s' => \$adfFile,
            'ADFLowEnergy|adfle=s' => \$adfLowEnergyFile)
or die("Error in command line arguments, use --help for assistance");

sub Print_Help_Sub {
     my $helptext = "\nThis script parses Gaussian 09 or ADF outfiles for molecular coordinates \n\n".
                    "Options available for the script :\n".
                    "--step=X ; The specify coordinates from geometric step X taken by the optimizer\n".
                    "         ; The default is the last step taken\n".
                    "--orca=X ; The ORCA file X containing molecular coordinates, prints last coords\n".
                    "--g09=X ; The G09 file X containing molecular coordinates, prints last coords\n".
                    "--adf=X ; The ADF file X containing molecular coordinates, prints last coords\n".
                    "--adfle=X ; The ADF file X containing molecular coordinates, prints lowest energy coords\n\n";
     my ($printHelpFlag) = @_;
     #print out help option is flag is specified
     if ( ($printHelpFlag) ){
         print $helptext."\n";
         exit 0;
     }
}

# global variables

# Hash containing atomic symbols.
# Note is is currently defined for elements 1-56,72-79,82
# if you need another elements please add to this hash using the same syntax
my %atomicNumberToSymbol = (
   1,"H",
   2,"He",
   3,"Li",
   4,"Be",
   5,"B",
   6,"C",
   7,"N",
   8,"O",
   9,"F",
  10,"Ne",
  11,"Na",
  12,"Mg",
  13,"Al",
  14,"Si",
  15,"P",
  16,"S",
  17,"Cl",
  18,"Ar",
  19,"K",
  20,"Ca",
  21,"Sc",
  22,"Ti",
  23,"V",
  24,"Cr",
  25,"Mn",
  26,"Fe",
  27,"Co",
  28,"Ni",
  29,"Cu",
  30,"Zn",
  31,"Ga",
  32,"Ge",
  33,"As",
  34,"Se",
  35,"Br",
  36,"Kr",
  37,"Rb",
  38,"Sr",
  39,"Y",
  40,"Zr",
  41,"Nb",
  42,"Mo",
  43,"Tc",
  44,"Ru",
  45,"Rh",
  46,"Pd",
  47,"Ag",
  48,"Cd",
  49,"In",
  50,"Sn",
  51,"Sb",
  52,"Te",
  53,"I",
  54,"Xe",
  55,"Cs",
  56,"Ba",
  72,"Hf",
  73,"Ta",
  74,"W",
  75,"Re",
  76,"Os",
  77,"Ir",
  78,"Pt",
  79,"Au",
  82,"Pb"
);

# SUBROUTINES 1
# read_G09_Outfile_Sub
# function: find last coordinates in Gaussian09 formated output file
sub read_G09_Outfile_Sub{

    my $standardflag=0;
    my $beginreadingflag=0;
    my @localOutputCoords;
    my $geometryStepper=0;
    open(INFILE, "<$g09File") or die;
    
    while (my $line = <INFILE>){
       if ($line =~ /Standard orientation/) { #search for Standard orientation line
           $standardflag=1;
           $geometryStepper++;
           undef(@localOutputCoords)
       }
       elsif ($line =~/-{20,69}/) {
           if (($standardflag eq 1) and ($beginreadingflag eq 0)) { #if line contains many ---, increase counter
              $beginreadingflag=1;
           }
           elsif (($standardflag eq 1) and ($beginreadingflag eq 1)) { # increase counter once more, prepare to being reading
              $beginreadingflag=2;
           }
           elsif (($standardflag eq 1) and ($beginreadingflag eq 2)) { # reset counter, coordinate block is over
              $beginreadingflag=0;
              $standardflag=0;
              if ($geometryStepper == $targetGeoStep) {
                  last;
              }
           }
        }
        elsif ( ( $standardflag eq 1) and ($beginreadingflag eq 2 )) { #begin reading coordinates in
              my @linesplit = split(/\s+/,$line);
              # print "@linesplit"."\n";;
              my $atomicNumber= $linesplit[2];
              my $atomicSymbol = $atomicNumberToSymbol{$atomicNumber};
              # print $atomicSymbol."\n";
              my $Xcoord= $linesplit[4];
              my $Ycoord= $linesplit[5];
              my $Zcoord= $linesplit[6];
              push(@localOutputCoords, [$atomicSymbol,$Xcoord,$Ycoord,$Zcoord]);
       }
       
    }
    close(INFILE);
    return @localOutputCoords;
}
# SUBROUTINE 2
# read_ORCA_Outfile_Sub
# function: find last coordinates in ORCA formated outputfile
sub read_ORCA_Outfile_Sub{

    my $standardflag=0;
    my $beginreadingflag=0;
    my @localOutputCoords;
    my $geometryStepper=0;
    open(INFILE, "<$orcaFile") or die;

    while (my $line = <INFILE>){
       if ($line =~ /CARTESIAN COORDINATES \(ANGSTROEM\)/) { #search for Cartesian Coordinates line
           $standardflag=1;
           $geometryStepper++;
           undef(@localOutputCoords)
       }
       elsif ($line =~/-{20,69}/) {
           if (($standardflag eq 1) and ($beginreadingflag eq 0)) { #if line contains many ---, increase counter
              $beginreadingflag=1;
           }
           elsif (($standardflag eq 1) and ($beginreadingflag eq 1)) { # reset counter, coordinate block is over
              $beginreadingflag=0;
              $standardflag=0;
              if ($geometryStepper == $targetGeoStep) {
                  last;
              }
           }
        }
        elsif ( ( $standardflag eq 1) and ($beginreadingflag eq 1 )) { #begin reading coordinates in
              if ($line =~ m/^\s*$/g) { next;}
              my @linesplit = split(/\s+/,$line);
              # print "@linesplit"."\n";;
              my $atomicSymbol = $linesplit[1];
              # print $atomicSymbol."\n";
              my $Xcoord= $linesplit[2];
              my $Ycoord= $linesplit[3];
              my $Zcoord= $linesplit[4];
              push(@localOutputCoords, [$atomicSymbol,$Xcoord,$Ycoord,$Zcoord]);
       }

    }
    close(INFILE);
    return @localOutputCoords;
}

# SUBROUTINE 3
# read_ADF_Outfile_sub
# function: read ADF formated output files, looking for last coordinate during optimziation
sub read_ADF_Outfile_Sub{

     my $standardflag=0;
     my $beginreadingflag=0;
     my @localOutputCoords;
     my $geometryStepper=0;

     open(INFILE, "<$adfFile") or die;
     while (my $line = <INFILE>){
        if ($line =~ /Coordinates \(Cartesian\)/) { #search for line beginning coordinate list
            $standardflag=1;
            $geometryStepper++;
            undef(@localOutputCoords)
        }
        elsif ($line =~/-{20,69}/) {
            if (($standardflag eq 1) and ($beginreadingflag eq 0)) { #if line contains many ---, increase counter
               $beginreadingflag=1;
            }
            elsif (($standardflag eq 1) and ($beginreadingflag eq 1)) { # reset counter, coordinate block is over
               $beginreadingflag=0;
               $standardflag=0;
               if ($geometryStepper == $targetGeoStep) {
                   last;
               }
            }
         }
         elsif ( ( $standardflag eq 1) and ($beginreadingflag eq 1 )) { #begin reading coordinates in
               my @linesplit = split(/\s+/,$line);
               # print "@linesplit"."\n";;
               my $atomicSymbol= $linesplit[2];
               # print $atomicSymbol."\n";
               my $Xcoord= $linesplit[6];
               my $Ycoord= $linesplit[7];
               my $Zcoord= $linesplit[8];
               push(@localOutputCoords, [$atomicSymbol,$Xcoord,$Ycoord,$Zcoord]);
        }
     }
     close(INFILE);
     return @localOutputCoords;

}

# SUBROUTINE 4 
# read_ADF_Outfile_Low_Energy_Sub
# function: read ADF formated output optimization jobs looking for lowest energy configuration
sub read_ADF_Outfile_Low_Energy_Sub {
    my $energyflag=0;
    my $standardflag=0;
    my $beginreadingflag=0;
    my $minimumEnergy=0;
    my @localOutputCoords;

    open(INFILE, "<$adfLowEnergyFile") or die;

     while (my $line = <INFILE>){
        if ($line =~ m/current energy\s*([\-\d\.]*)\s*Hartree/g ) { #look for minimal energy step/ search for line beginning coord list
             my $energy= $1;
             if ($energy < $minimumEnergy) {
                  $minimumEnergy= $energy;
                  $energyflag=1;
             }
        }
        elsif ($line =~ /Coordinates \(Cartesian\)/) { #search for line beginning coordinate list
            if ($energyflag eq 1) {
              $standardflag=1;
              undef(@localOutputCoords)
            }
        }
        elsif ($line =~/-{20,69}/) {
            if (($standardflag eq 1) and ($beginreadingflag eq 0)) { #if line contains many ---, increase counter
               $beginreadingflag=1;
            }
            elsif (($standardflag eq 1) and ($beginreadingflag eq 1)) { # reset counter, coordinate block is over
               $beginreadingflag=0;
               $standardflag=0;
               $energyflag=0;
            }
         }
         elsif ( ( $standardflag eq 1) and ($beginreadingflag eq 1 )) { #begin reading coordinates in
               my @linesplit = split(/\s+/,$line);
               # print "@linesplit"."\n";;
               my $atomicSymbol= $linesplit[2];
               # print $atomicSymbol."\n";
               my $Xcoord= $linesplit[6];
               my $Ycoord= $linesplit[7];
               my $Zcoord= $linesplit[8];
               push(@localOutputCoords, [$atomicSymbol,$Xcoord,$Ycoord,$Zcoord]);
        }
     }
     close(INFILE);
     return @localOutputCoords;
}

# SUBROUTINE 5
# print_XYZ_Sub
# Print the XYZ files to the passed output file
# Passed parameters are (string printfile, 2D array with atom &  XYZ coordinates)
sub print_XYZ_Sub{
    my ($printFile, @outputCoords) = @_;
    #begin printing the output coordinates
    my $atom_total = $#outputCoords + 1;
    if ( @outputCoords) {
        print $atom_total."\n";
        print $printFile."\n";
        for my $i (0 .. $#outputCoords){
            print $outputCoords[$i][0]."    ".$outputCoords[$i][1]."    ".$outputCoords[$i][2]."    ".$outputCoords[$i][3]."\n";
        }
    } 
    else { print "No coordinates found in ".$printFile."\n";}
}

# Main EXECUTABLE CALL

# Check if help is required
Print_Help_Sub($printHelp);
# If check on which flag is called
# G09 outfile
if ( ( $g09File ) and ( !$adfFile) and (!$adfLowEnergyFile) and ( !$orcaFile) ) {
    my @g09OutputCoords = read_G09_Outfile_Sub();
    print_XYZ_Sub($g09File, @g09OutputCoords);
}
# ADF, last geom
elsif ( ( !$g09File ) and ( $adfFile) and (!$adfLowEnergyFile) and (!$orcaFile) ) {
    my @adfOutputCoords = read_ADF_Outfile_Sub();
    print_XYZ_Sub($adfFile, @adfOutputCoords);
}
# ADF, low E
elsif ( ( !$g09File) and ( !$adfFile) and ($adfLowEnergyFile) and ( !$orcaFile) and ($targetGeoStep != -1 )) {
    my @adfleOutputCoords = read_ADF_Outfile_Low_Energy_Sub();
    print_XYZ_Sub($adfLowEnergyFile, @adfleOutputCoords);
}
# ORCA, last geom
elsif ( ( !$g09File ) and ( !$adfFile) and (!$adfLowEnergyFile) and ( $orcaFile) ) {
    my @orcaOutputCoords = read_ORCA_Outfile_Sub();
    print_XYZ_Sub($orcaFile, @orcaOutputCoords);
}
# print HELP if no valid flags found
else {
    print "Error in Input\n";
    Print_Help_Sub('true');
}


