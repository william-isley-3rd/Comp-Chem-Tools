#!/usr/bin/perl

#############################################################################################
##                                                                                          #
## VERSION :: 2.0.1                                                                         #
##                                                                                          #
## AUTHORS :: William Christian Isley III                                                   #
##                                                                                          #
##                                                                                          #
## SOFTWARE :: This script was created using Perl v 5.20                                    #
## This script has dependency on the pacakges, Getopt::Long,                                #
## This script is available via github.com/william-isley-3rd/Comp-Chem-Tools                #
##                                                                                          #
## INPUT :: This script has been tested to function with regularly formated Gaussian 09     #
## And Protein Databank (pdb) files                                                         #
##                                                                                          #       
##                                                                                          #
## DISCLAIMER :: This program is free software: you can redistribute it and/or modify       #
## it under the terms of the GNU General Public License as published by the Free Software   #
## Foundation, either version 3 of the License, or any later version                        #
## This program is distributed in the hope that it will be useful, but WIHOUT ANY WARRANTY  #
## without even the implied warranty of MERCANTABILITY or FITNESS FOR A PARTICULAR PURPOSE  #
## See the GNU General Public License for more details.                                     #
##                                                                                          #
#############################################################################################

use strict;
use warnings FATAL => 'uninitialized';
use Getopt::Long;
use Carp ();
$SIG{__DIE__} = \&Carp::confess;

sub trim { 
       my $s = shift; 
       $s =~ s/^\s+|\s+$//g; 
       return $s; 
}

my %resCharge = (
  'ALA' => '0',
  'VAL' => '0',
  'LEU' => '0',
  'ILE' => '0',
  'PRO' => '0',
  'TRP' => '0',
  'PHE' => '0',
  'MET' => '0',
  'GLY' => '0',
  'SER' => '0',
  'THR' => '0',
  'TYR' => '0',
  'CYS' => '0',
  'ASN' => '0',
  'GLN' => '0',
  'LYS' => '1',
  'ARG' => '1',
  'HIS' => '1',
  'ASP' => '-1',
  'GLU' => '-1',
  'WAT' => '0',
  'HOH' => '0',
);

################
# Help Prompt  #
################

my $help = "\nThis script generates protein clusters from pdb formatted gaussian input files

Options available for the script (include these before your the filename):
--cr=X         ; Change the cutoff radius from atoms on the residue of interest to be X Ang.
               ; 5 Ang is the default cutoff radius
--tr=X         ; Select the target residue to X. 
               ; This is required input. Default is residue 0.
--nh           ; no hydrogens on target residue (or test residues) used 
               ; for cluster generation determination

Examples:
Change cutoff radius to 15 Ang     ; perl protein_cluster_gen.pl --tr=97 --cr=15.0 filename.gjf

Any number of options can be combined as shown in the 2nd example";

############################
# Input Options, Variables #
############################

my $file=$ARGV[0]; #this reads the output file
my $printHelp=0; #flag to print help
my $cutoffRadius = 5.0; # default cutoff radius is 5 angstroms
my $targetResidue=0;       # the target residue
my $nohydrogenflag=0;    # flag to turn off protons determining kept residues

# the GetOpts (options) for input 
GetOptions ('help|h' => \$printHelp,
            'cutoffradius|cr=f' => \$cutoffRadius,
            'targetresidue|tr=i' => \$targetResidue,
            'nohydrogen|nh' => \$nohydrogenflag)
or die("Error in command line arguments");

# print out help option is flag is specified
if ( ($printHelp) or ($targetResidue == 0) ){
    print "$help\n";
    exit 0;
}

# global data parsed from PDB/GJF File
my @atomID;
my @PDBName;
my @ResName;
my @ResID;
my @Xcoord;
my @Ycoord;
my @Zcoord;
my @outputLine;

# Cluster Seeding Residue and Those Fitting Keep Criteria
my @targetRes;
my @keptRes;
my %residueTypes = ();

###################
# Read Input File #
###################

# Determine the format of the input file
# open up the inputfile and store data in temporary arrays

if ((index($file, ".gjf") != -1) or (index($file, ".com") != -1) or (index($file, ".inp") != -1))
{
    print "\n\nGaussian File Type found\n".
          "WARNING: If this file was not generated from a PDB with Residue information, ".
          "this script will not function\n\n";
    # File type is gaussian with ending gjf, inp or com
    # open up the inputfile and store data in temporary arrays
    open(INFILE, "<$file") or die;
    while (my $line = <INFILE>){
         # for a pdb formatted gaussian file, parse as below
         # sample input should match the regex below
         #  N(PDBName=N,ResName=SER,ResNum=42)                0   35.98200000    2.59600000    6.32900000 L
         #
         if ($line =~ m/ (\w+)\(PDBName\=(\w+)\,ResName\=(\w+)\,ResNum\=(\w+)\)\s+0\s+([\w\.\-]+)\s+([\w\.\-]+)\s+([\w\.\-]+)/ ) {
                push(@atomID,  $1);
                push(@PDBName, $2);
                push(@ResName, $3);
                push(@ResID,   $4);
                push(@Xcoord,  $5);
                push(@Ycoord,  $6);
                push(@Zcoord,  $7);
                push(@outputLine,$line);
                unless( $residueTypes{$4} ) { $residueTypes{$4} = $3;} 
         }
    }
    close(INFILE);
}
elsif (index($file, "pdb") != -1) {
    # File type is protein database with ending pdb
    # open up the inputfile and store data in temporary arrays
    print "\n\nProtein Databank File Type Found\n".
          "WARNING: If water solvents do not have WAT or HOH as the Residue Name, ".
          "this file will not function properly\n\n";
    open(INFILE, "<$file") or die;
    while (my $line = <INFILE>){
         if ($line =~ /^ATOM/ ) {
               my $subResID   = trim(substr($line,  22,  4));
               my $subResName = trim(substr($line,  17,  3)); 

               push(@atomID,   trim(substr($line,  76,  2)) ); # If error points here, you are missing element column
               push(@PDBName,  trim(substr($line,  12,  4)) );
               push(@ResName,  trim(substr($line,  17,  3)) );
               push(@ResID,    trim(substr($line,  22,  4)) );
               push(@Xcoord,   trim(substr($line,  30,  8)) );
               push(@Ycoord,   trim(substr($line,  38,  8)) );
               push(@Zcoord,   trim(substr($line,  46,  8)) );
               push(@outputLine, $line);

               unless ( $residueTypes{$subResID} ) { $residueTypes{$subResID} = $subResName;}
         }
         # get the water molecules as well, should match sample below
         # HETATM   45  O   HOH A   9      28.690  18.913   5.002  1.00  9.37           O
         elsif (($line =~/^HETATM/) and ( substr($line, 17, 3) =~ m/(HOH|WAT)/ ) ) {
               my $subResID   = trim(substr($line,  22,  4));
               my $subResName = trim(substr($line,  17,  3));

               push(@atomID,   trim(substr($line,  76,  2)) );
               push(@PDBName,  trim(substr($line,  12,  4)) );
               push(@ResName,  trim(substr($line,  17,  3)) );
               push(@ResID,    trim(substr($line,  22,  4)) );
               push(@Xcoord,   trim(substr($line,  30,  8)) );
               push(@Ycoord,   trim(substr($line,  38,  8)) );
               push(@Zcoord,   trim(substr($line,  46,  8)) );
               push(@outputLine, $line);

               unless ( $residueTypes{$subResID} ) { $residueTypes{$subResID} = $subResName; }
         }
    }
    close(INFILE);
}

print "\n".
      "FILE FOUND : $file \n".
      "READING    : SUCCESS\n".
      "Beginning Cluster Seeding.\n\n";

#########################################################
# Find All Residues Near Target Residue to Seed Cluster #
#########################################################

# find all atoms in the protein matching the target residue for seeding, put them in @targetRes
foreach my $resTemp (0..$#ResID) {
    if ($ResID[$resTemp] eq $targetResidue ) {
        push(@targetRes,$resTemp);
    }
}

print "TARGET RESIDUE NUMBER : $targetResidue \n".
      "TARGET RESIDUE TYPE   : $residueTypes{$targetResidue}. \n".
      "CUTOFF RADIUS (ANG)   : $cutoffRadius \n".
      "\nNow searching for residues and solvent within $cutoffRadius Angstroms of residue $targetResidue \n\n";

# find all residues in the protein within the cutoff radius of the target Residue
foreach my $probe (0..$#targetRes) {  # loop over target residue's atoms
     if ( ( $nohydrogenflag ) and ($atomID[$probe] eq "H") ) { 
        next;
        # if the atom is a hydrogen, and no hydrogen flag is on, skip it
     }
     else {
      
         foreach my $pairID (0..$#ResID) { # loop over all atoms
             if ( ( $nohydrogenflag ) and ($atomID[$pairID] eq "H") ) {
                 next;
             }
             else {

                 my $Xdist = $Xcoord[$targetRes[$probe]] - $Xcoord[$pairID]; # compute X distance
                 my $Ydist = $Ycoord[$targetRes[$probe]] - $Ycoord[$pairID]; # compute Y distance
                 my $Zdist = $Zcoord[$targetRes[$probe]] - $Zcoord[$pairID]; # compute Z distance
        
                 my $dist = sqrt( $Xdist**2 + $Ydist**2 + $Zdist**2); # compute radial distance
                 # if the atom is within the cutoff radius, add the residue to the keep list
                 if  ( $dist <= $cutoffRadius)  { 
                    push(@keptRes,$ResID[$pairID]);          # if less or equal, add Residue to kept residues
                 }  
            }
         }
     }
}

# Remove duplicates from @keptRes to obtain
# the unique residues that are within the cutoff radius
my %hashKeyRes = map { $_ => "kept"} @keptRes;         # map keptRes to a hash, removing replicates
my @uniqueKeptRes = keys %hashKeyRes;                  # put unique ResID's to an array
my @sortedKeptRes = sort { $a <=> $b } @uniqueKeptRes; # sort in numerical order the Kept Residues

##############################################
# Parse the Capping Instructions for Cluster #
##############################################


my @aminecaps;
my @carboxylcaps;
my @linkingResidues;
my @solvent;
# if the residue is one less, and residue isn't already kept, keep the C, O, and alpha_C for capping
# # if the residue is one more, and residue isn't already kept, kep the N, and alpha_C for capping
for my $i (0..$#sortedKeptRes) { # look for residues that need capped
    if ( $residueTypes{$sortedKeptRes[$i]} =~ m/(HOH|WAT)/ ) {
         push(@solvent, $sortedKeptRes[$i]);
    }
    elsif ( $residueTypes{$sortedKeptRes[$i]} !~ m/(HOH|WAT)/ ) {
         push(@linkingResidues, $sortedKeptRes[$i]);
         # added kept residue to list, now need to add the capping residues
         # check if residue below is in the kept list (if not, add carboxyl cap)
         # check if residue above is in the kept list (if not, add amine cap)
         if ( ( $i == 0) and 
             ( $residueTypes{$sortedKeptRes[$i]-1} !~ m/(HOH|WAT)/ ) ) { 
             push(@carboxylcaps, $sortedKeptRes[$i]-1);
         }
         elsif ( ( $sortedKeptRes[$i-1] != $sortedKeptRes[$i]-1) 
               and ( $residueTypes{$sortedKeptRes[$i]-1} !~ m/(HOH|WAT)/ ) ){
             push(@carboxylcaps, $sortedKeptRes[$i]-1);
         }
         if ( ( $i == $#sortedKeptRes) 
            and ( $residueTypes{$sortedKeptRes[$i]+1} !~ m/(WAT|HOH)/ ) ){ 
             push(@aminecaps, $sortedKeptRes[$i]+1);
         }
         elsif ( ( $sortedKeptRes[$i+1] != $sortedKeptRes[$i]+1) 
               and ( $residueTypes{$sortedKeptRes[$i]+1} !~ m/(WAT|HOH)/ ) ){ 
             push(@aminecaps, $sortedKeptRes[$i]+1);
        }
    }
}

# Instructions for which residues are linking vs capping
# Within capping residues, determine which ends are capped
my %linkerID = map { $_ => "linking"} @linkingResidues;
my %cappingInstructions = (%hashKeyRes, %linkerID);
foreach my $carboxylID (@carboxylcaps) {
     unless ( $cappingInstructions{$carboxylID}) { $cappingInstructions{$carboxylID} = "carboxyl"; }
}
foreach my $amineID (@aminecaps) {
     unless ($cappingInstructions{$amineID}) { $cappingInstructions{$amineID} = "amine";}
     if ($cappingInstructions{$amineID} =~ m/carboxyl/) { $cappingInstructions{$amineID} = "linking";} 
}
foreach my $all_Solvent (@solvent) {
     $cappingInstructions{$all_Solvent} = "solvent";
}


print "Done Cluster Seeding\n".
      "Statistics below for the seeded cluster\n".
      " ResNumber   ResType    Charge     Instruction  \n";
my @printInstructions;
my $totalCharge=0;
foreach my $resInstruction (sort keys %cappingInstructions) {
     my $printResNumber = $resInstruction;
     my $printResType = $residueTypes{$printResNumber};
     my $printResCharge=0;
     if ($resCharge{$printResType}) {$printResCharge = $resCharge{$printResType};} 

     my $tempPrint = sprintf('%10s', $printResNumber).sprintf('%10s', $printResType).sprintf('%10s', $printResCharge).sprintf('%15s', $cappingInstructions{$resInstruction});
     push(@printInstructions, $tempPrint);

     # adjust total charge to be printed based on presence of linking ligands
     if ( $cappingInstructions{$resInstruction} =~ m/linking/) {
         $totalCharge+=$printResCharge;
     }

}
my @sortedPrintInstructions = sort @printInstructions;
foreach my $printLine (@sortedPrintInstructions) {
    print "$printLine \n";
}

print "\nTOTAL CHARGE : $totalCharge  (Note, capping ligands do not count toward charge)\n".
      "WARNING: Ensure that accurate protonation state for capping ligands before performing further computations!\n\n".
      "Printing cluster coordinates to tr-$targetResidue-rad-$cutoffRadius.$file \n\n".

###########################################################
# Cluster Generation is Complete, Print Cluster to STDOUT #
###########################################################

open(OUTFILE, ">tr-$targetResidue-rad-$cutoffRadius.$file") or die;

# print header for gjf file and print the saved coordinates
print OUTFILE "# opt integral=ultrafine m06l/def2svp/auto

Title

$totalCharge 1\n";
# print out the kept atoms from each residue
# print out the linking, and capping residues
# protonation should be checked, especially for capped residues
foreach my $linkingID (sort keys %cappingInstructions) {
    for my $allResID (0..$#ResID) {
        # if a linking residue, print CA, N, C, O atoms  
        if ( ($linkingID eq $ResID[$allResID] ) and ($cappingInstructions{$linkingID} eq "linking")  ) {
             print OUTFILE $outputLine[$allResID];
        }
        # if the residue is an amine cap, take only CA and N atoms
        elsif  ( ($linkingID eq $ResID[$allResID] ) and ($cappingInstructions{$linkingID} eq "amine") and ( $PDBName[$allResID] =~ m/^(CA|N)$/ ) ) {
             print OUTFILE $outputLine[$allResID];
        }
        # if the residue is a carboxyl cap, take only CA, C and O atoms
        elsif( ($linkingID eq $ResID[$allResID] ) and ($cappingInstructions{$linkingID} eq "carboxyl")  and ( $PDBName[$allResID] =~ m/^(CA|C|O)$/) ) {
             print OUTFILE $outputLine[$allResID];
        }
        elsif( ($linkingID eq $ResID[$allResID] ) and ($cappingInstructions{$linkingID} eq "solvent") ) {
             print OUTFILE $outputLine[$allResID];
        }
     }
}

close(OUTFILE);
#######
# EOF #
#######

