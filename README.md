# Comp-Chem-Tools
####################

**Author**: William Isley III

**Tools for Computational Chemistry**

    Scripts used to assist data analysis of computational chemistry jobs
    These scripts are released under the GNU Public License, 
    and require Perl, with dependency on Getopt::Long


**Installation**

    When using the UNIX environment, place these scripts in your $HOME/bin directory
    In your .profile, add the line:
        export $PATH:$HOME/bin 

**Executing the Scripts**

    When executing the script, type 
        perl ~/bin/scripname.pl --help
    This will call the help text, and provide usage.

**G09-QuasiHarmonic-Freq-Corr.pl**

    This script performs the quasiharmonic oscillator frequency corrections for G09 Freq jobs
    The default is to replace frequencies lower than 50 cm^-1. 
    The script can also vary recompute thermodynamic quantities changing T (temperature),
    the cutoff frequency (50 cm^-1 to any number), and introduce a universal frequency 
    scaling factor. The script can handle geometric minima and saddle points (TS).
    usage: 
        perl ~/bin/G09-QuasiHarmonic-Freq-Corr.pl --help

**Get-Atomic-Coord.pl**

    This script parses Gaussian 09 or ADF outfiles for molecular coordinates.
    Can search for lowest energy structure, or extract coordinates from a specific step.
    Check comments of code for specific option instructions. 
    usage: 
        perl ~/bin/Get-Atomic-Coord.pl --help

**Protein-Cluster-Gen.pl**

    This script can generate molecular clusters from pdb or gaussian (pdb formated) files. 
    See help file for further details.
    usage: 
        perl ~/bin/Protein-Cluster-Gen.pl --help
