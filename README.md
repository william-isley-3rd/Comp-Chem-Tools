# Comp-Chem-Tools
Tools for Computational Chemistry
====================

Scripts used to assist data analysis of computational chemistry jobs

Scripts Contributed By William Isley III

Get-Atomic-Coord.pl 

    This script parses Gaussian 09 or ADF outfiles for molecular coordinates.
    Can search for lowest energy structure, or extract coordinates from a specific step.
    Check comments of code for specific option instructions. 
    usage: perl Get-Atomic-Coord.pl --help
    
 
G09-QuasiHarmonic-Freq-Corr.pl

    This script performs the quasiharmonic oscillator frequency corrections for G09 Freq jobs. 
    Default is 50 cm^-1. Script can vary temperature, cutoff frequency, 
    introduce a universal frequency scaling factor, and can handle a TS calculation. 
    usage: perl G09-QuasiHarmonic-Freq-Corr.pl --help

Protein-Cluster-Gen.pl   

    This script can generate molecular clusters from pdb or gaussian (pdb formated) files. 
    See help file for further details.
    usage: perl Protein-Cluster-Gen.pl --help
