#!/usr/bin/env python2

"""
<< parse_loci.py >>

Converts a .loci file from pyRAD into a FASTA file by taking
the first sequence from every locus. The locus names are
sequentially numbered.

Arguments
---------

	- infile <string>  : name of input loci file.
	- outfile <string> : name of output FASTA file.
"""

from __future__ import print_function
import argparse
from sys import argv, exit

################ OLD CODE ################
## read in data file (.loci)
#infile = open("Bann_Rapid2.loci", "r")

## create the output file
#outfile = open("Bann_Rapid2.fasta", "w")

## parse the loci in your file
#loci = infile.read().split("|\n")[:-1]

## write the first sequence from each locus to the output file in fasta format
#locus_counter = 0
#for loc in loci:
#	reads = loc.split("\n")
#	name, seq = reads[0].split()
#	locus_counter = locus_counter + 1
#	print(">locus_",locus_counter,"\n"+seq, file=outfile, sep='')
################ OLD CODE ################

if __name__ == "__main__":
	"""
	Runs the script.
	"""
	# Print docstring if only script name is given
	if len(argv) < 2:
		print(__doc__)
		exit(0)

	# Setup argument parser
	parser = argparse.ArgumentParser(description="Options for parse_loci.py",
	                                 add_help=True)
	required = parser.add_argument_group("required arguments")
	required.add_argument('-i', '--infile', action="store", type=str, metavar='\b',
	                      required=True, help="Input .loci file")
	required.add_argument('-o', '--outfile', action="store", type=str, metavar='\b',
	                      required=True, help="Name for output FASTA file")

	# Parse arguments and pass to variables
	args    = parser.parse_args()
	infile  = args.infile
	outfile = args.outfile

	# Open input and output files, then process .loci file
	# and print to the output FASTA file
	with open(infile, 'r') as f_in, open(outfile, 'a') as f_out:
		loci = f_in.read().split("|\n")[:-1]
		for loc in loci:
			reads = loc.split("\n")
			name, seq = reads[0].split()
			locus_counter = locus_counter + 1
			print(">locus_",locus_counter,"\n"+seq, file=f_out, sep='')
