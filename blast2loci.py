#!/usr/bin/env python2

"""
<< blast2loci.py >>

Take the results of BLASTing the pyRAD loci against the algal
genome to remove potential contaminants and write a new .loci
file.

Arguments
---------
  - infile <string>  : name of the input .loci file
  - blast <string>   : name of the BLAST results file
  - outfile <string> : name for the new .loci file
"""

from __future__ import print_function
import argparse
from sys import argv, exit

################ OLD CODE ################
## read in data file (.loci)
#infile = open("Bann_Rapid2.loci", "r")

## create the output file
#blastfile = open("Bann_Rapid2_aiptasia_out.txt", "r")

#outfile = open("Bann_Rapid2-test.loci", "w")

## parse the loci in your file
#loci = infile.read().split("|\n")[:-1]
#hits = blastfile.read().split("\n")[:-1]
#prev_hit = -9

## write the first sequence from each locus to the output file in fasta format
#for hit in hits:
#	#print(hit)
#	curr_hit = int(hit.split("\t")[0].split("_")[1])
#	if curr_hit == prev_hit:
#		pass
#	else:
#		for l in loci[curr_hit-1].split("\n")[:-1]:
#			print(l, file=outfile)
#		print(loci[curr_hit-1].split("\n")[-1], "|\n", end='', sep='', file=outfile)
#		#print(loci[curr_hit-1], file=outfile, sep='')
#		#prev_hit = curr_hit
#	#name, seq = reads[0].split()
#	#locus_counter = locus_counter + 1
#	#print(">locus_",locus_counter,"\n"+seq, file=outfile, sep='')
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
	parser = argparse.ArgumentParser(description="Options for blast2loci.py",
	                                 add_help=True)
	required = parser.add_argument_group("required arguments")
	required.add_argument('-i', '--infile', action="store", type=str, metavar='\b',
	                      required=True, help="Name of input .loci file")
	required.add_argument('-b', '--blast', action="store", type=str, metavar='\b',
	                      required=True, help="Name of BLAST results file")
	required.add_argument('-o', '--outfile', action="store", type=str, metavar='\b',
	                      required=True, help="Name for new .loci file")

	# Parse arguments and pass to variables
	args      = parser.parse_args()
	infile    = args.infile
	blastfile = args.blast
	outfile   = args.outfile

	with open(infile, 'r') as f_in:
		loci = f_in.read().split("|\n")[:-1]

	with open(blastfile, 'r') as b_in:
		hits = b_in.read().split("\n")[:-1]

	prev_hit = -9
	with open(outfile, 'a') as f_out:
		for hit in hits:
			curr_hit = int(hit.split("\t")[0].split("_")[1])
			if curr_hit == prev_hit:
				pass
			else:
				for l in loci[curr_hit-1].split("\n")[:-1]:
					print(l, file=f_out)
				print(loci[curr_hit-1].split("\n")[-1], "|\n", end='', sep='', file=f_out)
