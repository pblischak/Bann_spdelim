#!/usr/bin/env python2

"""
<< add_outgroup.py >>

Adds an outgroup sequence to a FASTA file combined with sequences
from the same locus locus in a pyRAD .loci file.
Each locus gets its own FASTA file that needs to be aligned.
For this we used Muscle (Edgar, 2004).

Arguments
---------
  - ref <string>   : name of outgroup reference sequences in FASTA format
  - loci <string>  : name of .loci file from pyRAD
  - blast <string> : name of file with BLAST results between .loci and
                     ref files
  - name <string>  : name of outgroup

Output
------
The output from the script is a new directory entitled 'FASTA/' with
a separate FASTA file for each locus. These can be aligned with Muscle
using the following command.

  cd FASTA/
  for f in *.fasta; do muscle -in $f -out $f.aln; done
"""

from __future__ import print_function
from Bio import SeqIO
import argparse
from sys import argv, exit

################ OLD CODE ################
# make an index for accessing the scaffolds of the reference genome
#idx1 = SeqIO.index("aiptasia_genome.scaffolds.fa", "fasta")
#loci = open("BA_FULL_snapp.loci").read().split("|\n")[:-1]
#headers = [i for i in idx.keys()]
#
# Read in the blast output
#blast_output = open("BA_FULL_snapp.txt").readlines()[:-1]
#
#for hit in blast_output:
#    info   = hit.split()
#    query  = info[0]
#    locus_number = int(query.split("_")[1]) - 1
#    outfasta = open("FASTA/locus_"+str(locus_number+1)+".fasta", 'wa')
#    ref    = info[1]
#    rstart = int(info[8]) - 1
#    rend   = int(info[9]) - 1
#    if rstart > rend:
#        print(">", "Exaiptasia", "\n", idx1[ref].seq[rend:rstart].reverse_complement(), sep='', file=outfasta)
#    else:
#        print(">", "Exaiptasia", "\n", idx1[ref].seq[rstart:rend], sep='', file=outfasta)
#
#    for locus in loci[locus_number].split("\n")[:-1]:
#        print(locus.split()[0], "\n", locus.split()[1], sep='', file=outfasta)
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
    parser = argparse.ArgumentParser(description="Options for add_outgroup.py",
                                     add_help=True)
    required = parser.add_argument_group("required arguments")
    required.add_argument('-r', '--ref', action="store", type=str, metavar='\b',
                          required=True, help="Reference genome in FASTA format")
    required.add_argument('-l', '--loci', action="store", type=str, metavar='\b',
                          required=True, help="PyRAD .loci file")
    required.add_argument('-b', '--blast', action="store", type=str, metavar='\b',
                          required=True, help="Blast results from ")
    required.add_argument('-n', '--name', action="store", type=str, metavar='\b',
                          required=True, help="Name of outgroup")

    args      = parser.parse_args()
    ref       = args.ref
    locfile   = args.loci
    blastfile = args.blast
    name      = args.name

    # make an index for accessing the scaffolds of the reference genome
    idx1 = SeqIO.index(ref, "fasta")

    with open(locfile, 'r') as f_in:
        loci = f_in.read().split("|\n")[:-1]

    with open(blastfile, 'r') as b_in:
        blast_output = b_in.readlines()[:-1]

    for hit in blast_output:
        info   = hit.split()
        query  = info[0]
        locus_number = int(query.split("_")[1]) - 1
        try:
            outfasta = open("FASTA/locus_"+str(locus_number+1)+".fasta", 'a')
        except:
            print("Couldn't make file: FASTA/locus_"+str(locus_number+1)+".fasta")
            continue
        ref    = info[1]
        rstart = int(info[8]) - 1
        rend   = int(info[9]) - 1
        if rstart > rend:
            print(">", name, "\n", idx1[ref].seq[rend:rstart].reverse_complement(), sep='', file=outfasta)
        else:
            print(">", name, "\n", idx1[ref].seq[rstart:rend], sep='', file=outfasta)

        for locus in loci[locus_number].split("\n")[:-1]:
            print(locus.split()[0], "\n", locus.split()[1], sep='', file=outfasta)
