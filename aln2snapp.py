#!/usr/bin/env python2

"""
<< aln2snapp.py >>


"""

from __future__ import print_function
from Bio import SeqIO, AlignIO
from sys import argv, exit
import numpy as np
import glob
from collections import Counter
import argparse

# Bases to reference other characters against
bases = ["A", "G", "C", "T"]

# Main dictionary for storing data as 0,1,2,?
# to run through SNAPP
SNAPP_Dict = {}

individuals = []

# set printing space between individual names and data
spacing = 12

def get_individuals(fasta_files):
    """
    loop through all fasta files and find all sequenced individuals.
    Create individual keys for SNAP_Dict.
    """
    for f in fasta_files:
        aln  = AlignIO.read(f, "fasta")
        inds = [ii.id for ii in aln]
        for i in inds:
            if i not in individuals:
                individuals.append(i)

    for i in individuals:
        SNAPP_Dict[i] = []

def get_biallelic_sites(alignment):
    """
    determine which sites are biallelic.
    """
    seq_dict = {}
    biallelic_sites = []
    curr_inds = [i.id for i in alignment]
    missing_ind = [i for i in individuals if i not in curr_inds]
    for a in alignment:
        seq_dict[a.id] = a.seq

    for s in range(alignment.get_alignment_length()):
        chars = list(set(alignment[:,s]))
        if "-" in chars:
            continue

        base_count = 0
        for c in chars:
            if c in bases:
                base_count += 1

        if base_count == 2 or (len(chars)==2 and base_count == 1 and ("N" not in chars)):
            #counts = Counter(list(alignment[:,s]))
            #keep = True
            #if len(counts.keys()) == 2:
            #    for k,v in count.items():
            #        if k in bases and v == 1:
            #            keep = False
            #if keep:
            biallelic_sites.append(s)

    if len(biallelic_sites) > 0:
        random_site = np.random.choice(biallelic_sites, 1)
        ref_base = seq_dict["Exaiptasia"][random_site[0]]
        for cind in curr_inds:
            if seq_dict[cind][random_site[0]] == "N":
                SNAPP_Dict[cind].append("?")
            elif seq_dict[cind][random_site[0]] == ref_base:
                SNAPP_Dict[cind].append("0")
            elif seq_dict[cind][random_site[0]] != ref_base and seq_dict[cind][random_site[0]] not in bases:
                SNAPP_Dict[cind].append("1")
            else:
                SNAPP_Dict[cind].append("2")
        for mind in missing_ind:
            SNAPP_Dict[mind].append("?")
    else:
        pass

if __name__ == '__main__':
    """
    Runs the script.
    """
    # Print docstring if only script name is given
    if len(argv) < 2:
        print(__doc__)
        exit(0)

    try:
        print("Processing *.aln files...")
        alignments = glob.glob("*.aln")
    except:
        print("ERROR:")
        print("  Couldn't find *.aln files.")
        exit(-1)
    print("Total number of alignments:", len(alignments))

    parser = argparse.ArgumentParser(description="Options for aln2snapp.py",
                                     add_help=True)
    required = parser.add_argument('-o', '--outfile', action="store", metavar='\b',
                                   required=True, help="Name of output SNAPP file")

    args = parser.parse_args()
    f_out = args.outfile
    try:
        outfile = open(f_out, 'wa')
    except:
        print("ERROR:")
        print("  Couldn't open output file: ", f_out, sep='')
        exit(-1)

    get_individuals(alignments)
    #idx = SeqIO.index(argv[1], "fasta")
    #aln = AlignIO.read(argv[1], "fasta")
    for a in alignments:
        try:
            aln = AlignIO.read(a, "fasta")
        except:
            print("Couldn't open alignment file: "+a)
            continue
        get_biallelic_sites(aln)

    print("#NEXUS", file=outfile)
    print("begin data;", file=outfile)
    print("\tdimensions ntax=",len(individuals), " nchar=", len(SNAPP_Dict.values()[0]), ";", sep='', file=outfile)
    print("\tformat datatype=integerdata missing=? gap=-;", file=outfile)
    print("matrix", file=outfile)
    for k,v in SNAPP_Dict.items():
        print(k, " " * (spacing - len(k)), sep='', end='', file=outfile)
        print(k, len(v))
        for i in v:
            print(i, sep='', end='', file=outfile)
        print("\n", end='', file=outfile)

    print(";", file=outfile)
    print("end;", file=outfile)

    """
    for k in idx.keys():
        print(">", k, ' ' * (spacing - len(k) + 1), sep='', end='')
        print(idx[k].seq)

    print("//", ' ' * 14, sep='', end='')
    for a in range(aln.get_alignment_length()):
        print(biallelic(aln[:,a]), end='')
    print("|\n", end='')
    """
