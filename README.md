## Python scripts for processing *B. annulata* ddRADseq data

Required python modules:

 - Biopython
 - numpy

To see the description for any of the scripts, type the name of the script
without any arguments.

**Example:**

```bash
# Print docstring for the parse_loci.py script
parse_loci.py
```

```bash
## << parse_loci.py >>
##
## Converts a .loci file from pyRAD into a FASTA file by taking
## the first sequence from every locus. The locus names are
## sequentially numbered.
##
## Arguments
## ---------
##   - infile <string>  : name of input loci file.
##   - outfile <string> : name of output FASTA file.
```

### `parse_loci.py`

Convert a `.loci` file from pyRAD to a FASTA file for BLASTing
against the *Exaiptasia* genome.

### `add_outgroup.py`

Add matching *Exaiptasia* sequence from BLAST and add to a
FASTA file for aligning.

### `blast2loci.py`

Take BLAST output

### `aln2snapp.py`

Convert aligned FASTA files (ending in `.aln`) for each locus back to `.loci`
format to complete pyRAD step 7 (generating output files).

## Data Processing Steps

```bash
# 1. Get loci to BLAST against outgroup genome
parse_loci.py -i Bann_Rapid2.loci -o Bann_Rapid2.fasta

# 2. Run BLAST by making database and
makeblastdb
blastn

# 3. Remove any loci that matched the algal symbiont
#    and reformat in .loci format to re-run step 7
#    in pyRAD
blast2loci.py

# 4. Add outgroup sequences using BLAST searches of
#    the RAD loci. This
add_outgroup.py -r aiptasia_genome.scaffolds.fa -l BA_FULL_snapp.loci \
                -b BA_FULL_snapp.txt

# 5.
cd FASTA/
for f in *.fasta; do muscle -in $f -out $f.aln; done

# 6. Convert these *.aln files into NEXUS format for running SNAPP.
fasta2snapp.py -o BANN_FULL_snapp.nex

```
