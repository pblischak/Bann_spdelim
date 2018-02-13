## Python scripts for processing *B. annulata* ddRADseq data

**Manuscript:**

Titus, B. M., P. D. Blischak, and M. Daly. Genomic signatures of
sympatric speciation with historical and contemporary gene flow in a tropical
anthozoan.

----

These scripts will process ddRADseq loci output by the program pyRAD in two ways:
(1) removing algal contaminants identified through BLAST searches and (2) adding
a homologous outgroup sequence to each locus (identified using BLAST).

All of the scripts were run using Python v2.7.

**Required Python modules:**

 - Biopython
 - numpy

**Additional dependencies:**

 - [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 - [Muscle](https://www.drive5.com/muscle/downloads.htm)

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
against the anemone/reference genome to identify cnidarian-only
loci.

### `blast2loci.py`

Take BLAST output for identifying anemone loci and make a new
`.loci` file that only contains loci that were cnidarian in
origin.

### `add_outgroup.py`

Add homologous outgroup sequence identified by BLAST to each RAD locus.

### `aln2snapp.py`

Convert aligned FASTA files (ending in `.aln`) for each locus into a matrix of
SNPs for input to the program SNAPP. If more than one SNP is present at a locus,
we randomly sample a single SNP.

## Data Processing Steps

```bash
##########################################
#### Extracting anemone specific loci ####
##########################################

# 1. Get loci to BLAST against outgroup genome
parse_loci.py -i Bann_Rapid2.loci -o Bann_Rapid2.fasta

# 2. Run BLAST by making a database and aligning ddRADseq
#    loci with BLASTN.
makeblastdb -dbtype nucl -in anemone-ref.fasta -out Anemone
blastn -db Anemone -query Bann_Rapid2.fasta -out Bann_Rapid2.txt -outfmt 6

# 3. Extract loci that matched the Anemone reference
#    and reformat in .loci format to re-run step 7
#    in pyRAD.
blast2loci.py -i Bann_Rapid2.loci -b Bann_Rapid2.txt -o Bann_Rapid2-decon.loci

#########################################
#### Adding an un-sequenced outgroup ####
#########################################

# 1. Make a database for the outgroup and BLAST the RAD loci
#    against this reference to add the homologous regions.
makeblastdb -dbtype nucl -in aiptasia_genome.scaffolds.fa -out Aiptasia
blastn -db Aiptasia -query Bann_Rapid2.fasta -out BA_FULL_snapp.txt -outfmt 6

# 2. Add outgroup sequences using BLAST searches of
#    the RAD loci. These are all written to new FASTA
#    files, one per locus, in a folder .
add_outgroup.py -r aiptasia_genome.scaffolds.fa -l BA_FULL_snapp.loci \
                -b BA_FULL_snapp.txt

# 3. After grabbing the corresponding outgroup sequence,
#    we realign all of the FASTA files using Muscle.
cd FASTA/
for f in *.fasta; do muscle -in $f -out $f.aln; done

# 4. Convert these *.aln files into NEXUS format for running SNAPP
#    by grabbing one randomly selected SNP per locus.
aln2snapp.py -o BANN_FULL_snapp.nex
```
