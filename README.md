# Script to identify peptide orthologs in human genome based on BLOSUM
# similarity and minimum number of mismatches


## Summary

*Inputs*: human proteome fasta file (downloaded by `src/run_PeptideOrtholog_finder.R`)

*Outputs*: csv file with human peptides that look like the target, 
with blosum scores and number of mismatches


## Process

Through R studio or on command line, run the `src/run_PeptideOrtholog_finder.R` script.

Number of mismatches is set to 4, but can be changed in the code. Note that
increasing number of mismatches has an non-linear impact on compute time. 

The script loads the human proteome reference fasta file. It generates wildcarded 
versions of the peptide (e.g., I*S*FL*LL, etc), then grep's (using linux command 
line) all of the wildcarded peptides against the proteome file simultaneously 
to find the lines that contain any peptide hits. It writes a new temporary 
proteome subset file that contains only those lines and then greps each WC 
peptide individually against the proteome subset. Then it adds the gene/protein
information back in from the reference file and compute the BLOSUM distances
(biochemical similarity) between the target and each mismatched hit from the 
proteome. Finally, it writes the matches to file 
(```outputs/${peptide}_(mm=${n_mismatches}).csv```). Some messages will print to
stdout, and log files in ```log/${peptide}.txt``` 
will track when analysis is completed.


## Requirements

### Run from linux OS - tested in Ubuntu 22.04.4 LTS and macOS 14.1.2

### R packages (CRAN unless specified) - version tested on in brackets
- Biostrings (from Bioconductor) (V2.64.1)
- fastmatch (V1.1.3)
- seqinr (V4.2.23)
- data.table (V1.14.6)
- parallel (V4.2.3)
- httr (V1.4.6)

## Runtime
3 minutes when tested, but will vary by peptide target and required
number of mismatches

## Authors
- [Andrew Walker](mailto:andysw90@gmail.com)
- [Agatha Treveil](mailto:Agatha.Treveil@immunocore.com)
