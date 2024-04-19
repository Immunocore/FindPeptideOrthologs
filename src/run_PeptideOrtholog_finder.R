# Script to identify peptide orthologs in human genome based on BLOSUM
# similarity and minimum number of mismatches
# Author: Andrew Walker, Agatha Treveil


# Load libraries and data -------------------------------------------------
library(Biostrings)
library(fastmatch)
library(data.table)
library(parallel)
library(seqinr)
#library(R.utils)
library(httr)

data(BLOSUM62)

dir.create("log")
dir.create("outputs")
dir.create("tmp")
dir.create("data")

source("src/PeptideOrthologs_fn.R")

# Set some parameters -----------------------------------------------------
# Put your peptide of interest here
peptide <- "RLPAKAPLL"

n_mismatches <- 4

# Source the input file --------------------------------------------------
# Put your genome of interest here
#uniprot_proteome_rest  <- "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29" #HBV
uniprot_proteome_rest  <- "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=proteome=UP000001584" #TB
fname_human            <- "data/UP000005640_9606.fasta"
fname_human_onePerLine <- "data/UP000005640_9606_onePerLine.fasta" 
if(!file.exists(fname_human)){
  response <- GET(uniprot_proteome_rest)
  writeLines(rawToChar(response$content),fname_human)
}

# Load data -------------------------------------------------
message("Loading data")

if(!file.exists(fname_human_onePerLine)){
  aln_raw <- read.alignment(fname_human,format="fasta")
  writeLines(paste0(">",aln_raw$nam,"\n",toupper(aln_raw$seq)),fname_human_onePerLine)
}

# Load alignments
HUMANseq <- read.alignment(fname_human_onePerLine,format="fasta")

# Process alignments as strings
genome_raw <- readLines(fname_human_onePerLine)
headers <- genome_raw[grep(">",genome_raw)]

blosum_names <- as.character(colnames(BLOSUM62))

# Run analysis -------------------------------------------------
writeLog("Running analysis...",peptide)

fname_out <- paste0("outputs/",peptide,"_(mm=",n_mismatches,").csv")

orthologs <- FindOrthologs(peptide = peptide,fname_genome = fname_human_onePerLine,n_mismatches = n_mismatches)
write.table(orthologs,file = fname_out)
