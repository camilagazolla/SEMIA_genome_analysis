## ...

...
The example files may be modified to be smaller.


## De novo Genome Assembly
...


## Download of Public Genome Sequence Data

Can consume 🐡

- [File with the GenBank release ID	accession numbers column](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/batchentrez_list.txt) from the file [assembly_result.txt](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/assembly_result_Rhizobiaceae.txt), generated generated via Send to: File / Format: ID Table after search on [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly)

Can provide 🍣

- Genome sequence files: Genomic FASTA (.fna), RNA from genomic FASTA (.fna), Genomic GenBank format (.gbff), etc.

**On web browser**

Access: https://www.ncbi.nlm.nih.gov/sites/batchentrez


## Collect RefSeq data

Use the FTP access to RefSeq to obtain additional data for our genome collection. 

Can provide 🍣
- [assembly_summary.txt](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/assembly_summary.txt), with metadata for all the genomes present at the FTP RefSeq


**On Unix/Linux terminal:**
``` 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
tail -n +2 assembly_summary.txt > assembly_summary_mod.txt
mv assembly_summary_mod.txt assembly_summary.txt
``` 


## 16S rRNA gene analysis using SILVA SSU r138_2019

Can consume 🐡

- [File](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/16S_genomic_sequences.fasta) with containing the 16S rRNA gene sequences extracted from RNA from genomic FASTA

- [TAB-delimited data](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/varnames.txt) with strain names and NCBI GenBank release ID

Can provide 🍣

- [Table](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/16S_SILVA.tab) with taxonomic assignments from order to genus

**On R:**

``` 
library(DECIPHER) 
library(dada2)

fas <- "~/16S_genomic_sequences.fasta" # path to FASTA file
fas <- readDNAStringSet(fas) # read
fas <- OrientNucleotides(fas) # reorient sequences
load("~/SILVA_SSU_r138_2019.RData") # path to SILVA_SSU_r138_2019.RData
ids <- IdTaxa(fas, trainingSet, strand="both", processors=NULL, verbose=TRUE)
ranks <- c("order", "family", "genus") # ranks of interest

# convert the output object of class "Taxa" to a matrix 
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

# insert additional information
taxid <- as.data.frame(taxid); colnames(taxid) <- ranks 
taxid$seq <- getSequences(fas) # sequence
taxid$seq_length <- nchar(taxid$seq, type = "chars", allowNA = FALSE, keepNA = NA)
taxid$access <- gsub("(GCF_\\d+\\.\\d).*", "\\1",rownames(taxid)) # id
taxid$location <-  gsub(".*locus_tag=(\\w+_\\w+).*", "\\1",rownames(taxid)) # locustag

# read archive with strain names and assembly_accession
varnames <- read.table(file = "~/varnames.txt", sep = "\t", header =T, stringsAsFactors=FALSE)

# merge varnames with taxid
taxid <- merge(taxid, varnames, by.x="access", by.y="assembly_accession") # pay attention! you might loose information here!

# remove sequences with less than 400 pb
taxid <- taxid[taxid$seq_length > 400, ]

# save
write.table(taxid, file= "16S_SILVA.tab")
``` 


## miComplete

Use [miComplete](https://doi.org/10.1093/bioinformatics/btz664) to retrieve basic statistics (size, GC-content, total number of CDS and contigs, as well as N- and L50, N- and L90), completeness and redundancy.

Can consume 🐡

- Genomic FASTA (.fna) files

Can provide 🍣

- [File with genome statistics](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/miComplete_results.tab)

**On Unix/Linux terminal:**

``` 
source micomplete/bin/activate # activate virtual environment
cd ./bins # directory containing the .fna files
find $(realpath .) -maxdepth 1 -type f -name "*.fna" | miCompletelist.sh > test_set.tab
miComplete test_set.tab --hmms Bact105 --weights Bact105 --threads 4 > miComplete_results.tab
``` 

## Genomic metrics with FastANI

Can consume 🐡

- Genomic FASTA (.fna) files

Can provide 🍣

- [File with FastANI results](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/fastANIout.txt)

**On Unix/Linux terminal:**
```
cd ./bins # directory containing the .fna files
mkdir out
ls *fna > list
mv list out
 conda activate ani
 for f in *fna; do fastANI -q "${f}" --rl out/list -o "${f}.fastANI" -t 7 --minFraction 0 ; mv "${f}.fastANI" out; done
cd out/
cat *ANI > fastANIout.txt
```

## Downstream analysis of the FastANI results with ProKlust and computation of genomic metrics with pyANI (ANIb)

Can consume 🐡

- [File with FastANI results](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/fastANIout.txt)

Can provide 🍣

- "Components files" containing the [isolated nodes and](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/components_all.csv)/or [groups formed of complete graphs](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/components_no_isolated.csv)
- pyANI results from the clusters obtained with FastANI

**On R:**
```
setwd("~/FASTANI/out/")
library(ProKlust)
library(tidyr)

# Importing FastANI results
identity <- (read.table(file = "fastANIout.txt", sep = "\t")) [1:3]
identity <- pivot_wider(identity, names_from =V1, values_from = V3)
identity <- as.data.frame(identity)
rownames(identity) <- identity$V2
identity.sorted <- identity[order(identity["V2"]),]
identity.sorted[,1] <- NULL

# Step 1. Analysis including all genomes
# Obtain clusters (this also saves the files at the current wd!)
basicResult<- prokluster(file = identity.sorted, cutoffs = 95, filterRemoveIsolated = FALSE)

# Step 2., creating list of genomes inputs for pyANI
dir.create("Input_for_pyANI")
setwd("~/FASTANI/out/Input_for_pyANI")

# Analysis excluding isolated nodes
filterIsolated <- prokluster(file = identity.sorted, cutoffs = 95, filterRemoveIsolated = TRUE)
filterIsolated$components$`components(g)$membership` <- paste0("comp",filterIsolated$components$`components(g)$membership`)

# Creatings lists to the PATH of the genomes to be analyzed together by pyANI
splitted_list <- split(filterIsolated$components, filterIsolated$components$`components(g)$membership`)

for (i in names(splitted_list)) {
 df <- splitted_list[i]
 dirarchives <- paste0("~/GenomeDB/",rownames(df[[i]])) #Full PATH to the folder containg the Genomic FASTA (.fna) files
 write(dirarchives, file = paste0("Input_for_pyANI/",i,".txt"))
}

```

**On Unix/Linux terminal:**
```
cd Input_for_pyANI

dos2unix *.txt # I work with a Windows Subsystem for Linux (WLS), so I had to correct the encoding.

mkdir genomes_comp

# compute pyANI only within components detected with FastANI/ProKlust 
for filename in *.txt; 
do mkdir genomes_comp/$(basename "$filename" .txt);
xargs -a $filename cp -t genomes_comp/$(basename "$filename" .txt);
average_nucleotide_identity.py -i genomes_comp/$(basename "$filename" .txt) -o genomes_comp/$(basename "$filename" .txt)/pyANI_out -v -m ANIb;
done
```

## Downstream analysis of the pyANI (ANIb) results

Can consume 🐡

- [assembly_summary.txt](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/assembly_summary.txt), generated at the "Collect RefSeq data" topic

Can provide 🍣

- assembly_summary_subset.tab, which is just the assembly_summary.txt including only our genome set
- ANIb_identity.csv
- 

**On R:**
```
# Call library
library(data.table)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Step 1 -  Remove the unwanted data from the assembly_summary.txt

# Read assembly_summary.txt
assembly_summary <- read.delim(file = "~PATH/assembly_summary.txt", header=T,
                               na.strings=".", stringsAsFactors=FALSE,
                               quote="")

# We need a vector containing the GenBank release ID accession numbers for the genomes included on the pyANI analysis
# I'm just extracting it for the previous data on the R env (ProKlust, Step 1, basicResult df) 

assembly_result <-  basicResult$components
assembly_result <- gsub("(GCF_...........).*", "\\1", row.names(assembly_result)) # simplify name

# Subset assembly_summary based on assembly_result
assembly_summary_subset <- setDT(assembly_summary, key = 'X..assembly_accession')[J(assembly_result)]

write.table(assembly_summary_subset, "assembly_summary_subset.tab")


# Step 2 -  Concatenate pyANI results (ANIb_percentage_identity and ANIb_percentage_identity) in a single file

setwd("~/Input_for_pyANI/genomes_comp") # PATH to the folder contaning all of the pyANI results

# Function to read ANI files
create_ANI_l <- function(files){
 ANIb_l <- list()
 for (i in files){
 name_archive <- gsub("/out/ANIb_percentage_identity.tab", "", i)
 df <-  read.table(file = i, sep = "\t", header = TRUE, row.names = 1)
 row.names(df) <- gsub("(GCF_...........).*", "\\1", row.names(df))
 names_full <- assembly_summary_subset[match(row.names(df), assembly_summary_subset$X..assembly_accession), c(1,8,9), drop=F]
 names_full$infraspecific_name <- gsub("strain=", "",names_full$infraspecific_name)
 names_full <- paste0(names_full$organism_name," ",names_full$infraspecific_name," (",names_full$X..assembly_accession, ")")
 rownames(df) <- names_full
 colnames(df) <- names_full
 ANIb_l[[name_archive]] <- df
 }
 return(ANIb_l)
}

files <- list.files(pattern = "ANIb_percentage_identity*", recursive = TRUE) # list files
ANIb_identity_l <- create_ANI_l(files = list.files(pattern = "ANIb_percentage_identity*", recursive = TRUE))
lapply(ANIb_identity_l, function(x) write.table(data.frame(x), 'ANIb_identity.csv', append= T, sep=',' ))
#Ignore the warning "appending column names to file"
  
files <- list.files(pattern = "ANIb_alignment_coverage*", recursive = TRUE) # list files
ANIb_coverage_l <- create_ANI_l(files = list.files(pattern = "ANIb_alignment_coverage*", recursive = TRUE))
lapply(ANIb_coverage_l, function(x) write.table(data.frame(x), 'ANIb_coverage.csv', append= T, sep=',' ))

```






## ...

Can consume 🐡

- [File](....) ...

Can provide 🍣

- [File](....) ...
**On R:**
```
```


## ...

Can consume 🐡

- [File](....) ...

Can provide 🍣

- [File](....) ...
**On R:**
```
```
