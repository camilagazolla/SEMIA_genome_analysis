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

mkdir genomes

# compute pyANI only within components detected with FastANI/ProKlust 
for filename in *.txt; 
do mkdir genomes/$(basename "$filename" .txt);
xargs -a $filename mv -t genomes/$(basename "$filename" .txt);
average_nucleotide_identity.py -i genomes/$(basename "$filename" .txt) -o genomes/$(basename "$filename" .txt)/out -v -m ANIb;
done
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


## ...

Can consume 🐡

- [File](....) ...

Can provide 🍣

- [File](....) ...
**On R:**
```
```
