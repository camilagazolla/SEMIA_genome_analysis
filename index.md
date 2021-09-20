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
- [ANIb_identity.csv](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/ANIb_identity.csv) and [ANIb_coverage.csv](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/ANIb_coverage.csv), ugly but usefull single archives for the pyANI results computed for all the folders on genomes_comp

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

## 16S rRNA gene phylogeny using

Can consume 🐡

- [FASTA with 16S gene unaligned_sequences](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/unaligned_sequences.fasta) ...

Can provide 🍣

- [File](....) ...

**On Unix/Linux terminal:**
```
sina -i unaligned_sequences.fasta -r ~/SILVA_138.1_SSURef_opt.arb -o aligned_sequences.fasta -p 30
```

**On R:**
```
library(DECIPHER)
library(phangorn)
library("ggtree")
library("tidyverse")
library(phytools)
library(phangorn)
library(treeio)

root_name <- "Mesorhizobium_loti_NBRC_14779_AB680660" # name of the root sequence

setwd("~/SEMIAS_Rhizobium_Agrobacterium") # PATH to the folder contaning the 16S sequences FASTA file 

alignment <- readRNAStringSet("aligned_sequences.fasta")

ntax = length(alignment)
b1 = floor(ntax*0.5)+1
b2 = floor(ntax*0.85)
b3 = 8
b4 = 10
b5 = "n"

# Run Gblock 
system(paste0("~/Gblocks_0.91b/Gblocks", "aligned_sequences.fasta -t=d", " -b1=", 
              b1, " -b2=", b2, " -b3=", b3, " -b4=", 
              b4, " -b5=", b5))
              
alignment_cleaned <- readRNAStringSet("aligned_sequences.fasta-gb") # Ignore warning message about ignoring invalid one-letter sequence codes

alignment_cleaned_nt <- alignment_cleaned@ranges@width[1] # final number of nucleotides remaining in the alignment

# view the alignment in a browser (optional)
BrowseSeqs(alignment_cleaned, highlight=0)

# construct a neighbor-joining tree
phangAlign <- phyDat(as(alignment_cleaned, "matrix"))
dm <- dist.ml(phangAlign) #compute pairwise distances
tree <- NJ(dm) # construct a neighbor-joining tree
mt <- modelTest(phangAlign, tree=tree, G = TRUE, I = TRUE, k = 4) # preper with modelTest
mt[order(mt$AIC),]
bestmodel <- mt$Model[which.min(mt$AIC)] # choose best model from the table according to AIC
write(bestmodel, file= "bestmodel.tab")
env = attr(mt, "env")
fitStart = eval(get(bestmodel, env), env) # let R search the table

# fit the maximum likelihood tree using the neighbor-joining tree as a starting point
mt.pml <- pml(tree, phangAlign, model=bestmodel, k=4)
mt.pml <- optim.pml(mt.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)

bs = bootstrap.pml(mt.pml, bs=500, optNni=TRUE, multicore = TRUE)

sink("MLparameters.txt")
print(mt.pml)
sink()

tree <- plotBS((reroot(fitStart$tree, (match(root_name,fitStart[["tree"]][["tip.label"]])))), bs, p = 50, type="p") # ROOT ON SELECTED NODE

write.tree(tree, file = "tree.newick", append = FALSE,
           digits = 10, tree.names = FALSE)

# Read tree
tree <- treeio::read.newick(file="tree.newick", node.label='support')
tree@phylo[["tip.label"]] <- gsub("_", " ", tree@phylo[["tip.label"]])
root <- rootnode(tree)  

ggtree(tree, ladderize=T, layout='rectangular', size=0.2)+ 
  geom_tiplab(size=2, fontface="italic")+
  theme(legend.position = c(1,-1)) + 
  geom_point2(aes(subset=!isTip & node != root & support > 70), 
              shape=21, size=0.5, color = "black", fill = "black") +
  xlim(-.1, 0.5) # + ylim(-10,35) # Change to improve the tree

ggsave("tree.pdf", plot = last_plot(), width = 21, height = 29.7, units = "cm", limitsize = FALSE)
ggsave("tree.png", plot = last_plot(), width = 21, height = 29.7, units = "cm", limitsize = FALSE)

# plot cladogram rooted on external group
ggtree(tree, ladderize=T, layout='rectangular', size=0.2, branch.length = "none")+ 
  geom_tiplab(size=2, fontface="italic")+
  theme(legend.position = c(1,-1)) + 
  geom_point2(aes(subset= !isTip & node != root & support > 70), 
              shape=21, size=0.5, color = "black", fill = "black")+
   xlim(0,150) # + ylim(-15,45)  # Change to improve the tree

ggsave("cladogram_bola.pdf", plot = last_plot(), width = 21/2.5, height = 29.7, units = "cm", limitsize = FALSE)
ggsave("cladogram_bola.png", plot = last_plot(), width = 21/2.5, height = 29.7, units = "cm", limitsize = FALSE)

```


## ...

Can consume 🐡

- [File](....) ...

Can provide 🍣

- [File](....) ...
**On R:**
```
```
