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

## 16S rRNA gene phylogeny using LPSN sequences

PRECISO OLHAR DE NOVO

Can consume 🐡

- [16S gene unaligned_sequences](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/unaligned_sequences.fasta) 

Can provide 🍣

- [Newick tree](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/tree.newick) 
- [Plotted tree](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/tree.pdf) and [cladogram](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/cladogram.pdf)

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

# Run Gs 
system(paste0("~/Gblocks_0.91b/Gblocks", "aligned_sequences.fasta -t=d", " -b1=",   ### -D ??????????????????????? 
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
############mt.pml <- optim.pml(mt.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)

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

# Plot tree rooted on external group
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

ggsave("cladogram.pdf", plot = last_plot(), width = 21/2.5, height = 29.7, units = "cm", limitsize = FALSE)
ggsave("cladogram.png", plot = last_plot(), width = 21/2.5, height = 29.7, units = "cm", limitsize = FALSE)

```

## Phylogenomics with GET_HOMOLOGUES and GET_PHYLOMARKERS

Can consume 🐡

- Genomic GenBank format (.gbff) files

Can provide 🍣

- [Tree plot with annotation](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/tree_annot.pdf) 

**On Unix/Linux terminal:**
```
cd get_homologues_X.Y # GET_HOMOLOGUES dir

# compute heuristic (fast) implementation of the bidirectional best-hit (BDBH) algorithm
./get_homologues.pl -b -e -n 4 -d <directory containing the .gbff files>

cd <get_homologues result directory *f0_alltaxa_algBDBHmin_e1_>

# run the GET_PHYLOMARKERS pipeline recomended for more divergent genomes
 ~/get_phylomarkers-1.3.2/run_get_phylomarkers_pipeline.sh -R 1 -t PROT

```
You well need the [IQ-TREE](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/IQT_MJRC_tree.nwk) (PhiPack/non_recomb_FAA_alns) and the [tree_labels.list](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/tree_labels.list) file generated by GET_PHYLOMARKERS. Access the [get_phylomarkers manual](https://vinuesa.github.io/get_phylomarkers/#get_homologues-get_phylomarkers-tutorials) for further information.

**On R:**
```
library("ggtree")
library(treeio)

# read tree
tree <- treeio::read.newick(file="IQT_MJRC_tree.nwk", node.label='support')

# read tree_labels
tree_labels <- read.table("tree_labels.list")
colnames(tree_labels) <- c("node", "name")
tree_labels$name <- gsub(".*\\]_(.*)","\\1",tree_labels$name) # change tip lab
tree_labels$name <- gsub("_"," ",tree_labels$name) # change tip lab

map = setNames(tree_labels$name,tree_labels$node) 

tree@phylo[["tip.label"]] <- map[tree@phylo[["tip.label"]]] #replace names

# annotate and save tree plot

  ggtree(tree, layout='daylight', size=0.1) +
   # geom_tiplab(size=1)+ # JUST TO HELP WITH ANNOT
   geom_point2(aes(subset=!isTip & support > 70), shape=21, size=1, color = "black", fill = "black") +
   # geom_text(aes(label=node), size=2, color= "red") + # JUST TO HELP WITH ANNOT
   geom_cladelabel (node=107, label= "G1", fontsize=5, hjust=0.5, vjust = 1.6, offset=-1) +
   geom_cladelabel (node=119, label= "G2", fontsize=5, hjust=0.5, vjust = 1.6) + 
   geom_cladelabel (node=121, label= "G3", fontsize=5, hjust=0.5, vjust = 1.6) +
   geom_cladelabel (node=126, label= "G4", fontsize=5, hjust=0.5, vjust = 1.6) +
   geom_cladelabel (node=131, label= "G5", fontsize=5, hjust= -0.1, vjust = 0.5) +
   geom_cladelabel (node=134, label= "G6", fontsize=5, hjust= -0.1, vjust = 0.5) +
   geom_cladelabel (node=138, label= "G7", fontsize=5, hjust= -0.1, vjust = 0.5) +
   geom_cladelabel (node=142, label= "G8", fontsize=5, hjust=0, vjust = 0, offset=-1) +
   geom_cladelabel (node=155, label= "G9", fontsize=5, hjust=0, vjust = 0, offset=-1) +
   geom_cladelabel (node=168, label= "G10", fontsize=4, hjust=0.3, vjust = -0.2, offset=-0.1) +
   geom_cladelabel (node=169, label= "G11", fontsize=4, hjust=0.5, vjust = 0, offset=-0.4) +
   geom_cladelabel (node=172, label= "G12", fontsize=5, hjust=0.5, vjust = 0, offset=-0.7) +
   geom_cladelabel (node=184, label= "G13", fontsize=5, hjust=0.5, vjust = 0, offset=-0.7) +
   geom_cladelabel (node=191, label= "G14", fontsize=5, hjust=0.5, vjust = 0) +
   geom_cladelabel (node=100, label= "G15", fontsize=5, hjust=1.1, vjust =0.5) +
   geom_tiplab(mapping=aes(subset = node %in% c(32,33,57,56,5,6,7)), offset = -0.3, size=3) +
   geom_tiplab(mapping=aes(subset = node %in% 36), offset = -0.3,angle=45, size=3)+
   geom_tiplab(mapping=aes(subset = node %in% 39), offset = -0.3,angle=55, size=3) +
   geom_tiplab(mapping=aes(subset = node %in% 40), offset = -0.3,angle=35, size=3) +
   geom_tiplab(mapping=aes(subset = node %in% 68), offset = -0.3,angle=55, size=3) +
   geom_tiplab(mapping=aes(subset = node %in% 1), offset = -0.3,angle=340, size=3) 

 ggsave("tree_annot.pdf", plot = last_plot(), width = 21, height = 21, units = "cm", limitsize = FALSE)
 ggsave("tree_annot.svg", plot = last_plot(), width = 21, height = 21, units = "cm", limitsize = FALSE)
```

## Housekeeping gene analysis

Performing phylogeny using housekeeping genes data to replace the lack of genomic sequences.

Can consume 🐡

- Genomic protein FASTA (.faa) files
- [List (.txt) data for type strain acession numbers of nucleotide sequences suitable for NCBI's Batch Entrez](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/atpD_batchentrez.txt)

Can provide 🍣

- [Tree plot with annotation]() 


First, use the NCBI's Batch Entrez to download the type strain sequences, which should be obtained from trustable sources. Check each sequence carefully. Second, to extract the genes from genomic protein FASTA sequences:

**On Unix/Linux terminal:**
```
cd [Genomic protein FASTA (.faa) files]
#add filenames to the sequences
for i in *.faa; do
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' $i > renamed_$i;
done

for i in renamed_*.faa; do
awk '/^>/ { p = ($0 ~ /recombinase RecA/)} p' $i >> recA_genomic.fasta;
done

for i in renamed_*.faa; do
awk '/^>/ { p = ($0 ~ /ATP synthase subunit beta/)} p' $i >> atpD_genomic.fasta;
done

for i in renamed_*.faa; do
awk '/^>/ { p = ($0 ~ /RNA polymerase subunit beta/)} p' $i >> rpoB_genomic.fasta;
done

for i in renamed_*.faa; do
awk '/^>/ { p = ($0 ~ /glutamine synthetase beta-grasp/)} p' $i >> glnII_genomic.fasta;
done

for i in renamed_*.faa; do
awk '/^>/ { p = ($0 ~ /RNA polymerase subunit alpha /)} p' $i >> rpoA_genomic.fasta;
done


```

Now, we should merge the files from SEMIA and the type strains in a unique FASTA file. 
We need to rename the files, to do it we firts simply them using wildcards

**On Notepad++:**
FIND: (>GCF_...........).*
REPLACE: \1

FIND: lcl\|
REPLACE:

FIND: .\d_prot.*
REPLACE:

After, we can use a dictionary to replace de names. First, write the following to a text archive named "foo.awk"

**On Notepad++:**
```
NR == FNR {
  rep[$1] = $2
  next
} 

{
  for (key in rep)
    gsub(key, rep[key])
  print
}
```

Then, we  rename the files, align with MUSCLE, concatenate the sequences and use Gblocks.

**On Unix/Linux terminal:**
```
# rename
awk -f foo.awk  dic.txt atpD.fasta > atpD_renamed.fasta
awk -f foo.awk  dic.txt recA.fasta > recA_renamed.fasta
awk -f foo.awk  dic.txt rpoB.fasta > rpoB_renamed.fasta
awk -f foo.awk  dic.txt glnII.fasta > glnII_renamed.fasta
awk -f foo.awk  dic.txt glnA.fasta > glnA_renamed.fasta

# align (if it doesnt work, delete the blanklines in the files)
muscle -in atpD_renamed.fasta -out atpD_aligned.fasta
muscle -in recA_renamed.fasta -out recA_aligned.fasta
muscle -in rpoB_renamed.fasta -out rpoB_aligned.fasta
muscle -in glnII_renamed.fasta -out glnII_aligned.fasta
muscle -in glnA_renamed.fasta -out glnA_aligned.fasta

# concatenate (maybe I can do it another way... but let it be for now)
seqkit concat rpoB_aligned.fasta atpD_aligned.fasta recA_aligned.fasta glnII_aligned.fasta glnA_aligned.fasta > concat_rpoB_atpD_recA_glnII_glnA.fasta

```

Inspect the aligment!

```
library(Biostrings)
library(phangorn)

root_name <- "Ensifer_alkalisoli_YIC4027"

alignment <- readAAStringSet("concat_rpoB_atpD_recA_glnII_glnA.fasta")

ntax = length(alignment)
b1 = floor(ntax*0.5)+1
b2 = floor(ntax*0.85)
b3 = 8
b4 = 10
b5 = "n"

# Run Gs 
system(paste0("~/Gblocks_0.91b/Gblocks", " concat_rpoB_atpD_recA_glnII_glnA.fasta ", " -b1= ", 
              b1, " -b2= ", b2, " -b3= ", b3, " -b4=", 
              b4, " -b5=", b5))	     
	     
alignment_cleaned <- readAAStringSet("concat_rpoB_atpD_recA_glnII_glnA.fasta-gb") # Ignore warning message about ignoring invalid one-letter sequence codes

alignment_cleaned_nt <- alignment_cleaned@ranges@width[1] # final number of nucleotides remaining in the alignment

# construct a neighbor-joining tree
phangAlign <- phyDat(as(alignment_cleaned, "matrix"), type = "AA") 
dm <- dist.ml(phangAlign) #compute pairwise distances
tree <- NJ(dm) # construct a neighbor-joining tree
mt <- modelTest(phangAlign, model=c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtArt", "MtZoa", "mtREV24", "VT", "RtREV", "HIVb", "FLU", "Blosum62", "Dayhoff_DCMut", "JTT_DCMut"), tree=tree, G = TRUE, I = TRUE) # preper with modelTest, for some reason "mtmam" and "HIVw" cannot be computed
mt[order(mt$AIC),]
bestmodel <- mt$Model[which.min(mt$AIC)] # choose best model from the table according to AIC
write(bestmodel, file= "bestmodel.tab")
env = attr(mt, "env")
fitStart = eval(get(bestmodel, env), env) # let R search the table

# fit the maximum likelihood tree using the neighbor-joining tree as a starting point
#check the best model
bestmodel

# CHANGE IF NEEDED!!!!!!!!!!!!
mt.pml <- pml(tree, phangAlign, model="LG", optGamma=TRUE, optInv=TRUE, optNni=TRUE, optQ=TRUE) # CHANGE IF NEEDED!!!!!!!!!!!!
# CHANGE IF NEEDED!!!!!!!!!!!!
mt.pml <- optim.pml(mt.pml)

bs = bootstrap.pml(mt.pml, bs=1000, optNni=TRUE, optInv=TRUE, multicore = TRUE)

sink("MLparameters.txt")
print(mt.pml)
sink()

library(phytools)

tree <- plotBS((reroot(fitStart$tree, (match(root_name,fitStart[["tree"]][["tip.label"]])))), bs, type="p") # ROOT ON SELECTED NODE

write.tree(tree, file = "tree.newick", append = FALSE,
           digits = 10, tree.names = FALSE)

library("ggtree")
library(treeio)
library(ggplot2)

# read tree
tree <- treeio::read.newick(file="tree.newick", node.label='support')

# annotate and save tree plot
root <- rootnode(tree) 

# plot cladogram rooted on external group
ggtree(tree, ladderize=T, layout='rectangular', size=0.2, branch.length = "none")+ 
  geom_tiplab(size=2, fontface="italic")+
  theme(legend.position = c(1,-1)) + 
  # geom_point2(aes(subset= !isTip & node != root & support > 70), 
  #             shape=21, size=1, color = "black", fill = "black")+
  xlim(0,150) + geom_label2(aes(label=support, 
                                subset = !is.na(as.numeric(support)) & as.numeric(support) > 70 &!isTip & node != root ), size =2, label.size=0,
                            nudge_x = 1,   label.padding = unit(0.0, "lines"))
```


##Pangenome analysis with anvi'o

The following section is based on the anvi’o version 2.1.0 workflow for microbial pangenomics, avaliable at https://merenlab.org/2016/11/08/pangenomics-v2/

Can consume 🐡

- Genomic FASTA (.fna) files
- [TAB-delimited (.txt) data for genomes (look at the example here!](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/external-genomes-all.txt)). Use ASCII letters, digits, and underscore... no spaces, and can't start with a digit.
- [TAB-delimited data layer (groups, species, etc.)](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/layer-additional-data-all.txt)

Can provide 🍣

- [ ]( )   ?????????????????????????

**On Unix/Linux terminal:**
```
# Activate Anvi'o 
conda activate anvio-7

# COPY THE GENOME DIR and cd to it
cd ./copy_bins # directory containing the .fna files

# Re-formatting genomes input 
for filename in *.fna; do anvi-script-reformat-fasta "${filename}" -o "${filename}.new" -l 0 --simplify-names; done

# Remove 
rm *.fna

# Rename genome .fna archives (simplify and change extension)
rename 's/\.1.*/.fna/' *

# Generate db (no external gene calls, uses prodigal, add functional annotation)
for i in `ls *fna | awk 'BEGIN{FS=".fna"}{print $1}'`
do
  anvi-gen-contigs-database -f $i.fna -o $i.db -T 200 -n "PANGENOME_SEMIA"
	anvi-run-hmms -T 40 -c $i.db
	anvi-run-kegg-kofams -T 40 -c $i.db
	anvi-run-ncbi-cogs -T 40 -c $i.db
	anvi-run-pfams -T 40 -c $i.db
	anvi-run-scg-taxonomy -T 40 -c $i.db
done
```

After this, I noticed how bad was the symbiotic genes annotations... Then, I've serched on https://www.uniprot.org/ for:

- NODULATION GENES: nodulation protein AND bacteria AND reviewed:yes
- NITROGEN FIXATION GENES: gene:nif* OR gene:vnf* OR gene:anf* "nitrogen fixation" AND reviewed:yes

I have downloaded the [Tab-separeted files](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/uniprot-gene_nif_%2BOR%2Bgene_vnf_%2BOR%2Bgene_anf_%2BAND%2Breviewed_yes%2BAND%2B_nitr--.tab) and curated the nodulation genes by hand. For the nitrogen fixation ones:

**On R:**
```
fixation_genes <- read.delim("uniprot-gene_nif_+OR+gene_vnf_+OR+gene_anf_+AND+reviewed_yes+AND+_nitr--.tab")

# extract rows with "nitrogen fixation" 
fixation_genes <- fixation_genes[grepl("nitrogen fixation", fixation_genes$Gene.ontology..biological.process.),]

# format gene names
fixation_genes$new_gene_names <- fixation_genes$Gene.names
fixation_genes$new_gene_names <- gsub("^(\\w+) .*","\\1", fixation_genes$new_gene_names)

# filter to remove Uncharacterized and blank
fixation_genes <- fixation_genes[!grepl("Uncharacterized", fixation_genes$Protein.names),]
fixation_genes <- fixation_genes[-which(fixation_genes$new_gene_names == ""), ]
 
# select the entry with max lenght for the same gene name
selected_entries <- data.frame()
for (i in unique(fixation_genes$new_gene_names)){
  df <- fixation_genes[grepl(i, fixation_genes$new_gene_names),]
  selected_entries <- rbind(selected_entries, df[which.max(df$Length),])
}

write.csv(selected_entries,"selected_genes.csv")

```
On https://www.uniprot.org/uploadlists/, the canonical and isoforms [.FASTA files](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/uniprot-yourlist_M202110234ABAA9BC7178C81CEBC9459510EDDEA32495DCL.fasta) were downloaded acording to the identifiers.

I have used wild cards to reformt the .FASTA headers quiclky on Notepad++:

FIND: ```>.*GN=(\w+).* ```

REPLACE: ```>\1[Nodulation] OR [Fixation]```

Then, the aa sequences from the .db anvi'o files were extracted.

**On Unix/Linux terminal:**
```
mkdir amino-acid-sequences
for i in `ls *db | awk 'BEGIN{FS=".db"}{print $1}'`
do
  anvi-get-sequences-for-gene-calls --get-aa-sequences -c $i.db -o amino-acid-sequences/$i.fa
done
```

And BLASTp was used to perform sequence aligment:

**On Unix/Linux terminal:**
```
# blastp
for i in *.fa; do for j in *.fasta; do blastp -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos positive" -query $i -subject $j  -evalue 0.0001 -qcov_hsp_perc 70  >> "${i%.fa}.txt" ; done; done

for i in *.txt; do sort -k1,1 -k12,12gr -k11,11g -k3,3gr $i | sort -u -k1,1 --merge > "bestHits_${i%.fa}"; done
#The first sort orders the blast output by query name then by the 12th column in descending order (bit score), then by 11th column ascending (e-value).
#The second sort picks the first line from each query. 
```

The [bestHist files](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/bestHits_GCF_014207095.txt) were further manipulated on R to match the columns needed for anvi'o.

**On R:**
```
setwd(path) # dir with bestHists

arch <- sort(list.files(path, pattern= ".txt", full.names = TRUE))

for (i in arch){
  name_arch <- gsub("~/amino-acid-sequences/bestHits_","",i)
  name_arch <- gsub(".txt","",name_arch)
  df_arch <- read.table(i)
  class(df_arch)
  df_arch <- as.data.frame(cbind(df_arch$V1, "blastp", df_arch$V2,df_arch$V2, df_arch$V11))
  df_arch$V4 <- gsub("_Avin_.*","",df_arch$V4)
  colnames(df_arch) <- c("gene_callers_id","source","accession","function","e_value")
  write.table(df_arch,paste0("bestHists_cleaned_",name_arch,".txt"), sep = "\t", quote = FALSE, row.names=FALSE)
}
```

Finally, the annot [bestHists_cleaned files](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/bestHists_cleaned_GCF_014207095.txt) were incorporated at the anvi'o .db files, along with the other annotation.

**On Unix/Linux terminal:**
```
for i in `ls *db | awk 'BEGIN{FS=".db"}{print $1}'`
do
  anvi-import-functions -c $i.db -i bestHists_cleaned_$i.txt
done

# Generate a genomes storage
anvi-gen-genomes-storage -e external-genomes-all.txt -o PANGENOME_SEMIA_ALL-GENOMES.db 
# had an issue here, solved with https://github.com/merenlab/anvio/issues/1199

# Running a pangenome analysis
anvi-pan-genome -g PANGENOME_SEMIA_ALL-GENOMES.db \
--project-name "PANGENOME_SEMIA" \
--output-dir PANGENOME_SEMIA \
--minbit 0.5 \
--mcl-inflation 2 \
--min-occurrence 1 \
--num-threads 8 \
--use-ncbi-blast \

anvi-import-misc-data layer-additional-data-all.txt \
                      -p PANGENOME_SEMIA/PANGENOME_SEMIA-PAN.db \
                      -D  default \
                      --target-data-table layers
                      
# Displaying the pan genome
anvi-display-pan -g PANGENOME_SEMIA_ALL-GENOMES.db -p PANGENOME_SEMIA/PANGENOME_SEMIA-PAN.db --server-only
```

Now, open the browser and paste http://localhost:8080/

After selecting the bins, we need to create a archive with anvi-summarize function for downstream analysis.

**On Unix/Linux terminal:**

```
anvi-summarize -p PANGENOME_SEMIA/PANGENOME_SEMIA-PAN.db \
           -g PANGENOME_SEMIA_ALL-GENOMES.db \
           -C Single_Double_Core_Etc \ ## Genome Collection Name
```


## PCA of COGs, KOs and Pfams

Can consume 🐡

- [TAB-delimited (.txt) generated with anvi-summarize](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/PANGENOME_SEMIA_gene_clusters_summary.txt)
- [TAB-delimited data layer (groups, species, etc.)](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/layer-additional-data-all.txt)

Can provide 🍣

- A series of [PCA plots for each bin/annotation source](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/Pfam_Core_plots.svg)
- [.tab with the contributions of variables](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/Pfam_Core_CONTRIB_VARIABLES.tab)
- Files with [PERMANOVA](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/Pfam_Core_PERMANOVA.TXT) AND [PERMDISP2](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/Pfam_Core_BETADISPER.TXT) statistics

**On R:**
```
# read data
gc_summary <- read.delim("PANGENOME_SEMIA_gene_clusters_summary.txt", sep = "\t")
dic <- read.table("layer-additional-data-all.txt", header = T) # group/strain data

# split bins
gc_summary_split <- split(gc_summary, gc_summary$bin_name)

# functional annotation rate
result_annot <- data.frame()
for (name in names(gc_summary_split)){
x <- gc_summary_split[[name]]

# sum results for each gene cluster with annot
result <- data.frame()
for (col in colnames(x)){
   x_tab <- as.data.frame(table(x[col])) # subset annot
   x_tab <-  subset(x_tab, Var1!="") # remove blank
   Freq <- sum(x_tab$Freq) # sum row with annot
   df <- cbind(col,Freq)
   result <- rbind(result, df)
}
result <- result[c(20,22,26,27),]
result$Group <- name
result_annot <- rbind(result_annot, result)
}

# just a quick fix
result_annot$Freq <- as.numeric(result_annot$Freq) 
result_annot$col <- gsub("_", " ", result_annot$col)
result_annot$col <- gsub("aa sequence", "Total GCs", result_annot$col)
result_annot$col <- gsub("COG20 FUNCTION", "COG20", result_annot$col)
result_annot$Group <- gsub("V10", "Variable ~10%", result_annot$Group)
result_annot$Group <- gsub("V25", "Variable ~25%", result_annot$Group)
result_annot$Group <- gsub("V50", "Variable ~50%", result_annot$Group)
result_annot$Group <- gsub("V90", "Variable ~90%", result_annot$Group)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- dim(table(result_annot$col))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

library(ggplot2)
pA <- ggplot(result_annot, aes(x = reorder(col, +Freq), y = Freq, fill = col)) + 
   geom_bar(stat = "identity") +
   facet_grid(~Group,  switch = "both") +
   labs(y = "Counts") + 
   theme_minimal(12) + scale_fill_manual(values = mycolors) +
   theme(axis.text.x = element_blank(), 
         legend.position = "bottom", legend.title = element_blank(), 
         panel.grid = element_blank(), 
         axis.title.x = element_blank(),
         plot.title = element_blank(),
         panel.spacing=unit(0.1, "lines"),
         strip.text.x = element_text(angle = 90, hjust = 1)) +
   guides(fill=guide_legend(nrow=2,byrow=TRUE))

# a look on cluster desconsidering bins
result <- list()
for (col in colnames(gc_summary)){
   x_tab <- as.data.frame(table(gc_summary[col])) # subset annot
   x_tab <-  subset(x_tab, Var1!="") # remove blank
   Freq <- sum(x_tab$Freq) # sum row with annot
   df <- cbind(col,Freq)
   result <- rbind(result, df)
}

ggsave(paste0("plot/Figure5.svg"), plot = last_plot(), width = 90, height = 120, units = "mm", limitsize = FALSE) 
ggsave(paste0("plot/Figure5.png"), plot = last_plot(), width = 90, height = 120, units = "mm", limitsize = FALSE) 
ggsave(paste0("plot/Figure5.pdf"), plot = last_plot(), width = 90, height = 120, units = "mm", limitsize = FALSE) 

# PCA 

library(dplyr)
library(tidyverse)
library(reshape2)

##! select bin/parameters to analyse !###

BIN_NAMES <- names(gc_summary_split)

annot_sources <- c("COG20_CATEGORY",
                  "COG20_FUNCTION",
                  "COG20_PATHWAY",
                  "KEGG_Class",
                  "KEGG_Module",
                  "KOfam",
                  "Pfam")
##!!! ##

pvalues <- NULL
for (BIN_NAME in BIN_NAMES){
x <- gc_summary_split[[BIN_NAME]]
for (annot_source in annot_sources){
# select
df <- NULL
df <- cbind(x$genome_name, x[[annot_source]])
colnames(df) <- c("genome_name", annot_source)
df <- as.data.frame(df)

# delete blank
df<-df[!(df[[annot_source]]==""),]

# obtain freq
freq_df <- as.data.frame(table(df))
freq_df$genome_name <- as.character(freq_df$genome_name)
freq_df$Group <- freq_df$genome_name

# add group/strain data
map = setNames(dic$group,dic$default) 
freq_df$Group <- map[freq_df$Group] #replace names
freq_df$Group <- gsub("_", " ", freq_df$Group)

levels_group <- c("G1","G2","G3","G4","SEMIA 4060","SEMIA 4080","G5","SEMIA 4064",
                  "G6","SEMIA 4085","SEMIA 4029","G7","G8","SEMIA 414","SEMIA 475","G9",
                  "SEMIA 4089","G10","G11","G12","G13","G14","SEMIA 4027","G15","34/80","SEMIA 442","SEMIA 4084")

levels_group <- intersect(levels_group,unique(freq_df$Group))
freq_df$Group <-  factor(freq_df$Group, levels=levels_group)
freq_df <- freq_df %>% group_by(genome_name) %>% mutate(relAbundByPath =Freq/sum(Freq))

# "unmelt" freq_df (need a count data set matrix or data.frame class)
freq_df_unmelt <- dcast(data = freq_df,formula = genome_name ~get(annot_source),fun.aggregate = sum,value.var = "Freq")
freq_df_unmelt <- freq_df_unmelt %>% remove_rownames %>% column_to_rownames(var="genome_name")

# some strains may be have been removed, arrange a dic cleaned
dic_cleaned <- filter(dic, default %in% rownames(freq_df_unmelt))

# cmultRepl 
library(zCompositions)

# Check if zeros are present
if(sum(colSums(freq_df_unmelt == 0)) < 1){
   
   Ab_temp_no0 <- freq_df_unmelt
   print("no zero found")
   
   } else {
      
      # Check columns with only 1 non-zero in the given data set
      checkNumZerosCol <- apply(freq_df_unmelt,2,function(x) sum(x==0))
      cases <- which(checkNumZerosCol == (nrow(freq_df_unmelt) - 1))
      
      if((length(cases) >=	1) == TRUE) {
         Ab_temp_no0 <- cmultRepl(freq_df_unmelt[,-cases], method = "GBM")
         } else {
            Ab_temp_no0 <- cmultRepl(freq_df_unmelt, method = "GBM")
         }
      }

# CLR 
library(mixOmics)
data.clr <- logratio.transfo(as.matrix(Ab_temp_no0), logratio = 'CLR', offset = 0) 
class(data.clr) <- "matrix"

# PCA
pca2c <- prcomp(data.clr, center = FALSE, scale = FALSE)
var_explained <- 100*(pca2c$sdev^2)/sum(pca2c$sdev^2) 
PC1 = format(var_explained[1], digits=2, nsmall=2) 
PC2 = format(var_explained[2], digits=2, nsmall=2) 

# graph of samples
df <- as.data.frame(pca2c$x)
df$genome_name <- rownames(df)

df$Group <- df$genome_name
df$Group <- map[df$Group]
df$Group <- gsub("_", " ", df$Group)

# Define the number of colors you want
nb.cols <- length(unique(dic$group))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

library(ggrepel)
pB <- NULL
pB <-  ggplot(df, aes(x=PC1,y=PC2)) +
   geom_point(aes_string(color="Group"), size=1) +
   theme_minimal(8) +
   labs(x = (paste0("PC1 (", PC1, "%)")), y = (paste0("PC2 (", PC2, "%)")))+
   geom_text_repel(aes_string(label = "Group", color = "Group"), size = 2, max.overlaps = 1000, segment.size=0.1) +
   theme(legend.position="none", panel.grid = element_blank()) +
   scale_color_manual(values = mycolors)

# SAVE FOR LATER
if(annot_source == "Pfam" & BIN_NAME == "Core"){
   pB_Pfam_Core <- pB 
}

# SAVE FOR LATER
if(annot_source == "Pfam" & BIN_NAME == "V50"){
   pB_Pfam_V50 <- pB 
}

# graph of variables
library(factoextra)
midpoint <- (max((facto_summarize(pca2c, element= "var", result = "contrib"))$contrib)/2)

options(ggrepel.max.overlaps = Inf)

contrib = (as.matrix(dim(pca2c$x))[2,1])
percent <- c(0.1, 0.25, 0.5, 0.75, 1)
p.list <- list() # create list to hold plots
for (n in percent){
      p <- fviz_pca_biplot(pca2c, repel = TRUE, geom.var = c("point", "text"), geom.ind = c("point"), col.ind = "black", alpha.ind = 0.1,
                        labelsize = 2, select.var = list(contrib =(as.integer(contrib*n))), col.var="contrib", axes.linetype = "blank") +
      scale_color_gradient2(low="#DDDDDD", mid="#009900", high="red", midpoint=midpoint)+
         # labs(title = (paste0("Top ", (n*100), "% contributing variables")), x = (paste0("PC1 (", PC1, "%)")), y = (paste0("PC2 (", PC2, "%)"))) +
         labs(x="", y="", title="")+
      theme_minimal(8) + theme(plot.title = element_text(size=8),panel.grid = element_blank() )
   
      p.list[[paste0(n)]] <- p
}
pD <- p.list[[1]] #10%

# SAVE FOR LATER
if(annot_source == "Pfam" & BIN_NAME == "Core"){
   pD_Pfam_Core <- pD 
}

# SAVE FOR LATER
if(annot_source == "Pfam" & BIN_NAME == "V50"){
   pD_Pfam_V50 <- pD 
}

contrib <- facto_summarize(pca2c, element= "var", result = "contrib")

write.table(contrib, file = paste0("plot/", annot_source, "_",BIN_NAME,"_CONTRIB_VARIABLES.tab"))

# multivariate Analysis
# calculate dist
library(vegan)
distab <- vegdist(data.clr, method = "euclidean")

D <- list(data = data.clr,
          D = as.matrix(distab),
          Coefficient = "Aitchison")

library(PERMANOVA)
permanova <- PERMANOVA(D, as.factor(dic_cleaned$group)) 
print(permanova)
# P-value
perMANOVA.p <- permanova[["pvalue"]]

sink(file = paste0("plot/",annot_source, "_",BIN_NAME,"_PERMANOVA.TXT"))
print(permanova)
sink()

# betadisper 
betadisperanova <- anova(betadisper(distab, as.factor(dic_cleaned$group)))  
print(betadisperanova)

sink(file = paste0("plot/",annot_source, "_", BIN_NAME,"_BETADISPER.TXT"))
print(betadisperanova)
sink()

# P-value
betadisper.p <- print(as.data.frame(betadisperanova)["Groups", "Pr(>F)"])

pvalues = rbind(pvalues, data.frame(perMANOVA.p, betadisper.p, annot_source, BIN_NAME)) # save p to a df

library(ggpubr)
pBpD <- ggarrange(pB,pD, ncol = 1)

pBpD <- annotate_figure(pBpD, bottom = text_grob(paste0("PERMANOVA p-value= ", round(perMANOVA.p, 2),
                                         " | PERMDISP2 p-value= ",  round(betadisper.p , 2)),
                                         size = 8))

annotate_figure(pBpD, top = text_grob(paste0(BIN_NAME, " (",annot_source,")"),
                                            size = 8))
					    
ggsave(paste0("plot/",annot_source,"_",BIN_NAME,"_plots.svg"), plot = last_plot(), width = 21, height = 21, units = "cm", limitsize = FALSE) #CHANGE THIS

}
}

write.table(pvalues, file = paste0("plot/pvalues_PERMANOVA_bdisper.tab"))

# "THE LATER"
pD_Pfam_V50 <- pD_Pfam_V50+theme(legend.position = "none")
pD_Pfam_Core <- pD_Pfam_Core + theme(legend.position = c(0.9, 0.95), legend.key.size = unit(0.3, "cm"))

ggarrange(pB_Pfam_Core, pB_Pfam_V50,
           pD_Pfam_Core, pD_Pfam_V50,
           labels = c("a", "b","c", "d"))

ggsave(paste0("plot/Figure6.svg"), plot = last_plot(), width = 24, height = 19, units = "cm", limitsize = FALSE)
ggsave(paste0("plot/Figure6.png"), plot = last_plot(), width = 24, height = 19, units = "cm", limitsize = FALSE)
ggsave(paste0("plot/Figure6.pdf"), plot = last_plot(), width = 24, height = 19, units = "cm", limitsize = FALSE)

```
## Symbiotic genes analysis

Can consume 🐡

- [TAB-delimited (.txt) generated with anvi-summarize](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/PANGENOME_SEMIA_gene_clusters_summary.txt)
- [TAB-delimited data layer (groups, species, etc.)](https://github.com/camilagazolla/SEMIA_genome_analysis/blob/gh-pages/layer-additional-data-all.txt)

Can provide 🍣

- ...

**On R:**
```
gc_summary <- read.delim("PANGENOME_SEMIA_gene_clusters_summary.txt", sep = "\t")

#Fixation
Fixation <- c("ccoO","cadA","ccoG","ccoS","dusB","fdxN","fixA","fixB","fixC","fixF","fixG","fixG","fixH","fixI",
"fixI1","fixI2","fixK","fixL","fixM","fixO","fixO","fixO1","fixO2","fixO3","fixQ","fixQ","fixQ1",
"fixQ2","fixQ3","fixR","fixS","fixS1","fixS2","fixT1","fixT2","fixT3","fixU","fixU","frxA",
"nifB","nifE","nifH","nifK","nifK1","nifK2","nifN","nifQ","nifR","nifS","nifT","nifU","nifV","nifW","nifW","nifX","nifZ","nifZ","sufS")
COG20_Fixation <- data.frame()
for (i in Fixation){
  df <- gc_summary[grep(i, gc_summary$COG20_FUNCTION, ignore.case=TRUE), ]
  COG20_Fixation <- rbind(COG20_Fixation,df)
}
COG20_Fixation$bin_name <- "Fixation"

#Nodulation
COG20_Nodulation <- data.frame()
Nodulation <- c("cysD","cysN","fcl","fabG","cysC","glmM","glmS","glmS2","gmd","manA","manC","maoC","nfeC","nfeD","nodA",
"nodB","nodC","nodD1","nodD2","nodD3","nodE","nodF","nodG","nodH","nodI","nodJ","nodL","nodM","nodN","nodN2",
"nodP","nodP1","nodP2","nodQ","nodQ1","nodQ2","nodS","nodU","nodV","nodY","nodZ","noeB","noeD","noeE","noeI",
"noeJ","noeK","noeL","nolA","nolB","nolE","nolG","nolK","nolL","nolM","nolN","nolO","nolR","nolU","nolV","nolY",
"nolZ","nopA","nopB","nopC","nopJ","nopL","nopM","nopP","nopT","nopX","nosF","nwsA","rhcC1","rhcL","y4nM")
for (i in Nodulation){
  df <- gc_summary[grep(i, gc_summary$COG20_FUNCTION, ignore.case=TRUE), ]
  COG20_Nodulation <- rbind(COG20_Nodulation,df)
}
COG20_Nodulation$bin_name <- "Nodulation"

#Symbiont fitness
COG20_Symbiont_fitness <- data.frame()
Symbiont_fitness <- c("betB","betB2","cyaO","chvA","attK","dinP","gabD3","hpaE","ndvA")
for (i in Symbiont_fitness){
  df <- gc_summary[grep(i, gc_summary$COG20_FUNCTION, ignore.case=TRUE), ]
  COG20_Symbiont_fitness <- rbind(COG20_Symbiont_fitness,df)
}
COG20_Symbiont_fitness$bin_name <- "Symbiont_fitness"

#Fixation,host benefit
COG20_Fixation_host_benefit <- data.frame()
Fixation_host_benefit <- c("ccoP","fixJ","fixN1","fixN","ccoN","fixN2","fixN3","fixP","fixP1","fixP2","fixP3","fixX","nifA","nifD","nifD1","nifD2","nodW","nwsB")
for (i in Fixation_host_benefit){
  df <- gc_summary[grep(i, gc_summary$COG20_FUNCTION, ignore.case=TRUE), ]
  COG20_Fixation_host_benefit <- rbind(COG20_Fixation_host_benefit,df)
}
COG20_Fixation_host_benefit$bin_name <- "Fixation_host_benefit"

library(reshape2)
library(dplyr)
library(tibble)

COG20 <- rbind(COG20_Fixation,COG20_Fixation_host_benefit,COG20_Nodulation,COG20_Symbiont_fitness)
COG20$COG20_FUNCTION <- paste0(COG20$COG20_FUNCTION," [",COG20$bin_name,"]")
COG20 <- cbind(COG20$genome_name, COG20$COG20_FUNCTION) #simplify
COG20 <- as.data.frame(COG20)
COG20 <- as.data.frame(table(COG20)) # subset

# "unmelt"
COG20 <- dcast(data = COG20,formula = V1 ~get("V2"),fun.aggregate = sum,value.var = "Freq")
COG20 <- COG20 %>% remove_rownames %>% column_to_rownames(var="V1")
COG20 <- as.matrix(COG20)
colSums(COG20)

#fix names, inspect results too!
#write.csv(colnames(COG20), file="new_colname.csv")
new_colname <- read.csv(file="new_colname.csv")
colnames(COG20) <- new_colname$new

COG20 <- as.data.frame(COG20)
# remove (i've write remove in the col)
COG20 <- COG20 %>% select(-contains("remove"))

library("pheatmap")
# creating annotation_col
dic <- read.table("layer-additional-data-all.txt", header = T) # group/strain data
rownames(dic) <- dic$default
dic <- dic[order(dic$default),]
dic$default <- NULL
colnames(dic) <- c("Cluster","Genus")
dic$Genus <- NULL

#set colors
ann_colors = list(
  Cluster = c(G1="#009999", G2="#00FFDD", G3="#7810D2", G4="#CC92FC", SEMIA_4060="#660066", SEMIA_4080="#9900CC",
              G5="#B177B1", SEMIA_4064="#6600CC", G6="#3333CC", SEMIA_4085="#FF00C4", SEMIA_4029="#C367F5",
              G7="#660066", G8="#FF0066", SEMIA_414="#CC00CC", SEMIA_475="#493DCC", G9="#B879B8", SEMIA_4089="#D9348F",
              G10="#5236CF", G11="#D10096", G12="#5F71B8", G13="#A612A3", G14="#924AFF", SEMIA_4027="#FF0000",
              G15="#8C8832", '34/80'="#000000", SEMIA_442="#4AAEFF", SEMIA_4084="#157D31"))

#  scale
COG20 <- scale(COG20)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

pheatmap(COG20, 
         annotation_row = dic, 
         annotation_colors = ann_colors,
         clustering_callback = sort_hclust,
         fontsize = 5,
         fontsize_col = 5,
         angle_col = "315",
         treeheight_col = 20,
         treeheight_row = 100,
         colorRampPalette(c("#BF0080", "#CE6EAE", "#dddddd", "#6EAE6E", "#008000"))(n = 8)
         )

dim(COG20)

# gonna check the presence of these genes in 1# bin
gc_summary_bin1 <- read.delim("/Users/Agrega/Desktop/Pangenoma/ANVIO/Anvio_db/PANGENOME_SEMIA/SUMMARY_Single_Double_Core_Etc/PANGENOME_SEMIA_gene_clusters_summary.txt", sep = "\t")
unique(gc_summary_bin1$bin_name)

#have to remove the annot between [] modifications i've added to the csv... 
new_colname$old <- gsub(" \\[Nodulation\\]","",new_colname$old)
new_colname$old <- gsub(" \\[Fixation\\]","",new_colname$old)
new_colname$old <- gsub(" \\[Fixation_host_benefit\\]","",new_colname$old)
new_colname$old <- gsub(" \\[Symbiont_fitness\\]","",new_colname$old)

gextracted <- gc_summary_bin1[gc_summary_bin1$COG20_FUNCTION %in% new_colname$old ,]

gextracted <- cbind(gextracted$bin_name,gextracted$COG20_FUNCTION)
gextracted <- as.data.frame(gextracted)
gextracted <- as.data.frame(table(gextracted))
gextracted <-  subset(gextracted, Freq!=0) # remove 0

# remove FabG [Nodulation] - NAD(P)-dependent dehydrogenase, short-chain alcohol dehydrogenase family (FabG) (PDB:6L1H) (MUITAS entradas!)
gextracted <-  subset(gextracted, V2!="NAD(P)-dependent dehydrogenase, short-chain alcohol dehydrogenase family (FabG) (PDB:6L1H)") # remove
View(gextracted)

gextracted$V2 <- as.character(gextracted$V2)

# replace using new_colname$old
for (i in 1:dim(new_colname)[1]){
for (j in 1:dim(gextracted)[1]){
  
 if(gextracted[j,2] == new_colname$old[i]){
   gextracted[j,2] <- new_colname$new[i]
 }
  }
}

# just a quick fix
gextracted$V1 <- gsub("V10", "Variable ~10%", gextracted$V1)
gextracted$V1 <- gsub("V25", "Variable ~25%", gextracted$V1)
gextracted$V1 <- gsub("V50", "Variable ~50%", gextracted$V1)
gextracted$V1 <- gsub("V90", "Variable ~90%", gextracted$V1)

library(ggplot2)
ggplot(gextracted, aes(x = reorder(V2, +Freq), y = Freq, fill = V2)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~V1,  strip.position = "right", scales = "free_y", ncol = 1) +
  labs(y = "") + 
  theme_minimal(8) + scale_y_continuous(position = "right") +
  theme(axis.text.x = element_text(size = 6, angle = 270, vjust = 0, hjust = 0),
        legend.position = "none", legend.title = element_blank(), 
        axis.title.x = element_blank(),
        plot.title = element_blank())

ggsave(paste0("Figure8.svg"), plot = last_plot(), width = 190, height = 220, units = "mm", limitsize = FALSE) 
ggsave(paste0("Figure8.png"), plot = last_plot(), width = 190, height = 220, units = "mm", limitsize = FALSE) 
ggsave(paste0("Figure8.pdf"), plot = last_plot(), width = 190, height = 220, units = "mm", limitsize = FALSE) 
```
