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

# Run Gblocks 
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

## Pangenome analysis with anvi'o

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
