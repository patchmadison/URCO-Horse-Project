# Code based on: https://www.nemabiome.ca/dada2_workflow



# libraries needed
library(DECIPHER)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr)
library(readr)
library(vegan)

# seed set for replication capabilities 
set.seed(106)

# Path to Miseq_Run files, adjust as needed
path <- "~/Miseq_Run"

# Define forward files vs reverse files
fwd_files <- sort(list.files(path, pattern = "R1", full.names = T))
rev_files <- sort(list.files(path, pattern = "R2", full.names = T))

# extract sample names from file names
samples = str_extract(basename(fwd_files), "^[^_]+")

# assign sample names to forward and reverse files
names(fwd_files) <- samples
names(rev_files) <- samples

# Define primer sequences and reverse complements
fwd_primer <- "ACGTCTGGTTCAGGGTTGTT"
rev_primer <- "TTAGTTTCTTTTCCTCCGCT"
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))

# confirm that primer is detected by counter function
count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = F)
  return(sum(num_hits> 0))
}
count_primers(fwd_primer, fwd_files[[1]])
count_primers(rev_primer, rev_files[[1]])



## Primer Removal 
# via cutadapt
cutadapt <- path.expand("~/cutadapt.exe")
#test that it is downloaded and works
system2(cutadapt, args = "--version")

# output directory to store clipped files
cut_dir <- file.path(path, "cutadapt")
if(!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(fwd_files))
rev_cut <- file.path(cut_dir, basename(rev_files))

names(fwd_cut) <- samples
names(rev_cut) <- samples

#log files 
cut_logs <- path.expand(file.path(cut_dir, paste0(samples, ".log")))

# cut adapt parameters
# -g: sequence to trim off the 5' end of the forward read (forward primer)
# -a: sequence to trim off the 3' end of the forward read (reverse complemented reverse primer)
# -G: sequence to trim off the 5' end of the reverse read (reverse primer)
# -A: sequence to trim off the 3' end of the reverse read (reverse complemented forward primer)
# -n: remove multiple primer hits
# --discard- untrimmed: only keep reads that contain a primer


cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
              "-G", rev_primer, "-A", fwd_primer_rev, 
              "-n", 2, "--discard-untrimmed")

# loop through files with cutadapt
for(i in seq_along(fwd_files)) {
  system2(cutadapt,
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i],
                   fwd_files[i], rev_files[i]),
  stdout = cut_logs[i])
}
head(list.files(cut_dir))


##Quality filtering

# check quality of first two reads
plotQualityProfile(fwd_cut[1:2]) + ggtitle("Forward")

plotQualityProfile(rev_cut[1:2]) + ggtitle("Reverse")

# Filter Reads
# output directory to store filtered files
filt_dir <- file.path(path, "filtered")
if(!dir.exists(filt_dir)) dir.create(filt_dir)

fwd_filt <- file.path(filt_dir, basename(fwd_files))
rev_filt <- file.path(filt_dir, basename(rev_files))

names(fwd_filt) <- samples
names(rev_filt) <- samples

filtered_out <- filterAndTrim(
  fwd = fwd_cut,
  filt = fwd_filt,
  rev = rev_cut,
  filt.rev = rev_filt,
  maxEE = c(2,5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = F
  )
head(filtered_out)


## Main workflow

# DADA2 learns error profile
err_fwd <- learnErrors(fwd_filt, multithread = T)
err_rev <- learnErrors(rev_filt, multithread = T)
plotErrors(err_fwd, nominalQ = T)

# denoise sequences
dada_fwd <- dada(fwd_filt, err = err_fwd, multithread = T)
dada_rev <- dada(rev_filt, err = err_rev, multithread = T)

# merge paired-end reads
mergers <- mergePairs(
  dadaF = dada_fwd,
  dadaR = dada_rev,
  derepF = fwd_filt,
  derepR = rev_filt,
  maxMismatch = 1,
  verbose = T
)

# create sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                    multithread = T, verbose = T)
dim(seqtab_nochim)

# saved after chimeras removed for ease of repeating steps afterwards
saveRDS(seqtab_nochim, file = "seqtab_nochim_28June24.rds")



## How did we do?

table(nchar(getSequences(seqtab_nochim)))

# funciton to get number of sequences
getN <- function(x) sum(getUniques(x))
track <- cbind(
  filtered_out,
  sapply(dada_fwd, getN),
  sapply(dada_rev, getN),
  sapply(mergers, getN),
  rowSums(seqtab_nochim)
)

# shows number of reads at each step
colnames(track) <- c("raw", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chim")
rownames(track) <- samples
track
saveRDS(track, file = "track_28June24.rds")

track <- readRDS("track_28June24.rds")



### The following code included additional samples unrelated to the project
### therefore indexing and certain sample names are not relevant to published research



## sum of total reads raw all samples and controls (except WAS -- unrelated samples) and duplicates: 790996
sum(track[c(-12,-25,-32,-48:-52),"raw"])

## percent of reads when finished filtered: 75.81%
(sum(track[c(-12,-25,-32,-48:-52),"no_chim"])/sum(track[c(-12,-25,-32,-48:-52),"raw"]))

## average # of reads for samples (no controls or WAS -- unrelated samples): 15,447.1
## standard deviation: 7,246.598
mean(track[c(-12,-25,-32,-43:-52),"no_chim"])
sd(track[c(-12,-25,-32,-43:-52),"no_chim"])


## Assign taxonomy with IDTAXA

# train and tax have been pulled from nemabiome website
train <- readDNAStringSet("~/Practice Dada2 pipeline/Nematode_ITS2_1.0.0_idtaxa.fasta")
tax <- read_tsv("~/Practice Dada2 pipeline/Nematode_ITS2_1.0.0_idtaxa.tax")

# train the classifier
trainingSet <- LearnTaxa(train, names(train), tax)
saveRDS(trainingSet, file = "trainingSet_28June24.rds")

seqtab_nochim <- readRDS("seqtab_nochim_28June24.rds")

trainingSet <- readRDS("trainingSet_28June24.rds")

dna <- DNAStringSet(getSequences(seqtab_nochim))

## tried with a lower threshold, see how it compares. load in training set and dna
## this threshold was used for subsequent analyses

idtaxa_lowerthreshold <- IdTaxa(dna,
                 trainingSet,
                 strand = "both",
                 threshold = 50,
                 bootstraps = 100,
                 processors = NULL,
                 verbose = T,
                 type = "extended")

ranks <- c("rootrank", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
taxid_lowerthreshold <- t(sapply(idtaxa_lowerthreshold, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid_lowerthreshold) <- ranks
rownames(taxid_lowerthreshold) <- getSequences(seqtab_nochim)

taxid_lowerthreshold <- taxid_lowerthreshold[,-1]

write.csv(taxid_lowerthreshold, file = "taxonomy_id_threshold_50_raw_18July24.csv")

# originally recommended parameters from nemabiome website, but threshold lowered for subsequent analyses
idtaxa <- IdTaxa(dna,
                 trainingSet,
                 strand = "both",
                 threshold = 60,
                 bootstraps = 100,
                 processors = NULL,
                 verbose = T,
                 type = "extended")


#from dada2 pipeline website
ranks <- c("rootrank", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
taxid <- t(sapply(idtaxa, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab_nochim)

taxid <- taxid[,-1]

write.csv(taxid, file = "taxonomy_id_raw_29June24.csv")


library(phyloseq)
samp_data <- data.frame(
  row.names = samples,
  sample = samples
)
write.csv(samp_data, "samp_metadata.csv")

### change samp_metadata.csv column to say sample name instead, then change code to match
samp_data <- read.csv("samp_metadata.csv")
samp_data <- samp_data[,-1]
row.names(samp_data) <- samples

asvs <- paste0("ASV_", 1:length(dna))
rownames(taxid) <- asvs
colnames(seqtab_nochim) <- asvs
names(dna) <- asvs

# creating phyloseq objects

# original physeq
physeq = phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = F),
  tax_table(taxid),
  sample_data(samp_data),
  dna
)
physeq


saveRDS(physeq, "physeq_raw_29June24.rds")

physeq <- readRDS("physeq_raw_29June24.rds")

# physeq with lower threshold taxid (ltt)
asvs <- paste0("ASV_", 1:length(dna))
rownames(taxid_lowerthreshold) <- asvs
colnames(seqtab_nochim) <- asvs
names(dna) <- asvs

physeq_ltt = phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = F),
  tax_table(taxid_lowerthreshold),
  sample_data(samp_data),
  dna
)
physeq_ltt

write_rds(physeq_ltt, file = "physeq_lower_taxonomy_threshold_19July24.rds")


physeq_ltt <- readRDS("physeq_lower_taxonomy_threshold_19July24.rds")

# tried but not used in analyses

library(phangorn)
alignment <- AlignSeqs(dna, anchor = NA)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data = phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, 
                    rearrangement = "stochastic", control = pml.control())
detach("package:phangorn", unload = TRUE)


# physeq with tree added
physeq_tree = phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = F),
  tax_table(taxid),
  sample_data(samp_data),
  phy_tree(fitGTR$tree),
  dna
)
physeq_tree



saveRDS(physeq_tree, file = "physeq_tree_5July24.rds")
physeq <- readRDS("physeq_tree_5July24.rds")



rank_names(physeq)
table(tax_table(physeq)[,"phylum"], exclude = NULL) # 451 Nematoda, 28 NA

## sums of asvs
asv_sums <- c()
for (asv in 1:479) {
  asv_sums <- c(asv_sums, sum(otu_table(physeq)[,asv]))
}
asv_sums


# remove WAS samples (not relevant samples)
WAS <- c("WAS03", "WAS06", "WAS07", "WAS11", "WAS14")
physeq_ltt_filt <- subset_samples(physeq_ltt, !Host %in% WAS)  


# remove all L3 duplicates
physeq_ltt_filt <- subset_samples(physeq_ltt_filt, Duplicate == "N")



## # of reads that are assigned to phylum, etc
#

# phylum
ps_NA_phylum <- subset_taxa(physeq_ltt_filt, is.na(phylum))
ps_nonNA_phylum <- subset_taxa(physeq_ltt_filt, !is.na(phylum))


reads_NA_phylum <- sum(taxa_sums(ps_NA_phylum)) # reads that phylum NA: 83
reads_nonNA_phylum <- sum(taxa_sums(ps_nonNA_class)) # reads that phylum not NA: 602,383
reads_nonNA_phylum/(reads_NA_phylum + reads_nonNA_phylum) # % assigned to phylum: 99.986%


# family
ps_NA_family <- subset_taxa(physeq_ltt_filt, is.na(family))
ps_nonNA_family <- subset_taxa(physeq_ltt_filt, !is.na(family))


reads_NA_family <- sum(taxa_sums(ps_NA_family)) # reads that family NA: 99
reads_nonNA_family <- sum(taxa_sums(ps_nonNA_family)) # reads that family not NA: 602,367
reads_nonNA_family/(reads_NA_family + reads_nonNA_family) # % assigned to family: 99.984%

# genus
ps_NA_genus <- subset_taxa(physeq_ltt_filt, is.na(genus))
ps_nonNA_genus <- subset_taxa(physeq_ltt_filt, !is.na(genus))


reads_NA_genus <- sum(taxa_sums(ps_NA_genus)) # reads that genus NA: 14652
reads_nonNA_genus <- sum(taxa_sums(ps_nonNA_genus)) # reads that genus not NA: 587,814
reads_nonNA_genus/(reads_NA_genus + reads_nonNA_genus) # % assigned to genus: 97.568%

# species
ps_NA_species <- subset_taxa(physeq_ltt_filt, is.na(species))
ps_nonNA_species <- subset_taxa(physeq_ltt_filt, !is.na(species))


reads_NA_species <- sum(taxa_sums(ps_NA_species)) # reads that species NA: 22,751
reads_nonNA_species <- sum(taxa_sums(ps_nonNA_species)) # reads that species not NA: 579,715
reads_nonNA_species/(reads_NA_species + reads_nonNA_species) # % assigned to species: 96.224%

# % asvs assigned to phylum
phylum_table <- table(tax_table(physeq_ltt_filt)[,"phylum"], exclude =NULL)
phylum_table[-length(phylum_table)]/sum(phylum_table) # 96.65% assigned phylum

#family
family_table <- table(tax_table(physeq_ltt_filt)[,"family"], exclude =NULL)
family_table[-length(family_table)]/sum(family_table) # 96.37% assigned family

# genus
genus_table <- table(tax_table(physeq_ltt_filt)[,"genus"], exclude =NULL)
sum(genus_table[-length(genus_table)])/sum(genus_table) # 84.36% assigned genus

# species
species_table <- table(tax_table(physeq_ltt_filt)[,"species"], exclude =NULL)
sum(species_table[-length(species_table)])/sum(species_table) # 77.37% assigned family


ctrls <- subset_samples(physeq_ltt, Host == 'CTRL')
ctrls_reads <- sample_sums(ctrls)

mean(ctrls_reads) ## mean 5.8 reads
sd(ctrls_reads) ## sd 7.95


## remove asvs with less than 15 reads in sample 
head(otu_table(physeq_ltt_filt), 10)

# if greater than 15, keep value, else change to 0
rm_low <- function(x){
  ifelse(x >= 16, x, 0)
}

# transform counts with function above to remove low counts in each sample
x1 <- transform_sample_counts(physeq_ltt_filt, rm_low)
head(otu_table(x1), 10)
# remove taxa with taxa with no reads
ps <- prune_taxa(taxa_sums(x1) > 0, x1)
head(otu_table(ps), 10)



# remove NA phylum taxa
psa <- subset_taxa(ps, !is.na(phylum))
psa <- subset_samples(psa, Host != "CTRL")

plot_bar(psa, x= "Sample_name")


meta_df1 <- data.frame(sample_data(psa))

# add column with # of reads
meta_df1$Reads <- sample_sums(psa)

meta_df1$Dry_weight_EPG <- meta_df1$Percent_dry * meta_df1$EPG


## adding # of reads to metadata
meta_df <- data.frame(sample_data(ps_final))

# add column with # of reads
meta_df$Reads <- sample_sums(ps_final)

meta_df$Dry_weight_EPG <- meta_df$Percent_dry * meta_df$EPG



saveRDS(psa, "psa_filtered_29Aug24.rds")

write.csv(meta_df, "metadata_29Aug24.csv")

psa <- read_rds("psa_filtered_29Aug24.rds")




