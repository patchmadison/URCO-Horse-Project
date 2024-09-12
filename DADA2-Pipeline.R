# libraries needed
library(DECIPHER)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr)
library(readr)
library(vegan)
set.seed(106)


path <- "C:/Users/madis/OneDrive/Documents/URCO-Horse-Project/Miseq_Run"
fwd_files <- sort(list.files(path, pattern = "R1", full.names = T))
rev_files <- sort(list.files(path, pattern = "R2", full.names = T))

samples = str_extract(basename(fwd_files), "^[^_]+")

names(fwd_files) <- samples
names(rev_files) <- samples

fwd_primer <- "ACGTCTGGTTCAGGGTTGTT"
rev_primer <- "TTAGTTTCTTTTCCTCCGCT"
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))


count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = F)
  return(sum(num_hits> 0))
}
count_primers(fwd_primer, fwd_files[[1]])
count_primers(rev_primer, rev_files[[1]])



## Primer Removal
cutadapt <- path.expand("C:/Users/madis/OneDrive/Desktop/cutadapt.exe")
system2(cutadapt, args = "--version")

cut_dir <- file.path(path, "cutadapt")
if(!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(fwd_files))
rev_cut <- file.path(cut_dir, basename(rev_files))

names(fwd_cut) <- samples
names(rev_cut) <- samples


cut_logs <- path.expand(file.path(cut_dir, paste0(samples, ".log")))
cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
              "-G", rev_primer, "-A", fwd_primer_rev, 
              "-n", 2, "--discard-untrimmed")
for(i in seq_along(fwd_files)) {
  system2(cutadapt,
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i],
                   fwd_files[i], rev_files[i]),
  stdout = cut_logs[i])
}
head(list.files(cut_dir))


##Quality filtering
plotQualityProfile(fwd_cut[1:2]) + ggtitle("Forward")

plotQualityProfile(rev_cut[1:2]) + ggtitle("Reverse")

# Filter Reads
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
err_fwd <- learnErrors(fwd_filt, multithread = T)
err_rev <- learnErrors(rev_filt, multithread = T)
plotErrors(err_fwd, nominalQ = T)

dada_fwd <- dada(fwd_filt, err = err_fwd, multithread = T)
dada_rev <- dada(rev_filt, err = err_rev, multithread = T)

mergers <- mergePairs(
  dadaF = dada_fwd,
  dadaR = dada_rev,
  derepF = fwd_filt,
  derepR = rev_filt,
  maxMismatch = 1,
  verbose = T
)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                    multithread = T, verbose = T)
dim(seqtab_nochim)

saveRDS(seqtab_nochim, file = "seqtab_nochim_28June24.rds")

## How did we do?
table(nchar(getSequences(seqtab_nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(
  filtered_out,
  sapply(dada_fwd, getN),
  sapply(dada_rev, getN),
  sapply(mergers, getN),
  rowSums(seqtab_nochim)
)
colnames(track) <- c("raw", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chim")
rownames(track) <- samples
track
saveRDS(track, file = "track_28June24.rds")

track <- readRDS("track_28June24.rds")

## sum of total reads raw all samples and controls except WAS and duplicates: 790996
sum(track[c(-12,-25,-32,-48:-52),"raw"])

## percent of reads when finished filtered: 75.81%
(sum(track[c(-12,-25,-32,-48:-52),"no_chim"])/sum(track[c(-12,-25,-32,-48:-52),"raw"]))

## average # of reads for samples (no controls or WAS): 15,447.1
## standard deviation: 7,246.598
mean(track[c(-12,-25,-32,-43:-52),"no_chim"])
sd(track[c(-12,-25,-32,-43:-52),"no_chim"])


## Assign taxonomy with IDTAXA
train <- readDNAStringSet("~/Practice Dada2 pipeline/Nematode_ITS2_1.0.0_idtaxa.fasta")
tax <- read_tsv("~/Practice Dada2 pipeline/Nematode_ITS2_1.0.0_idtaxa.tax")

trainingSet <- LearnTaxa(train, names(train), tax)
saveRDS(trainingSet, file = "trainingSet_28June24.rds")

seqtab_nochim <- readRDS("seqtab_nochim_28June24.rds")

trainingSet <- readRDS("trainingSet_28June24.rds")

dna <- DNAStringSet(getSequences(seqtab_nochim))

## try with a lowerer threshold, see how it compares. load in training set and dna

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


# remove was samples
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


## remove asvs with less than 15 reads in sample 
head(otu_table(physeq_ltt_filt), 10)

# if greater than 16, keep value, else change to 0
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




########Reads stats


plot_bar(subset_samples(psa, Sample_type == "Fecal"), x= "Sample_name")

plot_bar(subset_samples(psa, Sample_type == "Swab"), x= "Sample_name")

plot_bar(subset_samples(psa, Sample_type == "Larva"), x= "Sample_name")


ps_final <- subset_samples(psa, sample_sums(psa) > 1000)




#meta_df1 <- data.frame(sample_data(psa))
#meta_df1$Reads <- sample_sums(psa)
#meta_df1$Dry_weight_EPG <- meta_df1$Percent_dry * meta_df1$EPG


# boxplot of reads by method
readsxmethod <- ggplot(data = meta_df, mapping = aes(x = Sample_type, y = Reads)) 
readsxmethod + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = "red")+
  labs(title = "# of Reads by Processing Method",
       x = "Processing Method", y = "# of Reads")






readsxmethod

readsxmethod_aov <- aov(Reads ~ Sample_type, data = meta_df)
summary(readsxmethod_aov)

TukeyHSD(readsxmethod_aov)


readsxepg <- ggplot(data = meta_df1, mapping = aes(x = EPG, y = Reads)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "# of Reads by Eggs Per Gram")


readsxepg

readsxdry_epg <- ggplot(data = meta_df1, mapping = aes(x = Dry_weight_EPG, y = Reads)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "# of Reads by Eggs Per Gram")

readsxdry_epg


readsxepgxmethod <- ggplot(data = meta_df1, mapping = aes(x = EPG, y = Reads, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "# of Reads by Eggs Per Gram Grouped by Method")
readsxepgxmethod


readsxdry_epgxmethod <- ggplot(data = meta_df1, mapping = aes(x = Dry_weight_EPG, y = Reads, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "# of Reads by Eggs Per Gram Grouped by Method")
readsxdry_epgxmethod

interaction.plot(x.factor = meta_df$EPG, trace.factor = meta_df$Sample_type,
                 response = meta_df$Reads, ylab = "# of Reads", xlab = "Eggs per Gram",
                 trace.label = "Method", xpd = T, xtick = T)

interaction_aov <- aov(Reads ~ EPG * Sample_type, data = meta_df1)

readsxdry_epgxmethod_aov <- aov(Reads ~ Dry_weight_EPG * Sample_type, data = meta_df1)

summary(readsxdry_epgxmethod_aov)
summary(interaction_aov)


count_sample_type <- data.frame(count(meta_df, Sample_type))
count_sample_type

aov(Sample_type ~ n, data = count_sample_type)


########## Subsampling and diversity stats


min_sample.size <- min(sample_sums(ps_final))


ps_rarified <- rarefy_even_depth(ps_final, min_sample.size, replace = F)

ps_rarified_replaced <- rarefy_even_depth(ps_final, min_sample.size, replace = T)



plot_bar(ps_rarified, x="Sample_name")



# glom species, na **NOT** removed ** 
ps_rarified_sp_glom = tax_glom(ps_rarified, "species", NArm = F)
ps_rarified_sp_glom

ps_rarified_replace_sp_glom = tax_glom(ps_rarified_replaced, "species", NArm = F)
ps_rarified_replace_sp_glom


plot_bar(ps_rarified_sp_glom, x="Sample_name", fill="species")
plot_bar(physeq_ltt_filt_species_glom, x="Sample_name", fill="species")



richness_rarified <- estimate_richness(ps_rarified_sp_glom)

ps_final_rarified_df <- data.frame(sample_data(ps_rarified_replace_sp_glom))

ps_final_rarified_df$Richness <- richness_rarified[,"Observed"]



richnessxmethod <- ggplot(data = ps_final_rarified_df, mapping = aes(x = Sample_type, y = Richness)) 
richnessxmethod + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = "red")+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness")

richnessxhost <- ggplot(data = ps_final_rarified_df, mapping = aes(x = Host, y = Richness)) 
richnessxhost + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = "Sample_type")+
  labs(title = "Richness by Host",
       x = "Host", y = "Richness")





richnessxmethodxhost <- ggplot(data = ps_final_rarified_df, mapping = aes(x = Host, y = Richness, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Richness of Host Grouped by Method")
richnessxmethodxhost



summary(aov(Richness ~ Sample_type,data = ps_final_rarified_df))
TukeyHSD(aov(Richness ~ Sample_type,data = ps_final_rarified_df))




####################things that were tried but not used####################


## initial processing that not used now
# number of reads with NA for phylum
psNA <- subset_taxa(physeq, is.na(phylum))
plot_bar(psNA, x="Sample_name")
# plot of reads associated with genus labelled
plot_bar(physeq, x="Sample_name", fill="genus")

# number of reads with NA for phylum ** lower threshold
psNA_ltt <- subset_taxa(physeq_ltt, is.na(phylum))
plot_bar(psNA_ltt, x="Sample_name")
# plot of reads associated with genus labelled ** lower threshold
plot_bar(physeq_ltt, x="Sample_name", fill="genus")
plot_bar(physeq_ltt, x="Sample_name", fill="species")

physeq_ltt_filt <- subset_taxa(physeq_ltt, !is.na(phylum))
plot_bar(physeq_ltt_filt, x = "Sample_name")

spsCont <- subset_samples(physeq_ltt_filt, Host == "CTRL")
psC<-prune_taxa(taxa_sums(psCont)>0, psCont)
plot_bar(psC, x = "Sample_name", fill="genus")


head(sample_data(psNA))
unique(sample_data(physeq)$Sample_type)

prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x) {sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
plyr::ddply(prevdf, "family", function(df1) {cbind(mean(df1$Prevalence), sum(df1$Prevalence))})
ps1 = transform_sample_counts(ps0, function(x){x / sum(x)})
plot_bar(ps1, x = "Sample_name", fill = "species")


ptree = plot_tree(ps1, method = "treeonly", ladderize = "left")
ptree





## pruning taxa

## not what was intended, so not using
#physeq_ltt_filt <- subset_samples(physeq, Host =)
#physeq_ltt_filt<-prune_taxa(taxa_sums(physeq_ltt)>15, physeq_ltt)
#plot_bar(physeq_ltt_filt, x = "Sample_name", fill="genus")
#taxa_sums(physeq_ltt_filt)

### same here
#physeq_ltt_filt <- subset_taxa(physeq_ltt, taxa_sums(physeq_ltt) > 15)
#plot_bar(subset_samples(physeq_ltt_filt, Host =="CTRL"), x = "Sample_name", fill = "genus")


# remove samples that had 15 reads or less AKA controls removed
physeq_ltt_filt <- subset_samples(physeq_ltt_filt, sample_sums(physeq_ltt_filt) > 15)


### * remove percentage of reads or certain number of reads from each sample

#removed total reads with 15 or less
physeq_ltt_filt<-prune_taxa(taxa_sums(physeq_ltt)>15, physeq_ltt)



# glom species, na removed
physeq_ltt_filt_species_glom_narm = tax_glom(physeq_ltt_filt, "species", NArm = TRUE)
physeq_ltt_filt_species_glom_narm

# glom species, na not removed
physeq_ltt_filt_species_glom = tax_glom(physeq_ltt_filt, "species", NArm = F)
physeq_ltt_filt_species_glom

plot_bar(physeq_ltt_filt_species_glom, x="Sample_name", fill="genus")
plot_bar(physeq_ltt_filt_species_glom, x="Sample_name", fill="species")


### remove the LB samples no duplicates









### ATTEMPTING ALPHA DIVERSITY AND ANALYSIS




# richness for species 
##* also need to do genus glom and asvs
estimate_richness(physeq_ltt_filt)
richness <- estimate_richness(physeq_ltt_filt_species_glom)
richness['Sample_type'] <- sample_data(physeq_ltt_filt_species_glom)[,"Sample_type"]

write.csv(richness, file = "Sample_richness.csv")


# not needed for now
larva_samples <- subset_samples(physeq_ltt_filt_species_glom, Sample_type == "Larva")
fecal_samples <- subset_samples(physeq_ltt_filt_species_glom, Sample_type == "Fecal")
swab_samples <- subset_samples(physeq_ltt_filt_species_glom, Sample_type == "Swab")


# richness plot
boxplot(Observed ~ Sample_type,
        data = richness,
        main = "Richness by Sample Type",
        xlab = "Sample Type",
        ylab = "# of Species"
)

### do in ggplot, overlay points, etc
## find metric for just evenness

# shannon index box plot
boxplot(Shannon ~ Sample_type,
        data = richness,
        main = "Shannon Index by Sample Type",
        xlab = "Sample Type",
        ylab = "Shannon Index"
)

# fisher index box plot
boxplot(Fisher ~ Sample_type,
        data = richness,
        main = "Fisher Index by Sample Type",
        xlab = "Sample Type",
        ylab = "Fisher Index"
)




## attempting pcoa plot



## change meta data to include eggs per gram 
## wet weights vs dry weights





