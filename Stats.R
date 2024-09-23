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


# loading in RDS objects
ps_final <- readRDS("ps_final_12Sept24.rds")
metadata <- read.csv("metadata_29Aug24.csv")
psa <- readRDS("psa_filtered_29Aug24.rds")

# updating metadata
metadata <- metadata[,-1]
row.names(metadata) <- metadata$Sample_name
metadata$X._of_L3_counts

colnames(metadata)[colnames(metadata) == "X._of_L3_counts"] <- "Number_of_L3_counts"

## redefine metadata in ps final and psa objects

sample_data(ps_final) <- metadata

sample_data(psa) <- metadata




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

## keep all for asv richness
## remove na species richness
## remove na genera richness


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


