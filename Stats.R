# libraries needed
library(DECIPHER)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(phyloseq)
library(emmeans)
library(lme4)
library(tidyverse)
library(vegan)
library(betareg)
library(statmod)
library(multcomp)
set.seed(106)



# loading in RDS objects
ps_final <- readRDS("ps_final_12Sept24.rds")
metadata <- read.csv("metadata_29Aug24.csv")
psa <- readRDS("psa_filtered_29Aug24.rds")
epg_data <- read.csv("EPG_data_2Oct24.csv")


# loading in new metadata
metadata_raw_reads <- read.csv("metadata_raw_reads_25Oct24.csv")
metadata_rarefied <- read.csv("metadata_rarefied_16Oct24.csv")

# updating metadata
metadata <- metadata[,-1]
row.names(metadata) <- metadata$Sample_name
metadata$X._of_L3_counts
colnames(metadata)[colnames(metadata) == "X._of_L3_counts"] <- "Number_of_L3_counts"

metadata$Avg_L3_count <- metadata$L3_count_sum / metadata$Number_of_L3_counts
metadata$DC_wet_EPG <- metadata$DC_egg_count / metadata$Wet_weight_dc_egg_count
metadata$Dry_dc_EPG <- metadata$DC_wet_EPG * metadata$Percent_dry
  
## redefine metadata in ps final and psa objects

sample_data(ps_final) <- metadata
sample_data(psa) <- metadata

## relabeling dataframe columns for consistency




metadata_final <- metadata %>%
  select(Sample_name, Sample_type, Host, Wet_weight, Dry_weight, Percent_dry,
         EPG, Dry_weight_EPG, Wet_weight_dc_egg_count, DC_egg_count, DC_wet_EPG, 
         Dry_dc_EPG, L3_count_sum, Number_of_L3_counts, Avg_L3_count, Qubit, 
         Sub_conc, Reads)

metadata_final <- metadata_final %>%
  rename(
    EPG_McM_wet = EPG,
    EPG_McM_dry = Dry_weight_EPG,
    Weight_dc_count = Wet_weight_dc_egg_count,
    Egg_count_DC = DC_egg_count,
    EPG_DC_wet = DC_wet_EPG,
    EPG_DC_dry = Dry_dc_EPG
  )

## create new data frame with changed style of rows etc
## such a pain to do please clean up for later ***
epg_data <- NA
epg_data_McM_wet <- metadata_final %>%
  mutate(
    Method = "McMaster",
    Weight_type = "Wet",
    EPG = EPG_McM_wet
  ) %>%
  select(Host, EPG, Method, Weight_type)

epg_data_McM_dry <- metadata_final %>%
  mutate(
    Method = "McMaster",
    Weight_type = "Dry",
    
    EPG = EPG_McM_dry
  ) %>%
  select(Host, EPG, Method, Weight_type)
epg_data_DC_dry <- metadata_final %>%
  mutate(
    Method = "Double Centrifugation",
    Weight_type = "Dry",
    
    EPG = EPG_DC_dry
  ) %>%
  select(Host, EPG, Method, Weight_type)
epg_data_DC_wet <- metadata_final %>%
  mutate(
    Method = "Double Centrifugation",
    Weight_type = "Wet",
    
    EPG = EPG_DC_wet
  ) %>%
  select(Host, EPG, Method, Weight_type)

epg_data_not_finished <- rbind(epg_data_DC_dry, epg_data_DC_wet, epg_data_McM_dry, epg_data_McM_wet)



epg_data <- epg_data_not_finished %>%
  distinct()


write.csv(epg_data, "EPG_data_2Oct24.csv")
write.csv(metadata_final, "metadata_2Oct24.csv")

metadata <- read.csv("metadata_2Oct24.csv")
epg_data <- read.csv("EPG_data_2Oct24.csv")

## goodness gracious that was such a pain

## time to redo the plots and stats based on the new and improved dataframe


## reads by method and dry DC epg
readsxepgxmethod <- ggplot(data = metadata_final, mapping = aes(x = EPG_DC_dry, y = Reads, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Number of Reads by Eggs Per Gram Grouped by Method", x = 'Eggs Per Gram')+
  theme_bw()
readsxepgxmethod

readsxmethod_aov <- aov(Reads ~  Sample_type, data = metadata_final)
summary(readsxmethod_aov)
TukeyHSD(readsxmethod_aov)
 ## need to change the plot to a box plot instead 

readsxepgxmethodmdl <- aov(Reads ~ Sample_type * EPG_DC_dry, data = metadata_final)
summary(readsxepgxmethodmdl)

readsxmethodmdl <- aov(Reads ~ Sample_type, data = metadata_final)
summary(readsxmethodmdl)

anova(readsxmethodmdl, readsxepgxmethodmdl) ## not significant difference when EPG added p = 0.057






## might do more graphs later 
##but now onto diversity!


min_sample.size <- min(sample_sums(ps_final))
ps_rarified <- rarefy_even_depth(ps_final, min_sample.size)  ## seed has been set to 106

# combine asvs or multiples to on group of each level 
# with default settings NA remove 
ps_rarified_sp = tax_glom(ps_rarified, "species")
ps_rarified_g = tax_glom(ps_rarified, "genus")

# plot each level
plot_bar(ps_rarified, x = "Sample_name")
plot_bar(ps_rarified_sp, x="Sample_name", fill="species")
plot_bar(ps_rarified_g, x="Sample_name", fill="genus")

# get richness stats
asv_richness_rarified <- estimate_richness(ps_rarified)
sp_richness_rarified <- estimate_richness(ps_rarified_sp)
g_richness_rarified <- estimate_richness(ps_rarified_g)

# create new data frame with richness added
## look at doing this a different way
# assign sample(phyloseq) to meta data
metadata_rarefied <- metadata[-c(3,36,37),] 
metadata_rarefied$Richness_asv <- asv_richness_rarified[,"Observed"]
metadata_rarefied$Richness_sp <- sp_richness_rarified[,"Observed"]
metadata_rarefied$Richness_g <- g_richness_rarified[,"Observed"]

metadata_rarefied$Shannon_asv <- asv_richness_rarified[,"Shannon"]
metadata_rarefied$Shannon_sp <- sp_richness_rarified[,"Shannon"]
metadata_rarefied$Shannon_g <- g_richness_rarified[,"Shannon"]

metadata_rarefied$Evenness_asv <- metadata_rarefied$Shannon_asv / log(metadata_rarefied$Richness_asv)
metadata_rarefied$Evenness_sp <- metadata_rarefied$Shannon_sp / log(metadata_rarefied$Richness_sp)
metadata_rarefied$Evenness_g <- metadata_rarefied$Shannon_g / log(metadata_rarefied$Richness_g)
## Evenness is on scale of 0 to 1

mean(metadata_rarefied$Richness_asv) # 29
mean(metadata_rarefied$Richness_sp) # 7.22
mean(metadata_rarefied$Richness_g) # 4.92

sd(metadata_rarefied$Richness_asv) # 16.58
sd(metadata_rarefied$Richness_sp) # 3.53 
sd(metadata_rarefied$Richness_g) # 1.83


# start plots
# allllll richness
ggplot(metadata_rarefied, aes(x = Sample_type, y = Richness_asv, color = Host)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  labs(title = "Richness (ASVs) vs. Method")

richnessxmethod_asv <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Richness_asv)) 
richnessxmethod_asv + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness (# of ASVs)")
richnessxmethod_sp <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Richness_sp)) 
richnessxmethod_sp + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness (# of species)")
richnessxmethod_g <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Richness_g)) 
richnessxmethod_g + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness (# of genera)")



## try a plot?
asv_poisson <- glm(Richness_asv ~ Sample_type + EPG_DC_dry, data = metadata_rarefied, family = poisson())
summary(asv_poisson)
asv_emmeans <- emmeans(asv_poisson, pairwise ~ Sample_type)
summary(asv_emmeans)
sp_poisson <- glm(Richness_sp ~ Sample_type + EPG_DC_dry, data = metadata_rarefied, family = poisson())
summary(sp_poisson)
sp_emmeans <- emmeans(sp_poisson, pairwise ~ Sample_type)
summary(sp_emmeans)
g_poisson <- glm(Richness_g ~ Sample_type + EPG_DC_dry, data = metadata_rarefied, family = poisson())
summary(g_poisson) 
g_emmeans <- emmeans(g_poisson, pairwise ~ Sample_type)
summary(g_emmeans)

richness_asvxepgxmethod <- ggplot(data = metadata_rarefied, mapping = aes(x = EPG_DC_dry, y = Richness_asv, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Number of Unique ASVs by Eggs Per Gram Grouped by Method", x = 'Eggs Per Gram')+
  theme_bw()
richness_asvxepgxmethod

richness_spxepgxmethod <- ggplot(data = metadata_rarefied, mapping = aes(x = EPG_DC_dry, y = Richness_sp, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Number of Unique ASVs by Eggs Per Gram Grouped by Method", x = 'Eggs Per Gram')+
  theme_bw()
richness_spxepgxmethod

richness_gxepgxmethod <- ggplot(data = metadata_rarefied, mapping = aes(x = EPG_DC_dry, y = Richness_g, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Number of Unique ASVs by Eggs Per Gram Grouped by Method", x = 'Eggs Per Gram')+
  theme_bw()
richness_gxepgxmethod

### Richness ~ sample_type * EPG vs Richness ~ sample_type + EPG
## no significant difference 


# testing model significance asv
asvxepgxmethodmdl <- aov(Richness_asv ~ Sample_type + EPG_DC_dry, data = metadata_rarefied)
summary(asvxepgxmethodmdl)

asvxmethodmdl <- aov(Richness_asv ~ Sample_type, data = metadata_rarefied)
summary(asvxmethodmdl)

anova(asvxmethodmdl, asvxepgxmethodmdl) # significant difference p = 0.0001887

# testing model significance sp
spxepgxmethodmdl <- aov(Richness_sp ~ Sample_type + EPG_DC_dry, data = metadata_rarefied)
summary(spxepgxmethodmdl)

spxmethodmdl <- aov(Richness_sp ~ Sample_type, data = metadata_rarefied)
summary(spxmethodmdl)

anova(spxmethodmdl, spxepgxmethodmdl) # significant difference p = 0.001988

# testing model significance g
gxepgxmethodmdl <- aov(Richness_g ~ Sample_type + EPG_DC_dry, data = metadata_rarefied)
summary(gxepgxmethodmdl)

gxmethodmdl <- aov(Richness_g ~ Sample_type, data = metadata_rarefied)
summary(gxmethodmdl)

anova(gxmethodmdl, gxepgxmethodmdl) # significant difference p = 0.00524





# alllllll evenness Pielou
evennessxmethod_asv <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Evenness_asv)) 
evennessxmethod_asv + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Evenness by Processing Method",
       x = "Processing Method", y = "Richness (# of ASVs)")
evennessxmethod_sp <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Evenness_sp)) 
evennessxmethod_sp + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Evenness by Processing Method",
       x = "Processing Method", y = "Richness (# of species)")
evennessxmethod_g <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Evenness_g)) 
evennessxmethod_g + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Evenness by Processing Method",
       x = "Processing Method", y = "Richness (# of genera)")


df_long <- pivot_longer(metadata_rarefied, cols = c(Richness_asv, Richness_sp, Richness_g,
                                                    Evenness_asv, Evenness_sp, Evenness_g),
                        names_to = "Diversity Metric", values_to = "Diversity Values")
df_long <- separate(df_long, 'Diversity Metric', into = c("Diversity_metric", "Taxonomy_level"), sep = '_')

df_long$Taxonomy_level <- factor(df_long$Taxonomy_level, levels = c('asv', 'sp', 'g'))

## facet grid box plot richness and evenness
ggplot(df_long, mapping = aes(x= Sample_type, y = `Diversity Values`)) +
  geom_boxplot()+
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5, fill = 'red')+
  facet_grid(Diversity_metric ~ Taxonomy_level, scales = 'free_y')+
  theme_bw()





# redoing the above but for raw reads

# combine asvs or multiples to on group of each level 
# with default settings NA remove 
ps_final_asv <- ps_final
ps_final_sp = tax_glom(ps_final, "species")
ps_final_g = tax_glom(ps_final, "genus")

# plot each level
plot_bar(ps_final_asv, x = "Sample_name")
plot_bar(ps_final_sp, x="Sample_name", fill="species")
plot_bar(ps_final_g, x="Sample_name", fill="genus")

# get richness stats
asv_richness <- estimate_richness(ps_final_asv)
sp_richness <- estimate_richness(ps_final_sp)
g_richness <- estimate_richness(ps_final_g)

# create new data frame with richness added
## look at doing this a different way
# assign sample(phyloseq) to meta data
metadata_raw_reads <- metadata[-c(3,36,37),] 
metadata_raw_reads$Richness_asv <- asv_richness[,"Observed"]
metadata_raw_reads$Richness_sp <- sp_richness[,"Observed"]
metadata_raw_reads$Richness_g <- g_richness[,"Observed"]

metadata_raw_reads$Shannon_asv <- asv_richness[,"Shannon"]
metadata_raw_reads$Shannon_sp <- sp_richness[,"Shannon"]
metadata_raw_reads$Shannon_g <- g_richness[,"Shannon"]

metadata_raw_reads$Evenness_asv <- metadata_raw_reads$Shannon_asv / log(metadata_raw_reads$Richness_asv)
metadata_raw_reads$Evenness_sp <- metadata_raw_reads$Shannon_sp / log(metadata_raw_reads$Richness_sp)
metadata_raw_reads$Evenness_g <- metadata_raw_reads$Shannon_g / log(metadata_raw_reads$Richness_g)
## Evenness is on scale of 0 to 1

mean(metadata_raw_reads$Richness_asv) # 29
mean(metadata_raw_reads$Richness_sp) # 7.22
mean(metadata_raw_reads$Richness_g) # 4.92

sd(metadata_raw_reads$Richness_asv) # 16.58
sd(metadata_raw_reads$Richness_sp) # 3.53 
sd(metadata_raw_reads$Richness_g) # 1.83
#confused on how these numbers are exactly the same


# start plots
# allllll richness
ggplot(metadata_rarefied, aes(x = Sample_type, y = Richness_asv, color = Host)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  labs(title = "Richness (ASVs) vs. Method")

richnessxmethod_asv <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Richness_asv)) 
richnessxmethod_asv + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness (# of ASVs)")
richnessxmethod_sp <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Richness_sp)) 
richnessxmethod_sp + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness (# of species)")
richnessxmethod_g <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Richness_g)) 
richnessxmethod_g + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Richness by Processing Method",
       x = "Processing Method", y = "Richness (# of genera)")


asv_poisson <- glm(Richness_asv ~ Sample_type + EPG_DC_dry, data = metadata_rarefied, family = poisson())
summary(asv_poisson)
asv_emmeans <- emmeans(asv_poisson, pairwise ~ Sample_type)
summary(asv_emmeans)
sp_poisson <- glm(Richness_sp ~ Sample_type + EPG_DC_dry, data = metadata_rarefied, family = poisson())
summary(sp_poisson)
sp_emmeans <- emmeans(sp_poisson, pairwise ~ Sample_type)
summary(sp_emmeans)
g_poisson <- glm(Richness_g ~ Sample_type + EPG_DC_dry, data = metadata_rarefied, family = poisson())
summary(g_poisson) 
g_emmeans <- emmeans(g_poisson, pairwise ~ Sample_type)
summary(g_emmeans)

# alllllll evenness Pielou
evennessxmethod_asv <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Evenness_asv)) 
evennessxmethod_asv + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Evenness by Processing Method",
       x = "Processing Method", y = "Evenness of ASVs")
evennessxmethod_sp <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Evenness_sp)) 
evennessxmethod_sp + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Evenness by Processing Method",
       x = "Processing Method", y = "Evenness of species")
evennessxmethod_g <- ggplot(data = metadata_rarefied, mapping = aes(x = Sample_type, y = Evenness_g)) 
evennessxmethod_g + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Evenness by Processing Method",
       x = "Processing Method", y = "Evenness of genera")

summary(betareg(Evenness_sp ~ Sample_type,data = metadata_rarefied))



df_long_raw <- pivot_longer(metadata_raw_reads, cols = c(Richness_asv, Richness_sp, Richness_g,
                                                    Evenness_asv, Evenness_sp, Evenness_g),
                        names_to = "Diversity Metric", values_to = "Diversity Values")
df_long_raw <- separate(df_long_raw, 'Diversity Metric', into = c("Diversity_metric", "Taxonomy_level"), sep = '_')

df_long_raw$Taxonomy_level <- factor(df_long_raw$Taxonomy_level, levels = c('asv', 'sp', 'g'))

## facet grid box plot richness and evenness
ggplot(df_long_raw, mapping = aes(x= Sample_type, y = `Diversity Values`)) +
  geom_boxplot()+
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5, fill = 'red')+
  facet_grid(Diversity_metric ~ Taxonomy_level, scales = 'free_y', labeller = labeller(Taxonomy_level = c("asv" = "ASV", "sp" = "Species", "g" = "Genus")))+
  labs(title = 'Richness and Evenness by Taxonomy Level and Method', x = "Method", y= 'Diversity Values')+
  theme_bw()


# new concepts
# methods overlap or have differing groups of detectable species
# plots of percentage of samples that have specific genus, species, or asv
# ^ all methods would be combined into one to have just host

## trying something
ps_sp <- ps_final_sp
otu_table(ps_sp) <- otu_table(ps_sp) > 0
otu_table(ps_sp) <- otu_table(ps_sp) * 1

ps_method_sp <- merge_samples(ps_sp, 'Sample_type')


df_sp <- as.data.frame(otu_table(ps_method_sp))

df_merged_sp <- as.data.frame(otu_table(ps_merged_sp))

df_merged_sp['Host_total',] <- colSums(df_merged_sp)

df_sp['Host_total',] <- df_merged_sp['Host_total',]

df_sp_long_perc <- c()

for (i in 1:length(df_sp)){
  fecal_percent <- df_sp['Fecal',i]/df_sp['Host_total',i]
  larva_percent <- df_sp['Larva',i]/df_sp['Host_total',i]
  swab_percent <- df_sp['Swab',i]/df_sp['Host_total',i]
  df_sp_long_perc$Method <- c(df_sp_long_perc$Method, 'fecal', 'larva', 'swab')
  df_sp_long_perc$Host_perc <- c(df_sp_long_perc$Host_perc, fecal_percent, larva_percent, swab_percent)
}
asvs <- NULL
for (asv in colnames(df_sp)) {asvs <- c(asvs, rep(asv,3))}
df_sp_long_perc$asv <- asvs
df_sp_long_perc <- as.data.frame(df_sp_long_perc)

df_sp_long_perc$Method <- factor(df_sp_long_perc$Method, c("fecal", "swab", "larva"))

detectionxmethod_sp <- ggplot(data = df_sp_long_perc, mapping = aes(x = Method, y = Host_perc)) 
detectionxmethod_sp + geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red')+
  labs(title = "Detection by Processing Method",
       x = "Processing Method", y = "Percentage of hosts species detected in ")

detectionxmethod_sp <- ggplot(data = df_sp_long_perc, mapping = aes(x = Method, y = Host_perc, color = asv))  

detectionxmethod_sp +  
  geom_boxplot() + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5, fill = 'red') + 
  geom_line(aes(group = asv), alpha = 0.7) +  # Slight transparency for readability 
  labs(title = "Detection by Processing Method", 
       x = "Processing Method", y = "Percentage of hosts species detected") + 
  scale_color_manual(values = rainbow(length(unique(df_sp_long_perc$asv)))) # Assigns a different color to each asv 


perc_stats <- betareg(Host_perc ~ Method,data = df_sp_long_perc)
#glmmTMBmodel <- glmmTMB(Host_perc ~ Method + (1|asv), family = beta_family(), data = df_sp_long_perc)

summary(perc_stats)

emmeans(perc_stats, pairwise ~ Method, adjust = 'Tukey')
#glht(perc_stats)




# each final phyloseq glommed at different levels


ps_merged_asv <- merge_samples(ps_final_asv, 'Host')
otu_table(ps_merged_asv) <- otu_table(ps_merged_asv) > 0
otu_table(ps_merged_asv) <- otu_table(ps_merged_asv) * 1


ps_merged_sp <- merge_samples(ps_final_sp, 'Host')
otu_table(ps_merged_sp) <- otu_table(ps_merged_sp) > 0
otu_table(ps_merged_sp) <- otu_table(ps_merged_sp) * 1


ps_merged_g <- merge_samples(ps_final_g, 'Host')
otu_table(ps_merged_g) <- otu_table(ps_merged_g) > 0
otu_table(ps_merged_g) <- otu_table(ps_merged_g) * 1

df_asv <- as.data.frame(otu_table(ps_merged_asv))
df_asv$Host <- rownames(df_asv)

df_asv_long <- pivot_longer(df_asv, cols = -Host, names_to = 'ASV', values_to = 'Presence')
df_asv_presence <- df_asv_long %>%
  group_by(ASV) %>%
  mutate(Host_Count = sum(Presence))
df_asv_presence <- unique(df_asv_presence[,c(-1,-3)])

ggplot(df_asv_presence, aes(x = reorder(ASV, -Host_Count), y = Host_Count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



df_sp <- as.data.frame(otu_table(ps_merged_sp))
df_sp$Host <- rownames(df_sp)

df_sp_long <- pivot_longer(df_sp, cols = -Host, names_to = 'ASV', values_to = 'Presence')
df_sp_presence <- df_sp_long %>%
  group_by(ASV) %>%
  mutate(Host_Count = sum(Presence))
df_sp_presence <- unique(df_sp_presence[,c(-1,-3)])

tax_table_sp <- as.data.frame(tax_table(ps_final_sp))%>%
  rownames_to_column("ASV")
df_sp_presence <- as.data.frame(df_sp_presence)%>%
  right_join(tax_table_sp, by = 'ASV')

ggplot(df_sp_presence, aes(x = reorder(species, -Host_Count), y = Host_Count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(df_sp_presence, aes(x = reorder(species, -Host_Count), y = (Host_Count/13)*100, fill = genus)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(. ~genus, scales = "free", space = 'free') +
  labs(title = "Percentage of Samples Species is Present In", x = 'Species', y = 'Percentage of Hosts')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 10), 
        legend.position = 'None',  
        strip.text.x = element_text(angle = 55, hjust = 0.5, size = 12),
        strip.background = element_blank())

ggplot(df_sp_presence, aes(x = (Host_Count/13)*100, y = reorder(species, -Host_Count), fill = genus)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(genus ~ ., scales = "free", space='free') +
  labs(title = "Percentage of Samples Species is Present In", x = 'Percentage of Hosts', y = 'Species') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 11),
        legend.position = "none",
        strip.text.y = element_text(angle = 0, hjust = 0.5, size = 12))


ggplot(df_sp_presence, aes(x = reorder(species, -Host_Count), y = (Host_Count/13)*100, fill = genus)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(. ~ genus, scales = "free_x", space = "free") +
  labs(title = "Percentage of Samples Species is Present In", x = 'Species', y = 'Percentage of Hosts') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        strip.text.x = element_text(angle = 45, hjust = 1, size = 12, vjust = 1))


df_g <- as.data.frame(otu_table(ps_merged_g))
df_g$Host <- rownames(df_g)

df_g_long <- pivot_longer(df_g, cols = -Host, names_to = 'ASV', values_to = 'Presence')
df_g_presence <- df_g_long %>%
  group_by(ASV) %>%
  mutate(Host_Count = sum(Presence))
df_g_presence <- unique(df_g_presence[,c(-1,-3)])

tax_table_g <- as.data.frame(tax_table(ps_final_g))%>%
  rownames_to_column("ASV")
df_g_presence <- as.data.frame(df_g_presence)%>%
  right_join(tax_table_g, by = 'ASV')

ggplot(df_g_presence, aes(x = reorder(genus, -Host_Count), y = Host_Count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# trying to make venn diagram
ps_binary_sp <- ps_final_sp
otu_table(ps_binary_sp) <- otu_table(ps_binary_sp) > 0
otu_table(ps_binary_sp) <- otu_table(ps_binary_sp) * 1

df_binary_sp <- as.data.frame(otu_table(ps_binary_sp))
df_binary_sp$Host <- sample_data(ps_binary_sp)$Host
df_binary_sp$Method <- sample_data(ps_binary_sp)$Sample_type


df_binary_sp <- pivot_longer(df_binary_sp, cols = -c(Host, Method), names_to = 'ASV', values_to = 'Presence')

tax_table_sp <- as.data.frame(tax_table(ps_binary_sp))%>%
  rownames_to_column("ASV")
df_binary_sp <-left_join(df_binary_sp, tax_table_sp, by = 'ASV')

method_host_list <- split(df_binary_sp, list(df_binary_sp$Method, df_binary_sp$Host))
species_lists <- lapply(method_host_list, function(x) unique(x$species))

venn.diagram(species_lists,
             filename = NULL, # To display in R plot window
             imagetype = "png",
             fill = c("blue", "red", "green"),
             alpha = 0.5,
             cat.cex = 1.5,
             cex = 1.2)
## venn diagram is not working

## stats on if species detected by method for each horse

ps_presence_sp <- ps_final_sp
otu_table(ps_presence_sp) <- otu_table(ps_presence_sp) > 0
otu_table(ps_presence_sp) <- otu_table(ps_presence_sp) * 1


df_epg_host <- unique(metadata_rarefied[,c('Host','EPG_DC_dry')])


df_presence_sp <- as.data.frame(otu_table(ps_presence_sp))
df_presence_sp$Host <- sample_data(ps_presence_sp)$Host
df_presence_sp$Method <- sample_data(ps_presence_sp)$Sample_type

df_presence_sp_long <- pivot_longer(df_presence_sp, cols = -c(Host,Method), names_to = 'ASV', values_to = 'Detection')

df_long <- df_presence_sp %>%
  pivot_longer(cols = -c(Host, Method), names_to = "Species", values_to = "Detection") %>%
  arrange(Species, Method)

tax_table_sp <- as.data.frame(tax_table(ps_final_sp))%>%
  rownames_to_column("ASV")
tax_table_sp <- tax_table_sp[c(-2:-8)]
df_presence_sp_long <- as.data.frame(df_presence_sp_long)%>%
  left_join(tax_table_sp, by = 'ASV')

df_presence_sp_long <- rename(df_presence_sp_long, 'Species'= 'species')

df_presence_sp_long <- df_presence_sp_long %>%
  arrange(Species, Method)
present_species_vector <- c()

for (species in species_vector){
  df_presence_per_sp <- df_presence_sp_long[df_presence_sp_long$Species == species, ]
  df_presence_per_sp <- df_presence_per_sp[c(-2)]
  df_presence_per_sp <- distinct(df_presence_per_sp)
  num_of_hosts <- sum(df_presence_per_sp$Detection)
  print(num_of_hosts)
  if (num_of_hosts > 3){
    present_species_vector <- c(present_species_vector, species)
  }
}


df_presence_sp_long <- as.data.frame(df_presence_sp_long)%>%
  left_join(df_epg_host, by= 'Host')
df_presence_sp_long <- rename(df_presence_sp_long, 'EPG'='EPG_DC_dry')


df_Strongylus_equinus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == "Strongylus_equinus", ])
df_Cyathostomum_catinatum <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cyathostomum_catinatum',])
df_Cylicocyclus_elongatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicocyclus_elongatus',])
df_Triodontophorus_brevicauda <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Triodontophorus_brevicauda',])
df_Cylicostephanus_longibursatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicostephanus_longibursatus',])
df_Cylicocyclus_nassatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicocyclus_nassatus',])
df_Coronocyclus_coronatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Coronocyclus_coronatus',])
df_Cylicostephanus_goldi <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicostephanus_goldi',])
df_Petrovinema_poculatum <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Petrovinema_poculatum',])
df_Poteriostomum_imparidentatum <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Poteriostomum_imparidentatum',])
df_Strongylus_edentatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Strongylus_edentatus',])
df_Coronocyclus_labratus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Coronocyclus_labratus',])
df_Triodontophorus_serratus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Triodontophorus_serratus',])
df_Cylicostephanus_minutus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicostephanus_minutus',])
df_Poteriostomum_ratzii <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Poteriostomum_ratzii',])
df_Cylicocyclus_ultrajectinus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicocyclus_ultrajectinus',])
df_Cylicostephanus_calicatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicostephanus_calicatus',])
df_Strongylus_vulgaris <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Strongylus_vulgaris',])
df_Gyalocephalus_capitatus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Gyalocephalus_capitatus',])
df_Cylicocyclus_insigne <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == 'Cylicocyclus_insigne',])
df_Parapoteriostomum_mettami <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == "Parapoteriostomum_mettami", ])
df_Parapoteriostomum_euproctus <- data.frame(df_presence_sp_long[df_presence_sp_long$Species == "Parapoteriostomum_euproctus", ])


species_df_vector <- list(df_Strongylus_equinus, df_Cyathostomum_catinatum, 
                       df_Cylicocyclus_elongatus, df_Triodontophorus_brevicauda,
                       df_Cylicostephanus_longibursatus, df_Cylicocyclus_nassatus,
                       df_Coronocyclus_coronatus, df_Cylicostephanus_goldi,
                       df_Petrovinema_poculatum, df_Poteriostomum_imparidentatum, 
                       df_Strongylus_edentatus, df_Coronocyclus_labratus, 
                       df_Triodontophorus_serratus, df_Cylicostephanus_minutus,
                       df_Poteriostomum_ratzii, df_Cylicocyclus_ultrajectinus, 
                       df_Cylicostephanus_calicatus, df_Strongylus_vulgaris,
                       df_Gyalocephalus_capitatus, df_Cylicocyclus_insigne,
                       df_Parapoteriostomum_mettami,df_Parapoteriostomum_euproctus)

post_hoc_regression_species <- list()
for (species in 1:length(species_df_vector)){
  df <- species_df_vector[[species]]
  regression_model <- glm(Detection ~ Method, data = df, family = binomial())
  post_hoc <- summary(emmeans(regression_model, pairwise ~ Method))
  post_hoc_regression_species[[species]] <- post_hoc
}

for (species in 1:length(species_df_vector)){
  print(species_df_vector[[species]][1,'Species'])
  print(post_hoc_regression_species[[species]])
}



#for strongyluys_equinus, no difference in models, not significant difference by method
interaction <- glm(Detection ~ Method * EPG, data = df_Strongylus_equinus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Strongylus_equinus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Strongylus_equinus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Strongylus_equinus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))

#for Cyathostomum_catinatum, no interaction but epg is significant in model, but no significant diffference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Cyathostomum_catinatum, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Cyathostomum_catinatum, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Cyathostomum_catinatum, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method + EPG, data = df_Cyathostomum_catinatum, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Cylicocyclus_elongatus, no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Cylicocyclus_elongatus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Cylicocyclus_elongatus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Cylicocyclus_elongatus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Cylicocyclus_elongatus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Triodontophorus_brevicauda, no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Triodontophorus_brevicauda, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Triodontophorus_brevicauda, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Triodontophorus_brevicauda, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Triodontophorus_brevicauda, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Cylicostephanus_longibursatus, no interaction but epg is significant in model, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Cylicostephanus_longibursatus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Cylicostephanus_longibursatus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Cylicostephanus_longibursatus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method + EPG, data = df_Cylicostephanus_longibursatus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Cylicocyclus_nassatus, no interaction but epg is significant in model
#### #Fecal and Larva (p=0.0376)
#### #Larva and Swab (p = 0.0459)

interaction <- glm(Detection ~ Method * EPG, data = df_Cylicocyclus_nassatus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Cylicocyclus_nassatus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Cylicocyclus_nassatus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method + EPG, data = df_Cylicocyclus_nassatus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))



#for Coronocyclus_coronatus, no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Coronocyclus_coronatus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Coronocyclus_coronatus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Coronocyclus_coronatus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Coronocyclus_coronatus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))



#for Cylicostephanus_goldi, no interaction but epg is significant in model, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Cylicostephanus_goldi, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Cylicostephanus_goldi, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Cylicostephanus_goldi, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method + EPG, data = df_Cylicostephanus_goldi, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Petrovinema_poculatum,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Petrovinema_poculatum, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Petrovinema_poculatum, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Petrovinema_poculatum, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Petrovinema_poculatum, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))



#for Poteriostomum_imparidentatum,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Poteriostomum_imparidentatum, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Poteriostomum_imparidentatum, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Poteriostomum_imparidentatum, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Poteriostomum_imparidentatum, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Strongylus_edentatus,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Strongylus_edentatus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Strongylus_edentatus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Strongylus_edentatus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Strongylus_edentatus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Coronocyclus_labratus,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Coronocyclus_labratus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Coronocyclus_labratus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Coronocyclus_labratus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Coronocyclus_labratus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))


#for Triodontophorus_serratus,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Triodontophorus_serratus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Triodontophorus_serratus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Triodontophorus_serratus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method, data = df_Triodontophorus_serratus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))



#for Cylicostephanus_minutus,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Cylicostephanus_minutus, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Cylicostephanus_minutus, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Cylicostephanus_minutus, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method + EPG, data = df_Cylicostephanus_minutus, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))



#for Poteriostomum_ratzii,  no difference in models, not significant difference by method

interaction <- glm(Detection ~ Method * EPG, data = df_Poteriostomum_ratzii, family = binomial())
epg_model <- glm(Detection ~ Method + EPG, data = df_Poteriostomum_ratzii, family = binomial())
no_epg <- glm(Detection ~ Method, data = df_Poteriostomum_ratzii, family = binomial())

anova(interaction, epg_model)
anova(epg_model, no_epg)

summary(epg_interact_model <- glm(Detection ~ Method * EPG, data = df_Poteriostomum_ratzii, family = binomial()))
summary(emmeans(epg_interact_model, pairwise ~ Method))




## ignore the regression for now
# try relative abundance on the rarefied samples for reads
# use poisson model
# check the fit first with a few high abundance species or genus
# compare the abundance by method type
# don't worry about NAs since its part of relative abundance of ~3,000


ps_rarefied <- read_rds("ps_rarefied_16Oct24.rds")
ps_rarefied_sp <- read_rds("ps_rarefied_sp_16Oct24.rds")
ps_rarefied_g <- read_rds("ps_rarefied_g_16Oct24.rds")

plot_bar(ps_rarefied, x="Sample_name", fill = 'species')
plot_bar(ps_rarefied_sp, x="Sample_name", fill = 'species')



# rarefied to 3,413 reads
asv_rarefied <- sample_data(ps_rarefied)
asv_rarefied$Reads <- sample_sums(ps_rarefied)

# choose a species in particular
# create a data frame with sample name, sample type, and number of reads for that species
# then run poisson on count data technically relative abundance cause all sampled to same level



df_rel_abund_sp <- as.data.frame(otu_table(ps_rarefied_sp))
df_rel_abund_sp$Host <- sample_data(ps_rarefied_sp)$Host
df_rel_abund_sp$Method <- sample_data(ps_rarefied_sp)$Sample_type
df_rel_abund_sp$EPG <- metadata_rarefied$EPG_DC_dry

# just to try with one species
df_rel_abund_sp <- df_rel_abund_sp[c(1, 23, 24,25)]
df_rel_abund_sp_long <- pivot_longer(df_rel_abund_sp, cols = -c(Host,Method, EPG), names_to = 'ASV', values_to = 'Abundance')

tax_table_sp <- as.data.frame(tax_table(ps_rarefied_sp))%>%
  rownames_to_column("ASV")
tax_table_sp <- tax_table_sp[c(-2:-8)]
df_rel_abund_sp_long <- as.data.frame(df_rel_abund_sp_long)%>%
  left_join(tax_table_sp, by = 'ASV')

df_rel_abund_sp_long <- rename(df_rel_abund_sp_long, 'Species'= 'species')

rel_abund_model <- glm(Abundance ~ Method, data = df_rel_abund_sp_long, family = poisson())
summary(rel_abund_model)
emmeans(rel_abund_model, pairwise ~ Method)



# trying all now
df_rel_abund_sp <- as.data.frame(otu_table(ps_rarefied_sp))
df_rel_abund_sp$Host <- sample_data(ps_rarefied_sp)$Host
df_rel_abund_sp$Method <- sample_data(ps_rarefied_sp)$Sample_type
df_rel_abund_sp$EPG <- metadata_rarefied$EPG_DC_dry

df_rel_abund_sp_long <- pivot_longer(df_rel_abund_sp, cols = -c(Host,Method,EPG), names_to = 'ASV', values_to = 'Abundance')

tax_table_sp <- as.data.frame(tax_table(ps_rarefied_sp))%>%
  rownames_to_column("ASV")
tax_table_sp <- tax_table_sp[c(-2:-8)]
df_rel_abund_sp_long <- as.data.frame(df_rel_abund_sp_long)%>%
  left_join(tax_table_sp, by = 'ASV')

df_rel_abund_sp_long <- rename(df_rel_abund_sp_long, 'Species'= 'species')

species_vector <- tax_table_sp[,'species']

for (species in species_vector){
  df_rel_abund_per_sp <- df_rel_abund_sp_long[df_rel_abund_sp_long$Species == species, ]
  print(species)
  print(summary(glm(Abundance ~ Method, data = df_rel_abund_per_sp, family = poisson())))
}

## poisson was really bad model, trying negative binomial instead


for (species in species_vector){
  df_rel_abund_per_sp <- df_rel_abund_sp_long[df_rel_abund_sp_long$Species == species, ]
  print(species)
  print(summary(glmer.nb(Abundance ~ Method + (1|Host), data = df_rel_abund_per_sp)))
}


## see if epg interaction effect

df_strong_equin_rel_abund <- df_rel_abund_sp_long[df_rel_abund_sp_long$Species == 'Strongylus_equinus',]

df_strong_equin_rel_abund <- as.data.frame(df_strong_equin_rel_abund)%>%
  left_join(asv_rarefied, by = 'Host')

nb_model <- glmer.nb(Abundance ~ Method * EPG + (1|Host), data = df_strong_equin_rel_abund)
summary(nb_model)

nb_model_1 <- glmer.nb(Abundance ~ Method + (1|Host), data = df_strong_equin_rel_abund)
summary(nb_model_1)


anova(nb_model, nb_model_1)


# try interaction effect with EPG -> Abundance ~ Method * EPG + (1|Host)

species<-present_species_vector[2]

for (species in present_species_vector){
  df_rel_abund_per_sp <- df_rel_abund_sp_long[df_rel_abund_sp_long$Species == species, ]
  print(species)
  nb_model <- glmer.nb(Abundance ~  (1|Host), data = df_rel_abund_per_sp)
  
  nb_model_1 <- glmer.nb(Abundance ~ Method + (1|Host), data = df_rel_abund_per_sp)
  
  anova_result <- anova(nb_model, nb_model_1)
  p_value <- anova_result$`Pr(>Chisq)`[2]
  print(p_value)
  
  if (p_value < 0.05) {
    print(summary(nb_model))
    print(summary(nb_model_1))
    print(anova_result)
  }
  else {
    print('Not signficant difference between models')
  }
}


fecal_df <- metadata_final[metadata_final[,'Sample_type'] == 'Fecal', ]
swab_df <- metadata_final[metadata_final[,'Sample_type'] == 'Swab', ]
larva_df <- metadata_final[metadata_final[,'Sample_type'] == 'Larva', ]


# mean and sd of wet/dry weights

dry_weight_df <- unique(metadata_rarefied[c('Host', 'Wet_weight', 'Dry_weight', 'Percent_dry')])

mean(dry_weight_df[,'Wet_weight']) # mean 7.675
sd(dry_weight_df[,'Wet_weight'])  # sd 1.185

mean(dry_weight_df[,'Dry_weight']) # mean 2.2
sd(dry_weight_df[,'Dry_weight'])  # sd 0.524

mean(1-dry_weight_df[,'Percent_dry']) # mean 0.7137
sd(1-dry_weight_df[,'Percent_dry'])  # sd 0.054

max(1-dry_weight_df[,'Percent_dry']) # max 79.88 %
min(1-dry_weight_df[,'Percent_dry']) # min 0.5926

# mean, sd and range of epg

epg_df <- unique(metadata_rarefied[c('Host', 'EPG_McM_wet', 'EPG_DC_wet')])

min(epg_df[,'EPG_McM_wet']) # min 0
max(epg_df[, 'EPG_McM_wet']) #max 200

mean(epg_df[,'EPG_McM_wet']) # mean 57.84615
sd(epg_df[,'EPG_McM_wet']) # sd 62.75063

min(epg_df[,'EPG_DC_wet']) # min  0.6472
max(epg_df[, 'EPG_DC_wet']) #max  107.667

mean(epg_df[,'EPG_DC_wet']) # mean  29.896
sd(epg_df[,'EPG_DC_wet']) # sd  33.63558


# mean, sd and range of reads by method

min(fecal_df[,'Reads']) # min 0
max(fecal_df[, 'Reads']) # max 20914

mean(fecal_df[, 'Reads']) # mean 10795.54
sd(fecal_df[, 'Reads']) # sd 6656.37

min(swab_df[,'Reads']) # min 201
max(swab_df[, 'Reads']) # max 26864

mean(swab_df[, 'Reads']) # mean 16914.31
sd(swab_df[, 'Reads']) # sd 8447.56

min(larva_df[,'Reads']) # min 12168
max(larva_df[, 'Reads']) # max 23419

mean(larva_df[, 'Reads']) # mean 18572.69
sd(larva_df[, 'Reads']) # sd 3730.49


# mean and sd of richness by method

fecal_df <- metadata_rarefied[metadata_rarefied[,'Sample_type'] == 'Fecal', ]
swab_df <- metadata_rarefied[metadata_rarefied[,'Sample_type'] == 'Swab', ]
larva_df <- metadata_rarefied[metadata_rarefied[,'Sample_type'] == 'Larva', ]

mean(fecal_df[,'Richness_asv'])  # fecal ASVs  15.5833
sd(fecal_df[,'Richness_asv'])    # sd 10.526

mean(fecal_df[,'Richness_sp'])   # fecal sp  3.667
sd(fecal_df[,'Richness_sp'])     # sd     1.923  

mean(fecal_df[,'Richness_g'])    # fecal g 3.0833
sd(fecal_df[,'Richness_g'])      # sd  1.311


mean(swab_df[,'Richness_asv'])   # swab asv  31.636
sd(swab_df[,'Richness_asv'])     # sd     14.059

mean(swab_df[,'Richness_sp'])    # swab sp   7.636
sd(swab_df[,'Richness_sp'])      # sd     1.963

mean(swab_df[,'Richness_g'])     # swab g 5.454
sd(swab_df[,'Richness_g'])       # sd     1.214


mean(larva_df[,'Richness_asv'])  # larva asv  39.154
sd(larva_df[,'Richness_asv'])    # sd    15.410

mean(larva_df[,'Richness_sp'])   # larva sp  10.154
sd(larva_df[,'Richness_sp'])     # sd    2.764

mean(larva_df[,'Richness_g'])    # larva g   6.154
sd(larva_df[,'Richness_g'])      # sd   1.281



fecal_df <- df_sp_long_perc[df_sp_long_perc[,'Method'] == 'fecal', ]
swab_df <- df_sp_long_perc[df_sp_long_perc[,'Method'] == 'swab', ]
larva_df <- df_sp_long_perc[df_sp_long_perc[,'Method'] == 'larva', ]


min(fecal_df[,'Host_perc']) # min 0
max(fecal_df[, 'Host_perc']) # max 100

mean(fecal_df[, 'Host_perc']) # mean 22.85
sd(fecal_df[, 'Host_perc']) # sd 26.87

min(swab_df[,'Host_perc']) # min 0
max(swab_df[, 'Host_perc']) # max 100

mean(swab_df[, 'Host_perc']) # mean 57.08
sd(swab_df[, 'Host_perc']) # sd 26.67

min(larva_df[,'Host_perc']) # min 60.00
max(larva_df[, 'Host_perc']) # max 100

mean(larva_df[, 'Host_perc']) # mean 88.75
sd(larva_df[, 'Host_perc']) # sd 13.24





# sample richness by host 
ps_host_asv <- merge_samples(ps_rarified, 'Host')
otu_table(ps_host_asv) <- otu_table(ps_host_asv) > 0
otu_table(ps_host_asv) <- otu_table(ps_host_asv) * 1
asvxhost<- data.frame(sample_sums(ps_host_asv))

ps_host_sp <- merge_samples(ps_rarified_sp, 'Host')
otu_table(ps_host_sp) <- otu_table(ps_host_sp) > 0
otu_table(ps_host_sp) <- otu_table(ps_host_sp) * 1
spxhost <- data.frame(sample_sums(ps_host_sp))

ps_host_g <- merge_samples(ps_rarified_g, 'Host')
otu_table(ps_host_g) <- otu_table(ps_host_g) > 0
otu_table(ps_host_g) <- otu_table(ps_host_g) * 1
gxhost <- data.frame(sample_sums(ps_host_g))


mean(asvxhost[,1]) # asv by host mean 58.30769
sd(asvxhost[,1]) # asv by host sd 22.01282

mean(spxhost[,1]) # sp by host mean 11.38462
sd(spxhost[,1]) # sp by host sd 2.56705

mean(gxhost[,1]) # g by host mean 6.769231
sd(gxhost[,1]) # g by host sd 1.165751


tax_table_g <- as.data.frame(tax_table(ps_host_g))%>%
  rownames_to_column("ASV")
df_g_presence <- as.data.frame(otu_table(ps_host_g))
colnames(df_g_presence) <- c(tax_table_g[,'genus'])



## for each level, compare by host, method, epg, possibly all at same time??
## not super needed
asv_richnessxmethodxhost <- ggplot(data = metadata_rarefied, mapping = aes(x = Host, y = Richness_asv, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Richness of Host Grouped by Method", y = "Richness (# of ASVs)")
asv_richnessxmethodxhost

sp_richnessxmethodxhost <- ggplot(data = metadata_rarefied, mapping = aes(x = Host, y = Richness_sp, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Richness of Host Grouped by Method", y = "Richness (# of Species)")
sp_richnessxmethodxhost

g_richnessxmethodxhost <- ggplot(data = metadata_rarefied, mapping = aes(x = Host, y = Richness_g, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Richness of Host Grouped by Method", y = "Richness (# of Genera)")
g_richnessxmethodxhost


ggplot(metadata_rarefied, aes(x = EPG_DC_dry, y = Richness_asv, color = Host)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Sample_type ~ .) +
  labs(title = "Richness (ASVs) vs. EPG by Host and Method")

ggplot(metadata_rarefied, aes(x = EPG_DC_dry, y = Richness_sp, color = Host)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Sample_type ~ .) +
  labs(title = "Richness (species) vs. EPG by Host and Method")



# have no idea what this is doing, how it works, and what its telling me

#delete
lmer(Richness_asv ~ EPG_DC_dry * Host + (1 | Sample_type), data = metadata_rarefied)



## TRYING PCoA & NMDS!!!!!
pcoa_ord <- ordinate(ps_rarified_sp, "PCoA", "bray")
plot_ordination(ps_rarified_sp, pcoa_ord, color = "Host", shape = "Sample_type", title = "PCoA Plot")

nmds_ord <- ordinate(ps_rarified_sp, "NMDS", "bray")
plot_ordination(ps_rarified_sp, nmds_ord, color = "Host", shape = "Sample_type", title = "NMDS Plot")





## rename/redo all this crap








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










### considering below this extra code and reference for supplemental material
# wet mcm and dc
dc_vs_mcm_epg <- ggplot(data = metadata, mapping = aes(x = EPG_DC_wet, y = EPG_McM_wet, label = Host)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_text(nudge_x = 5, nudge_y = 10, size = 3)+
  labs(title = "Double centrifugation EPG vs McMaster EPG", y = "McMaster EPG", x = "Double Centrifugation EPG" )
dc_vs_mcm_epg



ggplot(epg_data, aes(x = Method, y = EPG, color = Weight_type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "EPG vs. Method by Weight_type")

lm_epg <- lm(EPG ~ Method + Weight_type, data = epg_data)
summary(lm_epg)


t.test(metadata$EPG_McM_wet, metadata$EPG_DC_wet, paired = TRUE)
t.test(metadata$EPG_McM_dry, metadata$EPG_DC_dry, paired = TRUE)

summary(lm(EPG_McM_dry ~ EPG_DC_dry, data = metadata))
summary(lm(EPG_McM_wet ~ EPG_DC_wet, data = metadata))






## redo larval




## creating combined data set metadata and epg

combined_metadata <- inner_join(metadata, epg_data, by = "Host")

summary(lm(Avg_L3_count ~ Method + Weight_type + EPG, combined_metadata))

ggplot(combined_metadata, aes(x = EPG, y = Avg_L3_count)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Method ~ Weight_type, scales = "free") +
  labs(title = "Avg_L3_count vs. EPG by Method and Weight_type")


L3_count_vs_EPG <- ggplot(data = metadata, mapping = aes(x = Avg_L3_count, y = EPG_, label = Host)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_text(nudge_x = 0.3, nudge_y = 10, size = 3)+
  labs(title = "Average L3 count vs McMaster EPG", y = "McMaster EPG", x = "Avg L3 count" )
L3_count_vs_EPG



## larval counts vs egg counts
L3_count_vs_EPG <- ggplot(data = sample_data(psa), mapping = aes(x = L3_count_sum / Number_of_L3_counts, y = EPG, label = Host)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_text(nudge_x = 0.3, nudge_y = 10, size = 3)+
  labs(title = "Average L3 count vs McMaster EPG", y = "McMaster EPG", x = "Avg L3 count" )
L3_count_vs_EPG

L3_count_vs_dry_EPG <- ggplot(data = sample_data(psa), mapping = aes(x = L3_count_sum / Number_of_L3_counts, y = EPG * Percent_dry, label = Host)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_text(nudge_x = 0.3, nudge_y = 2, size = 3)+
  labs(title = "Average L3 count vs McMaster Dry EPG", y = "McMaster Dry EPG", x = "Avg L3 count" )
L3_count_vs_dry_EPG

L3_count_vs_dc_egg_count <- ggplot(data = sample_data(psa), mapping = aes(x = L3_count_sum / Number_of_L3_counts, y = DC_egg_count * Percent_dry, label = Host)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_text(nudge_x = 0.3, nudge_y = 2, size = 3)+
  labs(title = "Average L3 count vs Double Centrifugation Dry EPG", y = "McMaster Dry EPG", x = "Avg L3 count" )
L3_count_vs_dc_egg_count

## TODO: dry weights vs wet weights

dry_vs_wet <- ggplot(data = sample_data(psa), mapping = aes(x = Dry_weight, y = Wet_weight)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "Dry weight vs wet weight")

dry_vs_wet

dry_hist <- ggplot(data = sample_data(psa), mapping = aes(x = Dry_weight)) +
  geom_histogram(bins = 10)+
  labs(title = "Dry weight distribution")
dry_hist


dry_perc_hist <- ggplot(data = sample_data(psa), mapping = aes(x = Percent_dry)) +
  geom_histogram(bins = 15)+
  labs(title = "Percent dry weight distribution")
dry_perc_hist

dry_perc_vs_epg <- ggplot(data = sample_data(psa), mapping = aes(x = Percent_dry, y = EPG)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = "Percent Dry weight vs EPG")
dry_perc_vs_epg




## TODO: Qubit vs epg by methods

qubitxepgxmethod <- ggplot(data = metadata, mapping = aes(x = EPG_DC_dry, y = Qubit, color = Sample_type)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(title = "Qubit by Eggs Per Gram Grouped by Method")
qubitxepgxmethod

qubitxepgxmethod_aov <- aov(Qubit ~ EPG_DC_dry * Sample_type, data = metadata)
summary(qubitxepgxmethod_aov)
TukeyHSD(qubitxepgxmethod_aov)





# boxplot of reads by method
readsxmethod <- ggplot(data = metadata_final, mapping = aes(x = Sample_type, y = Reads)) 
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


