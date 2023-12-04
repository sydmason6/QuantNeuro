# Load libraries

```{r, warning=F, message=F}
library(tidyverse)
library(cowplot)
library(ggpubr)
library(phyloseq)
library(vegan)
library(ape)
library(DESeq2)
theme_set(theme_bw(base_size=15))
```

# Import Metadata

```{r}
meta.df <- read.table("/Users/sydmason/Library/CloudStorage/Box-Box/2023-05-15_cohousing_separation/qiime_output/sample_metadata.txt", sep="\t", header=T) 

meta.df <- meta.df %>%
  mutate(X.sample.id=str_c("X", sample.id)) %>%
  mutate(age = case_when(
  group %in% c("Y", "CY", "exCY") ~ "Y",
  group %in% c("O", "CO", "exCO") ~ "O")) %>%
  # create metadata that is only relevant to the current timepoint
  # (e.g. at baseline, exCY doesn't mean anything)
  mutate(current.group = case_when(
    timepoint == "baseline" & age == "Y" ~ "Y",
    timepoint == "baseline" & age == "O" ~ "O",
    timepoint == "1st day separation" & group == "exCY" ~ "CY",
    timepoint == "1st day separation" & group == "exCO" ~ "CO",
    TRUE ~ group))
  ## makes Y precede old in ggplot
  ## mutate(age=factor(age, levels=c("Y", "O"))) %>%
  
```

# Import Data

```{r}
## Delete the # before OTU.ID in the feature table
counts.df <- read.table("/Users/sydmason/separation/exported/feature-table.tsv", 
                        sep="\t", header=T)
counts.df <- column_to_rownames(counts.df, var="OTU.ID")
dim(counts.df)
```

# Import Taxonomy Table

```{r}
tax.df <- read.table("/Users/sydmason/separation/exported/taxonomy.tsv", sep="\t", header=T)
tax.df <- tax.df %>% 
  separate(Taxon, c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep="; [[:alpha:]]__", remove=T, fill="right")
```

# Filtration

## Visualize filtration thresholds

```{r}
p1 <- data.frame(sample.sum=colSums(counts.df)) %>%
  ggplot(aes(sample.sum)) +
  geom_histogram() +
  geom_vline(xintercept=2e4, lty=2) +
  labs(title="Sample Filtration")

p2 <- data.frame(frac.of.samples=rowSums(counts.df > 0)/ncol(counts.df)) %>%
  ggplot(aes(frac.of.samples)) +
  geom_histogram() +
  geom_vline(xintercept=0.1, lty=2) +
  labs(title="ASV filtration")

plot_grid(p1, p2, ncol=2)
```
```{r}
samples.to.keep <- colnames(counts.df)[colSums(counts.df) > 2e4] #may lower cut off if monocolonized or looking for rare bacteria then may be less reads
ASVs.to.keep <- rownames(counts.df)[rowSums(counts.df > 0)/ncol(counts.df) > 0.05]
counts.filt.df <- counts.df[ASVs.to.keep, samples.to.keep]
```

```{r}
dim(counts.df)
dim(counts.filt.df)
```

# Create phyloseq object

```{r}
physeq.meta.df <- meta.df %>% 
  filter(X.sample.id %in% colnames(counts.filt.df)) %>%
  column_to_rownames("X.sample.id")

physeq.data.df <- counts.filt.df

physeq.tax.df <- tax.df %>%
  column_to_rownames("Feature.ID")

physeq <- phyloseq(
  physeq.data.df %>% as.matrix %>% otu_table(taxa_are_rows=T),
  physeq.meta.df %>% sample_data,
  physeq.tax.df %>% as.matrix %>% tax_table
)
```

# PCoA Separation

pcoa.separation <- physeq %>% 
  
  # subset to FMT samples
  #subset_samples(experiment == "co-housing separation") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.separation.df <- plot_ordination(physeq, pcoa.separation, type="samples", justDF=T)
# write_csv(pcoa.FMT.df,"/Users/timcox/Library/CloudStorage/Box-Box/2021-07_karthi_and_tim_16S/qiime_output/FMT.csv")

# Perform PCoA including PC3 and PC4
pcoa.separation.pc34 <- ordinate(physeq, method = "MDS", distance = "bray", number = c(3, 4))

# Plot PC3 and PC4
pcoa.separation.df.pc34 <- plot_ordination(physeq, pcoa.separation.pc34, type = "samples", justDF=T)

#PCoA Young

pcoa.young <- physeq %>% 
  
  # subset to FMT samples
  subset_samples(age == "Y") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.young.df <- plot_ordination(physeq, pcoa.young, type="samples", justDF=T)


#PCoA Old

pcoa.old <- physeq %>% 
  
  # subset to FMT samples
  subset_samples(age == "O") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.old.df <- plot_ordination(physeq, pcoa.old, type="samples", justDF=T)

#PCoA age + timepoint

pcoa.age.timepoint <- physeq %>% 
  
  # subset to FMT samples
  #subset_samples(age == "O") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.age.timepoint.df <- plot_ordination(physeq, pcoa.age.timepoint, type="samples", justDF=T)

#PCoA ex-co

pcoa.ex.co <- physeq %>% 
  
  # subset to FMT samples
  subset_samples(group == "exCO" | group == "exCY") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.ex.co.df <- plot_ordination(physeq, pcoa.ex.co, type="samples", justDF=T)

#PCoA co-house

pcoa.cohouse <- physeq %>% 
  
# subset to FMT samples
subset_samples(group == "CO" | group == "CY") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.cohouse.df <- plot_ordination(physeq, pcoa.cohouse, type="samples", justDF=T)

#PCoA co-house + ex

pcoa.cohouse.ex <- physeq %>% 
  
  # subset to FMT samples
  subset_samples(group == "CY" | group == "exCY") %>%
  
  # compute relative abundance
  transform_sample_counts(function(OTU) OTU/sum(OTU)) %>%
  
  # perform PCoA using Bray-Curtis distance
  ordinate(method="MDS", distance="bray")

# create coordinates
pcoa.cohouse.ex.df <- plot_ordination(physeq, pcoa.cohouse.ex, type="samples", justDF=T)

```
#Plots
```{r}
pcoa.separation.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=timepoint)) +
  geom_point(size=3) +
  stat_ellipse()

pcoa.separation.df.pc34 %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=timepoint)) +
  geom_point(size=3) +
  labs(x = "PC3", y = "PC4")

pcoa.young.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=group, shape=timepoint)) +
  geom_point(size=3) +
  stat_ellipse()

pcoa.old.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=group, shape=timepoint)) +
  geom_point(size=4)

pcoa.age.timepoint.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=timepoint, shape=group)) +
  geom_point(size=3)

pcoa.ex.co.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=group, shape=timepoint)) +
  geom_point(size=3) +
  stat_ellipse()

pcoa.cohouse.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=group, shape=timepoint)) +
  geom_point(size=3)

pcoa.cohouse.ex.df %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=group, shape=timepoint)) +
  geom_point(size=3) +
  stat_ellipse()

```
## Contribution of variables

# Extract eigenvalues (percent contribution) for each dimension
eigenvalues <- pcoa.separation$values$Eigenvalues
percent_contribution <- eigenvalues / sum(eigenvalues) * 100

# Create a bar graph
bar_data <- data.frame(
  Dimension = seq_along(percent_contribution),
  Percent_Contribution = percent_contribution
)

# Create the bar plot
ggplot(bar_data, aes(x = factor(Dimension), y = Percent_Contribution)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(
    title = "Percent Contribution of Dimensions",
    x = "Dimension",
    y = "Percent Contribution"
  ) +
  scale_x_discrete(labels = paste("PC", 1:length(percent_contribution))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)) # Label the dimensions as PC1, PC2, ...

# Display the plot





