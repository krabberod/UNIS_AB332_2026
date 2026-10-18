# Taxonomic Analysis of Isfjorden Microbial Eukaryote Communities - PHYLOSEQ VERSION
# AB332 - Community Ecology Course

# Load required libraries
library(phyloseq)
library(tidyverse)
library(RColorBrewer)

# Load data from previous labs
load("AB332_lab_III.RData")


### READ THE TAXONOMY ####
# The taxonomy is from PR2 (version 5.1.1), taken from the best BLAST hit of each OTU.
# Id is the % identity to the best hit. The taxonomic ranks in PR2 are:
# Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species
taxonomy_full <- read_tsv("data/AB332_taxtab_pr2.tsv")
taxonomy_full

# How similar are the OTUs to their best hit?
ggplot(taxonomy_full, aes(x = Id)) +
  geom_histogram(binwidth = 1, fill = "steelblue") +
  labs(x = "% identity to best hit", y = "Number of OTUs") +
  theme_minimal()

# Keep only the OTU name and the taxonomic ranks
taxonomy <- taxonomy_full %>%
  select(OTUname, Domain:Species)

# Which supergroups and divisions are in the data?
taxonomy %>% count(Supergroup, sort = TRUE)
taxonomy %>% count(Division, sort = TRUE)


### CREATE PHYLOSEQ OBJECT ####
# Create OTU table from the rarefied table in lab I (samples as rows, OTUs as columns)
otu_table_ps <- otu_table(otu.vegan.ss.nozero, taxa_are_rows = FALSE)

# Create taxonomy table
tax_table_ps <- taxonomy %>%
  column_to_rownames("OTUname") %>%
  as.matrix() %>%
  tax_table()

# Create the sample data from the metadata
sample_data_ps <- metadata %>%
  column_to_rownames("Sample_Name") %>%
  sample_data()

# Create phyloseq object of otu abundances, taxonomy and sample metadata:
# NB! Only OTUs that are found in both the OTU table and the taxonomy are kept
physeq <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
physeq

ntaxa(otu_table_ps) # Number of OTUs in the OTU table
ntaxa(physeq)       # Number of OTUs with taxonomy


### PHYLOSEQ-BASED VISUALIZATIONS ####
# Calculate relative abundances
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

### Simple plotting:
plot_bar(physeq_rel, fill = "Supergroup") +
  geom_bar(stat = "identity", position = "stack")

plot_bar(physeq_rel, fill = "Division", facet_grid = ~ seasons) +
  geom_bar(stat = "identity", position = "stack")


#### Supergroup-level composition ####
supergroup_plot_ps <- plot_bar(physeq_rel, fill = "Supergroup") +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative Abundance") +
  ggtitle("Supergroup-Level Composition")

print(supergroup_plot_ps)

#### Division-level composition ####
physeq_division <- tax_glom(physeq_rel, "Division")

top_divisions_ps <- taxa_sums(physeq_division) %>%
  sort(decreasing = TRUE) %>%
  head(10) %>%
  names()

physeq_division_top <- prune_taxa(top_divisions_ps, physeq_division)

division_plot_ps <- plot_bar(physeq_division_top, fill = "Division") +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative Abundance") +
  ggtitle("Top 10 Divisions Composition")

division_plot_ps

#### Composition within the most abundant division ####
top_division_name <- tax_table(physeq_division)[top_divisions_ps[1], "Division"] %>% as.character()
top_division_name

# Keep only the taxa belonging to this division and agglomerate at class level
top_division_only <- subset_taxa(physeq_rel, Division == top_division_name)
physeq_class <- tax_glom(top_division_only, "Class")

class_plot_ps <- plot_bar(physeq_class, fill = "Class") +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative Abundance") +
  ggtitle(paste("Classes within", top_division_name))
class_plot_ps

#### A selected group: Diatoms ####
# Try with other groups from the taxonomy, e.g. Subdivision == "Ciliophora" or Class == "Syndiniales"
physeq_diatoms <- subset_taxa(physeq_rel, Class %in% c("Bacillariophyceae", "Mediophyceae", "Coscinodiscophyceae"))
physeq_diatoms_order <- tax_glom(physeq_diatoms, "Order")

plot_bar(physeq_diatoms_order, fill = "Order") +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative Abundance") +
  ggtitle("Diatoms")


#### Alpha diversity by season ####
alpha <- estimate_richness(physeq) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(metadata, by = "Sample_Name")

ggplot(alpha, aes(x = seasons, y = Shannon, fill = seasons)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1) +
  theme_bw() +
  theme(legend.position = "none")


#### NMDS in phyloseq ####
set.seed(666)
ord <- ordinate(physeq_rel, method = "NMDS", distance = "bray")
plot_ordination(physeq_rel, ord, color = "seasons") +
  scale_colour_brewer(type = "qual", palette = "Set2") +
  theme_bw()

# Network of samples, samples are connected if the distance between them is less than maxdist
ig <- make_network(physeq, distance = "bray", maxdist = 0.4)
plot_network(ig, physeq, color = "seasons", label = NULL)

# Following the samples through time:
# plot_ordination(justDF = TRUE) returns the ordination scores and the metadata as a data frame
df_ord <- plot_ordination(physeq_rel, ord, justDF = TRUE) %>%
  arrange(date)

# Plot with lines connecting the samples in the order they were taken
ggplot(df_ord, aes(x = NMDS1, y = NMDS2)) +
  geom_path(aes(group = year), linewidth = 0.4, alpha = 0.6, colour = "grey50") +
  geom_point(aes(colour = seasons), size = 2) +
  facet_wrap(~ year) +
  coord_equal() +
  scale_colour_brewer(type = "qual", palette = "Set2") +
  labs(title = "NMDS with samples connected in time") +
  theme_bw()

### SAVE RESULTS ####
# Save phyloseq objects
save(physeq, physeq_rel, file = "AB332_phyloseq_results.RData")
