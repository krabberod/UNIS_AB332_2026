# Taxonomic Analysis of Isfjorden Microbial Eukaryote Communities - TIDYVERSE VERSION
# AB332 - Community Ecology Course

# Load required libraries
library(tidyverse)
library(vegan)

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
# OTUs that were not found in PR2 are called "Unassigned"
taxonomy <- taxonomy_full %>%
  select(OTUname, Domain:Species) %>%
  mutate(across(Domain:Species, ~ replace_na(.x, "Unassigned")))

# Which supergroups and divisions are in the data?
taxonomy %>% count(Supergroup, sort = TRUE)
taxonomy %>% count(Division, sort = TRUE)


### COMBINE OTU TABLE, TAXONOMY AND METADATA ####
# We make one long table with one row per OTU per sample (this is what phyloseq does for you)
otu.long <- otu.vegan.ss.nozero %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  pivot_longer(-Sample_Name, names_to = "OTUname", values_to = "reads") %>%
  left_join(taxonomy, by = "OTUname") %>%
  left_join(metadata, by = "Sample_Name") %>%
  group_by(Sample_Name) %>%
  mutate(rel_abund = reads / sum(reads)) %>% # Relative abundance within each sample
  ungroup()
otu.long


### VISUALIZATIONS ####
#### Supergroup-level composition ####
# Sum the relative abundance of all OTUs in the same supergroup in each sample
supergroup <- otu.long %>%
  group_by(Sample_Name, Supergroup) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

supergroup_plot <- ggplot(supergroup, aes(x = Sample_Name, y = rel_abund, fill = Supergroup)) +
  geom_col() +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Relative Abundance", title = "Supergroup-Level Composition")

supergroup_plot

#### Division-level composition by season ####
division_season <- otu.long %>%
  group_by(Sample_Name, seasons, Division) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

ggplot(division_season, aes(x = Sample_Name, y = rel_abund, fill = Division)) +
  geom_col() +
  facet_grid(~ seasons, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Relative Abundance")

#### Top 10 divisions ####
# fct_lump_n() keeps the 10 most abundant divisions and puts the rest in "Other"
division_top <- otu.long %>%
  mutate(Division = fct_lump_n(Division, n = 10, w = rel_abund)) %>%
  group_by(Sample_Name, Division) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

division_plot <- ggplot(division_top, aes(x = Sample_Name, y = rel_abund, fill = Division)) +
  geom_col() +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Relative Abundance", title = "Top 10 Divisions Composition")

division_plot

#### Composition within the most abundant division ####
top_division_name <- otu.long %>%
  count(Division, wt = rel_abund, sort = TRUE) %>%
  slice(1) %>%
  pull(Division)
top_division_name

# Keep only the OTUs belonging to this division and sum at class level
class_top <- otu.long %>%
  filter(Division == top_division_name) %>%
  group_by(Sample_Name, Class) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

class_plot <- ggplot(class_top, aes(x = Sample_Name, y = rel_abund, fill = Class)) +
  geom_col() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Relative Abundance", title = paste("Classes within", top_division_name))

class_plot

#### A selected group: Diatoms ####
# Try with other groups from the taxonomy, e.g. Subdivision == "Ciliophora" or Class == "Syndiniales"
diatoms <- otu.long %>%
  filter(Class %in% c("Bacillariophyceae", "Mediophyceae", "Coscinodiscophyceae")) %>%
  group_by(Sample_Name, Order) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

ggplot(diatoms, aes(x = Sample_Name, y = rel_abund, fill = Order)) +
  geom_col() +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Relative Abundance", title = "Diatoms")


#### Alpha diversity by season ####
alpha <- tibble(
  Sample_Name = rownames(otu.vegan.ss.nozero),
  Richness = specnumber(otu.vegan.ss.nozero),
  Shannon = diversity(otu.vegan.ss.nozero, index = "shannon")
) %>%
  left_join(metadata, by = "Sample_Name")

ggplot(alpha, aes(x = seasons, y = Shannon, fill = seasons)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1) +
  theme_bw() +
  theme(legend.position = "none")


#### NMDS ####
# We use the NMDS scores from lab II (already joined with the metadata)
ggplot(otu.vegan.ss.nozero.bray.nmds.scores, aes(x = NMDS1, y = NMDS2, colour = seasons)) +
  geom_point(size = 2) +
  scale_colour_brewer(type = "qual", palette = "Set2") +
  coord_equal() +
  theme_bw()

# Following the samples through time:
df_ord <- otu.vegan.ss.nozero.bray.nmds.scores %>%
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

# A network of samples can be made with the phyloseq package, see Computer_lab_IV-taxonomic_analysis_phyloseq.R

### SAVE RESULTS ####
save(otu.long, taxonomy, file = "AB332_lab_IV.RData")
