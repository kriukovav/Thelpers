library(Seurat)
library(tidyverse)
library(here)
library(pals)

set.seed(42)

#path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
intgr <- read_rds(path_to_intgr_seurat)


# extract UMAP coordinates for cells --------------------------------------

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  rownames_to_column(var = "barcode")



# clonality plot for supplements ------------------------------------------

# 100% - all TRB clones

umap_tx %>%
  filter(donor_id %in% c("D01", "D04", "D05")) %>%
  select(donor_id, Subset, TRB_clone_id) %>%
  drop_na() %>% # we take only TRB containing cells %>%
  group_by(donor_id, Subset) %>%
  mutate(n_cells_in_cluster = n()) %>%
  distinct() %>% # keep only unique clonotypes
  group_by(donor_id, Subset,  n_cells_in_cluster) %>%
  summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
  group_by(donor_id) %>%
  mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
         fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
         n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

clonality_data %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  ggplot(aes(x = n_cells_in_cluster, y = n_clones_in_cluster, shape = donor_id, color = colors)) +
  geom_abline(slope = 1, color = "grey80") +
  geom_point(size = 4) +
  coord_fixed() +
  theme_minimal() +
  ggtitle("TCR\u03b2 clonality") +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  scale_x_log10() +
  scale_y_log10() +
  ylab("Number of TCR\u03b2 clonotypes in the cluster") +
  xlab("Cells with a known TCR\u03b2 in cluster")



# clonality all donors
umap_tx %>%
  #filter(donor_id %in% c("D01", "D04", "D05")) %>%
  select(donor_id, Subset, TRB_clone_id) %>%
  drop_na() %>% # we take only TRB containing cells %>%
  group_by(donor_id, Subset) %>%
  mutate(n_cells_in_cluster = n()) %>%
  distinct() %>% # keep only unique clonotypes
  group_by(donor_id, Subset,  n_cells_in_cluster) %>%
  summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
  group_by(donor_id) %>%
  mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
         fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
         n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

clonality_data %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  ggplot(aes(x = n_cells_in_cluster, y = n_clones_in_cluster, color = colors)) +
  geom_abline(slope = 1, color = "grey80") +
  geom_point(size = 2) +
  coord_fixed() +
  theme_minimal() +
  facet_wrap(facets = vars(colors)) +
  ggtitle("TCR\u03b2 clonality") +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  #scale_x_log10() +
  #scale_y_log10() +
  ylab("Number of TCR\u03b2 clonotypes in the cluster") +
  xlab("Cells with a known TCR\u03b2 in cluster")


clonality_data %>%
  mutate(ratio = n_clones_in_cluster / n_cells_in_cluster) %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  ggplot(aes(y = Subset, x = ratio, color = colors)) +
  #geom_violin() +
  geom_jitter(size = 1, height = 0.2) +
  theme_minimal() +
  ggtitle("TCR\u03b2 clonality") +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  scale_color_identity(guide = "legend", labels = annotation_table$Subset, breaks = color_code$colors) +
  #scale_x_log10() +
  #scale_y_log10() +
  ylab("Subsets") +
  xlab("N TCRbeta clonotypes in cluster / N cells in cluster")

ggsave(here("outs", "scRNAseq", "figures", "clonality.pdf"), height = 6, width = 6)