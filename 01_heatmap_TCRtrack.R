library(Seurat)
library(tidyverse)
library(here)
library(patchwork)
library(pals)

set.seed(42)

dir.create(here("outs", "scRNAseq", "figures"), recursive = T)

# import data --------------------------------------------------

intgr <- read_rds(path_to_intgr_seurat) # read seurat object (fully prepared by Daniil Lukyanov)


# downsample data for heatmap rna ---------------------------------------------


# take equal numbers of cells from public and our own datasets 
intgrKriukova <- subset(intgr, subset = orig.ident %in% c("DL01_rep1", "DL01_rep2", "DL04_rep1", "DL05_rep1"))
intgrPublic <- subset(intgr, subset = orig.ident %in% c("Ncl_EMTAB10026", "Sanger_EMTAB10026", "Cambridge_EMTAB10026"))
intgrPublic <- subset(intgrPublic, cells = sample(1:length(rownames(intgrPublic[[]])), (rownames(intgrKriukova[[]]) %>% length())))

intgr_for_rna <- merge(intgrKriukova, intgrPublic)

rm(intgrKriukova)
rm(intgrPublic)

# subset 10000 cells for easier visualization
intgr_for_rna <- subset(intgr_for_rna, cells = sample(1:length(rownames(intgr_for_rna[[]])), 10000))


# heatmap rna -------------------------------------------------------------


intgr_for_rna[[]] %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  mutate(number = as.factor(as.numeric(.$number))) -> plot_data

plot_data %>%
  arrange(number) %>%
  select(barcode) %>%
  mutate(tile_order = 1:length(.$barcode)) -> barcode_tile_order #arrange cells in order for heatmap (x axis for geom_tile)

plot_data %>%
  arrange(number) %>%
  select(barcode, number) %>%
  mutate(tile_order = 1:length(.$barcode)) %>%
  group_by(number) %>%
  summarise(cluster_border = max(tile_order)) %>%
  mutate(text_position = c(0, .$cluster_border[1:length(.$cluster_border)-1]) + 
           (.$cluster_border - c(0, .$cluster_border[1:length(.$cluster_border)-1]))/2) -> cluster_borders # calculate x for cluster borders on heatmap and x for cluster label positions on heatmap


plot_data %>%
  left_join(barcode_tile_order, by = "barcode") -> plot_data


plot_data %>%
  ggplot() +
  geom_tile(aes(x = tile_order, y = 1, fill = colors), color = NA) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "white") +
  geom_text(data = cluster_borders, aes(x = text_position, y = 2, label = number), size = 3) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme_void() +
  scale_fill_identity(guide = "legend", labels = color_code$number, breaks = color_code$colors) +
  ylim(0.50, 2.50) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) -> plot_cluster_RNA

ggplot(plot_data) +
  geom_tile(aes(x = tile_order, y = 1, fill = orig.ident), color = NA) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "white") +
  theme_void() +
  ylim(0.50, 2.50) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) -> plot_cells_RNA


intgr_for_rna@assays$RNA@data %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  filter(feature %in% c("CXCR3", "CCR6", "CCR4", "CXCR5", "IL2RA", "PTGDR2", "TBX21", "RORC", "GATA3", "FOXP3")) %>%
  mutate(feature = factor(feature, levels = c("CXCR5", "FOXP3", "IL2RA", "PTGDR2", "GATA3", "CCR4", "RORC", "CCR6", "TBX21", "CXCR3"))) %>%
  pivot_longer(cols = (-feature), names_to = "barcode", values_to = "expression") %>%
  left_join(barcode_tile_order, by = "barcode") -> plot_data_expression_RNA # extract normalized RNA counts for heatmap 

ggplot(plot_data_expression_RNA) +
  geom_tile(aes(x = tile_order, y = feature, fill = expression), color = NA) +
  theme_bw() +
  scale_fill_continuous(low = "white", high = "black")+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "grey70", size = 0.5, linetype = "dashed") -> plot_RNA

plot_cluster_RNA + plot_RNA + plot_cells_RNA + plot_layout(heights = c(1, 7, 1)) -> plot_rna_combined
ggsave(here("outs", "scRNAseq", "figures", "plot_RNA_combined.pdf"), plot = plot_rna_combined)

# downsample data for heatmap tcrtrack -------------------------------------

# subset 10000 cells for easier visualization
# for tcrtrack we only include cells from this study
intgr_for_tcrtrack <- subset(intgr, subset = orig.ident %in% c("DL01_rep1", "DL01_rep2", "DL04_rep1", "DL05_rep1"))
intgr_for_tcrtrack <- subset(intgr_for_tcrtrack, cells = sample(1:length(rownames(intgr_for_tcrtrack[[]])), 10000))

# heatmap tcrtrack -------------------------------------------------------------

intgr_for_tcrtrack[[]] %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  mutate(number = as.factor(as.numeric(.$number))) -> plot_data

plot_data %>%
  arrange(number) %>%
  select(barcode) %>%
  mutate(tile_order = 1:length(.$barcode)) -> barcode_tile_order  #arrange cells in order for heatmap (x axis for geom_tile)

plot_data %>%
  arrange(number) %>%
  select(barcode, number) %>%
  mutate(tile_order = 1:length(.$barcode)) %>%
  group_by(number) %>%
  summarise(cluster_border = max(tile_order)) %>%
  mutate(text_position = c(0, .$cluster_border[1:length(.$cluster_border)-1]) + 
           (.$cluster_border - c(0, .$cluster_border[1:length(.$cluster_border)-1]))/2) -> cluster_borders # calculate x for cluster borders on heatmap and x for cluster label positions on heatmap


plot_data %>%
  pivot_longer(cols = starts_with("tcrtrack") & ends_with("trb"), names_to = "TCRtrack_phenotype", values_to = "freq") -> plot_data_long

plot_data_long %>% 
  mutate(freq = case_when(is.na(.$freq) ~ 0,
                          TRUE ~ 1)) -> plot_data_long

plot_data_long %>%
  left_join(barcode_tile_order, by = "barcode") -> plot_data_long


plot_data_long %>%
  ggplot() +
  geom_tile(aes(x = tile_order, y = 1, fill = colors), color = NA) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "white") +
  geom_text(data = cluster_borders, aes(x = text_position, y = 2, label = number), size = 3) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme_void() +
  scale_fill_identity(guide = "legend", labels = color_code$number, breaks = color_code$colors) +
  ylim(0.50, 2.50) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) -> plot_cluster_tcrtrack

ggplot(plot_data_long) +
  geom_tile(aes(x = tile_order, y = 1, fill = orig.ident), color = NA) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "white") +
  theme_void() +
  ylim(0.50, 2.50) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) -> plot_cells_tcrtrack

plot_data_long %>%
  mutate(TCRtrack_phenotype = str_remove_all(.$TCRtrack_phenotype, "tcrtrack\\_|\\_trb")) %>%
  mutate(TCRtrack_phenotype = str_replace(.$TCRtrack_phenotype, "\\.", "\\-")) %>%
  mutate(TCRtrack_phenotype = factor(TCRtrack_phenotype, levels = c("Tfh", "TREG", "Th2a", "Th2", "Th22", "Th17", "Th1-17", "Th1"))) %>%
  ggplot() +
  geom_tile(aes(x = tile_order, y = TCRtrack_phenotype, fill = as.character(freq)), color = NA) +
  theme_bw() +
  scale_fill_manual(values = c("white", "black")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "grey70", size = 0.5, linetype = "dashed") -> plot_tcrtrack

plot_cluster_tcrtrack + plot_tcrtrack + plot_cells_tcrtrack + plot_layout(heights = c(1, 7, 1))
ggsave(here("outs", "scRNAseq", "figures", "plot_tcrtrack_combined.pdf"), height = 8, width = 11, dpi = 600)


# downsample data for heatmap cite-seq -------------------------------------

# subset 10000 cells for easier visualization
# for cite-seq heatmap we only take cells from public data
intgr_for_adt <- subset(intgr, subset = orig.ident %in% c("Ncl_EMTAB10026", "Sanger_EMTAB10026", "Cambridge_EMTAB10026"))
intgr_for_adt <- subset(intgr_for_adt, cells = sample(1:length(rownames(intgr_for_adt[[]])), 10000))

# heatmap cite-seq --------------------------------------------------------

intgr_for_adt[[]] %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  mutate(number = as.factor(as.numeric(.$number))) -> plot_data

plot_data %>%
  arrange(number) %>%
  select(barcode) %>%
  mutate(tile_order = 1:length(.$barcode)) -> barcode_tile_order #arrange cells in order for heatmap (x axis for geom_tile)

plot_data %>%
  arrange(number) %>%
  select(barcode, number) %>%
  mutate(tile_order = 1:length(.$barcode)) %>%
  group_by(number) %>%
  summarise(cluster_border = max(tile_order)) %>%
  mutate(text_position = c(0, .$cluster_border[1:length(.$cluster_border)-1]) + 
           (.$cluster_border - c(0, .$cluster_border[1:length(.$cluster_border)-1]))/2) -> cluster_borders # calculate x for cluster borders on heatmap and x for cluster label positions on heatmap


plot_data %>%
  left_join(barcode_tile_order, by = "barcode") -> plot_data


plot_data %>%
  ggplot() +
  geom_tile(aes(x = tile_order, y = 1, fill = colors), color = NA) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "white") +
  geom_text(data = cluster_borders, aes(x = text_position, y = 2, label = number), size = 3) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme_void() +
  scale_fill_identity(guide = "legend", labels = color_code$number, breaks = color_code$colors) +
  ylim(0.50, 2.50) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) -> plot_cluster_protein

ggplot(plot_data) +
  geom_tile(aes(x = tile_order, y = 1, fill = orig.ident), color = NA) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "white") +
  theme_void() +
  ylim(0.50, 2.50) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) -> plot_cells_protein


intgr_for_adt@assays$ADT@data %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  filter(feature %in% c("AB-CXCR3", "AB-CCR6", "AB-CCR4", "AB-CXCR5", "AB-CD25", "AB-PTGDR2")) %>%
  mutate(feature = factor(feature, levels = c("AB-CXCR5", "AB-CD25", "AB-PTGDR2", "AB-CCR4", "AB-CCR6", "AB-CXCR3"))) %>%
  pivot_longer(cols = (-feature), names_to = "barcode", values_to = "expression") %>%
  left_join(barcode_tile_order, by = "barcode") -> expression_PROTEIN # extract normalized ADT counts for heatmap

# add Z-scoring within each ADT feature
# after z-scoring I trim all extreme z-score values higher than 99.5 quantile. This is needed for easier visualization in color scale
expression_PROTEIN %>%
  group_by(feature) %>%
  mutate(normalized_expression = (expression - mean(expression))/sd(expression)) %>%
  mutate(normalized_expression_trimmed = case_when(normalized_expression > quantile(normalized_expression, seq(0, 1, 0.005))[200] ~ quantile(normalized_expression, seq(0, 1, 0.005))[200], 
                                                   TRUE ~ normalized_expression)) -> plot_data_expression_PROTEIN

plot_data_expression_PROTEIN %>%
  ggplot() +
  geom_tile(aes(x = tile_order, y = feature, fill = normalized_expression_trimmed), color = NA) +
  theme_bw() +
  scale_fill_gradient2(low = "green", mid = "lightyellow", high = "magenta", na.value = "grey")+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = cluster_borders$cluster_border, color = "grey70", size = 0.5, linetype = "dashed") +
  labs(fill = "Trimmed normalized expression")-> plot_PROTEIN

plot_cluster_protein + plot_PROTEIN + plot_cells_protein + plot_layout(ncol = 1, heights = c(1, 7, 1)) 
ggsave(here("outs", "scRNAseq", "figures", "plot_CITEseq_combined.pdf"), height = 8, width = 12, dpi = 600)



plot_cluster_tcrtrack + plot_tcrtrack + plot_cluster_RNA + plot_RNA + plot_cluster_protein + plot_PROTEIN + plot_layout(ncol = 1, heights = c(1, 8, 1, 10, 1, 6)) 
ggsave(here("outs", "scRNAseq", "figures", "plot_RNA_SORT_CITE_combined.png"), height = 7, width = 13, bg = "white", dpi = 600)

