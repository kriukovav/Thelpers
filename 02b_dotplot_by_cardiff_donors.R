library(Seurat)
library(tidyverse)
library(here)
library(patchwork)
library(pals)

set.seed(42)

# path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
# folder_with_gating_results <- here("data", "gating_V2")"

dir.create(here("outs", "scRNAseq", "figures"), recursive = T)

# import data --------------------------------------------------

#intgr <- read_rds(path_to_intgr_seurat) # read seurat object (fully prepared by Daniil Lukyanov)


# downsample data for tcrtrack -------------------------------------

intgr_for_tcrtrack <- subset(intgr, subset = donor_id %in% donor) # subset the sample according to the orig.ident, selected in the main script

intgr_for_tcrtrack[[]] %>%
  rownames_to_column(var = "barcode") %>%
  left_join(annotation_table) %>%
  left_join(color_code) %>%
  mutate(number = as.factor(as.numeric(.$number))) -> plot_data

plot_data %>%
  pivot_longer(cols = starts_with("tcrtrack") & ends_with("trb"), names_to = "TCRtrack_phenotype", values_to = "freq") %>%
  mutate(TCRtrack_phenotype = str_remove_all(.$TCRtrack_phenotype, "tcrtrack\\_|\\_trb")) %>%
  mutate(TCRtrack_phenotype = str_replace(.$TCRtrack_phenotype, "\\.", "\\-")) -> plot_data_long

plot_data_long %>% 
  mutate(freq = case_when(is.na(.$freq) ~ 0,
                          TRUE ~ 1)) -> plot_data_long


plot_data_long %>%
  select(barcode, TCRtrack_phenotype, freq, number) %>%
  split(f = .$TCRtrack_phenotype) %>%
  map(function(x) {
    x %>%
      group_by(number) %>%
      summarise(n_cells_in_cluster = n(),
                Th_positive_cells = sum(freq),
                `% of cluster stained by TCRtrack` = mean(freq))
  }) %>%
  bind_rows(.id = "TCRtrack_phenotype") -> simple_heatmap_v1


plot_data_long %>%
  select(barcode, TCRtrack_phenotype, freq) %>%
  distinct() %>%
  group_by(TCRtrack_phenotype) %>%
  summarize(size_of_TCRtrack_subset = sum(freq)) -> size_of_TCRtrack_subset_table

plot_data_long %>%
  select(barcode, number) %>%
  distinct() %>%
  group_by(number) %>%
  summarize(size_of_seurat_cluster = n()) -> size_of_seurat_cluster

plot_data_long %>%
  select(barcode, number, TCRtrack_phenotype, freq) %>%
  left_join(size_of_TCRtrack_subset_table) %>% 
  left_join(size_of_seurat_cluster) %>%
  group_by(number, TCRtrack_phenotype, size_of_TCRtrack_subset, size_of_seurat_cluster) %>%
  summarize(n_cell_Th_positive_within_cluster = sum(freq)) %>% 
  mutate(`% of TCRtrack clonotypes mapped to cluster` = n_cell_Th_positive_within_cluster / size_of_TCRtrack_subset) %>%
  mutate(x_axis = paste0(number, " (", size_of_seurat_cluster, " cells)"),
         y_axis = paste0(TCRtrack_phenotype, " (", size_of_TCRtrack_subset, " cells)")) -> simple_heatmap_v2

ink_stain_plot_data <- simple_heatmap_v1 %>%
  full_join(simple_heatmap_v2, by = c("TCRtrack_phenotype", "number"), suffix = c("_v1", "_v2")) %>%
  mutate(x_axis = factor(x_axis, levels = paste0(.$number, " (", .$size_of_seurat_cluster, " cells)") %>% unique()))

y_order <- data.frame("TCRtrack_phenotype" = c("Tfh", "TREG", "Th2a", "Th2", "Th22", "Th17", "Th1-17", "Th1"),
                      "y_order" = 1:8)

# color_max <- ink_stain_plot_data %>% 
#   left_join(y_order) %>%
#   select(`% of cluster stained by TCRtrack`) %>%
#   pull %>%
#   max()

color_max = 0.75

dot_plot_tcrtrack <- ink_stain_plot_data %>%
  left_join(y_order) %>%
  ggplot() +
  geom_point(aes(x = x_axis, y = fct_reorder(y_axis, y_order), color = `% of cluster stained by TCRtrack`, size = `% of TCRtrack clonotypes mapped to cluster`)) +
  scale_color_gradient(low = "white", high = "black", limits = c(0, color_max)) +
  scale_size_continuous(limits = c(0.000001, 1)) + #I introduced limits to remove dots corresponding to 0 from the plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("clusters") +
  ylab("TCRtrack phenotypes") +
  coord_fixed() +
  ggtitle(paste(donor))

ggsave(plot = dot_plot_tcrtrack, filename = here("outs", "scRNAseq", "figures", paste0(donor, "dot_plot.pdf")), height = 6, width = 10, bg = "white")


