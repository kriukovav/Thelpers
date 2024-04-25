library(Seurat)
library(tidyverse)
library(here)
library(patchwork)
library(pals)

set.seed(42)

#path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
#folder_with_gating_results <- here("data", "gating_V2")
intgr <- read_rds(path_to_intgr_seurat)


# extract UMAP coordinates for cells --------------------------------------

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  rownames_to_column(var = "barcode") %>%
  slice(x = sample(1:length(.$barcode)), size = length(.$barcode), replace = FALSE) #mix rows for more clear plotting


alpha_min <- umap_tx %>%
  select(ends_with("trb")) %>%
  as.matrix() %>%
  min(na.rm = T)

alpha_max <- umap_tx %>%
  select(ends_with("trb")) %>%
  as.matrix() %>%
  max(na.rm = T)



# TCRtrack beta chain -----------------------------------------------------


p_Tfh <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Tfh_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Tfh_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Tfh") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Treg <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_TREG_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_TREG_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_TREG_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_TREG_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Treg") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th117 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th1.17_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th1.17_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th1.17_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th1.17_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th1-17") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th2a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2a_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th2a_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2a_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th2a_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th2a") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th17 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th17_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th17_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th17_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th17_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th17") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")


p_Th2 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th2_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th2_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th2") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")


p_Th22 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th22_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th22_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th22_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th22_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th22") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th1 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th1_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th1_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th1_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th1_trb)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th1") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")


wrap_plots(list(p_Th1, p_Th17, p_Th117, p_Th22, p_Th2, p_Th2a, p_Tfh, p_Treg)) + plot_layout(ncol = 2, byrow = F, guides = "collect")


# TCRtrack alpha chain ----------------------------------------------------


alpha_min <- umap_tx %>%
  select(ends_with("tra")) %>%
  as.matrix() %>%
  min(na.rm = T)

alpha_max <- umap_tx %>%
  select(ends_with("tra")) %>%
  as.matrix() %>%
  max(na.rm = T)

p_Tfh_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Tfh_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Tfh_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Tfh") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Treg_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_TREG_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_TREG_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_TREG_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_TREG_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Treg") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th117_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(`tcrtrack_Th1.17_tra`) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = `tcrtrack_Th1.17_tra`, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(`tcrtrack_Th1.17_tra`) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = `tcrtrack_Th1.17_tra`)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th1-17") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th2a_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2a_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th2a_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2a_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th2a_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th2a") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th17_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th17_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th17_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th17_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th17_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th17") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")


p_Th2_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th2_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th2_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th2_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th2") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")


p_Th22_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th22_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th22_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th22_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th22_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th22") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

p_Th1_a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th1_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Th1_tra, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th1_tra) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, color = "darkorchid1", aes(alpha = tcrtrack_Th1_tra)) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th1") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")


wrap_plots(list(p_Th1_a, p_Th17_a, p_Th117_a, p_Th22_a, p_Th2_a, p_Th2a_a, p_Tfh_a, p_Treg_a)) + plot_layout(ncol = 2, byrow = F)

ggsave(here("outs", "scRNAseq", "figures", "tcrtrack_alpha.png"), height = 6, width = 6, dpi = 600)



# map from citeseq --------------------------------------------------------


filenames <- list.files(folder_with_gating_results, pattern = "^t.*citeseq.txt", full.names = T)
filenames_short <- list.files(folder_with_gating_results, pattern = "^t.*citeseq.txt", full.names = F) %>% str_remove("\\_citeseq\\.txt")

citeseq_gated_th <- map(filenames, read_tsv, col_names = "barcode")
names(citeseq_gated_th) <- paste0("citeseq_", filenames_short)

citeseq_gated_th <- bind_rows(citeseq_gated_th, .id = "phenotype")

citeseq_gated_th %>%
  mutate(accessory_column = 1) %>%
  pivot_wider(names_from = phenotype, values_from = accessory_column, values_fill = 0) %>%
  right_join(umap_tx) -> umap_tx

cite_p_Tfh <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_tfh == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_tfh == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Tfh") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

cite_p_Treg <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_treg == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_treg == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Treg") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

cite_p_Th117 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_th1_17 == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_th1_17 == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Th1-17") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


cite_p_Th17 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_th17 == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_th17 == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Th17") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


cite_p_Th2 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_th2 == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_th2 == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Th2") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

cite_p_Th2a <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_th2a == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_th2a == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Th2a") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


cite_p_Th22 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_th22 == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_th22 == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Th22") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

cite_p_Th1 <- ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.01, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(citeseq_th1 == 1)), size = 0.3, shape = 16, aes(color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(citeseq_th1 == 1)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "#9DF20C")) +
  ggtitle("Th1") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


wrap_plots(list(p_Th1, p_Th17, p_Th117, p_Th22, p_Th2, p_Th2a, p_Tfh, p_Treg,
                plot_spacer(), plot_spacer(), plot_spacer(), plot_spacer(),
                cite_p_Th1, cite_p_Th17, cite_p_Th117, cite_p_Th22, cite_p_Th2, cite_p_Th2a, cite_p_Tfh, cite_p_Treg)) + plot_layout(ncol = 5, byrow = F, guides = "collect", widths = c(1, 1, 0.2, 1, 1))
ggsave(here("outs", "scRNAseq", "figures", "Fig2_tcrtrack_citeseq.png"), height = 6, width = 6, dpi = 600)

