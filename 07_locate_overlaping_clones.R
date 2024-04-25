library(tidyverse)
library(here)
library(Seurat)


# load data ---------------------------------------------------------------

filenames <- list.files(path_to_Kasatskaya_data, pattern = "Donor1|Donor4|Donor5", full.names=TRUE)
filenames_short <- list.files(path_to_Kasatskaya_data, pattern = "Donor1|Donor4|Donor5", full.names=F) %>% str_remove(".txt")

clonesets <- map(filenames, read_tsv) 
names(clonesets) <- filenames_short

clonesets <- map(clonesets, dplyr::slice, 1:2000) #filter top2000 clones
clonesets2 <- bind_rows(clonesets, .id = "sample_id") 

clonesets2 %>%
  separate(sample_id, into = c("donor_id", "replicate", "chain", "subset"), sep = "_") %>%
  filter(chain == "TRB") %>%
  group_by(donor_id, chain, subset, cdr3nt, cdr3aa, v, d, j) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(donor_id, chain, subset) %>%
  mutate(freq = count/sum(count)) %>%
  arrange(desc(freq), .by_group = T) %>%
  mutate(clone_id = paste(cdr3aa, cdr3nt, v, j, sep = "_")) %>%
  split(f = .$donor_id) -> data

map(data, function(x) {
  x %>%
    filter(subset %in% c("Th17", "Th22")) %>%
    group_by(cdr3nt,cdr3aa, v, d, j) %>%
    summarize(n = n()) %>%
    filter(n > 1) %>%
    mutate(clone_id = paste(cdr3aa, cdr3nt, v, j, sep = "_"))
}) -> shared_clonotypes


std <- function(x){
  sqrt(sum((x - mean(x))^2) / length(x))
}


map2(data, shared_clonotypes, function(x, y){
  x %>%
    filter(clone_id %in% y$clone_id) %>%
    group_by(cdr3nt, cdr3aa, v, j, clone_id) %>%
    summarise(std = std(freq))
}) -> std_clonotypes

# I check the histogram of the std, to select a threshold
std_clonotypes$Donor4 %>%
  ggplot(aes(x = std)) +
  geom_density()+
  scale_x_log10()

map(std_clonotypes, function(x) x %>% filter(std < 0.001)) -> std_low_clonotypes


map2(data, std_low_clonotypes, function(x, y){
  x %>%
    filter(clone_id %in% y$clone_id) %>%
    filter(subset %in% c("Th22", "Th17"))
}) %>% bind_rows(.id = "donor_id") -> final_shared_clonotypes



# load scRNA-Seq object -------------------------------------------------------------

intgr <- read_rds(path_to_intgr_seurat)



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


ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Tfh_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(TRB_clone_id %in% final_shared_clonotypes$clone_id)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("overlap Th22 Th17") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

ggsave(here("outs", "scRNAseq", "figures", "overlap_th22_th17.png"), height = 2, width = 2, dpi = 600)


ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Tfh_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th17_trb) & !TRB_clone_id %in% final_shared_clonotypes$clone_id)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th17 without overlapped") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

ggsave(here("outs", "scRNAseq", "figures", "Th17 without overlapped.png"), height = 2, width = 2, dpi = 600)



ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size = 0.1, shape = 16, color = "grey96", alpha = 1) +
  theme_minimal() +
  #geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Tfh_trb) & donor_id %in% c("D01", "D04", "D05"))), size = 0.3, shape = 16, aes(alpha = tcrtrack_Tfh_trb, color  = donor_id)) +
  geom_point(data = (umap_tx %>% filter(!is.na(tcrtrack_Th22_trb) & !TRB_clone_id %in% final_shared_clonotypes$clone_id)), size = 0.3, shape = 16, color = "darkorchid1", alpha = 0.2) +
  #scale_color_manual(values = c("#0838A6", "#F25C3D", "yellowgreen")) +
  ggtitle("Th22 without overlapped") +
  scale_alpha(limits = c(alpha_min, alpha_max)) +
  theme_void() +
  coord_fixed() +
  #guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  theme(legend.position = "none")

ggsave(here("outs", "scRNAseq", "figures", "Th22 without overlapped.png"), height = 2, width = 2, dpi = 600)





#### Th17 and Th22 in the same style


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

ggsave(here("outs", "scRNAseq", "figures", "th17.png"), plot = p_Th17, height = 2, width = 2, dpi = 600)
ggsave(here("outs", "scRNAseq", "figures", "th22.png"), plot = p_Th22, height = 2, width = 2, dpi = 600)