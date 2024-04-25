# This is a pipeline script which calls the analysis scripts
library(here)
library(pals)
library(tidyverse)


#01 heatmap figS6 
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")

color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))



source("01_heatmap_TCRtrack.R")
rm(list = ls(all.names = TRUE)) #will clear all objects

#02 dotplot comparison TCRtrack CITE-Seq
folder_with_gating_results <- here("data", "gating_V2")
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")

color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))



source("02_dot_plot_comparison_TCRtrack_cite.R")
rm(list = ls(all.names = TRUE)) #will clear all objects


#02b dotplots by cardiff donors
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
folder_with_gating_results <- here("data", "gating_V2")

color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))

intgr <- read_rds(path_to_intgr_seurat) # read seurat object (fully prepared by Daniil Lukyanov)

donor <- c("D01")
source("02b_dotplot_by_cardiff_donors.R")
donor <- c("D04")
source("02b_dotplot_by_cardiff_donors.R")
donor <- c("D05")
source("02b_dotplot_by_cardiff_donors.R")
rm(list = ls(all.names = TRUE)) #will clear all objects


#02c dotplots by status on day collection
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
folder_with_gating_results <- here("data", "gating_V2")

color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))

intgr <- read_rds(path_to_intgr_seurat) # read seurat object (fully prepared by Daniil Lukyanov)
intgr_for_adt <- subset(intgr, subset = orig.ident %in% c("Ncl_EMTAB10026", "Sanger_EMTAB10026", "Cambridge_EMTAB10026"))

rm(intgr)

status <- c("Moderate")
source("02c_dotplot_by_covid_status.R")
status <- c("Healthy")
source("02c_dotplot_by_covid_status.R")
status <- c("Severe")
source("02c_dotplot_by_covid_status.R")
status <- c("Mild")
source("02c_dotplot_by_covid_status.R")
status <- c("Non_covid")
source("02c_dotplot_by_covid_status.R")
status <- c("Critical")
source("02c_dotplot_by_covid_status.R")
status <- c("LPS_10hours")
source("02c_dotplot_by_covid_status.R")
status <- c("LPS_90mins")
source("02c_dotplot_by_covid_status.R")
status <- c("Asymptomatic")
source("02c_dotplot_by_covid_status.R")

rm(list = ls(all.names = TRUE)) #will clear all objects

#03 mapping of clonotypes from Kasatskaya et al to our UMAPs
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
folder_with_gating_results <- here("data", "gating_V2")

source("03_cardiff_coloring.R")
rm(list = ls(all.names = TRUE)) #will clear all objects


#04
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))


source("04_plot_clonality.R")
rm(list = ls(all.names = TRUE)) #will clear all objects



#05 annotation of the scRNA-Seq dataset
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))

source("05_annotation.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects



#06 IFN induced cells in COVID
path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")

color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))

source("06_IFNresponse_supplementaryFigure.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects


#07 Overlap between Th17 and Th22

path_to_intgr_seurat <- here("outs", "scRNAseq", "full_reference_return_model.rds")
path_to_Kasatskaya_data <- here("data", "Kasatskaya")


color_code <- data.frame("number" = as.character(0:15), "colors" = kelly()[c(22, 9, 11, 6, 8, 12, 21, 10, 5, 7, 13, 14, 19, 3, 17, 18)]) %>%
  mutate(colors = case_when(number == "9" ~ "azure2",
                            TRUE ~ colors))

annotation_table <-  data.frame("number" = as.character(0:15),
                                "Subset" = c("EffMem_Th2", "Treg", "Naive_RTE", "EffMem_Th1",
                                             "EffMem_Th17", "Naive", "EffMem_Th22", "CentMem1",
                                             "CentMem2", "Tfh", "Naive_IFN_induced", "PD1high",
                                             "EffMem_Th2a", "EffMem_IFN_induced", "Temra_cytotoxic_Th1", "Cycling"))


source("07_locate_overlaping_clones.R") 
rm(list = ls(all.names = TRUE)) #will clear all objects
