
library(ggplot2)
library(Cairo)
library(scattermore)
library(pals)
library(ggrepel)
library(patchwork)
library(purrr)
library(pool)
library(RSQLite)
library(dplyr)
library(magick)
library(stringr)
library(shinyalert)
library(fst)
library(dbplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))

#scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname ="/Volumes/Thunder//data/scEiaD_2022_02//MOARTABLES__anthology_limmaFALSE___4000-counts-universe-study_accession-scANVIprojection-15-5-0.1-50-20__pointRelease01.sqlite", idleTimeout = 3600000)
scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname ="~/data/scEiaD_2022_02/MOARTABLES__anthology_limmaFALSE___4000-counts-universe-study_accession-scANVIprojection-15-5-0.1-50-20__pointRelease01.sqlite", idleTimeout = 3600000)

# # testing
# load('~/data/scEiaD_CTP/xgboost_predictions/n_features-2000__transform-counts__partition-PR__covariate-batch__method-scVI__dims-20__epochs-50__dist-0.1__neighbors-50__knn-20__umapPredictions.Rdata')

x_dir <- -1
y_dir <- 1

#meta_filter <- read_fst('/Volumes/Thunder/data/scEiaD_2022_02//meta_filter.fst') %>%
meta_filter <- read_fst('~/data/scEiaD_2022_02//meta_filter.fst') %>%
  as_tibble() %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) && !is.na(CellType_predict) ~ 'Tabula Muris',
                                      is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) %>%
  # filter(study_accession != 'Bharti_Nguyen_iRPE_2D_3D') %>%
  mutate(UMAP_a = UMAP_1 * x_dir,
         UMAP_b = UMAP_2 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)


box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
  filter(organism == 'Homo sapiens', Source == 'Tissue') %>%
  collect()

box_data2 <- box_data %>%
  filter(Source == 'Tissue') %>% 
  mutate(Organ = case_when(Organ == 'Eye' ~ 'Eye', TRUE ~ 'Body')) %>% 
  filter(!Platform %in% c('SCRBSeq')) %>%
  mutate(Stage = factor(Stage, levels = c('Early Dev.', 'Late Dev.', 'Maturing', 'Mature'))) %>%
  #filter(!is.na(!!as.symbol(grouping_features))) %>%
  group_by_at(vars(one_of(c('Gene', 'study_accession',  'Organ',  'CellType_predict', 'Stage')))) %>%
  summarise(counts = sum(counts * cell_exp_ct) / sum(cell_exp_ct),
            cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
  full_join(., meta_filter %>%
              filter(organism == 'Homo sapiens',
                     Source == 'Tissue') %>%
              mutate(Organ = case_when(Organ == 'Eye' ~ 'Eye', TRUE ~ 'Body')) %>% 
              group_by_at(vars(one_of(c( 'study_accession',  'Organ',  'CellType_predict',  'Stage')))) %>%
              summarise(Count = n())) %>%
  mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
  mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
         Expression = round(counts * (`%` / 100), 2)) %>%
  select_at(vars(one_of(c('Gene', 'study_accession',  'Organ',  'CellType_predict',  'Stage',
                          'cell_exp_ct', 'Count', '%', 'Expression'
  )))) %>%
  arrange(-Expression) %>%
  dplyr::rename(`Cell # Detected` = cell_exp_ct,
                `Total Cells` = Count,
                `Mean log2(Counts + 1)` = Expression,
                `% of Cells Detected` = `%`) %>%
  mutate(`Mean log2(Counts + 1)` = case_when(is.na(`Mean log2(Counts + 1)`) ~ 0, TRUE ~ `Mean log2(Counts + 1)`))
# expand missing genes (nothing detected) to all genes used
# box_data_NA <- box_data2 %>% filter(is.na(Gene))
# if (nrow(box_data_NA > 0)){
#   box_data_NA_list <- list()
#   for (i in gene){
#     box_data_NA_list <- box_data_NA %>% mutate(Gene = i)
#   }
#   box_data2 <- bind_rows(box_data2 %>% filter(!is.na(Gene)),
#                          box_data_NA_list %>% bind_rows())
# }
scEiaD_CT_table <- box_data2 %>% filter(!is.na(Gene), `Total Cells` > 50 )






save(scEiaD_CT_table, file = 'data/scEiaD_CT_table_2023_03_01.Rdata')
