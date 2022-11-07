library(tidyverse)
library(tibble)
library(recount)
library(recount3)
library(dtplyr)
library(data.table)

source("~/git/EiaD_build/scripts/create_count_data_function.R")

# Importing Metadata
#eyeIntegration22 <- read_csv("https://hpc.nih.gov/~parikhpp/EiaD/2022_metadata.csv", col_types = cols(...1 = col_skip()))
emeta <- data.table::fread('data/eyeIntegration22_meta_2022_10_27.03.csv.gz') %>% as_tibble()
create_count_data_frames("http://duffel.rail.bio/recount3/",
                         c("SRP002881", "SRP011895", "SRP012585", "SRP015336", "SRP016140", "SRP034875",
                           "SRP035641", "SRP045639", "SRP053034", "SRP055101", "SRP061670", "SRP062870",
                           "SRP064956", "SRP070148", "SRP075990", "SRP079002", "SRP080002", "SRP090226",
                           "SRP091605", "SRP091675", "SRP092413", "SRP093877", "SRP094572", "SRP097696",
                           "SRP098761", "SRP105756", "SRP106457", "SRP108292", "SRP110135", "SRP111145",
                           "SRP115908", "SRP117613", "SRP119291", "SRP119766", "SRP151763", "SRP159246",
                           "ERP022243", "SRP186009", "SRP156453", "SRP136195"),
                         "recount3_TPM",
                         "long_recount3_TPM",
                         "recount3_count",
                         "long_recount3_count",
                         "recount3_mapping_information",
                         emeta, 
                         empty_cache = FALSE,
                         GTEX = FALSE)


create_count_data_frames("http://duffel.rail.bio/recount3/",
                              c("ADIPOSE_TISSUE", "MUSCLE", "BLOOD_VESSEL", "HEART", "OVARY", "UTERUS",
                                "VAGINA", "BREAST", "SKIN", "SALIVARY_GLAND", "BRAIN", "ADRENAL_GLAND",
                                "THYROID", "LUNG", "SPLEEN", "PANCREAS", "ESOPHAGUS", "STOMACH", "COLON",
                                "SMALL_INTESTINE", "PROSTATE", "TESTIS", "NERVE", "PITUITARY", "BLOOD",
                                "LIVER", "KIDNEY",   "CERVIX_UTERI", "FALLOPIAN_TUBE", "BLADDER", "BONE_MARROW"),
                              "GTEX_TPM",
                              "long_GTEX_TPM",
                              "GTEX_count",
                              "long_GTEX_count",
                              "gtex_mapping_information",
                              emeta, 
                              empty_cache = FALSE,
                              GTEX = TRUE)



create_count_data_frames("~/data/eiad_rse/rse/",
                         emeta %>% filter(data_location == 'local') %>% pull(study_accession) %>% unique(),
                         "local_TPM",
                         "long_local_TPM",
                         "local_count",
                         "long_local_count",
                         "local_mapping_information",
                         emeta, 
                         empty_cache = FALSE, 
                         GTEX = FALSE)




# input long TPM
recount3_TPM <- vroom::vroom("gene_counts/long_recount3_TPM.csv.gz") 
gtex_TPM <- vroom::vroom("gene_counts/long_GTEX_TPM.csv.gz")
local_TPM <- vroom::vroom("gene_counts/long_local_TPM.csv.gz") 

# input count matrices
recount3_counts <- vroom::vroom("gene_counts/recount3_count.csv.gz") %>% data.frame()
gtex_counts <- vroom::vroom("gene_counts/GTEX_count.csv.gz") %>% data.frame()
local_counts <- vroom::vroom("gene_counts/local_count.csv.gz")  %>% data.frame()

# Calling aggregated metadata
recount3_metadata <- vroom::vroom("mapping_data/recount3_mapping_information.csv.gz")
gtex_metadata <- vroom::vroom("mapping_data/gtex_mapping_information.csv.gz")
local_metadata <- vroom::vroom("mapping_data/local_mapping_information.csv.gz")
r3metadata <- bind_rows(recount3_metadata, gtex_metadata, local_metadata )

# Combine long TPM
gene_TPM <- bind_rows(recount3_TPM,
                      gtex_TPM,
                      local_TPM) %>% 
  unique()

# make matrix TPM
mat_TPM <- gene_TPM %>% unique() %>% pivot_wider(values_from = value, names_from = sample_accession)

# make matrix counts
mat <- cbind((recount3_counts %>% data.frame())[,2:ncol(recount3_counts)],
      (local_counts %>% data.frame())[,2:ncol(local_counts)],
      (gtex_counts %>% data.frame())[,2:ncol(gtex_counts)])
row.names(mat) <- recount3_counts[,1]

### Write files ---------
write_csv(gene_TPM, "gene_counts/gene_TPM.csv.gz", progress = TRUE)
write_csv(mat_TPM, "gene_counts/gene_TPM_matrix.csv.gz", progress = TRUE)
write_csv(mat, "gene_counts/gene_counts_matrix.csv.gz", progress = TRUE)
write_csv(r3metadata, "mapping_data/recount_metadata.csv.gz")
