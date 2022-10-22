library(tidyverse)
library(tibble)
library(recount)
library(recount3)
library(dtplyr)
library(data.table)

source("~/git/EiaD_build/scripts/create_count_data_function.R")

# Importing Metadata
eyeIntegration22 <- read_csv("https://hpc.nih.gov/~parikhpp/EiaD/2022_metadata.csv", col_types = cols(...1 = col_skip()))

create_count_data_frames("http://duffel.rail.bio/recount3/",
                         c("SRP002881", "SRP011895", "SRP012585", "SRP015336", "SRP016140", "SRP034875",
                           "SRP035641", "SRP045639", "SRP053034", "SRP055101", "SRP061670", "SRP062870",
                           "SRP064956", "SRP070148", "SRP075990", "SRP079002", "SRP080002", "SRP090226",
                           "SRP091605", "SRP091675", "SRP092413", "SRP093877", "SRP094572", "SRP097696",
                           "SRP098761", "SRP105756", "SRP106457", "SRP108292", "SRP110135", "SRP111145",
                           "SRP115908", "SRP117613", "SRP119291", "SRP119766", "SRP151763", "SRP159246",
                           "ERP022243", "SRP090027", "SRP090040", "SRP186009", "SRP156453", "SRP107937"),
                         "recount3_transformed_counts",
                         "aggregated_recount3_transformed_counts",
                         "recount3_mapping_information",
                         eyeIntegration22)


create_gtex_count_data_frames("http://duffel.rail.bio/recount3/",
                              c("ADIPOSE_TISSUE", "MUSCLE", "BLOOD_VESSEL", "HEART", "OVARY", "UTERUS",
                                "VAGINA", "BREAST", "SKIN", "SALIVARY_GLAND", "BRAIN", "ADRENAL_GLAND",
                                "THYROID", "LUNG", "SPLEEN", "PANCREAS", "ESOPHAGUS", "STOMACH", "COLON",
                                "SMALL_INTESTINE", "PROSTATE", "TESTIS", "NERVE", "PITUITARY", "BLOOD",
                                "LIVER", "KIDNEY",   "CERVIX_UTERI", "FALLOPIAN_TUBE", "BLADDER", "BONE_MARROW"),
                              "gtex_transformed_counts",
                              "aggregated_gtex_transformed_counts",
                              eyeIntegration22)


create_count_data_frames("~/data/eiad_rse/",
                         c("SRP018405", "SRP326606", "SRP323408", "SRP329409",
                           "SRP331221", "SRP090027", "SRP287234", "SRP070938",
                           "SRP090040", "SRP310948", "SRP287152", "SRP056957",
                           "SRP288670", "ERP126780", "SRP257684", "SRP080886",
                           "SRP300190", "SRP255891", "SRP273695"),
                         "local_transformed_counts",
                         "aggregated_local_transformed_counts",
                         "local_mapping_information",
                         eyeIntegration22)



# Calling aggregated data
aggregated_recount3_transformed_counts <- vroom::vroom("gene_counts/aggregated_recount3_transformed_counts.csv")
aggregated_gtex_transformed_counts <- vroom::vroom("gene_counts/aggregated_gtex_transformed_counts.csv")
aggregated_local_data_additions_transformed_counts <- vroom::vroom("gene_counts/aggregated_local_transformed_counts.csv")

# Calling aggregated metadata
aggregated_recount3_metadata <- vroom::vroom("mapping_data/recount3_mapping_information.csv")
aggregated_gtex_transformed_counts <- vroom::vroom("mapping_data/gtex_mapping_information.csv")
aggregated_local_data_additions_transformed_counts <- vroom::vroom("mapping_data/local_mapping_information.csv")
r3metadata <- bind_rows(aggregated_recount3_metadata, aggregated_gtex_transformed_counts, aggregated_local_data_additions_transformed_counts )

# Combine aggregated data
gene_counts <- bind_rows(aggregated_recount3_transformed_counts,
                         aggregated_gtex_transformed_counts,
                         aggregated_local_data_additions_transformed_counts)

mat <- gene_counts %>% select(-run_accession) %>% unique() %>% pivot_wider(values_from = value, names_from = sample_accession)

### Write final gene_counts csv file -----
write_csv(gene_counts, "gene_counts/gene_counts.csv.gz", progress = TRUE)
write_csv(mat, "gene_counts/gene_counts_matrix.csv.gz", progress = TRUE)
