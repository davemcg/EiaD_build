library(tidyverse)

# List of sample names
sample_names <- c("SRP389394", "DRP006678", "ERP126780", "SRP012682", 
                  "SRP016140", "SRP064956", "SRP080886", "SRP094572", 
                  "SRP098761", "SRP107937", "SRP136195", "SRP151763", 
                  "SRP186882", "SRP287234", "SRP331221", "SRP356976")
# List to store subset dataframes
subset_df_list <- list()
# List of column names you want to keep
columns_to_keep <- c("Age", "Sex", "Ethnicity", "AGE", "SEX", "ETHNICITY", "age", "sex", "ethnicity", "Run")
# Find common column names
common_columns <- character(0)
for (sample_name in sample_names) {
  df <- read.csv(paste("~/git/EiaD_build_pp/data/metadata_manual/", sample_name, ".txt", sep = ""))
  common_columns <- union(common_columns, intersect(columns_to_keep, colnames(df)))
}

# Loop through each sample name
for (sample_name in sample_names) {
  # Read the dataframe for the current sample
  df <- read.csv(paste("~/git/EiaD_build_pp/data/metadata_manual/", sample_name, ".txt", sep = ""))
  # Filter the list to keep only the common column names that exist in the current dataframe
  existing_columns_to_keep <- intersect(common_columns, colnames(df))
  # Create a new dataframe with unified structure and missing columns filled with NA
  subset_df <- data.frame(matrix(NA, nrow = nrow(df), ncol = length(common_columns)))
  colnames(subset_df) <- common_columns
  subset_df[existing_columns_to_keep] <- df[existing_columns_to_keep]
  # Save the subset dataframe to the list
  subset_df_list[[sample_name]] <- subset_df
}

# Combine all dataframes
identifier_metadata <- do.call(rbind, subset_df_list) %>% 
  rename(c(Ethnicity = ETHNICITY, Sex = sex, run_accession = Run, age_number = Age))
# Ensure uniform age numbers, ethnicity, and sex values
identifier_metadata$age_number <- gsub("y", "", identifier_metadata$age_number)
identifier_metadata$age_number <- gsub("gestation age 16-18wk", "16-18 Weeks (Gestation Age)", 
                                       identifier_metadata$age_number)
identifier_metadata$Sex <- gsub("not applicable", NA, 
                                identifier_metadata$Sex)
identifier_metadata$Ethnicity <- gsub("", NA, 
                                      identifier_metadata$Ethnicity)

# Combine with our eyeIntegration data
ei_23_metadata <- read_csv("data/eyeIntegration22_meta_2023_03_03.csv.gz")
# SRP064956 are all adult samples - Identified using the SRA run selector
ei_23_metadata <- within(ei_23_metadata, Age[study_accession == "SRP064956"] <- 'Adult')

# Join
ei_23_metadata <- ei_23_metadata %>% left_join(identifier_metadata) %>% select(-Ethnicity) # No observations
# Remove adult classification from the age_number variable
ei_23_metadata <- within(ei_23_metadata, age_number[study_accession == "SRP064956"] <- NA)

# Write CSV
write.csv(ei_23_metadata, file = "~/git/EiaD_build_pp/data/eyeIntegration23_meta_2023_08_28.csv",
          row.names = FALSE)
# Compress the CSV file using the gzip command-line tool
system("gzip ~/git/EiaD_build_pp/data/eyeIntegration23_meta_2023_08_28.csv")
