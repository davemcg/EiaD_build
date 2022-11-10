### All data files necessary for analysis can be found here (files must be uncompressed prior to use):
# https://hpc.nih.gov/~parikhpp/sex_inference/insert_file_name_here.csv.gz

# Set working directory
setwd("/Users/parikhpp/git/sex_inference")
# Import libraries and data
library(tidyverse)
regression_data <- read.csv("data/inference_regression_data.csv")

# Plot data

# X alignment
ggplot(regression_data, aes(x = sex, y = recount_qc.aligned_reads..chrx, color = sex)) +
  geom_boxplot() + theme_bw()
# Y alignment
ggplot(regression_data, aes(x = sex, y = recount_qc.aligned_reads..chry, color = sex)) +
  geom_boxplot() + theme_bw()
# X/Y ratio
ggplot(regression_data, 
       aes(x = sex, y = (recount_qc.aligned_reads..chrx/recount_qc.aligned_reads..chry), color = sex)) +
  geom_boxplot() + theme_bw()

# Which samples are overlapping?

regression_data_xy <- mutate(regression_data, 
                                x_y_ratio = (recount_qc.aligned_reads..chrx/recount_qc.aligned_reads..chry))
regression_data_xy %>% filter(sex == "female") %>% filter(x_y_ratio < 100)

# SRR2895372 - SRP064956
# SRR2895378 - SRP064956
# SRR2895380 - SRP064956
# SRR2895381 - SRP064956
# SRR4252548 - SRP090027
# SRR5535471 - SRP090027
# SRR5535778 - SRP090027
# SRR5535934 - SRP090027

# The above female samples are in the same tail range as 
# the tail for our male x to y ratio plot for percentage alignment

# Investigating whether these are the same outliers as are depicted in the female plot within our
# y alignment data

ggplot(regression_data, aes(x = sex, y = recount_qc.aligned_reads..chry, color = sex)) +
  geom_boxplot() + theme_bw()
# Find median
regression_data %>% filter(sex == "female") %>% pull(recount_qc.aligned_reads..chry) %>% median()
# Samples in tail
regression_data %>% filter(sex == "female") %>% filter(recount_qc.aligned_reads..chry > 0.01)

# The same 8 samples, including an additional ninth overlap between male and female samples

# SRR7461319 - SRP151763
# SRR2895372
# SRR2895378
# SRR2895380
# SRR2895381
# SRR4252548
# SRR5535471
# SRR5535778
# SRR5535934

# All samples were confirmed labeled as female from the Sequence Read Archive

# Through prediction models, we will aim to understand whether or not the provided data
# verifies the labeled sex