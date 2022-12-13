# Eye in a Disk (EiaD)

EiaD is the sqlite database at the core of eyeIntegration.nei.nih.gov

For 2022 (2023?), we have updated the "backend" of EiaD by switching from our homebuid quantification platform (Salmon based) to the recount3 / monorail pipeline. 


# Rough workflow
1. Snakerail (http://github.com/davemcg/Snakerail) wraps the pump (output) and unify (RSE) steps in monorail (https://github.com/langmead-lab/monorail-external)
2. In this repo, scripts/unifier_to_counts.R makes the count, TPM, and metadata matrices
3. scripts/build_eiad.R creates the 2022 sqlite files for eyeIntegration 2022
4. scripts/pca_workup_data_prep.R creates lightly processed TPM files and metadata for PCA / sample dist analysis
5. scripts/pca_workup.Rmd is a overview of QC steps to identify outliers
