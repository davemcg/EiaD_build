# Eye in a Disk (EiaD)

EiaD is the sqlite database at the core of eyeIntegration.nei.nih.gov

For 2023 (versions >= 2.0), we have updated the "backend" of EiaD in several significant ways:

1. Simplify the salmon-based quantification to better enable integration of our dataset with outside resources
2. Added many new samples and studies
3. More granular metadata schema that has five major categories:
  - Tissue (e.g. Retina)
  - Sub_Tissue (e.g. Macula)
  - Source (e.g. tissue or iPSC)
  - Age (e.g. fetal or adult)
  - Perturbation (e.g. None or AMD)
4. Added ML-based sex labels
5. Built a recount3 based quantification pipeline (http://github.com/davemcg/Snakerail) to enable base pair level coverage information
6. Used ML based approach to identify sample outliers for QC
7. Summarized cell type level gene tables imported from our plae.nei.nih.gov resource


# Workflow
1. Snakerail (http://github.com/davemcg/Snakerail) wraps the pump (output) and unify (RSE) steps in monorail (https://github.com/langmead-lab/monorail-external)
2. Snakefile runs a salmon-based quant to generate gene and transcript level counts
3. There are four "hand" steps that generate the 2023 EiaD datasets that are run in this order:
  - scripts/pull_scEiaD.R and scripts/build_eiad_2023_plae.R 
  - scripts/metamoRph_label.R and scripts/identify_outlier_samples.Rmd
  - scripts/pca_workup_data_prep.R
  - scripts/build_eiad_2023_bulk.R
