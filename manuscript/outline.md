# Abstract
  - write this last

# Introduction
  - brief discussion of how the eye works
  - role of genetics in eye development and human disease
  - what RNA-seq does
  - discussion of how differential expression used to identify 'marker' genes
    - crx?
    - nrl?
    - rpe (miller bharti papers)
  - eyeIntegration contains the broadest set of human eye RNA-seq data and dozens of other human tissues to compare against
    - also some scRNA-seq
  - other resources that aren't as useful
    - GTEx
      - no eye
    - http://ascot.cs.jhu.edu
      - few human samples
      - rod focused
    - http://retina.tigem.it
      - retina data only
    - at least one more I can't remember (had a bunch of mouse lens data)
  - brief mention of data reproducilblity (snakemake)
  - eyeIntegration.nei.nih.gov is a performant web app to serve the data to all users in a variety of ways

# Methods

  - search SRA / GEO for human eye-tissue RNA-seq
  - discussion of Snakemake pipeline
    - alignment free quantification
    - QC
    - etc. 
  - web app (written in R and Shiny, can be used online or deployed on local computer)
  
# Results
   - largest collection of processed and QC'ed human eye RNA-seq data
   - discuss improvements made since HMG paper
     - more samples
     - more tissues
     - snakemake
     - not just protein-coding
  - brief discussion of pan-tissue boxplot / fold change / heatmaps
  - scRNA-seq data
  - differential expression

# Conclusions

  
