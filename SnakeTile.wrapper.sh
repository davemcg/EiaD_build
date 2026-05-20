#!/bin/bash

# run in the data folder for this project
# on biowulf2:
# /data/mcgaugheyd/projects/nei/brooks/oca_rna-seq

mkdir -p 00log

#module load snakemake || exit 1


source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
mamba activate tiledb_py

# Execute using the Biowulf-optimized Snakemake 8+ profile
snakemake -s /home/mcgaugheyd/git/EiaD_build/SnakeTile \
  --profile biowulf_v8 \
  --jobs 500 \
  --local-cores 2 \
  --default-resources mem_mb=8000 runtime=60 \
  --latency-wait 120 \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --keep-going \
  --restart-times 0 \
  --groups ingest_sample=ingest_batch \
  --group-components ingest_batch=10

