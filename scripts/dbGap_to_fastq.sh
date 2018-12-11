#!/bin/bash

# Exports fastq from dbGaP accession SRA

cd /data/mcgaugheyd/dbGaP/11588

module load sratoolkit/2.9.2

SRR=$1
output_path=$2

fastq-dump $SRR --gzip --split-3 --readids --skip-technical

mv $SRR_\*.fastq.gz $output_path
