#!/bin/bash

# Exports fastq from dbGaP accession SRA

cd /data/mcgaugheyd/dbGaP/11588/sra

module load sratoolkit/2.9.1

SRR=$1
output_path=$2

fastq-dump $SRR --gzip --split-spot --readids

mv $SRR.fastq.gz $output_path
