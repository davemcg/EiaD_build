#!/bin/bash
# Submit one SLURM job per species x assay.
# Run from the project dir (e.g. EiaD_2026/) so relative paths resolve.
#
# usage:
#   run_build.sh --meta META --out_dir OUT [options]
#
# required:
#   --meta       FILE   EiaD metadata tsv / tsv.gz
#   --out_dir    DIR    output directory
# optional:
#   --config     FILE   pipeline config.yaml            [config.yaml]
#   --species    LIST   space-separated, quoted         ["human rat"]
#   --assays     LIST   any of: gene exon sj, quoted    ["gene exon sj"]
#   --version    STR    version stamp in output names   [today, YYYYMMDD]
#   --salmon_dir DIR    salmon quant dir (gene assay)   [salmon_quant]
#   --fc_dir     DIR    featureCounts dir (exon assay)  [featurecounts]
#   --star_dir   DIR    STAR align dir (sj assay)       [STAR_align]
#   --dry_run           print the sbatch commands, submit nothing
#
# examples:
#   run_build.sh --meta META.tsv.gz --out_dir counts
#   run_build.sh --meta META.tsv.gz --out_dir counts --species rat --assays "gene exon"
#   run_build.sh --meta META.tsv.gz --out_dir counts --config /path/config.yaml \
#       --star_dir /data/OGVFB_BG/EiaD_2026/STAR_align --version 20261120

set -euo pipefail

META=""; OUT=""; CONFIG="config.yaml"; SPECIES="human rat"; ASSAYS="gene exon sj"
VER=$(date +%Y%m%d)
SALMON_DIR="salmon_quant"; FC_DIR="featurecounts"; STAR_DIR="STAR_align"
DRY=0

usage() { sed -n '2,/^set -euo/p' "${BASH_SOURCE[0]}" | sed '$d' | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --meta)       META=$2; shift 2 ;;
    --out_dir)    OUT=$2; shift 2 ;;
    --config)     CONFIG=$2; shift 2 ;;
    --species)    SPECIES=$2; shift 2 ;;
    --assays)     ASSAYS=$2; shift 2 ;;
    --version)    VER=$2; shift 2 ;;
    --salmon_dir) SALMON_DIR=$2; shift 2 ;;
    --fc_dir)     FC_DIR=$2; shift 2 ;;
    --star_dir)   STAR_DIR=$2; shift 2 ;;
    --dry_run)    DRY=1; shift ;;
    -h|--help)    usage 0 ;;
    *) echo "unknown argument: $1" >&2; usage 1 ;;
  esac
done

[ -n "$META" ] && [ -n "$OUT" ] || { echo "--meta and --out_dir are required" >&2; usage 1; }
[ -f "$CONFIG" ] || { echo "config not found: $CONFIG" >&2; exit 1; }
[ -f "$META" ]   || { echo "metadata not found: $META" >&2; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p 00log "$OUT"

declare -A SCRIPT=( [gene]=build_tx_and_gene_h5_se.R [exon]=build_exon_h5_se.R [sj]=build_sj_h5_se.R )
declare -A INDIR=(  [gene]=$SALMON_DIR [exon]=$FC_DIR [sj]=$STAR_DIR )
# placeholder resources -- tune from `sacct -j <id> --format=MaxRSS,Elapsed` after a first run
declare -A RES=(
  [gene]="--cpus-per-task=8  --mem=64g  --time=4:00:00"
  [exon]="--cpus-per-task=8  --mem=96g  --time=4:00:00"
  [sj]="--cpus-per-task=16 --mem=300g --time=12:00:00" )

for a in $ASSAYS; do
  [ -n "${SCRIPT[$a]:-}" ] || { echo "unknown assay '$a' (gene|exon|sj)" >&2; exit 1; }
  [ -d "${INDIR[$a]}" ]   || echo "WARNING: input dir for '$a' not found from here: ${INDIR[$a]}" >&2
done

for sp in $SPECIES; do
  for a in $ASSAYS; do
    cmd=(sbatch ${RES[$a]}
      --job-name="eiad_${a}_${sp}"
      --output="00log/${a}_${sp}_%j.out"
      --wrap="module load R; Rscript ${SCRIPT_DIR}/${SCRIPT[$a]} --species ${sp} --config ${CONFIG} --meta ${META} --input_dir ${INDIR[$a]} --out_dir ${OUT} --version ${VER}")
    if [ "$DRY" = 1 ]; then printf '%q ' "${cmd[@]}"; echo; else "${cmd[@]}"; fi
  done
done
