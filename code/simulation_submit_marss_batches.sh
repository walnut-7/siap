#!/bin/bash

#SBATCH --job-name=submit_marss
#SBATCH --account=your_account_name ######!!!! Please fill in your slurm account name
#SBATCH --partition=standard
#SBATCH --output=/home/%u/solar/log/simulation_submit_marss.out   
#SBATCH --error=/home/%u/solar/log/simulation_submit_marss.err 
#SBATCH --mail-type=END 
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
module load R/4.3.2-mkl

# Control whether to restrict to jobs with diff == T (default: 1 / true).
DIFF_FLAG=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --diff)
      DIFF_FLAG="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: sbatch code/simulation_submit_marss_batches.sh [--diff 0|1]"
      exit 1
      ;;
  esac
done

if [[ "$DIFF_FLAG" == "1" ]]; then
  USE_DIFF=TRUE
  MODE_LABEL="diff-only"
elif [[ "$DIFF_FLAG" == "0" ]]; then
  USE_DIFF=FALSE
  MODE_LABEL="all-pending"
else
  echo "Invalid value for --diff: $DIFF_FLAG (expected 0 or 1)"
  exit 1
fi

SNAPSHOT_DIR="code/tmp"
SNAPSHOT_PATH="$SNAPSHOT_DIR/simulation_marss_submit_ids_${MODE_LABEL}_$(date '+%Y%m%d_%H%M%S')_$$.txt"

mkdir -p "$SNAPSHOT_DIR"

cleanup() {
  if [[ -f "$SNAPSHOT_PATH" ]]; then
    rm -f "$SNAPSHOT_PATH"
    echo "[submit-marss][driver] cleanup | removed snapshot $SNAPSHOT_PATH"
  fi
}

trap cleanup EXIT

echo "[submit-marss][driver] start | time=$(date '+%Y-%m-%d %H:%M:%S') | mode=$MODE_LABEL"
echo "[submit-marss][driver] snapshot | path=$SNAPSHOT_PATH"
for i in {1..16}  
do
  echo "[submit-marss][driver] batch $i/16 | launching at $(date '+%Y-%m-%d %H:%M:%S')"
  Rscript code/simulation_submit_marss_batches.R "$i" "$USE_DIFF" "$SNAPSHOT_PATH"
  echo "[submit-marss][driver] batch $i/16 | launch command finished at $(date '+%Y-%m-%d %H:%M:%S')"
  sleep 2h
done
echo "[submit-marss][driver] done | time=$(date '+%Y-%m-%d %H:%M:%S')"
