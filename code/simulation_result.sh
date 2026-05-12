#!/bin/bash
#SBATCH --job-name=simu_result
#SBATCH --account=your_account_name ######!!!! Please fill in your slurm account name
#SBATCH --partition=standard
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mail-type=END 
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
module load R/4.3.2-mkl

usage() {
  cat <<EOF
Usage: sbatch simulation_result.sh --switches "<parts>" [--add] [--model <name>] [--siap] [--preprocess <T|F>] [--beta-update] [--gp] [--diff <T|F>] [--marss] [--diff <T|F>]
Examples:
  sbatch simulation_result.sh --switches "1 2 3 4 5"
  sbatch simulation_result.sh --switches "1 2" --add --model sia
  sbatch simulation_result.sh --switches "1 2 3" --siap --preprocess T # default
  sbatch simulation_result.sh --switches "1 2 3" --siap --preprocess F --beta-update
  sbatch simulation_result.sh --switches "1 2 3 4 5" --gp --diff TRUE # default
  sbatch simulation_result.sh --switches "1 2 3 4 5" --gp --diff FALSE --marss --diff TRUE
  sbatch simulation_result.sh --switches "1 2 3" --add --model marss  # run only MARSS append
EOF
}

normalize_bool() {
  local value_upper
  value_upper="$(echo "$1" | tr '[:lower:]' '[:upper:]')"
  case "$value_upper" in
    T|TRUE|1) echo "TRUE" ;;
    F|FALSE|0) echo "FALSE" ;;
    *) return 1 ;;
  esac
}

SWITCHES=""
ADD="0"
MODEL="0"
MARSS_ONLY="0"
BETA_ONESHOT="T"
SIAP_PREPROCESS="with"
GP_DIFF="TRUE"
MARSS_DIFF="TRUE"
CURRENT_METHOD=""

while [ $# -gt 0 ]; do
  case "$1" in
    --switches)
      SWITCHES="$2"; shift 2;;
    --add)
      ADD="1"; shift 1;;
    --model)
      MODEL="$2"; shift 2;;
    --siap)
      CURRENT_METHOD="siap"; shift 1;;
    --gp)
      CURRENT_METHOD="gp"; shift 1;;
    --marss)
      CURRENT_METHOD="marss"; shift 1;;
    --preprocess)
      if [ "$CURRENT_METHOD" != "siap" ]; then
        echo "Error: --preprocess must appear after --siap."; usage; exit 1
      fi
      PREPROCESS_NORMALIZED="$(normalize_bool "$2")" || {
        echo "Error: --preprocess must be one of T/F/TRUE/FALSE/1/0."; usage; exit 1;
      }
      if [ "$PREPROCESS_NORMALIZED" = "TRUE" ]; then
        SIAP_PREPROCESS="with"
      else
        SIAP_PREPROCESS="without"
      fi
      shift 2;;
    --beta-update)
      if [ "$CURRENT_METHOD" != "siap" ]; then
        echo "Error: --beta-update must appear after --siap."; usage; exit 1
      fi
      BETA_ONESHOT="F"; shift 1;;
    --diff)
      DIFF_NORMALIZED="$(normalize_bool "$2")" || {
        echo "Error: --diff must be one of T/F/TRUE/FALSE/1/0."; usage; exit 1;
      }
      if [ "$CURRENT_METHOD" = "gp" ]; then
        GP_DIFF="$DIFF_NORMALIZED"
      elif [ "$CURRENT_METHOD" = "marss" ]; then
        MARSS_DIFF="$DIFF_NORMALIZED"
      else
        echo "Error: --diff must appear after --gp or --marss."; usage; exit 1
      fi
      shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      # Allow legacy positional first argument for switches
      if [ -z "$SWITCHES" ]; then
        SWITCHES="$1"; shift 1;
      else
        echo "Unknown option: $1"; usage; exit 1
      fi;;
  esac
done

if [ -z "$SWITCHES" ]; then
  echo "Error: --switches is required."; usage; exit 1
fi

# Route stdout/stderr to files that include switches
SWITCHES_TAG="$(echo "$SWITCHES" | tr ' ' '_' | tr -cd '[:alnum:]_.-')"
OUT_FILE="/home/${USER}/solar/log/simulation_result_${SWITCHES_TAG}.out"
ERR_FILE="/home/${USER}/solar/log/simulation_result_${SWITCHES_TAG}.err"
exec 1>"$OUT_FILE" 2>"$ERR_FILE"

# `--model marss` means run only the MARSS append path.
if [ "$MODEL" = "marss" ]; then
  MARSS_ONLY="1"
fi

echo "Running script with switches: $SWITCHES"
if [ "$MARSS_ONLY" = "1" ]; then
  echo "Mode: append MARSS outputs only"
else
  echo "Mode: run base results, then append MARSS outputs (default)"
  echo "Add a model's result (add: T; 0: F): $ADD"
  echo "The model to be added: $MODEL"
fi
echo "SIAP preprocess: $SIAP_PREPROCESS"
echo "SIAP beta-oneshot: $BETA_ONESHOT"
echo "GP diff: $GP_DIFF"
echo "MARSS diff: $MARSS_DIFF"

cd /home/yuxuank/solar

if [ "$MARSS_ONLY" = "1" ]; then
  Rscript code/simulation_result_add_marss.R "$SWITCHES" "$BETA_ONESHOT" "$SIAP_PREPROCESS" "$MARSS_DIFF" "$GP_DIFF"
else
  Rscript code/simulation_result.R "$ADD" "$MODEL" "$SWITCHES" "$BETA_ONESHOT" "$SIAP_PREPROCESS" "$GP_DIFF" "$MARSS_DIFF"
  Rscript code/simulation_result_add_marss.R "$SWITCHES" "$BETA_ONESHOT" "$SIAP_PREPROCESS" "$MARSS_DIFF" "$GP_DIFF"
fi
