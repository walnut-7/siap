#!/bin/bash

#SBATCH --job-name=simu_py
#SBATCH --account=your_account_name ######!!!! Please fill in your slurm account name
#SBATCH --partition=standard
#SBATCH --output=/home/%u/solar/log/%x/%A/%a.out
#SBATCH --error=/home/%u/solar/log/%x/%A/%a.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G

module load R/4.3.2-mkl
source "$HOME/miniconda3/etc/profile.d/conda.sh"  # adjust path to your Conda installation
conda activate siap

LOG_DIR="/home/${USER}/solar/log/${SLURM_JOB_NAME:-unknown}/${SLURM_JOB_ID:-unknown}"
mkdir -p "$LOG_DIR"

DEFAULT_ALGO="trmf"

ALGO="${ALGO:-$DEFAULT_ALGO}"

if [[ -n "${PARAM_FILE:-}" && ! -f "$PARAM_FILE" ]]; then
    echo "Parameter file ${PARAM_FILE} not found." >&2
    exit 1
fi

declare -a JOBS=()

if [[ -n "${PARAM_FILE:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then # array job mode
    CHUNK_INDEX="${CHUNK_INDEX:-0}"
    LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
    read -r PDT REPL <<< "$(
        awk -v chunk="$CHUNK_INDEX" -v line="$LINE_NUMBER" '
            (NF >= 3 && $1 == chunk) {
                count++
                if (count == line) { print $2, $3; exit }
            }
            (NF == 2 && chunk == 0) {
                legacy++
                if (legacy == line) { print $1, $2; exit }
            }
        ' "$PARAM_FILE"
    )"
    if [[ -z "$PDT" || -z "$REPL" ]]; then
        echo "No parameters found for chunk ${CHUNK_INDEX}, array index ${SLURM_ARRAY_TASK_ID}." >&2
        exit 1
    fi
    JOBS+=("$PDT $REPL")
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        JOB_NAME="chunk${CHUNK_INDEX}_task${SLURM_ARRAY_TASK_ID}"
        scontrol update jobid="${SLURM_JOB_ID}" jobname="$JOB_NAME" >/dev/null 2>&1 || true
    fi
elif [[ -n "${PARAM_FILE:-}" && -n "${CHUNK_INDEX:-}" ]]; then # chunked job mode
    mapfile -t JOBS < <(
        awk -v chunk="$CHUNK_INDEX" '
            (NF >= 3 && $1 == chunk) { print $2, $3 }
            (NF == 2 && chunk == 0) { print $1, $2 }
        ' "$PARAM_FILE"
    )
    if [[ ${#JOBS[@]} -eq 0 ]]; then
        echo "No parameters found for chunk ${CHUNK_INDEX} in ${PARAM_FILE}." >&2
        exit 1
    fi
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        JOB_NAME="chunk${CHUNK_INDEX}_seq"
        scontrol update jobid="${SLURM_JOB_ID}" jobname="$JOB_NAME" >/dev/null 2>&1 || true
    fi
elif [[ $# -ge 3 ]]; then # direct mode with algo specified
    ALGO="$1"
    PDT="$2"
    REPL="$3"
    JOBS+=("$PDT $REPL")
elif [[ $# -ge 2 ]]; then # direct mode with default algo
    PDT="$1"
    REPL="$2"
    JOBS+=("$PDT $REPL")
else # insufficient arguments
    echo "Usage: $0 [algo] PDT REPL" >&2
    exit 1
fi

if [[ ${#JOBS[@]} -eq 0 ]]; then
    echo "No jobs to run." >&2
    exit 1
fi

PY_SCRIPT="./code/simulation_${ALGO}_wrapper.py"
if [[ ! -f "$PY_SCRIPT" ]]; then
    echo "Algorithm script ${PY_SCRIPT} not found." >&2
    exit 1
fi

for JOB in "${JOBS[@]}"; do # process each (PDT, REPL) pair sequentially in each chunk
    read -r PDT REPL <<< "$JOB"
    JOB_PREFIX="simulation_py_${SLURM_JOB_ID:-local}_chunk${CHUNK_INDEX:-0}_p${PDT}_r${REPL}"
    JOB_OUT="${LOG_DIR}/${JOB_PREFIX}.out"
    JOB_ERR="${LOG_DIR}/${JOB_PREFIX}.err"

    {
        printf '[%s] Starting PDT=%s REPL=%s using %s\n' "$(date -Is)" "$PDT" "$REPL" "$PY_SCRIPT"
        Rscript ./code/simulation_repl.R "$PDT" "$REPL" 1 # prepare the missingness json file
        python "$PY_SCRIPT" "$PDT" "$REPL"
        printf '[%s] Finished PDT=%s REPL=%s\n' "$(date -Is)" "$PDT" "$REPL"
    } >"$JOB_OUT" 2>"$JOB_ERR"
done
