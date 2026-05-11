"""
Submit TRMF and LATC simulation jobs to Slurm via ``simulation_py.sh``.

This driver builds one or more chunk jobs so that the probabilities of dropout
(``pdt``), replicate indices (``repl``), and algorithm choice propagate
consistently from Python → shell → Slurm → the simulation scripts. Each chunk
runs sequentially on the same node to reuse the allocated resources.

Usage example:
python code/simulation_driver.py \
  --pdt 0.1 0.3 0.5\
  --repl $(seq 1 100) \
  --algo trmf \
  --chunk-size 20 \
  --dry-run

python code/simulation_driver.py \
  --pdt 0.1 0.3 0.5\
  --repl $(seq 1 100) \
  --algo latc \
  --chunk-size 20 \
  --array \
  --dry-run
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parent
SLURM_SCRIPT = REPO_ROOT / "simulation_py.sh"
ALGORITHMS = ("trmf", "latc")
PARAM_DIR = REPO_ROOT / ".slurm_params"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Submit jobs as one or more Slurm chunk jobs. "
            "Each chunk runs its (pdt, repl) combinations sequentially."
        )
    )
    parser.add_argument(
        "--pdt",
        type=float,
        nargs="+",
        required=True,
        help=("Missingness probabilities to test (written to the parameter file)."),
    )
    parser.add_argument(
        "--repl",
        type=int,
        nargs="+",
        required=True,
        help=("Replicate IDs to run for each probability (written to the parameter file)."),
    )
    parser.add_argument(
        "--algo",
        choices=ALGORITHMS,
        default="trmf",
        help="Algorithm to run inside the Slurm job.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=0,
        help=(
            "Maximum number of jobs per Slurm chunk job. "
            "Use 0 (default) to submit everything as one chunk."
        ),
    )
    parser.add_argument(
        "--array",
        action="store_true",
        help=(
            "Submit each chunk as a Slurm array job so the (pdt, repl) pairs "
            "run in parallel across array tasks."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands instead of submitting them to Slurm.",
    )
    return parser.parse_args(argv)


def chunk_pairs(
    pairs: Sequence[tuple[float, int]],
    chunk_size: int,
) -> list[list[tuple[float, int]]]:
    """Split the parameter combinations into sequential chunk jobs."""
    if chunk_size <= 0:
        chunk_size = len(pairs)
    if chunk_size <= 0:
        return []
    return [
        list(pairs[idx : idx + chunk_size])
        for idx in range(0, len(pairs), chunk_size)
    ]


def write_param_file(chunks: Sequence[Sequence[tuple[float, int]]]) -> Path:
    """Persist chunked (pdt, repl) rows for retrieval inside Slurm jobs."""
    PARAM_DIR.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = f"params-{timestamp}-{os.getpid()}.tsv"
    path = (PARAM_DIR / filename).resolve()
    with path.open("w", encoding="utf-8") as fh:
        for chunk_idx, chunk in enumerate(chunks):
            for pdt, repl in chunk:
                fh.write(f"{chunk_idx}\t{pdt}\t{repl}\n")
    return path


def submit_jobs(
    algo: str,
    pdts: Iterable[float],
    repls: Iterable[int],
    chunk_size: int,
    use_array: bool,
    dry_run: bool = False,
) -> None:
    pairs = list(product(pdts, repls))
    if not pairs:
        raise ValueError("No jobs to submit: provide at least one pdt and repl.")

    chunks = chunk_pairs(pairs, chunk_size)
    if not chunks:
        raise ValueError("Chunking produced no jobs; check --chunk-size.")
    param_file = write_param_file(chunks)

    for chunk_idx, chunk in enumerate(chunks):
        chunk_len = len(chunk)
        export_vars = (
            f"ALL,ALGO={algo},PARAM_FILE={param_file},"
            f"CHUNK_INDEX={chunk_idx},CHUNK_TASK_COUNT={chunk_len}"
        )
        cmd = ["sbatch"]
        if use_array:
            array_spec = "0" if chunk_len <= 1 else f"0-{chunk_len - 1}"
            cmd.append(f"--array={array_spec}")
        cmd.extend([
            f"--export={export_vars}",
            str(SLURM_SCRIPT),
        ])

        quoted = " ".join(shlex.quote(part) for part in cmd)
        if dry_run:
            print(
                "Dry run; would submit chunk job with:\n"
                f"{quoted}\nParameter file: {param_file}\nChunk size: {chunk_len}"
            )
            continue

        print(
            (
                f"Submitting chunk {chunk_idx} as an array with {chunk_len} tasks ..."
                if use_array
                else f"Submitting chunk {chunk_idx} with {chunk_len} sequential jobs ..."
            ),
            flush=True,
        )
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                "sbatch submission failed:\n"
                f"Command: {quoted}\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )
        print(result.stdout.strip() or f"Submitted {cmd}")


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv) # read arguments from command line
    if not SLURM_SCRIPT.exists():
        raise FileNotFoundError(
            f"Cannot find Slurm wrapper {SLURM_SCRIPT}. "
            "Confirm you are running the driver from the `code/` directory."
        )
    submit_jobs(
        args.algo,
        args.pdt,
        args.repl,
        args.chunk_size,
        args.array,
        args.dry_run,
    )


if __name__ == "__main__":
    main()
