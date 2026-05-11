"""
Utilities to persist and collect TRMF/LATC simulation results.
Usage example:
    python code/simulation_results.py \
        --algo trmf \
        --csv output/simulation/overall_rel_mrae_margin_trmf.csv
    python code/simulation_results.py \
        --algo latc \
        --csv output/simulation/overall_rel_mrae_margin_latc.csv
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULT_ROOT = REPO_ROOT / "output" / "simulation"


@dataclass(frozen=True)
class ResultIndex:
    algo: str
    pdt: float
    repl: int

    def filename(self) -> str:
        prob = f"{self.pdt:.3f}".rstrip("0").rstrip(".")
        return f"pdt{prob}_repl{self.repl}.json"

    def path(self) -> Path:
        return RESULT_ROOT / self.algo / self.filename()


def save_result(algo: str, pdt: float, repl: int, fit: dict) -> Path:
    result = ResultIndex(algo=algo, pdt=pdt, repl=repl)
    payload = {
        "algo": algo,
        "pdt": pdt,
        "repl": repl,
        "metrics": {
            "mrae": fit["mrae"],
            "mrae_w": fit["mrae_w"],
            "mrae_o": fit["mrae_o"],
            "rmse": fit["rmse"],
        },
    }
    path = result.path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path


def _iter_payloads(algo: str | None) -> Iterable[dict]:
    if algo:
        roots = [RESULT_ROOT / algo]
    else:
        roots = [path for path in RESULT_ROOT.glob("*") if path.is_dir()]
    for root in roots:
        for file in sorted(root.glob("*.json")):
            yield json.loads(file.read_text())


def collect_dataframe(algo: str | None = None):
    import pandas as pd  # Delay import so training jobs do not need pandas
    records = []
    for payload in _iter_payloads(algo):
        metrics = payload["metrics"]
        base = {
            "method": payload["algo"],
            "pdt": payload["pdt"],
            "repl": payload["repl"],
        }
        records.append({**base, "type": "overall", "mrae": metrics["mrae"]})
        records.append({**base, "type": "w", "mrae": metrics["mrae_w"]})
        records.append({**base, "type": "o", "mrae": metrics["mrae_o"]})
    return pd.DataFrame.from_records(records)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect simulation results into a DataFrame.")
    parser.add_argument("--algo", choices=["trmf", "latc"], default=None, help="Filter by algorithm.")
    parser.add_argument("--csv", type=Path, default=None, help="Optional path to write the DataFrame as CSV.")
    parser.add_argument("--head", type=int, default=0, help="Print the first N rows instead of the full frame.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    df = collect_dataframe(args.algo)
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.csv, index=False)
    if args.head > 0:
        print(df.head(args.head))
    else:
        print(df)


if __name__ == "__main__":
    main()
