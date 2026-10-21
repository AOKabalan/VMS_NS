"""Verification helper: run scripts/run_convergence_study.py and read its convergence table."""

import csv
import subprocess
import sys
from pathlib import Path

STUDY = Path(__file__).resolve().parents[2] / "scripts" / "run_convergence_study.py"


def run_study(folder: Path, base: str, name: str, *levels: str) -> list[dict[str, str]]:
    """Run a study of the base configuration text over the levels; return the table rows."""
    path = folder / f"{name}.toml"
    path.write_text(base, encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(STUDY), str(path), *levels],
        cwd=folder,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-4000:]
    case = base.split('"')[1]
    table = folder / "results" / case / f"{name}_study" / "convergence.csv"
    with table.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def finest_rate(table: list[dict[str, str]], name: str) -> float:
    """Rate of a diagnostic between the two finest levels."""
    return float(table[-1][f"{name}_rate"])
