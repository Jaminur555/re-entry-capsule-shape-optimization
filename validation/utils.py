"""Shared helpers for the validation scripts: reference-data paths, pass-band
evaluation and pass/fail reporting (RESULTS.md + JSON run log)."""

from __future__ import annotations

import csv
import datetime as _dt
import json
from pathlib import Path

import numpy as np

VALIDATION_DIR = Path(__file__).resolve().parent
REFERENCE_DIR = VALIDATION_DIR / "reference_data"
RESULTS_DIR = VALIDATION_DIR.parent / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
LOGS_DIR = RESULTS_DIR / "logs"
RESULTS_MD = VALIDATION_DIR / "RESULTS.md"
RUN_LOG = LOGS_DIR / "validation_runs.json"


def load_reference(name):
    """Read a reference_data CSV into a list of dict rows (strings kept)."""
    with open(REFERENCE_DIR / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_numeric(name, cols):
    """Read a reference_data CSV; return {col: float ndarray} for `cols`.

    The literal 'inf' (Mach-infinity Cp curves) converts to np.inf.
    """
    rows = load_reference(name)
    out = {}
    for c in cols:
        out[c] = np.array([np.inf if r[c] == "inf" else float(r[c]) for r in rows])
    return out


def report(check_id, description, passed, details=""):
    """Append one verdict row to RESULTS.md and the JSON run log; print it.

    Parameters
    ----------
    check_id : str      e.g. 'V2-CD-M10'
    description : str   one line saying what was compared
    passed : bool
    details : str       numbers worth keeping (band, max error, ...)
    """
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    LOGS_DIR.mkdir(exist_ok=True)

    verdict = "PASS" if passed else "FAIL"
    first = not RESULTS_MD.exists()
    with open(RESULTS_MD, "a", newline="", encoding="utf-8") as f:
        if first:
            f.write("# Validation results\n\n| id | check | verdict | details |\n"
                    "|---|---|---|---|\n")
        f.write(f"| {check_id} | {description} | {verdict} | {details} |\n")

    log = []
    if RUN_LOG.exists():
        try:
            log = json.loads(RUN_LOG.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            log = []
    log.append({
        "id": check_id,
        "description": description,
        "verdict": verdict,
        "details": details,
        "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
    })
    RUN_LOG.write_text(json.dumps(log, indent=1), encoding="utf-8")

    print(f"[{verdict}] {check_id}: {description}"
          + (f" -- {details}" if details else ""))
    return passed


def within_band(values, reference, band, label=""):
    """Return (passed, max_abs_error). NaN-safe elementwise comparison."""
    values, reference = np.asarray(values, float), np.asarray(reference, float)
    err = np.abs(values - reference)
    finite = np.isfinite(err)
    if not finite.any():
        return False, float("nan")
    max_err = float(np.nanmax(err))
    return bool(np.all(err[finite] <= band)), max_err


def set_paper_style():
    """Shared seaborn styling for all validation figures."""
    try:
        import seaborn as sns
    except ImportError:
        return
    import matplotlib.pyplot as plt
    sns.set_theme(style="ticks", context="paper", font_scale=1.15,
                  palette="colorblind")
    plt.rcParams.update({
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.linestyle": "--", "grid.alpha": 0.35,
        "legend.frameon": False, "axes.titlesize": 11,
        "axes.labelsize": 11, "lines.linewidth": 1.8,
        "lines.markersize": 4.5, "savefig.dpi": 300,
    })


def save_fig(fig, name):
    """Save a styled figure into results/figures (300 dpi)."""
    fig.tight_layout()
    out = FIGURES_DIR / name
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  figure -> {out}")
