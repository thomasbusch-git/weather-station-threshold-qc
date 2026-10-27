#!/usr/bin/env python3
"""
================================================================================
 Script:      Graph_Threshold_Errors.py
 Version:     1.0.0
 Description: Time-series plots for threshold (range) violations only.
              Reads the same station CSVs as Threshold_Errors.py and highlights
              out-of-range periods in red. No other check types (IQR, rate-of-
              change, etc.) are included, so you see only threshold flags.

 Author:      Thomas Busch
 Date:        2026-06-18

 Usage:
   # Save PNGs to ./threshold_graphs/
   python Graph_Threshold_Errors.py -i stationdata/ -o threshold_graphs/

   # Single station
   python Graph_Threshold_Errors.py -i stationdata/ -o threshold_graphs/ --station 0115

   # Single variable across all stations
   python Graph_Threshold_Errors.py -i stationdata/ -o threshold_graphs/ --variable tsoil_1

   # Interactive window (requires display)
   python Graph_Threshold_Errors.py -i stationdata/ --station 0115 --show

 Requirements:
   pip install polars matplotlib
================================================================================
"""

import argparse
import sys
from pathlib import Path

# ── resolve backend before pyplot import ─────────────────────────────────────
_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument("--show", action="store_true")
_show_mode = _pre.parse_known_args()[0].show

import matplotlib
if not _show_mode:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import polars as pl

# ── import shared config from sibling script ─────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
try:
    from Threshold_Errors import VARIABLES, THRESHOLDS, TIMESTAMP_ALIASES, STATION_ID_ALIASES
except ImportError:
    print("ERROR: Threshold_Errors.py must be in the same directory as this script.")
    sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────────
# STYLE
# ─────────────────────────────────────────────────────────────────────────────
LINE_COLOR   = "#3b2a8a"   # indigo — matches Threshold_Errors theme
FILL_COLOR   = "#ff6b6b"   # red shading for violations
THRESH_COLOR = "#c53030"   # darker red for threshold lines
GRID_COLOR   = "#e2e8f0"
SUBPLOT_BG   = "#fafbff"
COLS         = 3           # subplots per row


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def detect_delimiter(path: Path) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        line = f.readline()
    return "," if line.count(",") >= line.count("\t") else "\t"


def find_ts_col(columns):
    for alias in TIMESTAMP_ALIASES:
        if alias in columns:
            return alias
    for c in columns:
        if "timestamp" in c.lower():
            return c
    return None


def find_sid_col(columns):
    for alias in STATION_ID_ALIASES:
        if alias in columns:
            return alias
    for c in columns:
        if "station" in c.lower() and "id" in c.lower():
            return c
    return None


def find_var_col(columns, norm_name: str):
    """Return the first matching column name for a normalised variable key."""
    for alias in VARIABLES.get(norm_name, []):
        if alias in columns:
            return alias
    return None


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_station(station_id: str, input_dir: Path) -> pl.DataFrame | None:
    """
    Load and combine all CSV rows for station_id.
    Handles both named files (0115.csv) and multi-station files.
    Returns a DataFrame sorted by timestamp with a '_ts' Datetime column,
    or None if no data found.
    """
    frames = []

    # Fast path: if a per-station file exists, load only that file
    candidate = input_dir / f"{station_id}.csv"
    files_to_check = [candidate] if candidate.exists() else sorted(input_dir.glob("*.csv"))

    for f in files_to_check:
        sep = detect_delimiter(f)
        try:
            df = pl.read_csv(
                f, separator=sep,
                infer_schema_length=0,
                ignore_errors=True,
                truncate_ragged_lines=True,
            )
        except Exception:
            continue

        if df.is_empty():
            continue

        sid_col = find_sid_col(df.columns)
        if sid_col:
            df = df.filter(pl.col(sid_col).cast(pl.Utf8) == station_id)
        else:
            # No station-ID column — treat stem as station ID
            if f.stem != station_id:
                continue

        if df.is_empty():
            continue
        frames.append(df)

    if not frames:
        return None

    df = pl.concat(frames, how="diagonal_relaxed")

    ts_col = find_ts_col(df.columns)
    if ts_col is None:
        print(f"    WARNING: no timestamp column found — skipping station {station_id}")
        return None

    df = (
        df.with_columns(
            pl.col(ts_col)
            .str.strptime(pl.Datetime, "%Y-%m-%d %H:%M:%S", strict=False)
            .alias("_ts")
        )
        .filter(pl.col("_ts").is_not_null())
        .sort("_ts")
    )

    return df if not df.is_empty() else None


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────────────────────────────────────

def _violation_spans(ts, vals, lo, hi):
    """Return list of (start_ts, end_ts) for contiguous out-of-range runs."""
    spans = []
    in_viol = False
    vstart = None
    for t, v in zip(ts, vals):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            if in_viol:
                spans.append((vstart, t))
                in_viol = False
            continue
        bad = (lo is not None and v < lo) or (hi is not None and v > hi)
        if bad and not in_viol:
            vstart = t
            in_viol = True
        elif not bad and in_viol:
            spans.append((vstart, t))
            in_viol = False
    if in_viol and vstart is not None and ts:
        spans.append((vstart, ts[-1]))
    return spans


def plot_station(station_id: str, df: pl.DataFrame,
                 output_dir: Path, show: bool,
                 var_filter: str | None):
    """
    Build a subplot grid with one panel per monitored variable.
    Red shading = threshold violation.  Dashed lines = bounds.
    """
    # Determine which variables are present in this station's data
    target_vars = ([var_filter] if var_filter and var_filter in VARIABLES
                   else list(VARIABLES.keys()))
    found = []
    for norm_name in target_vars:
        col = find_var_col(df.columns, norm_name)
        if col is not None:
            found.append((norm_name, col))

    if not found:
        print(f"    No monitored variables found — skipping")
        return

    n      = len(found)
    ncols  = min(COLS, n)
    nrows  = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(7.5 * ncols, 3.5 * nrows),
        squeeze=False,
    )
    fig.suptitle(
        f"Station {station_id}  ·  Threshold Violations Only",
        fontsize=15, fontweight="bold", color=LINE_COLOR, y=1.01,
    )

    ts_list = df["_ts"].to_list()   # list[datetime]

    for idx, (norm_name, col) in enumerate(found):
        ax  = axes[idx // ncols][idx % ncols]
        lo, hi = THRESHOLDS.get(norm_name, (None, None))

        # ── values ────────────────────────────────────────────────────────
        try:
            raw = df[col].cast(pl.Float64, strict=False).to_list()
        except Exception:
            raw = [float("nan")] * len(df)
        vals = [v if v is not None else float("nan") for v in raw]

        # ── line ──────────────────────────────────────────────────────────
        ax.plot(ts_list, vals, color=LINE_COLOR, linewidth=0.6, alpha=0.85, zorder=2)

        # ── violation shading + threshold lines ───────────────────────────
        n_flags = 0
        if lo is not None or hi is not None:
            spans = _violation_spans(ts_list, vals, lo, hi)
            n_flags = sum(
                1 for v in vals
                if v is not None and not (isinstance(v, float) and np.isnan(v))
                and ((lo is not None and v < lo) or (hi is not None and v > hi))
            )
            for s, e in spans:
                ax.axvspan(s, e, color=FILL_COLOR, alpha=0.22, linewidth=0, zorder=1)

            if lo is not None:
                ax.axhline(lo, color=THRESH_COLOR, linewidth=0.9,
                           linestyle="--", alpha=0.85, label=f"Min: {lo}", zorder=3)
            if hi is not None:
                ax.axhline(hi, color=THRESH_COLOR, linewidth=0.9,
                           linestyle="--", alpha=0.85, label=f"Max: {hi}", zorder=3)

        # ── titles / labels ───────────────────────────────────────────────
        title_color = THRESH_COLOR if n_flags > 0 else "#2d3748"
        weight      = "bold"        if n_flags > 0 else "normal"
        flag_str    = f"  ⚠ {n_flags:,} flags" if n_flags > 0 else "  ✓ OK"
        ax.set_title(f"{norm_name}{flag_str}", fontsize=8.5,
                     color=title_color, fontweight=weight)

        ax.tick_params(axis="both", labelsize=7)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3, maxticks=7))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
        ax.grid(True, color=GRID_COLOR, linewidth=0.5, zorder=0)
        ax.set_facecolor(SUBPLOT_BG)
        ax.spines[["top", "right"]].set_visible(False)

        if (lo is not None or hi is not None) and n_flags > 0:
            ax.legend(fontsize=6.5, loc="upper right", framealpha=0.7)

    # ── hide unused subplot slots ─────────────────────────────────────────
    for idx in range(len(found), nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    plt.tight_layout()

    if show:
        plt.show()
    else:
        suffix = f"_{var_filter}" if var_filter else ""
        out_path = output_dir / f"threshold_graph_{station_id}{suffix}.png"
        fig.savefig(out_path, dpi=100, bbox_inches="tight")
        print(f"    → {out_path.name}")

    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Plot threshold violation time-series for Mesonet stations"
    )
    ap.add_argument("-i", "--input",    required=True,
                    help="Station data directory (same -i as Threshold_Errors.py)")
    ap.add_argument("-o", "--output",   default=None,
                    help="Output directory for PNGs  (default: <input>/../threshold_graphs)")
    ap.add_argument("--station",        default=None,
                    help="Graph a single station (e.g. 0115)")
    ap.add_argument("--variable",       default=None,
                    help="Graph only one variable (e.g. tsoil_1, sw_in)")
    ap.add_argument("--show",           action="store_true",
                    help="Open interactive matplotlib window instead of saving PNGs")
    args = ap.parse_args()

    input_dir  = Path(args.input)
    output_dir = (Path(args.output) if args.output
                  else input_dir.parent / "threshold_graphs")
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.variable and args.variable not in VARIABLES:
        print(f"ERROR: Unknown variable '{args.variable}'. "
              f"Valid options: {', '.join(VARIABLES)}")
        sys.exit(1)

    # ── discover stations ─────────────────────────────────────────────────
    all_csvs = sorted(input_dir.glob("*.csv"))
    if not all_csvs:
        print(f"No CSV files found in {input_dir}")
        sys.exit(1)

    if args.station:
        stations = [args.station]
    else:
        # Use file stems as station IDs (0115.csv → "0115")
        stations = [f.stem for f in all_csvs]

    print(f"{'='*60}")
    print(f"  Hawaii Mesonet  |  Threshold Graph  |  {len(stations)} station(s)")
    if args.variable:
        print(f"  Variable filter : {args.variable}")
    print(f"  Output dir      : {output_dir}")
    print(f"{'='*60}")

    for sid in stations:
        print(f"\nStation {sid}...")
        df = load_station(sid, input_dir)
        if df is None:
            print(f"  No data found — skipping")
            continue
        print(f"  {len(df):,} rows loaded")
        plot_station(sid, df, output_dir, args.show, args.variable)

    if not args.show:
        print(f"\nDone. PNGs saved to: {output_dir}\n")


if __name__ == "__main__":
    main()
