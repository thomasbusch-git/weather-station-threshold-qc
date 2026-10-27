#!/usr/bin/env python3
"""
================================================================================
 Script:      Threshold_Errors.py
 Version:     1.0.0
 Description: Simplified threshold-only error detection for Hawaii Mesonet
              stations. Flags values outside defined valid ranges (Critical).
              Consolidates consecutive out-of-range periods into single events
              using a 1 calendar day bridge tolerance.
              HTML output only, with sub-period dropdown per event.

 Author:      Thomas Busch
 Date:        2026-04-24

 Usage:
   python Threshold_Errors.py -i /path/to/csvs -o ./threshold_output
   python Threshold_Errors.py -i /path/to/csvs -o ./threshold_output --station 0235

 Requirements:
   pip install polars
================================================================================
"""

import polars as pl
import argparse
import sys
import json
import csv
import math
import ephem
from pathlib import Path
from datetime import datetime, timedelta, date as date_type, time as time_type
from dataclasses import dataclass, field
from typing import Optional, List
from functools import lru_cache

# ─────────────────────────────────────────────────────────────────────────────
# VARIABLES AND THRESHOLDS  (from project threshold table)
# ─────────────────────────────────────────────────────────────────────────────

# Maps internal name -> list of possible CSV column headers
VARIABLES = {
    "air_temp_1":     ["Air temperature, sensor 1 (°C)",  "Air temperature, sensor 1 (Â°C)",  "Tair_1_Avg"],
    "air_temp_2":     ["Air temperature, sensor 2 (°C)",  "Air temperature, sensor 2 (Â°C)",  "Tair_2_Avg"],
    "rh_1":           ["Relative humidity, sensor 1 (%)", "RH_1_Avg"],
    "rh_2":           ["Relative humidity, sensor 2 (%)", "RH_2_Avg"],
    "sw_in":          ["Incoming shortwave radiation (W/m²)", "Shortwave incoming radiation (W/m²)", "SWin_1_Avg"],
    "sw_out":         ["Outgoing shortwave radiation (W/m²)", "Shortwave outgoing radiation (W/m²)", "SWout_1_Avg"],
    "lw_in":          ["Incoming longwave radiation (W/m²)",  "Longwave incoming radiation (W/m²)",  "LWin_1_Avg"],
    "lw_out":         ["Outgoing longwave radiation (W/m²)",  "Longwave outgoing radiation (W/m²)",  "LWout_1_Avg"],
    "wind_speed":     ["Wind speed, scalar average (m/s)", "WS_1_Avg"],
    "wind_dir":       ["Wind direction, scalar average (°)", "Unit vector wind direction (°)", "WDuv_1_Avg"],
    "air_pressure":   ["Air pressure (kPa)", "P_1_Avg"],
    "surface_temp":   ["Surface temperature (°C)", "Tsrf_1_Avg"],
    "rainfall_1":     ["Rainfall, 5-minute total (mm)", "RF_1_Tot300s"],
    "rainfall_2":     ["Rainfall, 5-minute total, sensor 2 (mm)", "Rainfall total, sensor 2 (mm)", "RF_2_Tot300s"],
    "soil_vwc_1":     ["Soil volumetric water content (m³/m³)",           "Soil volumetric water content (mÂ³/mÂ³)",          "SM_1_Avg"],
    "soil_vwc_2":     ["Soil volumetric water content, sensor 2 (m³/m³)", "Soil volumetric water content, sensor 2 (mÂ³/mÂ³)", "SM_2_Avg"],
    "soil_vwc_3":     ["Soil volumetric water content, sensor 3 (m³/m³)", "Soil volumetric water content, sensor 3 (mÂ³/mÂ³)", "SM_3_Avg"],
    "soil_heat_flux": ["Soil heat flux at ground surface (W/m²)", "Soil heat flux at ground surface (W/mÂ²)", "SHFsrf_1_Avg"],
    "soil_hf_1":      ["Soil heat flux, sensor 1 (W/m²)", "Soil heat flux below ground, sensor 1 (W/m²)", "SHF_1_Avg"],
    "soil_hf_2":      ["Soil heat flux below ground, sensor 2 (W/m²)", "Soil heat flux, depth 2 (W/m²)", "SHF_2_Avg"],
    "tsoil_1":        ["Soil temperature (°C)",          "Soil temperature, depth 1 (°C)",  "Tsoil_1_Avg"],
    "tsoil_2":        ["Soil temperature, sensor 2 (°C)", "Soil temperature, depth 2 (°C)", "Tsoil_2_Avg"],
    "tsoil_3":        ["Soil temperature, sensor 3 (°C)", "Soil temperature, depth 3 (°C)", "Tsoil_3_Avg"],
    "tsoil_4":        ["Soil temperature, sensor 4 (°C)", "Soil temperature, depth 4 (°C)", "Tsoil_4_Avg"],
}

# Valid ranges — synced with MesonetAnomaly_v2.2 VALID_RANGES + TIERED_RANGES
THRESHOLDS = {
    "air_temp_1":     (-10,  40),
    "air_temp_2":     (-10,  40),
    "rh_1":           (0,    100),
    "rh_2":           (0,    100),
    "sw_in":          (0,    1400),
    "sw_out":         (0,    400),
    "lw_in":          (200,  500),
    "lw_out":         (200,  700),
    "wind_speed":     (0,    80),
    "wind_dir":       (0,    360),
    "air_pressure":   (60,   104),
    "surface_temp":   (0,    65),
    "rainfall_1":     (0,    150),
    "rainfall_2":     (0,    150),
    "soil_vwc_1":     (0,    0.9),
    "soil_vwc_2":     (0,    0.9),
    "soil_vwc_3":     (0,    0.9),
    "soil_heat_flux": (-400, 400),
    "soil_hf_1":      (-150, 150),
    "soil_hf_2":      (-150, 150),
    "tsoil_1":        (5,    45),
    "tsoil_2":        (5,    45),
    "tsoil_3":        (5,    45),
    "tsoil_4":        (5,    45),
}

# Standard display names — Campbell Scientific format (last alias in each VARIABLES entry)
DISPLAY_NAMES = {k: v[-1] for k, v in VARIABLES.items()}

# Derived variables produced by each measured variable (keyed by Campbell display name)
DERIVED_FROM = {
    "SWin_1_Avg":  ["SWnet_1_Avg", "Rin_1_Avg",  "Rnet_1_Avg", "Albedo_1_Avg"],
    "SWout_1_Avg": ["SWnet_1_Avg", "Rout_1_Avg", "Rnet_1_Avg", "Albedo_1_Avg"],
    "LWin_1_Avg":  ["LWnet_1_Avg", "Rin_1_Avg",  "Rnet_1_Avg", "Tsky_1_Avg"],
    "LWout_1_Avg": ["LWnet_1_Avg", "Rout_1_Avg", "Rnet_1_Avg", "Tsrf_1_Avg"],
    "Tair_1_Avg":  ["VPsat_1_Avg", "VP_1_Avg",   "VPD_1_Avg"],
    "RH_1_Avg":    ["VP_1_Avg",    "VPD_1_Avg"],
    "Tair_2_Avg":  ["VPsat_2_Avg", "VP_2_Avg",   "VPD_2_Avg"],
    "RH_2_Avg":    ["VP_2_Avg",    "VPD_2_Avg"],
    "WS_1_Avg":    ["WSrs_1_Avg",  "WDrs_1_Avg", "WD_1_StdCS"],
    "WDuv_1_Avg":  ["WSrs_1_Avg",  "WDrs_1_Avg", "WD_1_StdY",  "WD_1_StdCS"],
    "P_1_Avg":     ["Psl_1_Avg"],
    "SM_1_Avg":    ["SHFsrf_1_Avg"],
    "Tsoil_1_Avg": ["SHFsrf_1_Avg"],
    "SHF_1_Avg":   ["SHFsrf_1_Avg"],
    "SHF_2_Avg":   ["SHFsrf_1_Avg"],
}

TIMESTAMP_ALIASES  = ["Timestamp (Pacific/Honolulu)", "Timestamp"]
STATION_ID_ALIASES = ["Station ID"]

# Bridge tolerance: 1 calendar day
BRIDGE_TIMEDELTA = timedelta(days=1)

# sw_in / sw_out are only flagged during daylight hours.
# Sunrise/sunset computed per station per date using ephem + station_metadata.csv.
SW_DAYTIME_VARS = {"sw_in", "sw_out"}

def _load_station_coords():
    """Load lat/lon/elevation/timezone from station_metadata.csv."""
    coords = {}
    meta = Path(__file__).parent.parent / "data" / "station_metadata.csv"
    if not meta.exists():
        return coords
    tz_map = {"hawaii": timedelta(hours=-10), "american_samoa": timedelta(hours=-11)}
    with open(meta, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                lat  = float(row["lat"])
                lon  = float(row["lng"])
                elev = float(row["elevation"]) if row.get("elevation") else 0.0
            except (ValueError, TypeError):
                continue  # planned stations with missing coords
            coords[row["station_id"]] = {
                "lat": lat, "lon": lon, "elev": elev,
                "tz": tz_map.get(row.get("location", "hawaii"), timedelta(hours=-10)),
            }
    return coords

STATION_COORDS = _load_station_coords()

@lru_cache(maxsize=4096)
def _sw_window(station_id: str, date_ord: int):
    """Return (sunrise, sunset) as time objects in local time for the given station/date.

    Uses ephem with the station's lat, lon, and elevation so high-elevation
    stations (e.g. Haleakalā Summit at 2984 m) get an earlier sunrise and
    later sunset than sea-level stations.  Falls back to 07:30–17:30 if the
    station is not in the metadata or has missing coordinates.
    """
    c = STATION_COORDS.get(station_id)
    if c is None:
        return (time_type(7, 30), time_type(17, 30))

    d = date_type.fromordinal(date_ord)
    tz      = c["tz"]
    tz_h    = int(tz.total_seconds() / 3600)          # -10 or -11
    noon_utc_h = 12 - tz_h                            # local noon expressed in UTC hours

    _R   = 6371000.0
    elev = c["elev"]

    obs           = ephem.Observer()
    obs.lat       = str(c["lat"])
    obs.lon       = str(c["lon"])
    obs.elevation = elev
    # Depress the geometric horizon for elevated observers
    obs.horizon   = str(math.degrees(-math.acos(_R / (_R + elev)))) if elev > 0 else "0"
    obs.date      = f"{d.year}/{d.month}/{d.day} {noon_utc_h:02d}:00:00"

    try:
        sr_utc = obs.previous_rising(ephem.Sun()).datetime()
        ss_utc = obs.next_setting(ephem.Sun()).datetime()
    except (ephem.AlwaysUpError, ephem.NeverUpError):
        return (time_type(0, 0), time_type(23, 59))   # polar edge-case

    return ((sr_utc + tz).time(), (ss_utc + tz).time())


# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SubPeriod:
    start_time: str
    end_time: str
    duration_minutes: int
    bad_count: int


@dataclass
class ThresholdEvent:
    station_id: str
    variable: str
    vmin: float
    vmax: float
    start_time: Optional[str]
    end_time: Optional[str]
    duration_minutes: Optional[int]
    bad_reading_count: int
    sub_periods: List[SubPeriod] = field(default_factory=list)
    derived_affected: List[str] = field(default_factory=list)


@dataclass
class StationSummary:
    station_id: str
    rows: int = 0
    timestamp_min: str = ""
    timestamp_max: str = ""
    event_count: int = 0
    vars_found: List[str] = field(default_factory=list)


# ─────────────────────────────────────────────────────────────────────────────
# TIMESTAMP PARSING
# ─────────────────────────────────────────────────────────────────────────────

def parse_ts(ts_str):
    if ts_str is None:
        return None
    s = str(ts_str).strip()
    try:
        return datetime.strptime(s, "%Y-%m-%dT%H:%M:%S%z").replace(tzinfo=None)
    except Exception:
        pass
    if len(s) > 6 and s[-6] in ('+', '-'):
        try:
            return datetime.strptime(s[:-6], "%Y-%m-%dT%H:%M:%S")
        except Exception:
            pass
    for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except Exception:
            continue
    return None


# ─────────────────────────────────────────────────────────────────────────────
# CONSOLIDATION  (time-based bridge)
# ─────────────────────────────────────────────────────────────────────────────

def consolidate(flag_indices, timestamps_dt, bridge=BRIDGE_TIMEDELTA):
    """
    Merge flagged indices into events using a 1-day time bridge.
    Returns list of (start_idx, end_idx, sub_runs)
    where sub_runs are contiguous runs of bad rows within the event.
    """
    if not flag_indices:
        return []

    events   = []
    run_start = flag_indices[0]
    run_end   = flag_indices[0]
    # sub_runs tracks contiguous bad-row groups
    sub_runs  = [[flag_indices[0], flag_indices[0]]]

    for k in range(1, len(flag_indices)):
        idx      = flag_indices[k]
        prev_idx = flag_indices[k - 1]

        t_prev = timestamps_dt[prev_idx] if prev_idx < len(timestamps_dt) else None
        t_curr = timestamps_dt[idx]       if idx      < len(timestamps_dt) else None

        within_bridge = (t_prev and t_curr and (t_curr - t_prev) <= bridge)

        if within_bridge:
            run_end = idx
            # Extend last sub-run if contiguous, else start new one
            if idx == sub_runs[-1][1] + 1:
                sub_runs[-1][1] = idx
            else:
                sub_runs.append([idx, idx])
        else:
            events.append((run_start, run_end, [tuple(r) for r in sub_runs]))
            run_start = idx
            run_end   = idx
            sub_runs  = [[idx, idx]]

    events.append((run_start, run_end, [tuple(r) for r in sub_runs]))
    return events


def make_sub_periods(sub_runs, timestamps):
    result = []
    for (s, e) in sub_runs:
        ts_s = str(timestamps[s]) if s < len(timestamps) else ""
        ts_e = str(timestamps[e]) if e < len(timestamps) else ""
        dt_s, dt_e = parse_ts(ts_s), parse_ts(ts_e)
        if dt_s and dt_e:
            dur = int((dt_e - dt_s).total_seconds() / 60) + 5
        else:
            dur = (e - s + 1) * 5
        result.append(SubPeriod(
            start_time=ts_s, end_time=ts_e,
            duration_minutes=dur, bad_count=(e - s + 1),
        ))
    return result


# ─────────────────────────────────────────────────────────────────────────────
# LOADING
# ─────────────────────────────────────────────────────────────────────────────

def detect_delimiter(filepath):
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
    return "\t" if header.count("\t") > header.count(",") else ","


def resolve_col(df_columns, aliases):
    for a in aliases:
        if a in df_columns:
            return a
    return None


def load_csvs(input_path):
    path = Path(input_path)
    frames, n = [], 0
    files = [path] if path.is_file() else sorted(path.glob("*.csv"))
    if not files:
        print(f"[ERROR] No CSV files found."); sys.exit(1)

    for f in files:
        try:
            sep = detect_delimiter(f)
            df  = pl.read_csv(f, separator=sep, infer_schema_length=0,
                              try_parse_dates=False, truncate_ragged_lines=True,
                              ignore_errors=True)
            casts = [pl.col(c).cast(pl.Float64, strict=False).alias(c)
                     for c in df.columns
                     if "timestamp" not in c.lower() and "station" not in c.lower()]
            if casts:
                df = df.with_columns(casts)
            frames.append(df); n += 1
            print(f"  Loaded: {f.name} ({len(df):,} rows)")
        except Exception as e:
            print(f"  [WARN] {f.name}: {e}")

    combined = pl.concat(frames, how="diagonal_relaxed")
    print(f"  Combined: {len(combined):,} rows")
    return combined, n


# ─────────────────────────────────────────────────────────────────────────────
# PER-STATION ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_station(station_id, sdf, col_map, ts_col):
    if ts_col:
        sdf = sdf.sort(ts_col)

    timestamps    = sdf[ts_col].to_list() if ts_col else []
    timestamps_dt = [parse_ts(t) for t in timestamps]

    summary = StationSummary(station_id=station_id, rows=len(sdf))
    if timestamps:
        summary.timestamp_min = str(timestamps[0])
        summary.timestamp_max = str(timestamps[-1])

    events = []

    for norm_name, (vmin, vmax) in THRESHOLDS.items():
        actual_col = col_map.get(norm_name)
        if actual_col is None:
            continue
        summary.vars_found.append(DISPLAY_NAMES.get(norm_name, norm_name))
        try:
            series = sdf[actual_col].cast(pl.Float64, strict=False).to_list()
        except Exception:
            continue

        # Find flagged indices (sw_in/sw_out restricted to daytime only)
        is_sw = norm_name in SW_DAYTIME_VARS
        flag_indices = []
        for i, v in enumerate(series):
            if v is None:
                continue
            if is_sw and timestamps_dt[i] is not None:
                ts = timestamps_dt[i]
                sr, ss = _sw_window(station_id, ts.date().toordinal())
                if ts.time() < sr or ts.time() > ss:
                    continue
            if v < vmin or v > vmax:
                flag_indices.append(i)
        if not flag_indices:
            continue

        raw_events = consolidate(flag_indices, timestamps_dt)

        for (start_idx, end_idx, sub_runs) in raw_events:
            ts_s = str(timestamps[start_idx]) if start_idx < len(timestamps) else None
            ts_e = str(timestamps[end_idx])   if end_idx   < len(timestamps) else None
            dt_s, dt_e = parse_ts(ts_s), parse_ts(ts_e)
            duration = int((dt_e - dt_s).total_seconds() / 60) + 5 if dt_s and dt_e else None
            bad_count = sum(1 for i in range(start_idx, end_idx + 1)
                            if i < len(series) and series[i] is not None
                            and (series[i] < vmin or series[i] > vmax))

            sub_periods = make_sub_periods(sub_runs, timestamps)

            disp = DISPLAY_NAMES.get(norm_name, norm_name)
            events.append(ThresholdEvent(
                station_id=station_id, variable=disp,
                vmin=vmin, vmax=vmax,
                start_time=ts_s, end_time=ts_e,
                duration_minutes=duration, bad_reading_count=bad_count,
                sub_periods=sub_periods,
                derived_affected=DERIVED_FROM.get(disp, []),
            ))
            summary.event_count += 1

    return summary, events


# ─────────────────────────────────────────────────────────────────────────────
# HTML REPORT
# ─────────────────────────────────────────────────────────────────────────────

def fmt_dur(mins):
    if mins is None: return "—"
    if mins < 60:    return f"{mins} min"
    if mins < 1440:  return f"{mins//60}h {mins%60}m"
    return f"{mins//1440}d {(mins%1440)//60}h"


def generate_html(all_summaries, all_events, input_path, total_rows, total_files, output_dir):
    out = output_dir / "threshold_output.html"
    summaries_sorted = sorted(all_summaries, key=lambda s: -s.event_count)
    total_events = len(all_events)

    html = [f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Threshold Error Report</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'Segoe UI',Tahoma,sans-serif;max-width:1500px;margin:0 auto;
      padding:20px;background:#f0f4f8;color:#2d3748}}
h1{{background:linear-gradient(135deg,#3b2a8a,#6c4fd6);color:#fff;
    padding:18px 24px;border-radius:8px;margin-bottom:6px;font-size:1.5em}}
h2{{color:#3b2a8a;border-bottom:2px solid #3b2a8a;padding-bottom:6px;margin:28px 0 12px}}
h3{{margin:0;font-size:1em}}
.ts{{color:#718096;margin-bottom:18px;font-size:0.88em}}
.meta{{background:#fff;padding:16px;border-radius:8px;
       border:1px solid #cbd5e0;margin-bottom:18px}}
.meta table{{border-collapse:collapse;width:100%}}
.meta td{{padding:5px 12px;border-bottom:1px solid #eee;font-size:0.9em}}
.meta td:first-child{{font-weight:600;width:180px;color:#4a5568}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
       gap:12px;margin-bottom:20px}}
.card{{background:#fff;border-radius:8px;padding:14px;
       border:1px solid #cbd5e0;text-align:center;
       box-shadow:0 1px 3px rgba(0,0,0,.05)}}
.card .n{{font-size:2em;font-weight:700;line-height:1.1;color:#c53030}}
.card .l{{color:#718096;font-size:0.82em;margin-top:4px}}
.badge{{display:inline-block;padding:2px 10px;border-radius:12px;
        font-size:0.75em;font-weight:700;background:#c53030;color:#fff}}
table.dt{{border-collapse:collapse;width:100%;background:#fff;
          border-radius:8px;overflow:hidden;margin-bottom:14px;
          box-shadow:0 1px 3px rgba(0,0,0,.06)}}
table.dt th{{background:#3b2a8a;color:#fff;padding:9px 10px;
             text-align:left;font-size:0.78em;white-space:nowrap;
             cursor:pointer;user-select:none}}
table.dt th:after{{content:' ⇅';opacity:.5;font-size:.8em}}
table.dt th.asc:after{{content:' ▲';opacity:1}}
table.dt th.desc:after{{content:' ▼';opacity:1}}
table.dt td{{padding:7px 10px;border-bottom:1px solid #eee;
             font-size:0.8em;vertical-align:top}}
table.dt tr:hover{{background:#f7fafc}}
.ss-card{{background:#fff;border:1px solid #cbd5e0;border-radius:8px;
          padding:18px 20px;margin-bottom:16px;
          box-shadow:0 1px 3px rgba(0,0,0,.05)}}
.sh{{display:flex;justify-content:space-between;align-items:center;margin-bottom:12px}}
.nav{{background:#fff;border:1px solid #cbd5e0;border-radius:8px;
      padding:10px 14px;margin-bottom:20px;line-height:2}}
.nav a{{margin-right:6px;font-size:0.82em;color:#3b2a8a;text-decoration:none;
        padding:2px 7px;border:1px solid #e5deff;border-radius:4px}}
.nav a:hover{{background:#e5deff}}
/* Sub-period dropdown */
.sp-toggle{{cursor:pointer;color:#3b2a8a;font-size:0.78em;font-weight:600;
            border:1px solid #e5deff;background:none;padding:2px 8px;
            border-radius:4px;white-space:nowrap}}
.sp-toggle:hover{{background:#e5deff}}
.sp-panel{{display:none;margin-top:6px;padding:8px;
           background:#f7fafc;border-radius:6px;border:1px solid #e2e8f0}}
.sp-panel.open{{display:block}}
.sp-hdr{{font-weight:700;color:#4a5568;font-size:0.74em;
         display:grid;grid-template-columns:1fr 1fr auto auto;gap:8px;
         padding:2px 0 4px;border-bottom:2px solid #cbd5e0;margin-bottom:4px}}
.sp-row{{font-size:0.76em;padding:3px 0;border-bottom:1px solid #eee;
         display:grid;grid-template-columns:1fr 1fr auto auto;gap:8px}}
.sp-row:last-child{{border-bottom:none}}
.no-events{{color:#718096;font-style:italic;padding:8px 0}}
.filter-bar{{background:#f7fafc;border:1px solid #cbd5e0;border-radius:6px;
             padding:10px 14px;margin-bottom:10px;display:flex;flex-wrap:wrap;gap:8px;align-items:center}}
.filter-bar label{{font-size:0.8em;font-weight:600;color:#4a5568;margin-right:4px}}
.filter-btn{{cursor:pointer;padding:3px 10px;border-radius:12px;font-size:0.75em;
             font-weight:600;border:1px solid #cbd5e0;background:#fff;color:#4a5568}}
.filter-btn.derived{{border-color:#b7791f;color:#744210;background:#fffbeb}}
.filter-btn.derived.active{{background:#b7791f;color:#fff;border-color:#b7791f}}
.derived-tag{{font-size:0.72em;color:#b7791f;margin-top:3px;line-height:1.4}}
.filter-btn.active{{background:#3b2a8a;color:#fff;border-color:#3b2a8a}}
.filter-sep{{width:1px;background:#cbd5e0;height:20px;margin:0 4px}}
</style>
<script>
function toggleSP(btn) {{
  var panel = btn.nextElementSibling;
  panel.classList.toggle('open');
  btn.textContent = panel.classList.contains('open') ? '▲ Hide sub-periods' : '▼ Show sub-periods';
}}
function makeSortable(table) {{
  table.querySelectorAll('th').forEach((th, col) => {{
    th.addEventListener('click', () => {{
      var tbody = table.querySelector('tbody');
      var rows  = Array.from(tbody.querySelectorAll('tr'));
      var asc   = !th.classList.contains('asc');
      table.querySelectorAll('th').forEach(h => h.classList.remove('asc','desc'));
      th.classList.add(asc ? 'asc' : 'desc');
      rows.sort((a, b) => {{
        var av = a.cells[col]?.innerText.trim()||'';
        var bv = b.cells[col]?.innerText.trim()||'';
        var an = parseFloat(av), bn = parseFloat(bv);
        return (!isNaN(an)&&!isNaN(bn))?(asc?an-bn:bn-an):(asc?av.localeCompare(bv):bv.localeCompare(av));
      }});
      rows.forEach(r => tbody.appendChild(r));
    }});
  }});
}}
function initFilters(sid) {{
  var fb = document.getElementById('fb-' + sid);
  if (!fb) return;
  var card = document.getElementById('s-' + sid);
  var tbl  = card ? card.querySelector('table.dt') : null;
  if (!tbl) return;
  var activeVars    = new Set();
  var activeDerived = new Set();
  function applyFilters() {{
    tbl.querySelectorAll('tbody tr').forEach(function(row) {{
      var v   = row.getAttribute('data-var')     || '';
      var dv  = row.getAttribute('data-derived') || '';
      var dvSet = new Set(dv.split(' ').filter(Boolean));
      var showMeas = activeVars.size    === 0 || activeVars.has(v);
      var showDerv = activeDerived.size === 0 ||
                     [...activeDerived].some(function(d) {{ return dvSet.has(d); }});
      row.style.display = (showMeas && showDerv) ? '' : 'none';
    }});
  }}
  fb.querySelectorAll('.filter-btn[data-var]:not(.derived)').forEach(function(btn) {{
    btn.addEventListener('click', function() {{
      var v = btn.getAttribute('data-var');
      if (activeVars.has(v)) {{ activeVars.delete(v); btn.classList.remove('active'); }}
      else                   {{ activeVars.add(v);    btn.classList.add('active');    }}
      applyFilters();
    }});
  }});
  fb.querySelectorAll('.filter-btn.derived[data-derived]').forEach(function(btn) {{
    btn.addEventListener('click', function() {{
      var d = btn.getAttribute('data-derived');
      if (activeDerived.has(d)) {{ activeDerived.delete(d); btn.classList.remove('active'); }}
      else                      {{ activeDerived.add(d);    btn.classList.add('active');    }}
      applyFilters();
    }});
  }});
  var resetBtn = fb.querySelector('.filter-reset');
  if (resetBtn) {{
    resetBtn.addEventListener('click', function() {{
      activeVars.clear(); activeDerived.clear();
      fb.querySelectorAll('.filter-btn').forEach(function(b) {{ b.classList.remove('active'); }});
      applyFilters();
    }});
  }}
}}
document.addEventListener('DOMContentLoaded', () => {{
  document.querySelectorAll('table.dt').forEach(makeSortable);
  document.querySelectorAll('[id^="fb-"]').forEach(function(fb) {{
    var sid = fb.id.slice(3);
    initFilters(sid);
  }});
}});
</script>
</head>
<body>
<h1>🌺 Hawaii Mesonet — Threshold Error Report</h1>
<p class="ts">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} &nbsp;|&nbsp; Threshold_Errors v1.0.0</p>
<div class="meta"><table>
  <tr><td>Input Path</td><td>{input_path}</td></tr>
  <tr><td>Files Processed</td><td>{total_files:,}</td></tr>
  <tr><td>Total Rows</td><td>{total_rows:,}</td></tr>
  <tr><td>Stations</td><td>{len(all_summaries)}</td></tr>
  <tr><td>Total Threshold Events</td><td>{total_events:,}</td></tr>
</table></div>
<div class="grid">
  <div class="card"><div class="n" style="color:#2d3748">{len(all_summaries)}</div><div class="l">Stations</div></div>
  <div class="card"><div class="n" style="color:#2d3748">{total_rows:,}</div><div class="l">Total Rows</div></div>
  <div class="card"><div class="n">{total_events:,}</div><div class="l">Out of Range Events</div></div>
</div>
"""]

    # Overview table
    html.append("<h2>Station Overview</h2>")
    html.append("""<table class="dt"><thead><tr>
<th>Station</th><th>Rows</th><th>Date Range</th><th>Out of Range Events</th>
</tr></thead><tbody>""")
    for s in summaries_sorted:
        html.append(f"""<tr>
  <td><a href="#s-{s.station_id}" style="color:#3b2a8a;font-weight:600">{s.station_id}</a></td>
  <td>{s.rows:,}</td>
  <td style="font-size:0.78em;white-space:nowrap">
    {s.timestamp_min[:10] if s.timestamp_min else '?'} → 
    {s.timestamp_max[:10] if s.timestamp_max else '?'}
  </td>
  <td><span class="badge">{s.event_count}</span></td>
</tr>""")
    html.append("</tbody></table>")

    # Nav
    html.append('<div class="nav"><strong>Jump to station:</strong> ')
    for s in sorted(all_summaries, key=lambda x: x.station_id):
        html.append(f'<a href="#s-{s.station_id}">{s.station_id}</a>')
    html.append('</div>')

    # Per-station event tables
    html.append("<h2>Station Threshold Events</h2>")
    for s in summaries_sorted:
        sid = s.station_id
        s_events = sorted([e for e in all_events if e.station_id == sid],
                          key=lambda e: e.start_time or "")

        html.append(f"""<div class="ss-card" id="s-{sid}">
<div class="sh">
  <h3>Station {sid}
  &nbsp; <span style="font-weight:400;color:#718096;font-size:0.9em">
  {s.rows:,} rows &nbsp;|&nbsp; {s.event_count} out-of-range events
  </span></h3>
  <span style="color:#718096;font-size:0.8em">
    {s.timestamp_min[:10] if s.timestamp_min else ''} →
    {s.timestamp_max[:10] if s.timestamp_max else ''}
  </span>
</div>""")

        if s_events:
            _vars_with_events = set(e.variable for e in s_events)
            _all_vars = sorted(s.vars_found)
            # collect all derived vars affected at this station
            _all_derived: list = []
            _seen_derived: set = set()
            for _e in s_events:
                for _d in _e.derived_affected:
                    if _d not in _seen_derived:
                        _seen_derived.add(_d)
                        _all_derived.append(_d)
            _all_derived.sort()

            _fb = f'<div class="filter-bar" id="fb-{sid}"><label>Measured:</label>'
            for _v in _all_vars:
                if _v in _vars_with_events:
                    _fb += f'<button class="filter-btn" data-var="{_v}">{_v}</button>'
                else:
                    _fb += (f'<button class="filter-btn" data-var="{_v}" '
                            f'style="opacity:0.35;cursor:default" disabled>{_v}</button>')
            if _all_derived:
                _fb += '<div class="filter-sep"></div><label>Affected Derived:</label>'
                for _d in _all_derived:
                    _fb += (f'<button class="filter-btn derived" '
                            f'data-derived="{_d}">{_d}</button>')
            _fb += ('<div class="filter-sep"></div>'
                    '<button class="filter-btn filter-reset" '
                    'style="background:#e53e3e;color:#fff;border-color:#e53e3e">Reset</button></div>')
            html.append(_fb)
            html.append("""<table class="dt"><thead><tr>
<th>Variable</th><th>Valid Range</th><th>Start Time</th><th>End Time</th>
<th>Duration</th><th>Bad Readings</th><th>Sub-periods</th>
</tr></thead><tbody>""")
            for e in s_events:
                if len(e.sub_periods) > 1:
                    sp_html = (
                        '<button class="sp-toggle" onclick="toggleSP(this)">'
                        f'▼ Show sub-periods ({len(e.sub_periods)})</button>'
                        '<div class="sp-panel">'
                        '<div class="sp-hdr">'
                        '<span>Start</span><span>End</span>'
                        '<span>Duration</span><span>Bad Readings</span>'
                        '</div>'
                    )
                    for sp in e.sub_periods:
                        sp_html += (
                            f'<div class="sp-row">'
                            f'<span>{sp.start_time}</span>'
                            f'<span>{sp.end_time}</span>'
                            f'<span>{fmt_dur(sp.duration_minutes)}</span>'
                            f'<span>{sp.bad_count}</span>'
                            f'</div>'
                        )
                    sp_html += '</div>'
                else:
                    sp_html = '<span style="color:#a0aec0;font-size:0.78em">single period</span>'

                _derived_str = " ".join(e.derived_affected)
                _derived_td  = ""
                if e.derived_affected:
                    _derived_td = (
                        '<div class="derived-tag">⟶ Affects: '
                        + ", ".join(e.derived_affected)
                        + "</div>"
                    )
                html.append(f"""<tr data-var="{e.variable}" data-derived="{_derived_str}">
  <td style="font-weight:600">{e.variable}{_derived_td}</td>
  <td style="white-space:nowrap">[{e.vmin}, {e.vmax}]</td>
  <td style="font-size:0.78em;white-space:nowrap">{e.start_time or '—'}</td>
  <td style="font-size:0.78em;white-space:nowrap">{e.end_time   or '—'}</td>
  <td style="white-space:nowrap">{fmt_dur(e.duration_minutes)}</td>
  <td>{e.bad_reading_count:,}</td>
  <td style="min-width:180px">{sp_html}</td>
</tr>""")
            html.append("</tbody></table>")
        else:
            html.append('<p class="no-events">No threshold violations detected.</p>')

        html.append("</div>")

    html.append("""<hr style="margin-top:40px;border-color:#e2e8f0">
<p style="color:#a0aec0;font-size:0.78em;margin-top:10px">
  Threshold_Errors v1.0.0 — Thomas Busch</p>
</body></html>""")

    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(html))
    print(f"  HTML: {out}")


def generate_json(summaries, events, input_path, total_rows, total_files, output_dir):
    """Write threshold_report.json.

    Each measured-variable event is emitted with check_types='threshold' and a
    derived_affected list.  For every derived variable that depends on a flagged
    measured variable, an additional entry is emitted with check_types='derived'.
    """
    all_entries = []

    for e in events:
        # measured variable event
        all_entries.append({
            "station_id":        e.station_id,
            "variable":          e.variable,
            "check_types":       "threshold",
            "severity":          "Critical",
            "start_time":        e.start_time,
            "end_time":          e.end_time,
            "duration_minutes":  e.duration_minutes,
            "bad_reading_count": e.bad_reading_count,
            "detail":            f"Outside valid range [{e.vmin}, {e.vmax}]",
            "derived_affected":  e.derived_affected,
            "sub_periods": [
                {
                    "start_time":        sp.start_time,
                    "end_time":          sp.end_time,
                    "duration_minutes":  sp.duration_minutes,
                    "bad_reading_count": sp.bad_count,
                }
                for sp in e.sub_periods
            ],
        })
        # one entry per affected derived variable
        for derived_var in e.derived_affected:
            all_entries.append({
                "station_id":        e.station_id,
                "variable":          derived_var,
                "check_types":       "derived",
                "severity":          "Warning",
                "start_time":        e.start_time,
                "end_time":          e.end_time,
                "duration_minutes":  e.duration_minutes,
                "bad_reading_count": e.bad_reading_count,
                "detail":            f"Source variable {e.variable} had an error",
                "derived_affected":  [],
                "sub_periods":       [],
            })

    out = {
        "generated":     datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version":       "1.0",
        "input_path":    str(input_path),
        "total_rows":    total_rows,
        "total_files":   total_files,
        "station_count": len(summaries),
        "events":        all_entries,
    }
    json_path = Path(output_dir) / "threshold_report.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"  JSON : {json_path}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Threshold-only error detection for Hawaii Mesonet stations."
    )
    parser.add_argument("--input",   "-i", required=True)
    parser.add_argument("--output",  "-o", default="./threshold_output")
    parser.add_argument("--station", type=str, default=None)
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_path = Path(args.input)
    csv_files  = [input_path] if input_path.is_file() else sorted(input_path.glob("*.csv"))
    if not csv_files:
        print(f"[ERROR] No CSV files found."); sys.exit(1)

    total_rows    = 0
    total_files   = len(csv_files)
    all_summaries = []
    all_events    = []
    processed     = set()

    print(f"\n[PHASE 1-3] Loading and checking thresholds for {total_files} file(s)...")
    for fi, f in enumerate(csv_files, 1):
        try:
            sep = detect_delimiter(f)
            df  = pl.read_csv(f, separator=sep, infer_schema_length=0,
                              try_parse_dates=False, truncate_ragged_lines=True,
                              ignore_errors=True)
            casts = [pl.col(c).cast(pl.Float64, strict=False).alias(c)
                     for c in df.columns
                     if "timestamp" not in c.lower() and "station" not in c.lower()]
            if casts:
                df = df.with_columns(casts)
        except Exception as e:
            print(f"  [WARN] {f.name}: {e}"); continue

        total_rows += len(df)
        ts_col  = resolve_col(df.columns, TIMESTAMP_ALIASES)
        sid_col = resolve_col(df.columns, STATION_ID_ALIASES)
        col_map = {norm: resolve_col(df.columns, aliases)
                   for norm, aliases in VARIABLES.items()}

        if fi == 1:
            found   = [n for n,c in col_map.items() if c]
            missing = [n for n,c in col_map.items() if not c]
            print(f"  Variables found: {len(found)}  |  Not found: {len(missing)}")
            if missing:
                print(f"  Missing: {', '.join(missing)}")

        stations = sorted([str(s) for s in df[sid_col].drop_nulls().unique().to_list()]) \
                   if sid_col else ["unknown"]

        if args.station:
            if args.station not in stations:
                continue
            stations = [args.station]

        for sid in stations:
            if sid in processed:
                continue
            sdf = df.filter(pl.col(sid_col).cast(pl.Utf8) == sid) if sid_col else df
            print(f"  [{fi:2d}/{total_files}] {f.name} \u2192 Station {sid} ({len(sdf):,} rows)...",
                  end="", flush=True)
            summary, events = analyze_station(sid, sdf, col_map, ts_col)
            all_summaries.append(summary)
            all_events.extend(events)
            processed.add(sid)
            print(f" {summary.event_count} events")

    total_events = sum(s.event_count for s in all_summaries)
    print(f"\n[PHASE 4] Generating reports...")
    print(f"  Total threshold events: {total_events:,}")

    generate_html(
        all_summaries=all_summaries,
        all_events=all_events,
        input_path=str(input_path),
        total_rows=total_rows,
        total_files=total_files,
        output_dir=output_dir,
    )
    generate_json(
        summaries=all_summaries,
        events=all_events,
        input_path=input_path,
        total_rows=total_rows,
        total_files=total_files,
        output_dir=output_dir,
    )

    print(f"\n  Output dir : {output_dir}")
    print(f"  Done.\n")


if __name__ == "__main__":
    main()
