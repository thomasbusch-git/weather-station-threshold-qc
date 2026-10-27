#!/usr/bin/env python3
"""
generate_sample_data.py

Creates synthetic (fake) weather station data in the same CSV format the
threshold scripts expect, with known errors planted in it. This lets anyone
run the scripts and see a real error report without needing actual station data.

Writes:
    data/station_metadata.csv
    data/stationdata/9001.csv, 9002.csv, 9003.csv

Uses only the Python standard library. The random seed is fixed, so the
output is identical every time you run it.

Usage:
    python generate_sample_data.py
"""

import csv
import math
import random
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data"
STATION_DIR = DATA_DIR / "stationdata"

START = datetime(2026, 3, 1, 0, 0)
DAYS = 14
STEP = timedelta(minutes=5)
SEED = 42

# Fictional stations (IDs 9xxx so they can't be confused with real ones)
STATIONS = [
    # id,    name,             lat,     lng,       elev_m, base_temp_C
    ("9001", "Sample Upland",  20.80, -156.30,  650, 19.0),
    ("9002", "Sample Coastal", 21.30, -157.85,   10, 24.0),
    ("9003", "Sample Ridge",   19.80, -155.50, 1200, 15.0),
]

COLUMNS = [
    "Timestamp (Pacific/Honolulu)", "Station ID",
    "Incoming longwave radiation (W/m²)", "Outgoing longwave radiation (W/m²)",
    "Air pressure (kPa)",
    "Rainfall, 5-minute total (mm)", "Rainfall, 5-minute total, sensor 2 (mm)",
    "Relative humidity, sensor 1 (%)", "Relative humidity, sensor 2 (%)",
    "Net radiation (W/m²)",
    "Soil heat flux below ground, sensor 2 (W/m²)", "Soil heat flux at ground surface (W/m²)",
    "Soil volumetric water content (m³/m³)",
    "Soil volumetric water content, sensor 2 (m³/m³)",
    "Soil volumetric water content, sensor 3 (m³/m³)",
    "Incoming shortwave radiation (W/m²)", "Outgoing shortwave radiation (W/m²)",
    "Air temperature, sensor 1 (°C)", "Air temperature, sensor 2 (°C)",
    "Soil temperature (°C)", "Soil temperature, sensor 2 (°C)",
    "Soil temperature, sensor 3 (°C)", "Soil temperature, sensor 4 (°C)",
    "Surface temperature (°C)",
    "Vapor pressure, sensor 1 (kPa)", "Vapor pressure, sensor 2 (kPa)",
    "Vapor pressure deficit, sensor 1 (kPa)", "Vapor pressure deficit, sensor 2 (kPa)",
    "Wind direction, scalar average (°)", "Wind speed, scalar average (m/s)",
]

# Short names for the columns errors get planted in
COL = {
    "air_temp_2":   "Air temperature, sensor 2 (°C)",
    "rh_1":         "Relative humidity, sensor 1 (%)",
    "sw_in":        "Incoming shortwave radiation (W/m²)",
    "lw_in":        "Incoming longwave radiation (W/m²)",
    "air_pressure": "Air pressure (kPa)",
    "tsoil_3":      "Soil temperature, sensor 3 (°C)",
    "wind_speed":   "Wind speed, scalar average (m/s)",
    "wind_dir":     "Wind direction, scalar average (°)",
    "soil_vwc_2":   "Soil volumetric water content, sensor 2 (m³/m³)",
    "rainfall_1":   "Rainfall, 5-minute total (mm)",
}

# ── Planted errors ────────────────────────────────────────────────────────────
# (station, variable, start, end, bad value, what it simulates)
# Every timestamp from start to end (inclusive) gets the bad value.
def t(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M")

PLANTED_ERRORS = [
    ("9001", "air_temp_2", t("2026-03-10 14:35"), t("2026-03-10 14:35"), 55.3,
     "single-reading spike"),
    ("9001", "sw_in",      t("2026-03-04 10:00"), t("2026-03-04 12:00"), -25.0,
     "negative radiation during daylight (sensor fault)"),
    ("9001", "rh_1",       t("2026-03-07 02:00"), t("2026-03-07 03:30"), 103.0,
     "humidity over 100% in fog (burst 1)"),
    ("9001", "rh_1",       t("2026-03-07 05:00"), t("2026-03-07 05:45"), 101.5,
     "humidity over 100% in fog (burst 2, merges with burst 1)"),
    ("9001", "rh_1",       t("2026-03-12 01:00"), t("2026-03-12 01:30"), 102.0,
     "humidity over 100%, days later (separate event)"),
    ("9002", "tsoil_3",    t("2026-03-08 00:00"), t("2026-03-08 23:55"), -6999.0,
     "logger missing-value code for a full day"),
    ("9002", "air_pressure", t("2026-03-05 13:00"), t("2026-03-05 13:25"), 0.0,
     "pressure sensor dropout"),
    ("9002", "wind_speed", t("2026-03-03 15:00"), t("2026-03-03 15:00"), 99.0,
     "wind spike (part 1)"),
    ("9002", "wind_speed", t("2026-03-03 18:20"), t("2026-03-03 18:20"), 99.0,
     "wind spike (part 2, merges with part 1)"),
    ("9002", "wind_speed", t("2026-03-10 09:00"), t("2026-03-10 09:00"), 99.0,
     "wind spike, days later (separate event)"),
    ("9003", "soil_vwc_2", t("2026-03-06 06:00"), t("2026-03-06 12:00"), 1.2,
     "soil moisture above physical limit"),
    ("9003", "lw_in",      t("2026-03-09 20:00"), t("2026-03-09 23:00"), 150.0,
     "longwave radiation too low"),
    ("9003", "rainfall_1", t("2026-03-11 16:40"), t("2026-03-11 16:40"), 999.0,
     "impossible 5-minute rainfall"),
    ("9003", "wind_dir",   t("2026-03-02 08:15"), t("2026-03-02 08:30"), 400.0,
     "wind direction over 360°"),
]

# Fraction of cells left blank to mimic real gaps in logger data
MISSING_RATE = 0.002


# ── Helpers ───────────────────────────────────────────────────────────────────
def sat_vp(temp_c):
    """Saturation vapor pressure in kPa (Tetens equation)."""
    return 0.6108 * math.exp(17.27 * temp_c / (temp_c + 237.3))


def daylight(h):
    """0 at night, rising to 1 at solar noon (roughly 06:30–18:30 in March)."""
    if 6.5 < h < 18.5:
        return math.sin(math.pi * (h - 6.5) / 12.0)
    return 0.0


def generate_station(rng, sid, elev, base_temp):
    rows = []
    n = DAYS * 24 * 12
    pressure_base = 101.3 * math.exp(-elev / 8400.0)
    vwc = [0.35, 0.40, 0.38]
    wind_dir = 70.0
    rain_left = 0

    for i in range(n):
        ts = START + i * STEP
        h = ts.hour + ts.minute / 60.0
        day = daylight(h)
        diurnal = math.cos(2 * math.pi * (h - 14) / 24.0)   # peaks at 14:00
        clouds = 0.75 + 0.25 * math.sin(i / 37.0) + rng.uniform(-0.1, 0.1)

        air1 = base_temp + 4.0 * diurnal + rng.gauss(0, 0.15)
        air2 = air1 + rng.gauss(0, 0.05)
        rh1 = min(99.0, max(30.0, 78 - 14 * diurnal + rng.gauss(0, 1.0)))
        rh2 = min(99.0, max(30.0, rh1 + rng.gauss(0, 0.3)))

        if day > 0:
            sw_in = 1000 * day * clouds + rng.gauss(0, 5)
            sw_in = max(0.0, sw_in)
        elif h < 6 or h > 19:
            sw_in = rng.uniform(-2.0, 0.0)      # normal small negative at night
        else:
            sw_in = rng.uniform(0.0, 2.0)
        sw_out = max(0.0, 0.18 * sw_in + rng.gauss(0, 0.5))

        surf = air1 + 8 * day + rng.gauss(0, 0.3)
        lw_in = 345 + 10 * diurnal + rng.gauss(0, 4)
        lw_out = 400 + 4.5 * (surf - 18) + rng.gauss(0, 2)
        net = sw_in - sw_out + lw_in - lw_out

        # Rain: occasional showers lasting 15–60 minutes
        if rain_left == 0 and rng.random() < 0.004:
            rain_left = rng.randint(3, 12)
        if rain_left > 0:
            rain1 = round(rng.uniform(0.2, 3.0), 1)
            rain_left -= 1
            vwc = [min(0.50, v + r) for v, r in zip(vwc, (0.0020, 0.0014, 0.0008))]
        else:
            rain1 = 0.0
            vwc = [max(0.20, v - d) for v, d in zip(vwc, (0.00006, 0.00004, 0.00002))]
        rain2 = round(max(0.0, rain1 + rng.gauss(0, 0.1)), 1) if rain1 else 0.0

        shf_surf = 60 * diurnal * (0.5 + 0.5 * day) + rng.gauss(0, 3)
        shf_2 = 20 * math.cos(2 * math.pi * (h - 17) / 24.0) + rng.gauss(0, 1)

        ts1 = base_temp + 2 + 3.0 * math.cos(2 * math.pi * (h - 15) / 24.0) + rng.gauss(0, 0.05)
        ts2 = base_temp + 2 + 1.5 * math.cos(2 * math.pi * (h - 17) / 24.0)
        ts3 = base_temp + 2 + 0.6 * math.cos(2 * math.pi * (h - 20) / 24.0)
        ts4 = base_temp + 2 + 0.2 * math.cos(2 * math.pi * (h - 23) / 24.0)

        es1, es2 = sat_vp(air1), sat_vp(air2)
        vp1, vp2 = es1 * rh1 / 100, es2 * rh2 / 100

        wind_dir = (wind_dir + rng.gauss(0, 4)) % 360
        if abs(((wind_dir - 70 + 180) % 360) - 180) > 50:   # stay near trade winds
            wind_dir = 70 + rng.gauss(0, 10)
        wind_dir %= 360
        wind_speed = max(0.0, 3 + 1.5 * day + rng.gauss(0, 0.8))

        rows.append({
            COLUMNS[0]:  ts.strftime("%Y-%m-%d %H:%M:%S"),
            COLUMNS[1]:  sid,
            COLUMNS[2]:  round(lw_in, 3),
            COLUMNS[3]:  round(lw_out, 3),
            COLUMNS[4]:  round(pressure_base + 0.15 * math.sin(2 * math.pi * h / 12)
                               + rng.gauss(0, 0.01), 2),
            COLUMNS[5]:  rain1,
            COLUMNS[6]:  rain2,
            COLUMNS[7]:  round(rh1, 1),
            COLUMNS[8]:  round(rh2, 1),
            COLUMNS[9]:  round(net, 3),
            COLUMNS[10]: round(shf_2, 2),
            COLUMNS[11]: round(shf_surf, 2),
            COLUMNS[12]: round(vwc[0], 3),
            COLUMNS[13]: round(vwc[1], 3),
            COLUMNS[14]: round(vwc[2], 3),
            COLUMNS[15]: round(sw_in, 3),
            COLUMNS[16]: round(sw_out, 3),
            COLUMNS[17]: round(air1, 2),
            COLUMNS[18]: round(air2, 2),
            COLUMNS[19]: round(ts1, 3),
            COLUMNS[20]: round(ts2, 2),
            COLUMNS[21]: round(ts3, 2),
            COLUMNS[22]: round(ts4, 2),
            COLUMNS[23]: round(surf, 2),
            COLUMNS[24]: round(vp1, 3),
            COLUMNS[25]: round(vp2, 3),
            COLUMNS[26]: round(es1 - vp1, 3),
            COLUMNS[27]: round(es2 - vp2, 3),
            COLUMNS[28]: round(wind_dir, 1),
            COLUMNS[29]: round(wind_speed, 3),
        })

    # Random blank cells (never the timestamp or station ID)
    for row in rows:
        for c in COLUMNS[2:]:
            if rng.random() < MISSING_RATE:
                row[c] = ""

    # Plant the errors
    by_ts = {row[COLUMNS[0]]: row for row in rows}
    for (esid, var, start, end, value, _desc) in PLANTED_ERRORS:
        if esid != sid:
            continue
        cur = start
        while cur <= end:
            by_ts[cur.strftime("%Y-%m-%d %H:%M:%S")][COL[var]] = value
            cur += STEP

    rows.reverse()   # real Mesonet exports are newest-first
    return rows


def main():
    rng = random.Random(SEED)
    STATION_DIR.mkdir(parents=True, exist_ok=True)

    with open(DATA_DIR / "station_metadata.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["station_id", "location", "name", "full_name", "lat", "lng",
                    "elevation", "vegh", "skn", "nws_id", "status"])
        for sid, name, lat, lng, elev, _ in STATIONS:
            w.writerow([sid, "hawaii", name, f"{name} (synthetic)", lat, lng,
                        elev, 0.1, "", "", "active"])

    for sid, name, _lat, _lng, elev, base_temp in STATIONS:
        rows = generate_station(rng, sid, elev, base_temp)
        out = STATION_DIR / f"{sid}.csv"
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(rows)
        print(f"  Wrote {out.relative_to(ROOT)}  ({len(rows):,} rows)")

    print(f"\n  Planted {len(PLANTED_ERRORS)} error periods:")
    for sid, var, start, end, value, desc in PLANTED_ERRORS:
        span = start.strftime("%m-%d %H:%M")
        if end != start:
            span += " to " + end.strftime("%m-%d %H:%M")
        print(f"    {sid}  {var:<13} = {value:<8} {span:<28} {desc}")


if __name__ == "__main__":
    main()
