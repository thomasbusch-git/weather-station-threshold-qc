# Weather Station Threshold QC

Threshold-based error detection for weather station data.
The scripts flag readings outside each sensor's valid range, group nearby bad
readings into events, and produce an HTML/JSON error report plus time-series
graphs that highlight the bad periods.

The repo ships with **synthetic station data that has errors planted in it**,
so you can run everything and see a full report without any real data.

## Repository layout

```
weather-station-threshold-qc/
├── generate_sample_data.py      # rebuilds the synthetic data in data/
├── requirements.txt
├── data/
│   ├── station_metadata.csv     # station lat/lon/elevation (for sunrise/sunset)
│   └── stationdata/             # one CSV per station (9001, 9002, 9003)
└── threshold/
    ├── Threshold_Errors.py          # detection → HTML + JSON report
    ├── Graph_Threshold_Errors.py    # PNG graphs with bad periods shaded red
    ├── threshold_graph.ipynb        # interactive Plotly chart for one station
    ├── Run_Threshold_Errors.bat     # Windows double-click shortcuts
    └── Run_Graph_Threshold_Errors.bat
```

## Quick start

```bash
pip install -r requirements.txt

cd threshold
python Threshold_Errors.py -i ../data/stationdata -o output
python Graph_Threshold_Errors.py -i ../data/stationdata -o output/graphs
```

Open `threshold/output/threshold_output.html` to see the report.
On Windows you can double-click the two `.bat` files instead.

For the interactive chart, open `threshold/threshold_graph.ipynb` in Jupyter
**after** running `Threshold_Errors.py` (it reads the JSON report), set
`STATION_ID` and `VARIABLE` in the first code cell, and run all cells.

## How detection works

1. Each variable is checked against a valid range (the `THRESHOLDS` table in
   `Threshold_Errors.py`), e.g. air temperature −10 to 40 °C, RH 0–100 %.
2. Incoming/outgoing shortwave radiation is only checked between sunrise and
   sunset, computed per station from `station_metadata.csv`, since small
   negative readings at night are normal.
3. Bad readings within **1 day** of each other are merged into one event;
   each event keeps its separate bursts as sub-periods.
4. Derived variables that depend on a bad reading (e.g. vapor pressure from
   air temperature) are listed as affected.

## Sample data

`generate_sample_data.py` creates 14 days (March 1–14, 2026) of 5-minute data
for three fictional stations, using the same column headers as Mesonet CSV
exports. It uses a fixed random seed, so the output is the same every time.
These errors are planted:

| Station | Variable | Bad value | When | What it simulates |
|---|---|---|---|---|
| 9001 | Air temp, sensor 2 | 55.3 °C | Mar 10 14:35 | Single-reading spike |
| 9001 | Incoming shortwave | −25 W/m² | Mar 4 10:00–12:00 | Sensor fault in daylight |
| 9001 | RH, sensor 1 | 103 / 101.5 % | Mar 7 02:00–03:30 and 05:00–05:45 | Two fog bursts → **one** event |
| 9001 | RH, sensor 1 | 102 % | Mar 12 01:00–01:30 | Days later → separate event |
| 9002 | Soil temp, sensor 3 | −6999 °C | Mar 8, all day | Logger missing-value code |
| 9002 | Air pressure | 0 kPa | Mar 5 13:00–13:25 | Sensor dropout |
| 9002 | Wind speed | 99 m/s | Mar 3 15:00 and 18:20 | Two spikes → **one** event |
| 9002 | Wind speed | 99 m/s | Mar 10 09:00 | Days later → separate event |
| 9003 | Soil moisture, sensor 2 | 1.2 m³/m³ | Mar 6 06:00–12:00 | Above physical limit |
| 9003 | Incoming longwave | 150 W/m² | Mar 9 20:00–23:00 | Too low |
| 9003 | Rainfall, sensor 1 | 999 mm | Mar 11 16:40 | Impossible 5-min total |
| 9003 | Wind direction | 400° | Mar 2 08:15–08:30 | Over 360° |

**Expected result:** 12 threshold events (4 per station). Nighttime negative
shortwave values in the data are *not* flagged. `soil_hf_1` is reported as
not found because the sample files don't include that sensor.

## Using your own data

Point `-i` at a folder of station CSVs with the same column headers (see
`VARIABLES` in `Threshold_Errors.py` for accepted header names). Put your
station metadata at `data/station_metadata.csv`. Real data in
`data/stationdata/` is ignored by git; only the `900x` sample files are tracked.

## Requirements

Python 3.10+ and the packages in `requirements.txt`
(polars, ephem, matplotlib, numpy, pandas, plotly).
