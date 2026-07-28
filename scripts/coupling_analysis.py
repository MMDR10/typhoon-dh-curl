"""
天災(颱風 dH_curl) × 地災(USGS earthquakes) × 水災(ERA5 divergence proxy)
三域藕合探測 — MVP v1
"""
import csv, json, os, sys
import numpy as np
from datetime import datetime, timedelta
import urllib.request
import xarray as xr

# ========== 1. TYPHOON LANDFALL DATES ==========
# Known landfall dates (UTC)
typhoon_landfalls = {
    "meranti_2016": datetime(2016, 9, 14, 0, 0),   # Batanes/Taiwan
    "hato_2017":    datetime(2017, 8, 23, 0, 0),    # HK/Macau
    "mangkhut_2018":datetime(2018, 9, 16, 0, 0),    # HK
    "hagibis_2019": datetime(2019, 10, 12, 0, 0),   # Japan
    "goni_2020":    datetime(2020, 11, 1, 0, 0),    # Philippines
    "saola_2023":   datetime(2023, 9, 1, 0, 0),     # HK
}

DATA_DIR = "projects/typhoon-backtest/data/6typhoon_results"
results = {}

# ========== 2. LOAD dH_curl DATA ==========
for fname in sorted(os.listdir(DATA_DIR)):
    if not fname.endswith("_dh_curl.csv"):
        continue
    parts = fname.replace(".csv", "").split("_")
    dh_idx = next(i for i, p in enumerate(parts) if p == "dh")
    name = "_".join(parts[dh_idx-2:dh_idx])  # e.g. hato_2017
    path = os.path.join(DATA_DIR, fname)
    hours, dh_curl, dh_roll = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            hours.append(float(row["hours"]))
            dh_curl.append(float(row["dH_curl"]))
            dh_roll.append(float(row["dH_roll"]))
    results[name] = {
        "hours": np.array(hours),
        "dH_curl": np.array(dh_curl),
        "dH_roll": np.array(dh_roll),
    }

# ========== 3. QUERY USGS EARTHQUAKE CATALOG ==========
def query_usgs(start_dt, end_dt, min_mag=4.0, lat_range=(0, 50), lon_range=(100, 180)):
    """Query USGS earthquake catalog"""
    url = (
        f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
        f"&starttime={start_dt.strftime('%Y-%m-%d')}"
        f"&endtime={end_dt.strftime('%Y-%m-%d')}"
        f"&minlatitude={lat_range[0]}&maxlatitude={lat_range[1]}"
        f"&minlongitude={lon_range[0]}&maxlongitude={lon_range[1]}"
        f"&minmagnitude={min_mag}"
    )
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
        events = []
        for feat in data["features"]:
            p = feat["properties"]
            g = feat["geometry"]["coordinates"]
            t_ms = p["time"]
            if t_ms:
                dt = datetime.utcfromtimestamp(t_ms / 1000.0)
            else:
                continue
            events.append({
                "time": dt,
                "mag": p["mag"],
                "depth": g[2],
                "place": p["place"],
                "lat": g[1],
                "lon": g[0],
            })
        return events
    except Exception as e:
        print(f"  USGS query error: {e}", file=sys.stderr)
        return []

# ========== 4. BIN EARTHQUAKES RELATIVE TO LANDFALL ==========
def bin_seismic_relative(events, landfall_dt, window_h=168):
    """Bin earthquake counts by hour relative to landfall"""
    bins = np.zeros(window_h * 2 + 1)  # -window to +window
    hours_rel = np.arange(-window_h, window_h + 1)
    
    for e in events:
        delta = (e["time"] - landfall_dt).total_seconds() / 3600
        idx = int(round(delta)) + window_h
        if 0 <= idx < len(bins):
            bins[idx] += 1
    
    # Also compute cumulative energy (mag-weighted)
    energy_bins = np.zeros_like(bins)
    for e in events:
        delta = (e["time"] - landfall_dt).total_seconds() / 3600
        idx = int(round(delta)) + window_h
        if 0 <= idx < len(bins):
            energy_bins[idx] += 10 ** (1.5 * e["mag"])  # seismic moment proxy
    
    return hours_rel, bins, energy_bins

# ========== 5. CROSS-CORRELATION ANALYSIS ==========
def cross_correlate(sig1, sig2, max_lag=72):
    """Cross-correlation between two signals"""
    # Ensure same length
    n = min(len(sig1), len(sig2))
    s1 = sig1[:n] - np.mean(sig1[:n])
    s2 = sig2[:n] - np.mean(sig2[:n])
    
    lags = np.arange(-max_lag, max_lag + 1)
    cors = np.zeros(len(lags))
    for i, lag in enumerate(lags):
        if lag < 0:
            cors[i] = np.corrcoef(s1[-lag:], s2[:lag])[0, 1]
        elif lag > 0:
            cors[i] = np.corrcoef(s1[:-lag], s2[lag:])[0, 1]
        else:
            cors[i] = np.corrcoef(s1, s2)[0, 1]
    return lags, cors

# ========== 6. MAIN ANALYSIS ==========
print("=" * 75)
print("天災(dH_curl) × 地災(USGS M≥4) × 水災(divergence) 藕合探測")
print("=" * 75)

coupling_results = {}

for typhoon_key, landfall_dt in typhoon_landfalls.items():
    print(f"\n{'='*50}")
    print(f"🌀 {typhoon_key}  (landfall: {landfall_dt.strftime('%Y-%m-%d')})")
    print(f"{'='*50}")
    
    if typhoon_key not in results:
        print(f"  ❌ No dH_curl data for {typhoon_key}")
        continue
    
    d = results[typhoon_key]
    
    # ---- 3a. Query USGS ----
    window_days = 10
    start_dt = landfall_dt - timedelta(days=window_days)
    end_dt = landfall_dt + timedelta(days=window_days)
    
    events = query_usgs(start_dt, end_dt, min_mag=4.0)
    print(f"  📡 USGS: {len(events)} M≥4 events in ±{window_days}d window")
    
    # Bin seismic by hour
    h_seis, seis_count, seis_energy = bin_seismic_relative(events, landfall_dt, window_h=window_days*24)
    
    # ---- 3b. Align dH_curl with seismic ----
    # dH_curl is already relative to landfall (hours)
    # Need to resample to common grid
    common_hours = np.arange(-window_days*24, window_days*24 + 1)
    
    # Interpolate dH_curl to common grid
    dh_common = np.interp(common_hours, d["hours"], d["dH_roll"])
    
    # ---- 3c. Cross-correlation ----
    lags, cors_dh_seis = cross_correlate(dh_common, seis_count, max_lag=72)
    lags2, cors_dh_energy = cross_correlate(dh_common, seis_energy, max_lag=72)
    
    peak_idx = np.argmax(np.abs(cors_dh_seis))
    peak_lag = lags[peak_idx]
    peak_cor = cors_dh_seis[peak_idx]
    
    peak_idx_e = np.argmax(np.abs(cors_dh_energy))
    peak_lag_e = lags2[peak_idx_e]
    peak_cor_e = cors_dh_energy[peak_idx_e]
    
    print(f"  🔗 dH_roll × seismic_count:  peak r={peak_cor:+.4f} @ lag={peak_lag}h")
    print(f"  🔗 dH_roll × seismic_energy: peak r={peak_cor_e:+.4f} @ lag={peak_lag_e}h")
    
    # ---- 3d. Phase analysis ----
    # Pre-landfall vs post-landfall seismic activity
    pre_mask = common_hours < 0
    post_mask = common_hours >= 0
    pre_seis = seis_count[pre_mask].sum()
    post_seis = seis_count[post_mask].sum()
    pre_energy = seis_energy[pre_mask].sum()
    post_energy = seis_energy[post_mask].sum()
    
    seis_ratio = post_seis / (pre_seis + 0.001)
    energy_ratio = post_energy / (pre_energy + 0.001)
    
    print(f"  📊 Pre-landfall seismic: {pre_seis:.0f} events, energy={pre_energy:.1e}")
    print(f"  📊 Post-landfall seismic: {post_seis:.0f} events, energy={post_energy:.1e}")
    print(f"  📊 Post/Pre ratio: events={seis_ratio:.2f}, energy={energy_ratio:.2f}")
    
    coupling_results[typhoon_key] = {
        "n_events": len(events),
        "peak_cor_dh_seis": float(peak_cor),
        "peak_lag_dh_seis": int(peak_lag),
        "peak_cor_dh_energy": float(peak_cor_e),
        "peak_lag_dh_energy": int(peak_lag_e),
        "pre_seis_events": int(pre_seis),
        "post_seis_events": int(post_seis),
        "post_pre_ratio": float(seis_ratio),
        "post_pre_energy_ratio": float(energy_ratio),
        "common_hours": common_hours.tolist(),
        "seis_count": seis_count.tolist(),
        "seis_energy": seis_energy.tolist(),
        "dh_roll_common": dh_common.tolist(),
        "lags": lags.tolist(),
        "cors_dh_seis": cors_dh_seis.tolist(),
        "cors_dh_energy": cors_dh_energy.tolist(),
    }

# ========== 7. SUMMARY ==========
print("\n" + "=" * 75)
print("📊 CROSS-TYPHOON COUPLING SUMMARY")
print("=" * 75)
print(f"{'Typhoon':<20} {'n_eq':>6} {'peak_r':>8} {'lag':>6} {'pre_eq':>7} {'post_eq':>7} {'ratio':>7}")
print("-" * 65)
for key, cr in coupling_results.items():
    name = key.split("_")[0]
    print(f"{name:<20} {cr['n_events']:>6d} {cr['peak_cor_dh_seis']:>+8.4f} {cr['peak_lag_dh_seis']:>+5d}h {cr['pre_seis_events']:>7d} {cr['post_seis_events']:>7d} {cr['post_pre_ratio']:>7.2f}")

# Save results
output_path = "projects/typhoon-backtest/output/coupling_analysis_v1.json"
os.makedirs("projects/typhoon-backtest/output", exist_ok=True)
with open(output_path, "w") as f:
    json.dump(coupling_results, f, indent=2, default=str)
print(f"\n✅ Results saved to {output_path}")
