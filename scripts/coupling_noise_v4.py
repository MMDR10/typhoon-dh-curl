"""
天-地-水藕合 v4: Noise Topography 方法
- 用 inter-event interval (地震間時間差) 做 seismic noise proxy
- 用 dH_curl of precipitation sequence 做水文噪音 proxy
- 對比 typhoon dH_curl 同 seismic/precipitation noise 結構
- 信號用 dH_curl (noise structure)，唔係 raw event count
"""
import csv, json, os, sys, math
import numpy as np
from datetime import datetime, timedelta
import urllib.request

typhoon_info = {
    "meranti_2016": {"landfall": datetime(2016, 9, 14, 0), "track": [(15.0,128.0),(18.0,125.0),(20.0,122.0),(22.0,119.0)], "center": (20.0,122.0)},
    "hato_2017":    {"landfall": datetime(2017, 8, 23, 0), "track": [(18.0,118.0),(20.0,115.0),(22.0,114.0),(23.0,112.0)], "center": (22.0,114.0)},
    "mangkhut_2018":{"landfall": datetime(2018, 9, 16, 0), "track": [(14.0,130.0),(17.0,122.0),(19.0,116.0),(22.0,114.0)], "center": (22.0,114.0)},
    "hagibis_2019": {"landfall": datetime(2019, 10, 12, 0),"track": [(25.0,142.0),(30.0,140.0),(33.0,139.0),(35.0,140.0)], "center": (35.0,140.0)},
    "goni_2020":    {"landfall": datetime(2020, 11, 1, 0), "track": [(12.0,130.0),(13.0,126.0),(14.0,124.0),(15.0,120.0)], "center": (14.0,124.0)},
    "saola_2023":   {"landfall": datetime(2023, 9, 1, 0), "track": [(18.0,122.0),(20.0,118.0),(22.0,114.0),(23.0,113.0)], "center": (22.0,114.0)},
}

RADIUS_KM = 500
WINDOW = 15  # ±15 days for longer noise time series

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1))*math.cos(math.radians(lat2))*math.sin(dlon/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

def dist_to_track(lat, lon, track):
    min_dist = float('inf')
    for i in range(len(track)-1):
        lat1, lon1 = track[i]; lat2, lon2 = track[i+1]
        dx = lon2-lon1; dy = lat2-lat1
        if dx==0 and dy==0:
            d = haversine_km(lat, lon, lat1, lon1)
        else:
            t = max(0, min(1, ((lon-lon1)*dx+(lat-lat1)*dy)/(dx*dx+dy*dy)))
            d = haversine_km(lat, lon, lat1+t*dy, lon1+t*dx)
        min_dist = min(min_dist, d)
    return min_dist

def query_usgs(start_dt, end_dt, min_mag=4.0):
    url = (
        f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
        f"&starttime={start_dt.strftime('%Y-%m-%d')}"
        f"&endtime={end_dt.strftime('%Y-%m-%d')}"
        f"&minlatitude=0&maxlatitude=50&minlongitude=100&maxlongitude=180"
        f"&minmagnitude={min_mag}"
    )
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v4/1.0'})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read())
        events = []
        for feat in data["features"]:
            p = feat["properties"]; g = feat["geometry"]["coordinates"]
            t_ms = p["time"]
            if t_ms:
                events.append({
                    "time": datetime.utcfromtimestamp(t_ms/1000.0),
                    "mag": p["mag"], "depth": g[2], "lat": g[1], "lon": g[0],
                })
        return events
    except Exception as e:
        print(f"  USGS error: {e}", file=sys.stderr)
        return []

def query_openmeteo_precip(lat, lon, start_date, end_date):
    from urllib.parse import urlencode
    params = {
        "latitude": lat, "longitude": lon,
        "start_date": start_date.strftime('%Y-%m-%d'),
        "end_date": end_date.strftime('%Y-%m-%d'),
        "daily": "precipitation_sum", "timezone": "UTC",
    }
    url = "https://archive-api.open-meteo.com/v1/archive?" + urlencode(params)
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v4/1.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
        if "error" in data: return []
        return list(zip(data["daily"]["time"], data["daily"]["precipitation_sum"]))
    except Exception as e:
        print(f"  open-meteo error: {e}", file=sys.stderr)
        return []

def compute_dh_curl(seq, window=3):
    """Compute dH_curl for a 1D sequence"""
    seq = np.array(seq, dtype=float)
    n = len(seq)
    if n < 2*window+1:
        return np.full_like(seq, np.nan)
    curl = np.full(n, np.nan)
    for i in range(window, n-window):
        local = seq[i-window:i+window+1]
        local = local[~np.isnan(local)]
        if len(local) < 3:
            continue
        # dH = local complexity (permutation entropy-like)
        # curl = local directional change
        dH = np.std(np.diff(local)) / (np.mean(np.abs(local)) + 1e-8)
        curl[i] = dH
    return curl

print("="*75)
print("天-地-水藕合 v4: Noise Topology Approach")
print("用 inter-event interval + dH_curl of noise time series")
print("="*75)

all_results = {}

for typhoon_key, info in typhoon_info.items():
    landfall_dt = info["landfall"]
    track = info["track"]
    center = info["center"]
    short = typhoon_key.split("_")[0]
    
    print(f"\n{'='*55}")
    print(f"🌀 {short}  center={center}  landfall={landfall_dt.strftime('%Y-%m-%d')}")
    
    start_dt = landfall_dt - timedelta(days=WINDOW)
    end_dt = landfall_dt + timedelta(days=WINDOW)
    
    # --- 1. Load typhoon dH_curl ---
    DATA_DIR = "projects/typhoon-backtest/data/6typhoon_results"
    dh_data = None
    for fname in sorted(os.listdir(DATA_DIR)):
        if not fname.endswith("_dh_curl.csv"): continue
        parts = fname.replace(".csv","").split("_")
        dh_idx = next(i for i,p in enumerate(parts) if p=="dh")
        name = "_".join(parts[dh_idx-2:dh_idx])
        if name == typhoon_key.replace("-2016","_2016").replace("-2017","_2017").replace("-2018","_2018").replace("-2019","_2019").replace("-2020","_2020").replace("-2023","_2023"):
            hours,dh_curl,dh_roll=[],[],[]
            with open(os.path.join(DATA_DIR,fname)) as f:
                for row in csv.DictReader(f):
                    hours.append(float(row["hours"]))
                    dh_curl.append(float(row["dH_curl"]))
                    dh_roll.append(float(row["dH_roll"]))
            dh_data = {"hours": np.array(hours), "dH_curl": np.array(dh_curl), "dH_roll": np.array(dh_roll)}
            break
    
    if dh_data is None:
        print(f"  ❌ No dH_curl data, skipping")
        continue
    
    # --- 2. Seismic inter-event intervals ---
    all_events = query_usgs(start_dt, end_dt, min_mag=4.0)
    
    # Filter near events to track
    near_events = sorted([e for e in all_events if dist_to_track(e["lat"],e["lon"],track) <= RADIUS_KM], key=lambda x: x["time"])
    far_events = sorted([e for e in all_events if dist_to_track(e["lat"],e["lon"],track) > RADIUS_KM], key=lambda x: x["time"])
    
    print(f"  📡 M≥4 events: all={len(all_events)}  near({RADIUS_KM}km)={len(near_events)}")
    
    # Compute inter-event intervals (hours) - ensure sorted by time
    def compute_intervals(events):
        if len(events) < 3: return []
        events_sorted = sorted(events, key=lambda x: x["time"])
        intervals = []
        for i in range(1, len(events_sorted)):
            delta = (events_sorted[i]["time"] - events_sorted[i-1]["time"]).total_seconds() / 3600
            if delta > 0:  # guard against duplicate timestamps
                intervals.append(delta)
        return np.array(intervals)
    
    near_intervals = compute_intervals(near_events)
    far_intervals = compute_intervals(far_events)
    all_intervals = compute_intervals(all_events)
    
    print(f"     Inter-event intervals: near n={len(near_intervals)}, far n={len(far_intervals)}")
    
    # --- 3. Precipitation ---
    lat, lon = center
    precip_data = query_openmeteo_precip(lat, lon, start_dt, end_dt)
    if precip_data:
        precip_vals = np.array([v for _,v in precip_data])
        print(f"  🌧️  Precip: total={precip_vals.sum():.0f}mm  n_days={len(precip_vals)}")
    else:
        precip_vals = None
        print(f"  🌧️  No precip data")
    
    # --- 4. Noise Structure Comparison ---
    # Metric A: dH_curl of seismic intervals vs typhoon dH_curl
    results = {}
    
    for label, intervals in [("near", near_intervals), ("far", far_intervals), ("all", all_intervals)]:
        if len(intervals) < 5:
            results[f"{label}_ie_cv"] = None
            results[f"{label}_ie_mean"] = None
            results[f"{label}_ie_std"] = None
            continue
        
        # Noise metrics of inter-event intervals
        cv = np.std(intervals) / (np.mean(intervals) + 1e-8)
        skew = float(np.mean((intervals - np.mean(intervals))**3) / (np.std(intervals)**3 + 1e-8))
        kurt = float(np.mean((intervals - np.mean(intervals))**4) / (np.var(intervals)**2 + 1e-8) - 3)
        
        results[f"{label}_n_events"] = len(intervals) + 1
        results[f"{label}_ie_mean_h"] = float(np.mean(intervals))
        results[f"{label}_ie_std_h"] = float(np.std(intervals))
        results[f"{label}_ie_cv"] = float(cv)
        results[f"{label}_ie_skew"] = skew
        results[f"{label}_ie_kurt"] = kurt
        
        print(f"     {label}: CV={cv:.3f}  skew={skew:.3f}  kurt={kurt:.3f}  mean={np.mean(intervals):.1f}h")
    
    # Metric B: Seismic event rate dH_curl (noise structure of seismic time series)
    # Bin events by day
    n_days = 2*WINDOW + 1
    daily_bins_near = np.zeros(n_days)
    daily_bins_all = np.zeros(n_days)
    for e in near_events:
        day_rel = (e["time"] - landfall_dt).days + WINDOW
        if 0 <= day_rel < n_days: daily_bins_near[day_rel] += 1
    for e in all_events:
        day_rel = (e["time"] - landfall_dt).days + WINDOW
        if 0 <= day_rel < n_days: daily_bins_all[day_rel] += 1
    
    # Compute dH_curl of daily seismic rate sequence
    dh_seismic_near = compute_dh_curl(daily_bins_near, window=2)
    dh_seismic_all = compute_dh_curl(daily_bins_all, window=2)
    
    # Metric C: dH_curl of precipitation sequence
    dh_precip = compute_dh_curl(precip_vals, window=2) if precip_vals is not None else None
    
    # Cross-correlate: typhoon dH_roll vs seismic noise dH_curl (common time axis)
    # Typhoon dH data is 6-hourly, seismic/precip is daily
    # Downsample typhoon data to daily
    ty_days = np.floor(dh_data["hours"] / 24).astype(int)
    ty_unique_days = np.unique(ty_days)
    ty_daily_dh = np.array([np.mean(dh_data["dH_roll"][ty_days == d]) for d in ty_unique_days])
    
    # Align
    ty_start_day = int(np.floor(dh_data["hours"].min() / 24))
    shift = WINDOW + ty_start_day
    if shift < 0:
        pad_left = -shift
        ty_daily_dh = np.pad(ty_daily_dh, (pad_left, 0), constant_values=np.nan)[:n_days]
    elif shift > 0:
        ty_daily_dh = np.pad(ty_daily_dh, (0, shift), constant_values=np.nan)[shift:]
    
    # Truncate to same length
    min_len = min(len(ty_daily_dh), len(daily_bins_near), len(daily_bins_all))
    ty_aligned = ty_daily_dh[:min_len]
    seis_near = daily_bins_near[:min_len]
    seis_all = daily_bins_all[:min_len]
    dh_seis_near = dh_seismic_near[:min_len]
    dh_seis_all = dh_seismic_all[:min_len]
    
    # Correlations
    mask_ty_seis = ~(np.isnan(ty_aligned) | np.isnan(seis_near))
    if mask_ty_seis.sum() > 3:
        r_ty_vs_seis = float(np.corrcoef(ty_aligned[mask_ty_seis], seis_near[mask_ty_seis])[0,1])
    else:
        r_ty_vs_seis = None
    
    mask_ty_dhseis = ~(np.isnan(ty_aligned) | np.isnan(dh_seis_near))
    if mask_ty_dhseis.sum() > 3:
        r_ty_vs_dhseis = float(np.corrcoef(ty_aligned[mask_ty_dhseis], dh_seis_near[mask_ty_dhseis])[0,1])
    else:
        r_ty_vs_dhseis = None
    
    # Precip correlation
    if precip_vals is not None and dh_precip is not None:
        min_len_p = min(len(ty_aligned), len(precip_vals))
        ty_p = ty_aligned[:min_len_p]
        pr = np.array(precip_vals[:min_len_p], dtype=float)
        dp = np.array(dh_precip[:min_len_p], dtype=float)
        mask_tp = ~(np.isnan(ty_p) | np.isnan(pr))
        if mask_tp.sum() > 3:
            r_ty_vs_precip = float(np.corrcoef(ty_p[mask_tp], pr[mask_tp])[0,1])
        else:
            r_ty_vs_precip = None
        mask_tdp = ~(np.isnan(ty_p) | np.isnan(dp))
        if mask_tdp.sum() > 3:
            r_ty_vs_dhprecip = float(np.corrcoef(ty_p[mask_tdp], dp[mask_tdp])[0,1])
        else:
            r_ty_vs_dhprecip = None
    else:
        r_ty_vs_precip = r_ty_vs_dhprecip = None
    
    results.update({
        "r_ty_vs_seis_raw": r_ty_vs_seis,
        "r_ty_vs_seis_noise": r_ty_vs_dhseis,
        "r_ty_vs_precip_raw": r_ty_vs_precip,
        "r_ty_vs_precip_noise": r_ty_vs_dhprecip,
    })
    
    print(f"  🔗 颱風 dH_roll vs 地震頻率(raw): r={r_ty_vs_seis}")
    print(f"  🔗 颱風 dH_roll vs 地震雜訊(dH_curl): r={r_ty_vs_dhseis}")
    print(f"  🔗 颱風 dH_roll vs 降水(raw): r={r_ty_vs_precip}")
    print(f"  🔗 颱風 dH_roll vs 降水雜訊(dH_curl): r={r_ty_vs_dhprecip}")
    
    # Metric D: Cross-correlation between seismic noise dH_curl and precip dH_curl (地↔水 noise test)
    if dh_precip is not None and dh_seis_near is not None:
        min_len_np = min(len(dh_seis_near), len(dh_precip))
        ds = dh_seis_near[:min_len_np]
        dp = dh_precip[:min_len_np]
        mask_dsdp = ~(np.isnan(ds) | np.isnan(dp))
        if mask_dsdp.sum() > 3:
            r_seis_noise_vs_precip_noise = float(np.corrcoef(ds[mask_dsdp], dp[mask_dsdp])[0,1])
        else:
            r_seis_noise_vs_precip_noise = None
    else:
        r_seis_noise_vs_precip_noise = None
    
    results["r_seis_noise_vs_precip_noise"] = r_seis_noise_vs_precip_noise
    print(f"  🔗 地震雜訊(dH_curl) vs 降水雜訊(dH_curl): r={r_seis_noise_vs_precip_noise}")
    
    all_results[typhoon_key] = results

# ========== SUMMARY ==========
print(f"\n{'='*75}")
print("📊 Noise Topology 藕合彙總")
print(f"{'='*75}")

print(f"\n{'Typhoon':<12} {'ty_vs_seis':>12} {'ty_vs_dhseis':>14} {'ty_vs_precip':>14} {'ty_vs_dhprecip':>16} {'seis_dh_vs_pr_dh':>16}")
print("-"*90)
for key, cr in all_results.items():
    name = key.split("_")[0]
    r1 = f"{cr.get('r_ty_vs_seis_raw', 0):+.4f}" if cr.get('r_ty_vs_seis_raw') else "N/A"
    r2 = f"{cr.get('r_ty_vs_seis_noise', 0):+.4f}" if cr.get('r_ty_vs_seis_noise') else "N/A"
    r3 = f"{cr.get('r_ty_vs_precip_raw', 0):+.4f}" if cr.get('r_ty_vs_precip_raw') else "N/A"
    r4 = f"{cr.get('r_ty_vs_precip_noise', 0):+.4f}" if cr.get('r_ty_vs_precip_noise') else "N/A"
    r5 = f"{cr.get('r_seis_noise_vs_precip_noise', 0):+.4f}" if cr.get('r_seis_noise_vs_precip_noise') else "N/A"
    print(f"{name:<12} {r1:>12} {r2:>14} {r3:>14} {r4:>16} {r5:>16}")

# Aggregate
print(f"\n{'='*50}")
print("聚合統計")
noise_impacts = []
for key, cr in all_results.items():
    raw = cr.get('r_ty_vs_seis_raw')
    noise = cr.get('r_ty_vs_seis_noise')
    if raw is not None and noise is not None and not np.isnan(raw) and not np.isnan(noise):
        noise_impacts.append(noise - raw)
if noise_impacts:
    print(f"  Noise 比 raw 改善幅度: mean={np.mean(noise_impacts):+.4f}  "
          f"max={np.max(noise_impacts):+.4f}  min={np.min(noise_impacts):+.4f}")
    print(f"  個別: {[f'{d:+.4f}' for d in noise_impacts]}")

# Save
os.makedirs("projects/typhoon-backtest/output", exist_ok=True)
with open("projects/typhoon-backtest/output/coupling_analysis_v4.json","w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\n✅ Saved to coupling_analysis_v4.json")
