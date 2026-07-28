"""
天-地-水藕合 v2: Spatial filtering + precipitation
- Filter earthquakes to within 500km of typhoon track
- Add precipitation data via open-meteo historical API
"""
import csv, json, os, sys, math
import numpy as np
from datetime import datetime, timedelta
import urllib.request

# ========== TYPHOON TRACKS (approximate) ==========
# Format: (lat, lon) points defining the track near landfall
typhoon_info = {
    "meranti_2016": {
        "landfall": datetime(2016, 9, 14, 0),
        "track": [(15.0, 128.0), (18.0, 125.0), (20.0, 122.0), (22.0, 119.0)],  # Philippines→Taiwan
        "center": (20.0, 122.0),
    },
    "hato_2017": {
        "landfall": datetime(2017, 8, 23, 0),
        "track": [(18.0, 118.0), (20.0, 115.0), (22.0, 114.0), (23.0, 112.0)],  # SCS→HK
        "center": (22.0, 114.0),
    },
    "mangkhut_2018": {
        "landfall": datetime(2018, 9, 16, 0),
        "track": [(14.0, 130.0), (17.0, 122.0), (19.0, 116.0), (22.0, 114.0)],  # PH→HK
        "center": (22.0, 114.0),
    },
    "hagibis_2019": {
        "landfall": datetime(2019, 10, 12, 0),
        "track": [(25.0, 142.0), (30.0, 140.0), (33.0, 139.0), (35.0, 140.0)],  # Pacific→Japan
        "center": (35.0, 140.0),
    },
    "goni_2020": {
        "landfall": datetime(2020, 11, 1, 0),
        "track": [(12.0, 130.0), (13.0, 126.0), (14.0, 124.0), (15.0, 120.0)],  # Pacific→PH
        "center": (14.0, 124.0),
    },
    "saola_2023": {
        "landfall": datetime(2023, 9, 1, 0),
        "track": [(18.0, 122.0), (20.0, 118.0), (22.0, 114.0), (23.0, 113.0)],  # PH→HK
        "center": (22.0, 114.0),
    },
}

RADIUS_KM = 500.0

# ========== HELPERS ==========
def haversine_km(lat1, lon1, lat2, lon2):
    """Distance between two points in km"""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

def dist_to_track(lat, lon, track):
    """Minimum distance from point to typhoon track (line segments)"""
    min_dist = float('inf')
    for i in range(len(track) - 1):
        lat1, lon1 = track[i]
        lat2, lon2 = track[i+1]
        
        # Project point onto line segment
        dx = lon2 - lon1
        dy = lat2 - lat1
        if dx == 0 and dy == 0:
            d = haversine_km(lat, lon, lat1, lon1)
        else:
            t = ((lon - lon1)*dx + (lat - lat1)*dy) / (dx*dx + dy*dy)
            t = max(0, min(1, t))
            proj_lon = lon1 + t * dx
            proj_lat = lat1 + t * dy
            d = haversine_km(lat, lon, proj_lat, proj_lon)
        min_dist = min(min_dist, d)
    return min_dist

def query_usgs(start_dt, end_dt, min_mag=4.0):
    """Query USGS earthquake catalog for broad NW Pacific"""
    url = (
        f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
        f"&starttime={start_dt.strftime('%Y-%m-%d')}"
        f"&endtime={end_dt.strftime('%Y-%m-%d')}"
        f"&minlatitude=0&maxlatitude=50&minlongitude=100&maxlongitude=180"
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
                "time": dt, "mag": p["mag"], "depth": g[2],
                "place": p["place"], "lat": g[1], "lon": g[0],
            })
        return events
    except Exception as e:
        print(f"  USGS error: {e}", file=sys.stderr)
        return []

def query_openmeteo_precip(lat, lon, start_date, end_date):
    """Query open-meteo historical for daily precipitation"""
    url = (
        f"https://archive-api.open-meteo.com/v1/archive"
        f"?latitude={lat}&longitude={lon}"
        f"&start_date={start_date.strftime('%Y-%m-%d')}"
        f"&end_date={end_date.strftime('%Y-%m-%d')}"
        f"&daily=precipitation_sum"
        f"&timezone=UTC"
    )
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
        dates = data["daily"]["time"]
        precip = data["daily"]["precipitation_sum"]
        return list(zip(dates, precip))
    except Exception as e:
        print(f"  open-meteo error: {e}", file=sys.stderr)
        return []

# ========== LOAD dH_curl DATA ==========
DATA_DIR = "projects/typhoon-backtest/data/6typhoon_results"
dh_data = {}
for fname in sorted(os.listdir(DATA_DIR)):
    if not fname.endswith("_dh_curl.csv"):
        continue
    parts = fname.replace(".csv", "").split("_")
    dh_idx = next(i for i, p in enumerate(parts) if p == "dh")
    name = "_".join(parts[dh_idx-2:dh_idx])
    path = os.path.join(DATA_DIR, fname)
    hours, dh_curl, dh_roll = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            hours.append(float(row["hours"]))
            dh_curl.append(float(row["dH_curl"]))
            dh_roll.append(float(row["dH_roll"]))
    dh_data[name] = {"hours": np.array(hours), "dH_curl": np.array(dh_curl), "dH_roll": np.array(dh_roll)}

# ========== MAIN ==========
print("=" * 75)
print("天災(dH_curl) × 地災(spatial-filtered) × 水災(precip) 藕合 v2")
print("=" * 75)

all_results = {}

for typhoon_key, info in typhoon_info.items():
    landfall_dt = info["landfall"]
    track = info["track"]
    center = info["center"]
    short = typhoon_key.split("_")[0]
    
    print(f"\n{'='*55}")
    print(f"🌀 {short}  center={center}  landfall={landfall_dt.strftime('%Y-%m-%d')}")
    print(f"{'='*55}")
    
    if typhoon_key not in dh_data:
        print(f"  ❌ No dH_curl data")
        continue
    
    d = dh_data[typhoon_key]
    window_days = 10
    
    # --- SEISMIC: spatial-filtered ---
    start_dt = landfall_dt - timedelta(days=window_days)
    end_dt = landfall_dt + timedelta(days=window_days)
    
    all_events = query_usgs(start_dt, end_dt, min_mag=4.0)
    
    # Filter by distance to track
    near_events = []
    far_events = []
    for e in all_events:
        d_km = dist_to_track(e["lat"], e["lon"], track)
        if d_km <= RADIUS_KM:
            near_events.append(e)
        else:
            far_events.append(e)
    
    print(f"  📡 All USGS M≥4: {len(all_events)}  |  ≤{RADIUS_KM:.0f}km: {len(near_events)}  |  >{RADIUS_KM:.0f}km: {len(far_events)}")
    
    # Bin near events by hour
    bin_hours = np.arange(-window_days*24, window_days*24 + 1)
    near_bins = np.zeros(len(bin_hours))
    near_energy = np.zeros(len(bin_hours))
    for e in near_events:
        delta = (e["time"] - landfall_dt).total_seconds() / 3600
        idx = int(round(delta)) + window_days*24
        if 0 <= idx < len(bin_hours):
            near_bins[idx] += 1
            near_energy[idx] += 10 ** (1.5 * e["mag"])
    
    # Pre/post split
    pre_mask = bin_hours < 0
    post_mask = bin_hours >= 0
    pre_n = near_bins[pre_mask].sum()
    post_n = near_bins[post_mask].sum()
    pre_e = near_energy[pre_mask].sum()
    post_e = near_energy[post_mask].sum()
    
    # Cross-correlation
    dh_common = np.interp(bin_hours, d["hours"], d["dH_roll"])
    
    max_lag = 72
    lags = np.arange(-max_lag, max_lag + 1)
    cors = np.zeros(len(lags))
    n = len(bin_hours)
    s1 = dh_common - np.mean(dh_common)
    s2 = near_bins - np.mean(near_bins)
    for i, lag in enumerate(lags):
        if lag < 0:
            if -lag < n:
                cors[i] = np.corrcoef(s1[-lag:], s2[:lag])[0, 1]
        elif lag > 0:
            if lag < n:
                cors[i] = np.corrcoef(s1[:-lag], s2[lag:])[0, 1]
        else:
            cors[i] = np.corrcoef(s1, s2)[0, 1]
    
    peak_idx = np.nanargmax(np.abs(cors))
    peak_lag = lags[peak_idx]
    peak_cor = cors[peak_idx]
    
    ratio = post_n / (pre_n + 0.001)
    e_ratio = post_e / (pre_e + 0.001)
    
    print(f"  🔗 Spatial-filtered: peak r={peak_cor:+.4f} @ lag={peak_lag:+d}h")
    print(f"  📊 Pre={pre_n:.0f} events, Post={post_n:.0f} events, ratio={ratio:.2f}")
    print(f"  📊 Energy: pre={pre_e:.1e}, post={post_e:.1e}, ratio={e_ratio:.2f}")
    
    # --- PRECIPITATION ---
    lat, lon = center
    precip_data = query_openmeteo_precip(lat, lon, start_dt, end_dt)
    
    if precip_data:
        # Convert to relative days
        precip_rel = []
        for date_str, val in precip_data:
            dt = datetime.strptime(date_str, "%Y-%m-%d")
            day_rel = (dt - landfall_dt).days
            precip_rel.append((day_rel, val))
        
        total_precip = sum(v for _, v in precip_rel)
        max_precip = max(v for _, v in precip_rel)
        pre_precip = sum(v for d, v in precip_rel if d < 0)
        post_precip = sum(v for d, v in precip_rel if d >= 0)
        
        print(f"  🌧️  Precip @ ({lat},{lon}): total={total_precip:.1f}mm, max={max_precip:.1f}mm")
        print(f"  🌧️  Pre-landfall={pre_precip:.1f}mm, Post-landfall={post_precip:.1f}mm")
        
        # Correlate daily precip with daily dH_roll
        daily_dh = []
        daily_precip_vals = []
        for day_rel, pval in precip_rel:
            hour_start = day_rel * 24
            hour_end = hour_start + 24
            mask = (d["hours"] >= hour_start) & (d["hours"] < hour_end)
            if mask.sum() > 0:
                daily_dh.append(d["dH_roll"][mask].mean())
                daily_precip_vals.append(pval)
        
        if len(daily_dh) > 3:
            dh_precip_cor = np.corrcoef(daily_dh, daily_precip_vals)[0, 1]
            print(f"  🔗 dH_roll × daily_precip: r={dh_precip_cor:+.4f}")
        else:
            dh_precip_cor = np.nan
    else:
        precip_rel = []
        total_precip = max_precip = pre_precip = post_precip = dh_precip_cor = np.nan
        print(f"  🌧️  No precipitation data")
    
    all_results[typhoon_key] = {
        "n_all_eq": len(all_events),
        "n_near_eq": len(near_events),
        "n_far_eq": len(far_events),
        "near_ratio": float(ratio),
        "near_energy_ratio": float(e_ratio),
        "peak_cor_spatial": float(peak_cor),
        "peak_lag_spatial": int(peak_lag),
        "total_precip_mm": float(total_precip) if not np.isnan(total_precip) else None,
        "max_precip_mm": float(max_precip) if not np.isnan(max_precip) else None,
        "precip_pre_post_ratio": float(post_precip/pre_precip) if pre_precip > 0 else None,
        "cor_dh_precip": float(dh_precip_cor) if not np.isnan(dh_precip_cor) else None,
    }

# ========== SUMMARY ==========
print(f"\n{'='*75}")
print("📊 SPATIAL-FILTERED COUPLING SUMMARY")
print(f"{'='*75}")
print(f"{'Typhoon':<12} {'all_eq':>6} {'near':>5} {'ratio':>7} {'e_ratio':>8} {'peak_r':>8} {'lag':>5} {'precip':>8} {'dh×prcp':>8}")
print("-" * 75)
for key, cr in all_results.items():
    name = key.split("_")[0]
    prec_str = f"{cr['total_precip_mm']:.0f}mm" if cr['total_precip_mm'] else "N/A"
    dhpc_str = f"{cr['cor_dh_precip']:+.4f}" if cr['cor_dh_precip'] is not None else "N/A"
    print(f"{name:<12} {cr['n_all_eq']:>6d} {cr['n_near_eq']:>5d} {cr['near_ratio']:>7.2f} {cr['near_energy_ratio']:>8.2f} {cr['peak_cor_spatial']:>+8.4f} {cr['peak_lag_spatial']:>+4d}h {prec_str:>8} {dhpc_str:>8}")

# Save
os.makedirs("projects/typhoon-backtest/output", exist_ok=True)
with open("projects/typhoon-backtest/output/coupling_analysis_v2.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\n✅ Saved to coupling_analysis_v2.json")
