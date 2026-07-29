"""
地-水藕合 v3: Seismic × Precipitation cross-correlation
- Test whether precipitation leads seismic activity (水→地) or vice versa (地→水)
- Uses same 6 typhoon temporal windows, USGS earthquake catalog, open-meteo precipitation
- Two spatial scales: all NW Pacific + ≤500km from track
- Lagged cross-correlation: precip lead/lag seismic by up to ±10 days
"""
import csv, json, os, sys, math
import numpy as np
from datetime import datetime, timedelta
import urllib.request

# ========== TYPHOON INFO ==========
typhoon_info = {
    "meranti_2016": {
        "landfall": datetime(2016, 9, 14, 0),
        "track": [(15.0, 128.0), (18.0, 125.0), (20.0, 122.0), (22.0, 119.0)],
        "center": (20.0, 122.0),
    },
    "hato_2017": {
        "landfall": datetime(2017, 8, 23, 0),
        "track": [(18.0, 118.0), (20.0, 115.0), (22.0, 114.0), (23.0, 112.0)],
        "center": (22.0, 114.0),
    },
    "mangkhut_2018": {
        "landfall": datetime(2018, 9, 16, 0),
        "track": [(14.0, 130.0), (17.0, 122.0), (19.0, 116.0), (22.0, 114.0)],
        "center": (22.0, 114.0),
    },
    "hagibis_2019": {
        "landfall": datetime(2019, 10, 12, 0),
        "track": [(25.0, 142.0), (30.0, 140.0), (33.0, 139.0), (35.0, 140.0)],
        "center": (35.0, 140.0),
    },
    "goni_2020": {
        "landfall": datetime(2020, 11, 1, 0),
        "track": [(12.0, 130.0), (13.0, 126.0), (14.0, 124.0), (15.0, 120.0)],
        "center": (14.0, 124.0),
    },
    "saola_2023": {
        "landfall": datetime(2023, 9, 1, 0),
        "track": [(18.0, 122.0), (20.0, 118.0), (22.0, 114.0), (23.0, 113.0)],
        "center": (22.0, 114.0),
    },
}

RADIUS_KM = 500.0
WINDOW_DAYS = 10

# ========== HELPERS ==========
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

def dist_to_track(lat, lon, track):
    min_dist = float('inf')
    for i in range(len(track) - 1):
        lat1, lon1 = track[i]
        lat2, lon2 = track[i+1]
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
    url = (
        f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
        f"&starttime={start_dt.strftime('%Y-%m-%d')}"
        f"&endtime={end_dt.strftime('%Y-%m-%d')}"
        f"&minlatitude=0&maxlatitude=50&minlongitude=100&maxlongitude=180"
        f"&minmagnitude={min_mag}"
    )
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v3/1.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
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
    from urllib.parse import urlencode
    params = {
        "latitude": lat, "longitude": lon,
        "start_date": start_date.strftime('%Y-%m-%d'),
        "end_date": end_date.strftime('%Y-%m-%d'),
        "daily": "precipitation_sum", "timezone": "UTC",
    }
    url = "https://archive-api.open-meteo.com/v1/archive?" + urlencode(params)
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v3/1.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
        if "error" in data:
            print(f"  open-meteo API error: {data}", file=sys.stderr)
            return []
        dates = data["daily"]["time"]
        precip = data["daily"]["precipitation_sum"]
        return list(zip(dates, precip))
    except Exception as e:
        print(f"  open-meteo error: {e}", file=sys.stderr)
        return []

# ========== CROSS-CORRELATION ==========
def lagged_crosscorrelation(x, y, max_lag_days):
    """
    Compute lagged cross-correlation between two 1D arrays.
    x = precipitation (daily), y = earthquake count (daily)
    Positive lag = x leads y (precip → seismic, i.e. 水→地)
    Negative lag = y leads x (seismic → precip, i.e. 地→水)
    """
    lags = np.arange(-max_lag_days, max_lag_days + 1)
    cors = np.zeros(len(lags))
    n = len(x)
    for i, lag in enumerate(lags):
        if lag < 0:
            if -lag < n:
                cors[i] = np.corrcoef(x[-lag:], y[:lag])[0, 1] if len(x[-lag:]) > 2 else np.nan
        elif lag > 0:
            if lag < n:
                cors[i] = np.corrcoef(x[:-lag], y[lag:])[0, 1] if len(x[:-lag]) > 2 else np.nan
        else:
            cors[i] = np.corrcoef(x, y)[0, 1] if n > 2 else np.nan
    return lags, cors

# ========== MAIN ==========
print("=" * 75)
print("地災(地震) × 水災(降水) 藕合 v3")
print("測試：暴雨是否觸發地震？地震是否影響降水？")
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
    
    start_dt = landfall_dt - timedelta(days=WINDOW_DAYS)
    end_dt = landfall_dt + timedelta(days=WINDOW_DAYS)
    
    # --- SEISMIC ---
    all_events = query_usgs(start_dt, end_dt, min_mag=4.0)
    
    # Split by spatial filter
    near_events = []
    far_events = []
    for e in all_events:
        d_km = dist_to_track(e["lat"], e["lon"], track)
        if d_km <= RADIUS_KM:
            near_events.append(e)
        else:
            far_events.append(e)
    
    print(f"  📡 USGS M≥4: all={len(all_events)}  near(≤{RADIUS_KM:.0f}km)={len(near_events)}  far={len(far_events)}")
    
    # Bin seismic events by day (relative to landfall)
    n_days = 2 * WINDOW_DAYS + 1
    day_bins_all = np.zeros(n_days)
    day_bins_near = np.zeros(n_days)
    day_energy_all = np.zeros(n_days)
    day_energy_near = np.zeros(n_days)
    
    for e in all_events:
        day_rel = (e["time"] - landfall_dt).days
        idx = day_rel + WINDOW_DAYS
        if 0 <= idx < n_days:
            day_bins_all[idx] += 1
            day_energy_all[idx] += 10 ** (1.5 * e["mag"])
    
    for e in near_events:
        day_rel = (e["time"] - landfall_dt).days
        idx = day_rel + WINDOW_DAYS
        if 0 <= idx < n_days:
            day_bins_near[idx] += 1
            day_energy_near[idx] += 10 ** (1.5 * e["mag"])
    
    # --- PRECIPITATION ---
    lat, lon = center
    precip_data = query_openmeteo_precip(lat, lon, start_dt, end_dt)
    
    if not precip_data:
        print(f"  ❌ No precipitation data, skipping")
        continue
    
    # Align precipitation to daily bins
    precip_bins = np.zeros(n_days)
    for date_str, val in precip_data:
        dt = datetime.strptime(date_str, "%Y-%m-%d")
        day_rel = (dt - landfall_dt).days
        idx = day_rel + WINDOW_DAYS
        if 0 <= idx < n_days:
            precip_bins[idx] = val
    
    total_precip = precip_bins.sum()
    pre_precip = precip_bins[:WINDOW_DAYS].sum()
    post_precip = precip_bins[WINDOW_DAYS:].sum()
    max_precip_day = np.argmax(precip_bins) - WINDOW_DAYS
    
    print(f"  🌧️  Precip: total={total_precip:.1f}mm  pre={pre_precip:.1f}mm  post={post_precip:.1f}mm  max_day={max_precip_day:+d}")
    
    # --- CROSS-CORRELATION: precip × seismic ---
    max_lag = 10  # days
    
    # 4 combinations: (all vs near) × (count vs energy)
    combos = [
        ("all×count", precip_bins, day_bins_all),
        ("all×energy", precip_bins, day_energy_all),
        ("near×count", precip_bins, day_bins_near),
        ("near×energy", precip_bins, day_energy_near),
    ]
    
    results = {}
    
    for label, x, y in combos:
        lags, cors = lagged_crosscorrelation(x, y, max_lag)
        
        # Find peak
        valid = ~np.isnan(cors)
        if not valid.any():
            results[label] = {"peak_r": np.nan, "peak_lag": np.nan, "lag0_r": np.nan}
            continue
        
        peak_idx = np.nanargmax(np.abs(cors[valid]))
        peak_lag = lags[valid][peak_idx]
        peak_cor = cors[valid][peak_idx]
        lag0_cor = cors[lags == 0][0] if 0 in lags else np.nan
        
        # Find strongest positive (precip→seismic) and negative (seismic→precip)
        pos_mask = (lags >= 0) & valid
        neg_mask = (lags <= 0) & valid
        best_pos = (lags[pos_mask][np.nanargmax(cors[pos_mask])], np.nanmax(cors[pos_mask])) if pos_mask.any() else (np.nan, np.nan)
        best_neg = (lags[neg_mask][np.nanargmax(cors[neg_mask])], np.nanmax(cors[neg_mask])) if neg_mask.any() else (np.nan, np.nan)
        
        results[label] = {
            "peak_r": float(peak_cor),
            "peak_lag": int(peak_lag),
            "lag0_r": float(lag0_cor),
            "best_pos_lag": int(best_pos[0]) if not np.isnan(best_pos[0]) else None,
            "best_pos_r": float(best_pos[1]) if not np.isnan(best_pos[1]) else None,
            "best_neg_lag": int(best_neg[0]) if not np.isnan(best_neg[0]) else None,
            "best_neg_r": float(best_neg[1]) if not np.isnan(best_neg[1]) else None,
        }
        
        print(f"  🔗 {label:<15} peak r={peak_cor:+.4f} @ lag={peak_lag:+d}d  lag0 r={lag0_cor:+.4f}")
        print(f"     水→地(pos lag): r={best_pos[1]:+.4f} @ {best_pos[0]:+.0f}d  |  地→水(neg lag): r={best_neg[1]:+.4f} @ {best_neg[0]:+.0f}d")
    
    all_results[typhoon_key] = {
        "n_all_eq": len(all_events),
        "n_near_eq": len(near_events),
        "total_precip_mm": float(total_precip),
        "pre_precip_mm": float(pre_precip),
        "post_precip_mm": float(post_precip),
        "max_precip_day": int(max_precip_day),
        "cross_correlations": results,
    }

# ========== SUMMARY ==========
print(f"\n{'='*75}")
print("📊 地×水 藕合彙總")
print(f"{'='*75}")

print(f"\n{'Typhoon':<12} {'n_eq':>5} {'n_near':>7} {'precip':>8} {'all×ct r':>9} {'lag':>5} {'near×ct r':>9} {'lag':>5} {'all×E r':>9} {'near×E r':>9}")
print("-" * 90)
for key, cr in all_results.items():
    name = key.split("_")[0]
    c = cr["cross_correlations"]
    ac = c.get("all×count", {})
    nc = c.get("near×count", {})
    ae = c.get("all×energy", {})
    ne = c.get("near×energy", {})
    print(f"{name:<12} {cr['n_all_eq']:>5d} {cr['n_near_eq']:>7d} {cr['total_precip_mm']:>7.0f}mm "
          f"{ac.get('peak_r', 0):>+9.4f} {ac.get('peak_lag', 0):>+4d}d "
          f"{nc.get('peak_r', 0):>+9.4f} {nc.get('peak_lag', 0):>+4d}d "
          f"{ae.get('peak_r', 0):>+9.4f} {ne.get('peak_r', 0):>+9.4f}")

# --- Aggregate statistics ---
print(f"\n{'='*75}")
print("📊 聚合統計（跨 6 颱風）")
print(f"{'='*75}")

for label in ["all×count", "all×energy", "near×count", "near×energy"]:
    peak_rs = []
    lags = []
    for key, cr in all_results.items():
        if label in cr["cross_correlations"]:
            r = cr["cross_correlations"][label]["peak_r"]
            lag = cr["cross_correlations"][label]["peak_lag"]
            if not np.isnan(r):
                peak_rs.append(r)
                lags.append(lag)
    
    if peak_rs:
        mean_r = np.mean(np.abs(peak_rs))
        mean_lag = np.mean(lags)
        max_r = np.max(np.abs(peak_rs))
        # Test: are the signed correlations significantly different from zero?
        t_stat = np.mean(peak_rs) / (np.std(peak_rs, ddof=1) / np.sqrt(len(peak_rs))) if len(peak_rs) > 1 else 0
        
        print(f"  {label:<15} n={len(peak_rs)}  mean|r|={mean_r:.4f}  max|r|={max_r:.4f}  mean_lag={mean_lag:+.1f}d  t={t_stat:+.2f}")
        print(f"     individual: {[f'{r:+.3f}' for r in peak_rs]}")

# --- Directionality test ---
print(f"\n{'='*75}")
print("📊 方向性檢驗：水→地 (precip leads seismic) vs 地→水 (seismic leads precip)")
print(f"{'='*75}")

for label in ["all×count", "near×count"]:
    pos_rs = []
    neg_rs = []
    for key, cr in all_results.items():
        if label in cr["cross_correlations"]:
            c = cr["cross_correlations"][label]
            if c.get("best_pos_r") is not None and not np.isnan(c["best_pos_r"]):
                pos_rs.append(c["best_pos_r"])
            if c.get("best_neg_r") is not None and not np.isnan(c["best_neg_r"]):
                neg_rs.append(c["best_neg_r"])
    
    if pos_rs and neg_rs:
        mean_pos = np.mean(np.abs(pos_rs))
        mean_neg = np.mean(np.abs(neg_rs))
        print(f"  {label:<15} 水→地 mean|r|={mean_pos:.4f}  地→水 mean|r|={mean_neg:.4f}  Δ={mean_pos-mean_neg:+.4f}")
        print(f"     水→地: {[f'{r:+.3f}' for r in pos_rs]}")
        print(f"     地→水: {[f'{r:+.3f}' for r in neg_rs]}")

# Save
os.makedirs("projects/typhoon-backtest/output", exist_ok=True)
with open("projects/typhoon-backtest/output/coupling_analysis_v3.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\n✅ Saved to coupling_analysis_v3.json")
