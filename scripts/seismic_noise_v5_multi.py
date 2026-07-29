"""
v5 Multi-Typhoon: IRIS continuous waveform microseism noise
Fetch 1-hour BHZ windows at 6-hourly intervals, ±3 days around landfall
Match to typhoon dH_roll timeseries for coupling analysis
"""
import numpy as np
import json, os, csv, sys
from datetime import datetime, timedelta
import urllib.request

# ========== TYPHOON DEFINITIONS ==========
TYPHOONS = {
    "meranti": {
        "landfall": datetime(2016, 9, 14, 0, 0),
        "station": ("IU", "TATO", "00", "BHZ"),
        "center": (20.0, 122.0),
        "note": "Taiwan/Batanes"
    },
    "hato": {
        "landfall": datetime(2017, 8, 23, 0, 0),
        "station": ("HK", "HKPS", "00", "BHZ"),
        "center": (22.0, 114.0),
        "note": "HK/Macau"
    },
    "mangkhut": {
        "landfall": datetime(2018, 9, 16, 0, 0),
        "station": ("HK", "HKPS", "00", "BHZ"),
        "center": (22.0, 114.0),
        "note": "HK/Guangdong"
    },
    "hagibis": {
        "landfall": datetime(2019, 10, 12, 0, 0),
        "station": ("G", "INU", "00", "BHZ"),
        "center": (35.0, 140.0),
        "note": "Japan"
    },
    "goni": {
        "landfall": datetime(2020, 11, 1, 0, 0),
        "station": ("IU", "DAV", "00", "BHZ"),
        "center": (14.0, 124.0),
        "note": "Philippines (DAV station)"
    },
    "saola": {
        "landfall": datetime(2023, 9, 1, 0, 0),
        "station": ("HK", "HKPS", "00", "BHZ"),
        "center": (22.0, 114.0),
        "note": "HK/Guangdong"
    }
}

DATA_DIR = "projects/typhoon-backtest/data/6typhoon_results"
OUTPUT_DIR = "projects/typhoon-backtest/output"

# ========== HELPERS ==========
def fetch_1hour(net, sta, loc, cha, dt):
    """Fetch 1 hour of continuous BHZ data at given datetime"""
    start = dt.strftime("%Y-%m-%dT%H:%M:%S")
    end = (dt + timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%S")
    url = (
        f"http://service.iris.edu/fdsnws/dataselect/1/query"
        f"?net={net}&sta={sta}&loc={loc}&cha={cha}"
        f"&starttime={start}&endtime={end}"
    )
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v5/1.0'})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = resp.read()
        # miniSEED format - extract samples
        # FDSN dataselect returns miniSEED, we need to parse it
        # MiniSEED header is 56 bytes, then 4-byte integers
        # But actual parsing is complex. Let's try a simpler approach:
        # use the timeseries API which returns ASCII
        return None  # Will use timeseries API instead
    except Exception as e:
        return None

def fetch_1hour_ascii(net, sta, loc, cha, dt):
    """Fetch 1 hour of continuous data as ASCII via irisws/timeseries"""
    start = dt.strftime("%Y-%m-%dT%H:%M:%S")
    end = (dt + timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%S")
    url = (
        f"http://service.iris.edu/irisws/timeseries/1/query"
        f"?net={net}&sta={sta}&loc={loc}&cha={cha}"
        f"&starttime={start}&endtime={end}&output=ascii"
    )
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v5/1.0'})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = resp.read().decode()
        lines = data.strip().split("\n")
        values = []
        for line in lines[1:]:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        values.append(float(parts[-1]))
                    except (ValueError, IndexError):
                        continue
        return np.array(values) if len(values) > 100 else None
    except Exception as e:
        return None

def compute_rms(values):
    """Compute RMS noise amplitude (detrended)"""
    if values is None or len(values) < 100:
        return None
    detrended = values - np.mean(values)
    return float(np.sqrt(np.mean(detrended**2)))

def dh_curl(seq, window=2):
    """Compute dH_curl of a sequence (noise of noise)"""
    seq = np.array(seq, dtype=float)
    curl = np.full_like(seq, np.nan)
    for i in range(window, len(seq)-window):
        local = seq[i-window:i+window+1]
        curl[i] = float(np.std(np.diff(local)) / (np.mean(np.abs(local)) + 1e-8))
    return curl

def load_typhoon_dh(typhoon_key):
    """Load typhoon dH_curl data from CSV"""
    for fname in sorted(os.listdir(DATA_DIR)):
        if typhoon_key in fname.lower() and fname.endswith("_dh_curl.csv"):
            with open(os.path.join(DATA_DIR, fname)) as f:
                reader = csv.DictReader(f)
                hours = np.array([float(r["hours"]) for r in reader])
            with open(os.path.join(DATA_DIR, fname)) as f:
                reader = csv.DictReader(f)
                dh_curl_vals = np.array([float(r["dH_curl"]) for r in reader])
            with open(os.path.join(DATA_DIR, fname)) as f:
                reader = csv.DictReader(f)
                dh_roll_vals = np.array([float(r["dH_roll"]) for r in reader])
            return hours, dh_curl_vals, dh_roll_vals
    return None, None, None

# ========== MAIN ==========
os.makedirs(OUTPUT_DIR, exist_ok=True)

all_results = {}

for ty_name, ty_info in TYPHOONS.items():
    print(f"\n{'='*70}")
    print(f"🌀 {ty_name.upper()}  ({ty_info['note']})")
    print(f"   Landfall: {ty_info['landfall'].strftime('%Y-%m-%d %H:%M')} UTC")
    print(f"   Station:  {'.'.join(ty_info['station'])}")
    print(f"{'='*70}")
    
    net, sta, loc, cha = ty_info['station']
    landfall = ty_info['landfall']
    
    # Load typhoon dH_curl
    dh_hours, dh_curl_vals, dh_roll_vals = load_typhoon_dh(ty_name)
    if dh_hours is None:
        print(f"  ❌ No dH_curl data for {ty_name}")
        continue
    
    # Filter to ±3 days and subsample to 6-hourly
    WINDOW_DAYS = 3
    SUBSAMPLE = 6
    mask = np.abs(dh_hours / 24) <= WINDOW_DAYS
    filtered_hours = dh_hours[mask]
    filtered_roll = dh_roll_vals[mask]
    
    # Subsample
    sel_hours = filtered_hours[::SUBSAMPLE]
    sel_roll = filtered_roll[::SUBSAMPLE]
    n_points = len(sel_hours)
    
    print(f"   Time points: {n_points} (every ~6h, ±{WINDOW_DAYS}d)")
    print(f"   Range: {sel_hours[0]:.0f}h to {sel_hours[-1]:.0f}h from landfall")
    
    # Fetch seismic noise at each time point
    noise_rms = []
    timestamps = []
    
    for i, h in enumerate(sel_hours):
        dt = landfall + timedelta(hours=h)
        ts = dt.strftime("%Y-%m-%d %H:%M")
        
        print(f"  [{i+1}/{n_points}] t={h:+.0f}h  {ts} ...", end=" ", flush=True)
        
        values = fetch_1hour_ascii(net, sta, loc, cha, dt)
        rms = compute_rms(values)
        
        if rms is None:
            print("--")
            noise_rms.append(None)
        else:
            print(f"RMS={rms:.0f}")
            noise_rms.append(rms)
        
        timestamps.append(ts)
    
    # Analysis
    valid = [(i, r) for i, r in enumerate(noise_rms) if r is not None]
    n_valid = len(valid)
    
    print(f"\n   Valid noise measurements: {n_valid}/{n_points}")
    
    result = {
        "typhoon": f"{ty_name}_{landfall.year}",
        "station": f"{net}.{sta}.{loc}.{cha}",
        "landfall": str(landfall),
        "n_points": n_points,
        "n_valid": n_valid,
        "timestamps": timestamps,
        "noise_rms": noise_rms,
        "typhoon_dh_hours": [float(h) for h in sel_hours],
    }
    
    if n_valid < 5:
        print("   Too few valid points for analysis")
        all_results[ty_name] = result
        continue
    
    valid_rms = np.array([r for _, r in valid])
    valid_idx = np.array([i for i, _ in valid])
    valid_hours = sel_hours[valid_idx]
    valid_roll = sel_roll[valid_idx]
    
    # Basic stats
    result["rms_min"] = float(valid_rms.min())
    result["rms_max"] = float(valid_rms.max())
    result["rms_mean"] = float(valid_rms.mean())
    result["rms_std"] = float(valid_rms.std())
    
    print(f"   RMS range: {valid_rms.min():.0f} - {valid_rms.max():.0f}")
    print(f"   RMS mean:  {valid_rms.mean():.0f}")
    
    # Pre/post landfall
    lf_idx = np.argmin(np.abs(valid_hours))
    pre_mask = np.arange(len(valid_rms)) < lf_idx
    post_mask = np.arange(len(valid_rms)) >= lf_idx
    
    if pre_mask.sum() > 2 and post_mask.sum() > 2:
        pre_mean = float(valid_rms[pre_mask].mean())
        post_mean = float(valid_rms[post_mask].mean())
        ratio = post_mean / pre_mean if pre_mean > 0 else float('nan')
        
        result["pre_landfall_mean_rms"] = pre_mean
        result["post_landfall_mean_rms"] = post_mean
        result["pre_post_ratio"] = ratio
        
        # Baseline = first valid point
        result["baseline_rms"] = float(valid_rms[0])
        result["peak_rms"] = float(valid_rms.max())
        result["peak_hours"] = float(valid_hours[np.argmax(valid_rms)])
        
        print(f"   Pre-landfall:  {pre_mean:.0f} (n={int(pre_mask.sum())})")
        print(f"   Post-landfall: {post_mean:.0f} (n={int(post_mask.sum())})")
        print(f"   Ratio: {ratio:.2f}x")
    
    # Correlation: typhoon dH_roll × noise RMS
    r_ty_vs_noise = float(np.corrcoef(valid_roll, valid_rms)[0, 1])
    result["r_typhoon_dHroll_vs_RMS"] = r_ty_vs_noise
    print(f"   Typhoon dH_roll × Noise RMS: r = {r_ty_vs_noise:+.4f}")
    
    # Time correlation
    r_time_vs_noise = float(np.corrcoef(valid_hours, valid_rms)[0, 1])
    result["r_time_vs_RMS"] = r_time_vs_noise
    
    # dH_curl of noise (noise topology)
    noise_dh = dh_curl(valid_rms)
    valid_dh_mask = ~np.isnan(noise_dh)
    if valid_dh_mask.sum() > 3:
        vd = noise_dh[valid_dh_mask]
        vr = valid_roll[valid_dh_mask]
        r_dh = float(np.corrcoef(vr, vd)[0, 1])
        result["r_typhoon_dHroll_vs_noise_dHcurl"] = r_dh
        
        # Lagged cross-correlation
        max_lag = 5
        lags = np.arange(-max_lag, max_lag + 1)
        cors = []
        for lag in lags:
            if lag < 0:
                x, y = vr[-lag:], vd[:lag]
            elif lag > 0:
                x, y = vr[:-lag], vd[lag:]
            else:
                x, y = vr, vd
            if len(x) > 3:
                cors.append(float(np.corrcoef(x, y)[0, 1]))
            else:
                cors.append(float('nan'))
        
        cors = np.array(cors)
        peak_idx = int(np.nanargmax(np.abs(cors)))
        peak_lag = int(lags[peak_idx])
        peak_r = float(cors[peak_idx])
        
        result["lagged_cross_corr"] = {
            "lags": [int(l) for l in lags],
            "correlations": cors.tolist(),
            "peak_lag": peak_lag,
            "peak_r": peak_r
        }
        
        print(f"   Typhoon dH_roll × Noise dH_curl: r = {r_dh:+.4f}")
        print(f"   Peak lagged x-corr: lag={peak_lag:+d}, r={peak_r:+.4f}")
        
        # Print full lagged corr
        for l, c in zip(lags, cors):
            marker = "  <-- PEAK" if l == peak_lag else ""
            print(f"     lag={l:+d}: r={c:+.4f}{marker}")
    
    all_results[ty_name] = result

# ========== SAVE ==========
with open(os.path.join(OUTPUT_DIR, "seismic_noise_v5_multi.json"), "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\n{'='*70}")
print(f"✅ Saved to output/seismic_noise_v5_multi.json")
print(f"{'='*70}")

# ========== SUMMARY TABLE ==========
print(f"\n{'='*70}")
print("SUMMARY: Ty-phoon Seismic Noise Coupling v5")
print(f"{'='*70}")
print(f"{'Typhoon':<12} {'Station':<16} {'Pre->Post':<12} {'Peak/baseline':<16} {'r(dHroll,RMS)':<16} {'Peak lag r':<12}")
print(f"{'-'*12} {'-'*16} {'-'*12} {'-'*16} {'-'*16} {'-'*12}")

for ty_name in TYPHOONS:
    r = all_results.get(ty_name, {})
    station = r.get("station", "N/A")
    ratio = r.get("pre_post_ratio", float('nan'))
    peak_bl = r.get("peak_rms", 0) / r.get("baseline_rms", 1) if r.get("baseline_rms", 0) > 0 else float('nan')
    r_simple = r.get("r_typhoon_dHroll_vs_RMS", float('nan'))
    
    lags = r.get("lagged_cross_corr", {})
    peak_lag_r = lags.get("peak_r", float('nan')) if lags else float('nan')
    
    ratio_str = f"{ratio:.2f}x" if not np.isnan(ratio) else "N/A"
    peak_str = f"{peak_bl:.1f}x" if not np.isnan(peak_bl) else "N/A"
    r_str = f"{r_simple:+.3f}" if not np.isnan(r_simple) else "N/A"
    lag_str = f"{peak_lag_r:+.3f}" if not np.isnan(peak_lag_r) else "N/A"
    
    print(f"{ty_name:<12} {station:<16} {ratio_str:<12} {peak_str:<16} {r_str:<16} {lag_str:<12}")
