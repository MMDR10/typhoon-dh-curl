"""
v5 Parallel: Fetch IRIS continuous waveform noise for all 6 typhoons
Concurrent API calls for speed
"""
import numpy as np
import json, os, csv, sys, time
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request

TYPHOONS = {
    "meranti": {
        "landfall": datetime(2016, 9, 14, 0, 0),
        "station": ("IU", "TATO", "00", "BHZ"),
    },
    "hato": {
        "landfall": datetime(2017, 8, 23, 0, 0),
        "station": ("HK", "HKPS", "00", "BHZ"),
    },
    "mangkhut": {
        "landfall": datetime(2018, 9, 16, 0, 0),
        "station": ("HK", "HKPS", "00", "BHZ"),
    },
    "hagibis": {
        "landfall": datetime(2019, 10, 12, 0, 0),
        "station": ("G", "INU", "00", "BHZ"),
    },
    "goni": {
        "landfall": datetime(2020, 11, 1, 0, 0),
        "station": ("IU", "DAV", "00", "BHZ"),
    },
    "saola": {
        "landfall": datetime(2023, 9, 1, 0, 0),
        "station": ("HK", "HKPS", "00", "BHZ"),
    }
}

DATA_DIR = "projects/typhoon-backtest/data/6typhoon_results"
OUTPUT_DIR = "projects/typhoon-backtest/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def fetch_1h(net, sta, loc, cha, dt):
    """Fetch 1h continuous BHZ as ASCII, return RMS"""
    start = dt.strftime("%Y-%m-%dT%H:%M:%S")
    end = (dt + timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%S")
    url = f"http://service.iris.edu/irisws/timeseries/1/query?net={net}&sta={sta}&loc={loc}&cha={cha}&starttime={start}&endtime={end}&output=ascii"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-v5/1.0'})
        with urllib.request.urlopen(req, timeout=25) as r:
            data = r.read().decode()
        lines = data.strip().split("\n")
        vals = []
        for line in lines[1:]:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    try: vals.append(float(parts[-1]))
                    except: pass
        arr = np.array(vals)
        if len(arr) > 100:
            det = arr - np.mean(arr)
            return float(np.sqrt(np.mean(det**2)))
    except:
        pass
    return None

def process_typhoon(ty_name, ty_info):
    """Process one typhoon: fetch all time points and analyze"""
    net, sta, loc, cha = ty_info["station"]
    lf = ty_info["landfall"]
    
    # Load typhoon dH_curl CSV
    csv_file = None
    for fn in sorted(os.listdir(DATA_DIR)):
        if ty_name in fn.lower() and fn.endswith("_dh_curl.csv"):
            csv_file = os.path.join(DATA_DIR, fn)
            break
    if not csv_file:
        return {"error": f"No dH_curl data for {ty_name}"}
    
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        dh_hours = np.array([float(r["hours"]) for r in reader])
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        dh_roll = np.array([float(r["dH_roll"]) for r in reader])
    
    # Filter ±3 days, subsample 6-hourly
    mask = np.abs(dh_hours / 24) <= 3
    fh = dh_hours[mask][::6]
    fr = dh_roll[mask][::6]
    
    # Fetch in parallel
    tasks = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for h in fh:
            dt = lf + timedelta(hours=h)
            tasks.append(ex.submit(fetch_1h, net, sta, loc, cha, dt))
        
        rmss = []
        for i, ft in enumerate(as_completed(tasks)):
            rms = ft.result()
            rmss.append(rms)
            sys.stdout.write(f"\r  {ty_name}: {i+1}/{len(fh)}  ({rms if rms else '--'})    ")
            sys.stdout.flush()
    print()
    
    # Align results
    rmss_arr = np.array([r if r is not None else np.nan for r in rmss])
    valid = ~np.isnan(rmss_arr)
    v_rms = rmss_arr[valid]
    v_hours = fh[valid]
    v_roll = fr[valid]
    
    result = {
        "typhoon": f"{ty_name}_{lf.year}",
        "station": f"{net}.{sta}.{loc}.{cha}",
        "landfall": str(lf),
        "n_total": len(fh),
        "n_valid": int(valid.sum()),
        "hours": fh.tolist(),
        "rms": rmss_arr.tolist(),
    }
    
    if valid.sum() < 5:
        return result
    
    # Pre/post landfall
    lf_idx = np.argmin(np.abs(v_hours))
    pre = v_rms[:lf_idx]
    post = v_rms[lf_idx:]
    pre_m = float(pre.mean()) if len(pre) > 0 else float('nan')
    post_m = float(post.mean()) if len(post) > 0 else float('nan')
    result["pre_mean_rms"] = pre_m
    result["post_mean_rms"] = post_m
    result["ratio_pre_post"] = float(post_m / pre_m) if pre_m > 0 else float('nan')
    result["peak_rms"] = float(v_rms.max())
    result["peak_hour"] = float(v_hours[np.argmax(v_rms)])
    result["baseline_rms"] = float(v_rms[0])
    result["peak_baseline_ratio"] = float(v_rms.max() / v_rms[0]) if v_rms[0] > 0 else float('nan')
    
    # Correlations
    result["r_hours_vs_rms"] = float(np.corrcoef(v_hours, v_rms)[0, 1])
    result["r_dhroll_vs_rms"] = float(np.corrcoef(v_roll, v_rms)[0, 1])
    
    # dH_curl of noise
    def _dh_curl(seq, w=2):
        s = np.array(seq)
        c = np.full_like(s, np.nan)
        for i in range(w, len(s) - w):
            local = s[i-w:i+w+1]
            c[i] = float(np.std(np.diff(local)) / (np.mean(np.abs(local)) + 1e-8))
        return c
    
    ndh = _dh_curl(v_rms)
    vdh = ndh[~np.isnan(ndh)]
    vdr = v_roll[~np.isnan(ndh)]
    
    if len(vdh) > 3:
        result["r_dhcurl_noise_vs_dhroll"] = float(np.corrcoef(vdr, vdh)[0, 1])
        
        # Lagged cross-corr
        max_lag = 5
        lags = np.arange(-max_lag, max_lag + 1)
        cors = []
        for lag in lags:
            if lag < 0: x, y = vdr[-lag:], vdh[:lag]
            elif lag > 0: x, y = vdr[:-lag], vdh[lag:]
            else: x, y = vdr, vdh
            cors.append(float(np.corrcoef(x, y)[0, 1]) if len(x) > 3 else float('nan'))
        cors = np.array(cors)
        pk = int(np.nanargmax(np.abs(cors)))
        result["lagged_corr"] = {
            "lags": [int(l) for l in lags],
            "r": cors.tolist(),
            "peak_lag": int(lags[pk]),
            "peak_r": float(cors[pk]),
        }
    
    return result

# === MAIN ===
print(f"{'='*60}")
print(f"v5 Multi-Typhoon Seismic Noise (Parallel)")
print(f"{'='*60}")
t0 = time.time()

all_results = {}
for ty_name, ty_info in TYPHOONS.items():
    print(f"\n{'─'*50}")
    print(f"🌀 {ty_name.upper()} @ {'.'.join(ty_info['station'])}")
    print(f"{'─'*50}")
    result = process_typhoon(ty_name, ty_info)
    all_results[ty_name] = result
    
    if "error" not in result:
        r = result
        ratio = r.get("ratio_pre_post", float('nan'))
        peak = r.get("peak_baseline_ratio", float('nan'))
        r_simple = r.get("r_dhroll_vs_rms", float('nan'))
        lag = r.get("lagged_corr", {})
        lr = lag.get("peak_r", float('nan')) if lag else float('nan')
        ll = lag.get("peak_lag", float('nan')) if lag else float('nan')
        print(f"   Valid: {r['n_valid']}/{r['n_total']}")
        print(f"   Pre→Post: {r['pre_mean_rms']:.0f}→{r['post_mean_rms']:.0f} ({ratio:.2f}x)" if not np.isnan(ratio) else "   Pre→Post: N/A")
        print(f"   Peak/baseline: {peak:.1f}x @ t={r['peak_hour']:+.0f}h" if not np.isnan(peak) else "   Peak: N/A")
        print(f"   r(dHroll,RMS): {r_simple:+.4f}" if not np.isnan(r_simple) else "   r(dHroll,RMS): N/A")
        print(f"   Peak lag: lag={ll:+d}, r={lr:+.4f}" if not np.isnan(lr) else "   Peak lag: N/A")

elapsed = time.time() - t0
print(f"\n{'='*60}")
print(f"✅ Done in {elapsed:.0f}s")
print(f"{'='*60}")

# Save
with open(os.path.join(OUTPUT_DIR, "seismic_noise_v5_multi.json"), "w") as f:
    json.dump(all_results, f, indent=2, default=str)

# Summary table
print(f"\n{'='*80}")
print(f"{'Typhoon':<12} {'Station':<16} {'Valid':<8} {'Ratio':<10} {'Peak/bl':<10} {'r(RMS)':<10} {'Lag r':<10}")
print(f"{'-'*72}")
for ty in TYPHOONS:
    r = all_results.get(ty, {})
    if "error" in r:
        print(f"{ty:<12} ERROR")
        continue
    v = f"{r.get('n_valid','?'):>2}/{r.get('n_total','?'):>2}"
    rt = f"{r.get('ratio_pre_post',float('nan')):.2f}x" if not np.isnan(r.get('ratio_pre_post',float('nan'))) else "N/A"
    pk = f"{r.get('peak_baseline_ratio',float('nan')):.1f}x" if not np.isnan(r.get('peak_baseline_ratio',float('nan'))) else "N/A"
    rs = f"{r.get('r_dhroll_vs_rms',float('nan')):+.3f}" if not np.isnan(r.get('r_dhroll_vs_rms',float('nan'))) else "N/A"
    lg = r.get('lagged_corr',{})
    lr = f"{lg.get('peak_r',float('nan')):+.3f}" if lg and not np.isnan(lg.get('peak_r',float('nan'))) else "N/A"
    print(f"{ty:<12} {r.get('station','?'):<16} {v:<8} {rt:<10} {pk:<10} {rs:<10} {lr:<10}")
print(f"{'='*80}")
