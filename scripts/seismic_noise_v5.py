"""
v5: Smart seismic noise fetch - 1-hour windows at typhoon dH_curl time points
Match IRIS continuous waveform to typhoon 6-hourly dH_curl timestamps
Batched requests for efficiency
"""
import numpy as np
import json, os, csv
from datetime import datetime, timedelta
import urllib.request

# ========== CONFIG ==========
LANDFALL = datetime(2016, 9, 14, 0)  # Meranti
STATION = "TATO"
NET = "IU"
LOC = "00"
CHA = "BHZ"
WINDOW_DAYS = 3  # ±3 days around landfall
SUBSAMPLE = 6     # take every Nth point to reduce API calls

# ========== HELPERS ==========
def fetch_1hour(dt):
    """Fetch 1 hour of continuous BHZ data at given datetime"""
    start = dt.strftime("%Y-%m-%dT%H:%M:%S")
    end = (dt + timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%S")
    url = (
        f"https://service.iris.edu/irisws/timeseries/1/query"
        f"?net={NET}&sta={STATION}&loc={LOC}&cha={CHA}"
        f"&starttime={start}&endtime={end}&output=ascii"
    )
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v5/1.0'})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = resp.read().decode()
        lines = data.strip().split("\n")
        values = []
        for line in lines[1:]:
            line = line.strip()
            if line:
                try:
                    values.append(float(line.split()[-1]))
                except (ValueError, IndexError):
                    continue
        return np.array(values) if values else None
    except Exception as e:
        print(f"    Error: {e}")
        return None

def compute_rms(values):
    """Compute RMS noise amplitude (detrended)"""
    if values is None or len(values) < 100:
        return None
    detrended = values - np.mean(values)
    return float(np.sqrt(np.mean(detrended**2)))

def dh_curl(seq, window=2):
    """Compute dH_curl of a sequence"""
    seq = np.array(seq, dtype=float)
    curl = np.full_like(seq, np.nan)
    for i in range(window, len(seq)-window):
        local = seq[i-window:i+window+1]
        curl[i] = np.std(np.diff(local)) / (np.mean(np.abs(local)) + 1e-8)
    return curl

# ========== GET TYPHOON dH_curl TIME POINTS ==========
datadir = "projects/typhoon-backtest/data/6typhoon_results"
dh_hours = None
for fname in sorted(os.listdir(datadir)):
    if "meranti" in fname.lower() and fname.endswith("_dh_curl.csv"):
        with open(os.path.join(datadir, fname)) as f:
            reader = csv.DictReader(f)
            dh_hours = np.array([float(r["hours"]) for r in reader])
        break

if dh_hours is None:
    print("ERROR: Cannot find Meranti dH_curl data")
    exit(1)

# Filter to ±WINDOW_DAYS
mask = np.abs(dh_hours / 24) <= WINDOW_DAYS
dh_hours = dh_hours[mask]
# Subsample to reduce API calls
dh_hours = dh_hours[::SUBSAMPLE]
n_points = len(dh_hours)
print(f"Typhoon dH_curl: {n_points} time points (every ~6h, ±{WINDOW_DAYS}d)")
print(f"  Range: {dh_hours[0]:.0f}h to {dh_hours[-1]:.0f}h from landfall")

# ========== FETCH SEISMIC NOISE AT EACH TYPHOON TIME POINT ==========
print(f"\nFetching {n_points} IRIS 1-hour windows...")
print(f"{'='*60}")

noise_rms = []
timestamps = []

for i, h in enumerate(dh_hours):
    dt = LANDFALL + timedelta(hours=h)
    ts = dt.strftime("%Y-%m-%d %H:%M")
    
    print(f"  [{i+1}/{n_points}] t={h:+.0f}h  {ts} ...", end=" ", flush=True)
    
    values = fetch_1hour(dt)
    rms = compute_rms(values)
    
    if rms is None:
        print(f"❌ (no data)")
        noise_rms.append(None)
    else:
        noise_rms.append(rms)
        print(f"✅ RMS={rms:.0f} counts")
    
    timestamps.append(ts)

# ========== ANALYSIS ==========
print(f"\n{'='*60}")
print("Noise Analysis")
print(f"{'='*60}")

valid = [(i, r) for i, r in enumerate(noise_rms) if r is not None]
n_valid = len(valid)
print(f"Valid noise measurements: {n_valid}/{n_points}")

if n_valid < 5:
    print("Too few valid points for meaningful analysis")
else:
    valid_rms = np.array([r for _, r in valid])
    valid_idx = np.array([i for i, _ in valid])
    valid_hours = dh_hours[valid_idx]
    
    print(f"\nNoise RMS range: {valid_rms.min():.0f} - {valid_rms.max():.0f} counts")
    print(f"Noise RMS mean:  {valid_rms.mean():.0f} ± {valid_rms.std():.0f} counts")
    print(f"Noise RMS CV:    {valid_rms.std()/valid_rms.mean():.3f}")
    
    # Load typhoon dH_roll at these points
    with open(os.path.join(datadir, fname)) as f:
        reader = csv.DictReader(f)
        all_roll = [float(r["dH_roll"]) for r in reader]
    dh_roll_full = np.array(all_roll)
    
    # Align with filtered time points
    dh_roll_aligned = dh_roll_full[mask]
    
    # Cross-correlation: typhoon dH_roll × seismic noise RMS
    valid_roll = dh_roll_aligned[valid_idx]
    r_ty_vs_noise = float(np.corrcoef(valid_roll, valid_rms)[0,1])
    print(f"\n🔗 Typhoon dH_roll × Seismic Noise RMS: r = {r_ty_vs_noise:+.4f}")
    
    # Check pre/post landfall difference
    lf_idx = np.argmin(np.abs(dh_hours))  # closest to landfall
    pre_mask = valid_idx < lf_idx
    post_mask = valid_idx >= lf_idx
    
    if pre_mask.sum() > 2 and post_mask.sum() > 2:
        pre_noise = valid_rms[pre_mask]
        post_noise = valid_rms[post_mask]
        pre_mean = pre_noise.mean()
        post_mean = post_noise.mean()
        print(f"\nLandfall noise comparison:")
        print(f"  Pre-landfall  ({pre_mask.sum()} pts): {pre_mean:.0f} ± {pre_noise.std():.0f}")
        print(f"  Post-landfall ({post_mask.sum()} pts): {post_mean:.0f} ± {post_noise.std():.0f}")
        print(f"  Ratio: {post_mean/pre_mean:.3f}")
    
    # dH_curl of noise RMS (noise of the noise)
    noise_dh = dh_curl(valid_rms)
    valid_dh = noise_dh[~np.isnan(noise_dh)]
    if len(valid_dh) > 3:
        valid_roll_dh = valid_roll[~np.isnan(noise_dh)]
        r_ty_dh_vs_noise_dh = float(np.corrcoef(valid_roll_dh, valid_dh)[0,1])
        print(f"\n🔗 Typhoon dH_roll × Noise dH_curl: r = {r_ty_dh_vs_noise_dh:+.4f}")
        
        # Directionality: does typhoon noise lead seismic noise?
        max_lag = 5
        lags = np.arange(-max_lag, max_lag+1)
        cors = []
        for lag in lags:
            if lag < 0:
                x, y = valid_roll_dh[-lag:], valid_dh[:lag]
            elif lag > 0:
                x, y = valid_roll_dh[:-lag], valid_dh[lag:]
            else:
                x, y = valid_roll_dh, valid_dh
            if len(x) > 3:
                cors.append(np.corrcoef(x, y)[0,1])
            else:
                cors.append(np.nan)
        cors = np.array(cors)
        peak_idx = np.nanargmax(np.abs(cors))
        print(f"  Lagged cross-correlation:")
        for l, c in zip(lags, cors):
            marker = " ← PEAK" if l == lags[peak_idx] else ""
            print(f"    lag={l:+d}: r={c:+.4f}{marker}")

# ========== SAVE ==========
results = {
    "typhoon": "Meranti_2016",
    "station": f"{NET}.{STATION}.{LOC}.{CHA}",
    "landfall": str(LANDFALL),
    "n_points": n_points,
    "n_valid": n_valid,
    "timestamps": timestamps,
    "noise_rms": noise_rms,
    "typhoon_dh_hours": [float(h) for h in dh_hours],
}

os.makedirs("projects/typhoon-backtest/output", exist_ok=True)
with open("projects/typhoon-backtest/output/seismic_noise_v5.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\n✅ Saved to output/seismic_noise_v5.json")
