"""
Fetch IRIS continuous seismic waveform data for one station + one typhoon
Compute hourly RMS noise amplitude and dH_curl for cross-correlation
"""
import numpy as np
import subprocess, json, os, sys, math
from datetime import datetime, timedelta
import urllib.request

# ========== CONFIG ==========
# Meranti 2016 + IU.TATO (Taipei, 24.97N, 121.50E)
# TATO is ~250km from Meranti's landfall path
LANDFALL = datetime(2016, 9, 14, 0)
STATION = "TATO"
NET = "IU"
LOC = "00"
CHA = "BHZ"
WINDOW = 7  # ±7 days (14 days total, manageable)
OUTPUT = "projects/typhoon-backtest/output/tato_noise.json"

# TATO coordinates
TATO_LAT, TATO_LON = 24.9735, 121.4971

def fetch_24h(date):
    """Fetch 24 hours of BHZ data as ASCII"""
    start = date.strftime("%Y-%m-%dT00:00:00")
    end = (date + timedelta(days=1)).strftime("%Y-%m-%dT00:00:00")
    url = (
        f"https://service.iris.edu/irisws/timeseries/1/query"
        f"?net={NET}&sta={STATION}&loc={LOC}&cha={CHA}"
        f"&starttime={start}&endtime={end}&output=ascii"
    )
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'DR-coupling-v5/1.0'})
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read().decode()
        # Parse: skip header, extract values
        lines = data.strip().split("\n")
        if len(lines) < 2:
            return None
        # Parse header
        header = lines[0]
        parts = header.split(",")
        sps = 20  # known from earlier
        # Parse data lines: "YYYY-MM-DDTHH:MM:SS.ffffff  VALUE"
        values = []
        for line in lines[1:]:
            line = line.strip()
            if line:
                try:
                    val = float(line.split()[-1])
                    values.append(val)
                except (ValueError, IndexError):
                    continue
        return np.array(values), sps, header
    except Exception as e:
        print(f"  Fetch error for {date}: {e}", file=sys.stderr)
        return None

def compute_noise_metrics(amplitudes, sps):
    """Compute noise metrics from continuous amplitude data"""
    n = len(amplitudes)
    if n < 100:
        return None
    
    # 1. Basic stats
    mean = np.mean(amplitudes)
    std = np.std(amplitudes)
    
    # 2. Hourly RMS (detrended to remove DC offset)
    hourly_samples = sps * 3600  # samples per hour
    n_hours = n // hourly_samples
    hourly_rms = []
    for h in range(n_hours):
        seg = amplitudes[h*hourly_samples:(h+1)*hourly_samples]
        seg_detrended = seg - np.mean(seg)
        rms = np.sqrt(np.mean(seg_detrended**2))
        hourly_rms.append(rms)
    
    # 3. Frequency bands: power in different bands
    # Compute segments of ~10 min for dH_curl
    ten_min = sps * 600
    n_seg = n // ten_min
    seg_rms = np.array([np.sqrt(np.mean((amplitudes[i*ten_min:(i+1)*ten_min] - 
                                         np.mean(amplitudes[i*ten_min:(i+1)*ten_min]))**2)) 
                        for i in range(n_seg) if i*ten_min < n])
    
    # 4. dH_curl of RMS sequence (noise of the noise)
    def dh_curl(seq, w=2):
        seq = np.array(seq, dtype=float)
        curl = np.full_like(seq, np.nan)
        for i in range(w, len(seq)-w):
            local = seq[i-w:i+w+1]
            dH = np.std(np.diff(local)) / (np.mean(np.abs(local)) + 1e-8)
            curl[i] = dH
        return curl
    
    seg_dh = dh_curl(seg_rms, w=2)
    
    return {
        "n_samples": int(n),
        "mean": float(mean),
        "std": float(std),
        "cv": float(std/abs(mean)) if mean != 0 else None,
        "hourly_rms": [float(r) for r in hourly_rms],
        "n_hours": len(hourly_rms),
        "seg_rms": [float(r) for r in seg_rms],
        "n_seg": len(seg_rms),
        "seg_dh_curl": [float(x) for x in seg_dh if not np.isnan(x)],
        "mean_seg_rms": float(np.mean(seg_rms)),
        "std_seg_rms": float(np.std(seg_rms)),
        "cv_seg_rms": float(np.std(seg_rms) / np.mean(seg_rms)) if np.mean(seg_rms) > 0 else None,
    }

# ========== MAIN ==========
print("="*60)
print(f"Seismic Noise Fetch: {NET}.{STATION}.{LOC}.{CHA}")
print(f"Typhoon: Meranti landfall {LANDFALL.date()}")
print(f"Window: ±{WINDOW} days")
print("="*60)

start_date = LANDFALL - timedelta(days=WINDOW)
end_date = LANDFALL + timedelta(days=WINDOW)

all_metrics = {}
dates = []
current = start_date
while current <= end_date:
    dates.append(current)
    current += timedelta(days=1)

print(f"Fetching {len(dates)} days of data...")

for date in dates:
    date_str = date.strftime("%Y-%m-%d")
    print(f"  [{date_str}] fetching 24h...", end=" ", flush=True)
    result = fetch_24h(date)
    
    if result is None:
        print(f"❌")
        all_metrics[date_str] = None
        continue
    
    amplitudes, sps, header = result
    metrics = compute_noise_metrics(amplitudes, sps)
    
    if metrics is None:
        print(f"⚠️ (too few samples: {len(amplitudes)})")
        all_metrics[date_str] = {"n_samples": len(amplitudes), "error": "too_few"}
        continue
    
    print(f"✅ {metrics['n_hours']}h  mean_RMS={metrics['mean_seg_rms']:.1f}  CV={metrics['cv_seg_rms']:.3f}")
    all_metrics[date_str] = metrics

# ========== SUMMARY ==========
print(f"\n{'='*60}")
print("Noise Metrics Summary")
print(f"{'='*60}")

valid_days = {k: v for k, v in all_metrics.items() if v and "n_hours" in v}
print(f"Valid days: {len(valid_days)}/{len(dates)}")

if valid_days:
    # Aggregate statistics
    daily_mean_rms = [v["mean_seg_rms"] for v in valid_days.values()]
    daily_cv_rms = [v["cv_seg_rms"] for v in valid_days.values()]
    
    print(f"\nDaily segment RMS noise:")
    print(f"  Range: {min(daily_mean_rms):.1f} - {max(daily_mean_rms):.1f}")
    print(f"  Mean: {np.mean(daily_mean_rms):.1f} ± {np.std(daily_mean_rms):.1f}")
    print(f"  CV across days: {np.std(daily_mean_rms)/np.mean(daily_mean_rms):.3f}")
    
    print(f"\nNoise variability (CV of segment RMS):")
    print(f"  Range: {min(daily_cv_rms):.3f} - {max(daily_cv_rms):.3f}")
    
    # Check if noise changes around landfall
    days_sorted = sorted(valid_days.keys())
    landfall_str = LANDFALL.strftime("%Y-%m-%d")
    pre = [valid_days[d] for d in days_sorted if d < landfall_str]
    post = [valid_days[d] for d in days_sorted if d >= landfall_str]
    
    if pre and post:
        pre_mean = np.mean([v["mean_seg_rms"] for v in pre])
        post_mean = np.mean([v["mean_seg_rms"] for v in post])
        print(f"\nLandfall-relative noise change:")
        print(f"  Pre-landfall mean RMS:  {pre_mean:.1f} (n={len(pre)} days)")
        print(f"  Post-landfall mean RMS: {post_mean:.1f} (n={len(post)} days)")
        print(f"  Ratio: {post_mean/pre_mean:.3f}")

# Save
os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
with open(OUTPUT, "w") as f:
    json.dump(all_metrics, f, indent=2, default=str)
print(f"\n✅ Saved to {OUTPUT}")
