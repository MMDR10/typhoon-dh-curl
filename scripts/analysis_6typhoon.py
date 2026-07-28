import csv
import json
import numpy as np
from collections import defaultdict
import os

# Load all CSVs
media_dir = "media"
typhoons = {}

for fname in os.listdir(media_dir):
    if not fname.endswith("_dh_curl.csv"):
        continue
    name = fname.replace(".csv", "").split("_", 2)[-1]  # e.g. hato_2017_dh_curl
    path = os.path.join(media_dir, fname)
    hours, dh_curl, dh_roll = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            hours.append(float(row["hours"]))
            dh_curl.append(float(row["dH_curl"]))
            dh_roll.append(float(row["dH_roll"]))
    typhoons[name] = {"hours": np.array(hours), "dH_curl": np.array(dh_curl), "dH_roll": np.array(dh_roll)}

# Load summary
with open(os.path.join(media_dir, "a945787a310c4c53b3e8679972f8a433_6typhoon_summary.json")) as f:
    summary = json.load(f)

# ===== ANALYSIS 1: Statistical distributions =====
print("=" * 70)
print("ANALYSIS 1: dH_curl Statistical Distributions")
print("=" * 70)
print(f"{'Typhoon':<20} {'mean':>8} {'std':>8} {'skew':>8} {'kurt':>8} {'min':>8} {'max':>8} {'var_ratio':>10}")
print("-" * 78)

stats = {}
for key, data in typhoons.items():
    c = data["dH_curl"]
    r = data["dH_roll"]
    mean_c = np.mean(c)
    std_c = np.std(c)
    skew_c = ((c - mean_c)**3).mean() / (std_c**3) if std_c > 0 else 0
    kurt_c = ((c - mean_c)**4).mean() / (std_c**4) - 3  # excess kurtosis
    var_ratio = np.var(r) / np.var(c) if np.var(c) > 0 else 0
    
    short = key.split("_")[0]
    print(f"{short:<20} {mean_c:>8.4f} {std_c:>8.4f} {skew_c:>8.3f} {kurt_c:>8.3f} {np.min(c):>8.3f} {np.max(c):>8.3f} {var_ratio:>10.4f}")
    stats[key] = {"mean": mean_c, "std": std_c, "skew": skew_c, "kurt": kurt_c, "var_ratio": var_ratio}

print()
print("Interpretation:")
print("  - Negative skew = more extreme negative dH_curl spikes (landfall impact)")
print("  - High kurtosis = fat tails, extreme events dominate")
print("  - var_ratio > 1 = dH_roll amplifies signal; < 1 = dH_roll smooths")

# ===== ANALYSIS 2: Phase structure =====
print()
print("=" * 70)
print("ANALYSIS 2: Phase Decomposition (dH_roll)")
print("=" * 70)

phase_windows = {
    "pre_-120_-96": (-120, -96),
    "pre_-96_-72":  (-96, -72),
    "pre_-72_-48":  (-72, -48),
    "pre_-48_-24":  (-48, -24),
    "pre_-24_0":    (-24, 0),
    "post_0_24":    (0, 24),
    "post_24_72":   (24, 72),
    "post_72_120":  (72, 120),
    "post_120_168": (120, 168),
}

header = f"{'Typhoon':<20}"
for pw in phase_windows:
    header += f" {pw:>12}"
print(header)
print("-" * len(header))

phase_data = {}
for key, data in typhoons.items():
    h = data["hours"]
    r = data["dH_roll"]
    short = key.split("_")[0]
    line = f"{short:<20}"
    pd = {}
    for pw, (lo, hi) in phase_windows.items():
        mask = (h >= lo) & (h < hi)
        if mask.sum() > 0:
            val = r[mask].mean()
        else:
            val = np.nan
        pd[pw] = val
        line += f" {val:>12.4f}"
    print(line)
    phase_data[key] = pd

# ===== ANALYSIS 3: Autocorrelation / Persistence =====
print()
print("=" * 70)
print("ANALYSIS 3: Persistence (dH_roll autocorrelation)")
print("=" * 70)

for key, data in typhoons.items():
    r = data["dH_roll"]
    short = key.split("_")[0]
    n = len(r)
    if n > 24:
        ac1 = np.corrcoef(r[:-1], r[1:])[0, 1]
        ac24 = np.corrcoef(r[:-24], r[24:])[0, 1]
        decay = ac24 / ac1 if ac1 > 0 else 0
        print(f"  {short:<20} ac(1h)={ac1:>7.4f}  ac(24h)={ac24:>7.4f}  decay/day={decay:>7.4f}")

# ===== ANALYSIS 4: Extreme event detection =====
print()
print("=" * 70)
print("ANALYSIS 4: Extreme dH_curl Events (beyond +-2sigma)")
print("=" * 70)

for key, data in typhoons.items():
    c = data["dH_curl"]
    h = data["hours"]
    short = key.split("_")[0]
    mean_c = stats[key]["mean"]
    std_c = stats[key]["std"]
    neg_extreme = c < (mean_c - 2*std_c)
    pos_extreme = c > (mean_c + 2*std_c)
    n_neg = neg_extreme.sum()
    n_pos = pos_extreme.sum()
    
    neg_hours = h[neg_extreme]
    pos_hours = h[pos_extreme]
    
    print(f"  {short:<20} neg_extreme: {n_neg:>4d} ({n_neg/len(c)*100:>5.1f}%)  pos_extreme: {n_pos:>4d} ({n_pos/len(c)*100:>5.1f}%)")
    if n_neg > 0:
        print(f"    neg cluster: {min(neg_hours):+.0f}h to {max(neg_hours):+.0f}h  |  median: {np.median(neg_hours):+.0f}h")

# ===== ANALYSIS 5: Cross-typhoon correlation =====
print()
print("=" * 70)
print("ANALYSIS 5: dH_roll Cross-Correlation Matrix (aligned @ landfall)")
print("=" * 70)

names = sorted(typhoons.keys())
n = len(names)
corr_matrix = np.zeros((n, n))
for i, ki in enumerate(names):
    for j, kj in enumerate(names):
        if i == j:
            corr_matrix[i, j] = 1.0
        else:
            ri = typhoons[ki]["dH_roll"]
            rj = typhoons[kj]["dH_roll"]
            min_len = min(len(ri), len(rj))
            corr_matrix[i, j] = np.corrcoef(ri[:min_len], rj[:min_len])[0, 1]

short_names = [k.split("_")[0] for k in names]
print(f"{'':>12}", end="")
for sn in short_names:
    print(f" {sn:>8}", end="")
print()
for i, sn in enumerate(short_names):
    print(f"{sn:>12}", end="")
    for j in range(n):
        print(f" {corr_matrix[i,j]:>8.4f}", end="")
    print()

# ===== ANALYSIS 6: dH_curl power spectrum =====
print()
print("=" * 70)
print("ANALYSIS 6: dH_curl Dominant Frequencies (FFT)")
print("=" * 70)

for key, data in typhoons.items():
    c = data["dH_curl"]
    short = key.split("_")[0]
    n = len(c)
    fft = np.abs(np.fft.rfft(c - np.mean(c)))**2
    freqs = np.fft.rfftfreq(n, d=1.0)
    
    peak_idx = np.argsort(fft[1:])[-3:] + 1
    peak_power = fft[peak_idx]
    peak_periods = 1.0 / freqs[peak_idx]
    
    sort_idx = np.argsort(peak_periods)
    print(f"  {short:<20} top periods: ", end="")
    for idx in sort_idx:
        print(f"{peak_periods[idx]:.1f}h ({peak_power[idx]:.0f})  ", end="")
    print()

# ===== ANALYSIS 7: Cumulative dH_roll energy =====
print()
print("=" * 70)
print("ANALYSIS 7: Cumulative dH Energy Budget")
print("=" * 70)

for key, data in typhoons.items():
    r = data["dH_roll"]
    h = data["hours"]
    short = key.split("_")[0]
    
    pre_mask = h < 0
    post_mask = h >= 0
    
    pre_energy = abs(r[pre_mask].sum()) if pre_mask.sum() > 0 else 0
    post_energy = abs(r[post_mask].sum()) if post_mask.sum() > 0 else 0
    total_energy = abs(r.sum())
    
    pre_pct = pre_energy / total_energy * 100 if total_energy > 0 else 0
    
    print(f"  {short:<20} |pre_E|={pre_energy:>8.2f}  |post_E|={post_energy:>8.2f}  |total|={total_energy:>8.2f}  pre_frac={pre_pct:>5.1f}%")

# ===== ANALYSIS 8: Recovery time =====
print()
print("=" * 70)
print("ANALYSIS 8: Recovery -- Time to dH_roll return to baseline")
print("=" * 70)

for key, data in typhoons.items():
    r = data["dH_roll"]
    h = data["hours"]
    short = key.split("_")[0]
    
    min_idx = np.argmin(r)
    min_h = h[min_idx]
    min_val = r[min_idx]
    
    post_min_mask = (h > min_h) & (r > -0.1)
    if post_min_mask.sum() > 0:
        recovery_h = h[post_min_mask][0]
        recovery_time = recovery_h - min_h
    else:
        recovery_h = h[-1]
        recovery_time = recovery_h - min_h
    
    print(f"  {short:<20} dH_roll_min={min_val:>7.3f} @ {min_h:>+6.0f}h")
    print(f"    recovery to >-0.1 @ {recovery_h:>+6.0f}h  ({recovery_time:>5.0f}h = {recovery_time/24:.1f}d)")

# ===== ANALYSIS 9: Typhoon "fingerprint" classification =====
print()
print("=" * 70)
print("ANALYSIS 9: Typhoon dH Fingerprint Classification")
print("=" * 70)

# Classify by shape
classifications = {}
for key, data in typhoons.items():
    r = data["dH_roll"]
    h = data["hours"]
    short = key.split("_")[0]
    s = summary.get(short, {})
    
    mode = s.get("mode", "unknown")
    u_str = s.get("u_strength", 0)
    dH_range = s.get("dH_range", 0)
    
    # Deepen rate: how fast does dH_roll drop pre-landfall?
    pre_72_0 = h[(h >= -72) & (h < 0)]
    pre_72_0_r = r[(h >= -72) & (h < 0)]
    if len(pre_72_0) > 0:
        deepen_rate = (pre_72_0_r[-1] - pre_72_0_r[0]) / 72  # per hour
    else:
        deepen_rate = 0
    
    # Post-landfall deepen
    post_0_72 = h[(h >= 0) & (h < 72)]
    post_0_72_r = r[(h >= 0) & (h < 72)]
    if len(post_0_72) > 0:
        post_rate = (post_0_72_r[-1] - post_0_72_r[0]) / 72
    else:
        post_rate = 0
    
    classifications[key] = {
        "mode": mode,
        "u_str": u_str,
        "dH_range": dH_range,
        "deepen_rate": deepen_rate,
        "post_rate": post_rate
    }
    
    print(f"  {short:<20} mode={mode:<12} pre72_deepen={deepen_rate:>+8.4f}/h  post72_rate={post_rate:>+8.4f}/h")

# ===== ANALYSIS 10: Detect "U-shape completeness" =====
print()
print("=" * 70)
print("ANALYSIS 10: U-shape Completeness Metric")
print("=" * 70)
print("Completeness = 1 means perfect symmetric U; < 0.5 means asymmetric or flat")
print()

for key, data in typhoons.items():
    r = data["dH_roll"]
    h = data["hours"]
    short = key.split("_")[0]
    
    # Pre-landfall deepening (overall trend -120 to 0)
    pre_mask = h < 0
    post_mask = h >= 0
    
    pre_min = r[pre_mask].min() if pre_mask.sum() > 0 else 0
    post_min = r[post_mask].min() if post_mask.sum() > 0 else 0
    
    # Pre minimum depth vs post minimum depth
    if pre_min < 0 and post_min < 0:
        symmetry = min(abs(pre_min), abs(post_min)) / max(abs(pre_min), abs(post_min))
    else:
        symmetry = 0
    
    # Bottom persistence: how long dH stays near minimum
    overall_min = r.min()
    near_min_mask = r < (overall_min * 0.8)  # within 80% of min
    bottom_duration = near_min_mask.sum()  # hours
    
    print(f"  {short:<20} pre_min={pre_min:>7.3f}  post_min={post_min:>7.3f}  symmetry={symmetry:>5.3f}  bottom_dur={bottom_duration:>4d}h ({bottom_duration/24:.1f}d)")

print()
print("=" * 70)
print("ALL ANALYSES COMPLETE")
print("=" * 70)
