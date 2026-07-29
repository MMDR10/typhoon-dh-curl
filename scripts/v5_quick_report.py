#!/usr/bin/env python3
"""Quick V5 Multi-Typhoon Report — no f-string nesting issues"""
import json, numpy as np

with open('projects/typhoon-backtest/output/seismic_noise_v5_multi.json') as f:
    data = json.load(f)

typhoons = ["meranti", "hato", "mangkhut", "hagibis", "goni", "saola"]
stations = {
    "meranti": "IU.TATO (Taipei)",  "hato": "HK.HKPS (HK)",
    "mangkhut": "HK.HKPS (HK)",    "hagibis": "G.INU (Japan)",
    "goni": "IU.DAV (Philippines)", "saola": "HK.HKPS (HK)",
}

print("=" * 80)
print("V5 MULTI-TYPHOON ANALYSIS")
print("6 Typhoons x IRIS continuous seismic noise")
print("=" * 80)

# Extract data
rows = {}
for ty in typhoons:
    r = data.get(ty, {})
    pre = r.get("pre_landfall_mean_rms", float('nan'))
    post = r.get("post_landfall_mean_rms", float('nan'))
    ratio = r.get("pre_post_ratio", float('nan'))
    peak = r.get("peak_rms", float('nan'))
    base = r.get("baseline_rms", float('nan'))
    pk_hr = r.get("peak_hours", float('nan'))
    pk_base = peak / base if base and base > 0 else float('nan')
    rs = r.get("r_typhoon_dHroll_vs_RMS", float('nan'))
    lag = r.get("lagged_cross_corr", {})
    if isinstance(lag, dict):
        lr = lag.get("peak_r", float('nan'))
        ll = lag.get("peak_lag", float('nan'))
    else:
        lr = float('nan'); ll = float('nan')
    rows[ty] = (pre, post, ratio, peak, base, pk_hr, pk_base, rs, lr, ll)

# 1. Table
print("\n--- 1. PER-TYPHOON TABLE ---")
hdr = f"{'Typhoon':<12} {'PreRMS':>8} {'PostRMS':>9} {'Ratio':>7} {'PeakRMS':>9} {'Pk/Base':>8} {'Pk@hr':>6} {'r(dH,RM)':>8} {'Lag r':>6} {'Lag':>4}"
print(hdr)
print("-" * 85)

for ty in typhoons:
    pre, post, ratio, peak, base, pk_hr, pk_base, rs, lr, ll = rows[ty]
    def fmt(v, d=1):
        if np.isnan(v): return f"{'N/A':>7}"
        return f"{v:>7.{d}f}"
    def fmt2(v, d=1):
        if np.isnan(v): return f"{'N/A':>5}"
        return f"{v:>5.{d}f}"
    line = f"{ty:<12} {fmt(pre,0):>8} {fmt(post,0):>9} {fmt(ratio,2):>7} {fmt(peak,0):>9} {fmt(pk_base,1):>8} {fmt2(pk_hr,0):>6} {fmt(rs,3):>8} {fmt2(lr,3):>6} {fmt2(ll,0):>4}"
    print(line)

# 2. Stats
print("\n--- 2. STATISTICS ---")
ratios = np.array([rows[t][2] for t in typhoons])
pbs = np.array([rows[t][6] for t in typhoons])
pks = np.array([rows[t][5] for t in typhoons])
lrs = np.array([rows[t][8] for t in typhoons])
lls = np.array([rows[t][9] for t in typhoons])

print(f"  Pre/Post ratio: mean={ratios.mean():.2f}x  range=[{ratios.min():.2f}, {ratios.max():.2f}]")
print(f"  All > 1? {((ratios > 1).all())}")
print(f"  Peak/Baseline: mean={pbs.mean():.1f}x  max={pbs.max():.1f}x ({typhoons[np.argmax(pbs)]})")
print(f"  Peak timing: mean={np.mean(pks):+.0f}h from landfall")
print(f"  Lagged r(dH,noise): mean={np.mean(lrs):+.3f} @ mean lag={np.mean(lls):+.0f}h")
print(f"  Dominant sign: {'NEGATIVE' if np.mean(lrs) < 0 else 'POSITIVE'}")

# 3. Ranking
print("\n--- 3. RANKED ---")
ranked = sorted(typhoons, key=lambda t: rows[t][2], reverse=True)
for i, ty in enumerate(ranked):
    print(f"  {i+1}. {ty:<12} ratio={rows[ty][2]:.2f}x  peak/base={rows[ty][6]:.1f}x  station={stations[ty]}")

# 4. Station groups
hkps = ["hato", "mangkhut", "saola"]
other = ["meranti", "hagibis", "goni"]
for label, grp in [("HK.HKPS", hkps), ("Non-HKPS", other)]:
    gr = np.mean([rows[t][2] for t in grp])
    print(f"\n  {label}: mean ratio={gr:.2f}x")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
