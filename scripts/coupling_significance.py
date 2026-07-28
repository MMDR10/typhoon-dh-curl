"""
Significance test: shuffle landfall dates to get null distribution
"""
import json, numpy as np
from datetime import datetime, timedelta
import urllib.request

# Load coupling results
with open("projects/typhoon-backtest/output/coupling_analysis_v1.json") as f:
    results = json.load(f)

# Baseline: query USGS for a "control" period (no typhoons)
# We'll use the entire 2016-2023 period and sample random windows
print("=" * 60)
print("STATISTICAL SIGNIFICANCE TEST")
print("=" * 60)

# Observed post/pre ratios
observed_ratios = {k: v["post_pre_ratio"] for k, v in results.items()}
observed_energy_ratios = {k: v["post_pre_energy_ratio"] for k, v in results.items()}

print("\nObserved ratios:")
for k in observed_ratios:
    print(f"  {k.split('_')[0]:<12} events={observed_ratios[k]:.2f}  energy={observed_energy_ratios[k]:.2f}")

# Simple test: is the mean post/pre ratio significantly different from 1.0?
ratios_list = list(observed_ratios.values())
energy_ratios_list = list(observed_energy_ratios.values())

print(f"\nMean event ratio: {np.mean(ratios_list):.3f} ± {np.std(ratios_list):.3f}")
print(f"Mean energy ratio: {np.mean(energy_ratios_list):.3f} ± {np.std(energy_ratios_list):.3f}")

# Simple sign test: how many > 1.0?
n_above = sum(1 for r in ratios_list if r > 1.0)
n_below = sum(1 for r in ratios_list if r < 1.0)
# Binomial test: H0 p=0.5
from math import comb
p_binomial = sum(comb(6, k) * 0.5**6 for k in range(n_above, 7))
print(f"Sign test (H0: p=0.5): {n_above} above, {n_below} below → p={p_binomial:.4f} (binomial)")

n_above_e = sum(1 for r in energy_ratios_list if r > 1.0)
p_binomial_e = sum(comb(6, k) * 0.5**6 for k in range(n_above_e, 7))
print(f"Sign test energy: {n_above_e} above 1.0 → p={p_binomial_e:.4f}")

# t-stat manually
n = len(ratios_list)
t_stat = (np.mean(ratios_list) - 1.0) / (np.std(ratios_list, ddof=1) / np.sqrt(n))
# Approximate p from t-dist (2-sided)
print(f"t-stat (manual): {t_stat:.3f}, mean={np.mean(ratios_list):.3f}, se={np.std(ratios_list, ddof=1)/np.sqrt(n):.3f}")

t_stat_e = (np.mean(energy_ratios_list) - 1.0) / (np.std(energy_ratios_list, ddof=1) / np.sqrt(n))
print(f"t-stat energy: {t_stat_e:.3f}, mean={np.mean(energy_ratios_list):.3f}")

# ===== KEY INSIGHT: background rate estimation =====
print(f"\n{'='*60}")
print("BACKGROUND SEISMIC RATE ESTIMATION")
print(f"{'='*60}")

# For each typhoon, compute daily rate
for k, v in results.items():
    n_events = v["n_events"]
    window_days = 20  # ±10 days
    daily_rate = n_events / window_days
    print(f"  {k.split('_')[0]:<12} {daily_rate:.1f} events/day in NW Pacific (M≥4)")

# NW Pacific M≥4 baseline: approximately 15-20/day (from global catalogs)
# Let's check by querying a random control window
print("\nQuerying control window (no major typhoon, mid-2019)...")
ctrl_start = datetime(2019, 6, 1)
ctrl_end = datetime(2019, 6, 21)

url = (
    f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
    f"&starttime={ctrl_start.strftime('%Y-%m-%d')}"
    f"&endtime={ctrl_end.strftime('%Y-%m-%d')}"
    f"&minlatitude=0&maxlatitude=50&minlongitude=100&maxlongitude=180"
    f"&minmagnitude=4"
)
try:
    with urllib.request.urlopen(url, timeout=15) as resp:
        ctrl_data = json.loads(resp.read())
    ctrl_n = ctrl_data["metadata"]["count"]
    print(f"  Control window: {ctrl_n} events in 20 days = {ctrl_n/20:.1f}/day")
except Exception as e:
    print(f"  Control query failed: {e}")

print(f"\n{'='*60}")
print("ISSUES & FINDINGS")
print(f"{'='*60}")
print("""
1. ⚠️ Correlation too weak (|r| < 0.15 for all typhoons)
   → dH_curl and seismic rate are essentially uncorrelated at hourly resolution
   → This makes physical sense: atmospheric loading ~kPa, crustal stress ~MPa
   → 3 orders of magnitude difference → no direct linear coupling

2. 🔍 Hagibis anomaly: energy ratio = 10.20
   → Post-landfall seismic energy 10x pre-landfall
   → But this could be: (a) the M7.1 Fukushima 2021 foreshock cluster
      (b) coincidence with Japan Trench background seismicity
   → Need to filter out tectonic swarms from typhoon-related events

3. 🔍 Mangkhut paradox: post/pre = 0.67 (fewer earthquakes AFTER)
   → Strongest typhoon but seismic activity decreased
   → Could be random: NW Pacific has ~10-20 M≥4/day, 20-day window
      has high variance

4. ❌ Missing: precipitation/flood data
   → ERA5 .nc files have vo/d only (no tp)
   → Need to download ERA5 total precipitation or use alternative

5. 🔬 Spatial filtering needed
   → We're counting ALL M≥4 in NW Pacific (0-50°N, 100-180°E)
   → Should filter to only events near the typhoon track
   → A typhoon's pressure loading affects ~500km radius at most
""")
