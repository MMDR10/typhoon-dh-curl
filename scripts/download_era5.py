#!/usr/bin/env python3
"""
🌀 Typhoon Cyclogenesis Backtest — ERA5 Data Pipeline
=====================================================
Downloads ERA5 hourly vorticity data for 3 classic typhoons,
then runs Ô-HAT dH_curl early-warning detection.

Typhoons:
  1. 山竹 Mangkhut (2018) — long-track RI, JTWC TD: 2018-09-07 06Z
  2. 天鴿 Hato     (2017) — near-shore RI, JTWC TD: 2017-08-20 18Z
  3. 蘇拉 Saola    (2023) — recent, complex, JTWC TD: 2023-08-24 00Z

Dependencies: cdsapi, numpy, scipy, netCDF4, xarray
"""

import cdsapi
import numpy as np
import os, json, sys
from datetime import datetime, timedelta

# ─── CDS Config ───
CDS_CLIENT = cdsapi.Client()

# ─── Typhoon Definitions ───
TYPHOONS = {
    "mangkhut_2018": {
        "name": "Mangkhut (山竹)",
        "year": 2018,
        "jtwc_td_time": "2018-09-07T06:00",    # JTWC first warning as TD
        "genesis_lat": 12.0, "genesis_lon": 168.0,
        "bbox": [5, 120, 25, 180],              # N, W, S, E
        "days_before": 6,
        "days_after": 8,
        "peak_kts": 155,                         # Super Typhoon Cat 5
    },
    "hato_2017": {
        "name": "Hato (天鴿)",
        "year": 2017,
        "jtwc_td_time": "2017-08-20T18:00",
        "genesis_lat": 18.5, "genesis_lon": 129.0,
        "bbox": [10, 110, 30, 145],
        "days_before": 6,
        "days_after": 4,
        "peak_kts": 100,                         # Cat 3
    },
    "saola_2023": {
        "name": "Saola (蘇拉)",
        "year": 2023,
        "jtwc_td_time": "2023-08-24T00:00",
        "genesis_lat": 17.5, "genesis_lon": 126.0,
        "bbox": [10, 110, 30, 145],
        "days_before": 6,
        "days_after": 10,
        "peak_kts": 140,                         # Super Typhoon Cat 5
    },
}

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# ─── Download ERA5 ───
def download_era5(typhoon_id):
    """Download ERA5 hourly vorticity for a typhoon."""
    tc = TYPHOONS[typhoon_id]
    td_dt = datetime.fromisoformat(tc["jtwc_td_time"])
    
    start = (td_dt - timedelta(days=tc["days_before"])).strftime("%Y-%m-%d")
    end   = (td_dt + timedelta(days=tc["days_after"])).strftime("%Y-%m-%d")
    
    outfile = os.path.join(DATA_DIR, f"{typhoon_id}_era5.nc")
    
    if os.path.exists(outfile):
        print(f"  ✅ Already downloaded: {outfile}")
        return outfile
    
    bbox = tc["bbox"]  # [N, W, S, E]
    
    print(f"  Downloading {tc['name']}...")
    print(f"    Period: {start} → {end}")
    print(f"    BBox: {bbox}")
    print(f"    Variables: 850hPa vo, 200hPa d, MSLP")
    
    try:
        CDS_CLIENT.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': [
                    'relative_vorticity',       # vo at each pressure level
                    'divergence',                # d at each pressure level
                ],
                'pressure_level': ['200', '850'],
                'year': str(tc["year"]),
                'month': [f"{m:02d}" for m in set(
                    [td_dt.month, 
                     (td_dt - timedelta(days=tc["days_before"])).month,
                     (td_dt + timedelta(days=tc["days_after"])).month]
                )],
                'day': [f"{d:02d}" for d in range(1, 32)],
                'time': [f"{h:02d}:00" for h in range(24)],
                'area': bbox,
            },
            outfile
        )
        print(f"  ✅ Saved: {outfile}")
        return outfile
    
    except Exception as e:
        print(f"  ❌ Failed: {e}")
        return None

# ─── Test CDS connection ───
def test_cds():
    """Quick test that CDS API is working."""
    print("Testing CDS API connection...")
    try:
        # Tiny test download
        test_file = os.path.join(DATA_DIR, "_test_cds.nc")
        if os.path.exists(test_file):
            os.remove(test_file)
        CDS_CLIENT.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': ['relative_vorticity'],
                'pressure_level': ['850'],
                'year': '2023',
                'month': ['01'],
                'day': ['01'],
                'time': ['00:00'],
                'area': [25, 120, 20, 125],  # tiny box near Luzon
            },
            test_file
        )
        if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
            print(f"  ✅ CDS API working! Test file: {test_file} ({os.path.getsize(test_file)} bytes)")
            os.remove(test_file)
            return True
        else:
            print(f"  ⚠️  Test file too small or missing")
            return False
    except Exception as e:
        print(f"  ❌ CDS API error: {e}")
        return False

# ═══════════════ MAIN ═══════════════

def main():
    print("="*60)
    print("🌀 Typhoon Backtest — ERA5 Data Pipeline")
    print("="*60)
    
    # 1. Test CDS
    print("\n[0/4] Testing CDS connection...")
    if not test_cds():
        print("CDS API not working — aborting")
        sys.exit(1)
    
    # 2. Download each typhoon
    results = {}
    for i, (tid, tc) in enumerate(TYPHOONS.items(), 1):
        print(f"\n[{i}/4] {tc['name']} ({tid})")
        print(f"  JTWC TD: {tc['jtwc_td_time']}")
        print(f"  Genesis: {tc['genesis_lat']}N {tc['genesis_lon']}E")
        filepath = download_era5(tid)
        results[tid] = {
            "file": filepath,
            "success": filepath is not None,
        }
    
    # 3. Summary
    print("\n" + "="*60)
    print("Download Summary:")
    for tid, r in results.items():
        status = "✅" if r["success"] else "❌"
        print(f"  {status} {TYPHOONS[tid]['name']}: {r['file']}")
    
    # 4. Save metadata
    meta = {tid: {**tc, "era5_file": results[tid]["file"]} 
            for tid, tc in TYPHOONS.items()}
    meta_file = os.path.join(DATA_DIR, "typhoons_meta.json")
    with open(meta_file, "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print(f"\n  Metadata: {meta_file}")
    print("="*60)

if __name__ == "__main__":
    main()
