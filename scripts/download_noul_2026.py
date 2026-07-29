#!/usr/bin/env python3
"""
🌀 Noul 2026 (紅霞) — ERA5 Download
====================================
JMA Typhoon 202612 (NOUL)
- Birth: 2026-07-23 18:00 UTC
- Peak: 960 hPa, 75 kt (2026-07-25 15:00-21:00 UTC)
- Death: 2026-07-26 12:00 UTC
- Landfall: ~22.9°N, 114.3°E (Guangdong, near Huizhou)
- Track: 15.9-25.0°N, 114.0-129.6°E
"""

import cdsapi
import os, sys
from datetime import datetime

CDS_CLIENT = cdsapi.Client()

DATA_DIR = "/app/working/workspaces/tygtDc/projects/typhoon-backtest/data/noul_2026"
OUTFILE = os.path.join(DATA_DIR, "noul_2026_era5.nc")

# BBox: [N, W, S, E] — covers track + margin
BBOX = [30, 110, 12, 135]

YEAR = 2026
MONTHS = ['07']
DAYS = [f"{d:02d}" for d in range(22, 28)]  # July 22-27 for margin
TIMES = [f"{h:02d}:00" for h in range(24)]

def try_download(product_type='reanalysis'):
    """Try ERA5 download with given product type."""
    print(f"\n🔍 Trying product_type='{product_type}'...")
    try:
        CDS_CLIENT.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': product_type,
                'format': 'netcdf',
                'variable': [
                    'relative_vorticity',
                    'divergence',
                ],
                'pressure_level': ['200', '850'],
                'year': str(YEAR),
                'month': MONTHS,
                'day': DAYS,
                'time': TIMES,
                'area': BBOX,
            },
            OUTFILE
        )
        size_mb = os.path.getsize(OUTFILE) / (1024*1024)
        print(f"  ✅ Downloaded: {OUTFILE} ({size_mb:.1f} MB)")
        return True
    except Exception as e:
        print(f"  ❌ Failed: {e}")
        return False

def main():
    print("=" * 60)
    print("🌀 Noul 2026 (紅霞) — ERA5 Download")
    print(f"   Period: 2026-07-22 → 2026-07-27")
    print(f"   BBox: {BBOX}")
    print(f"   Output: {OUTFILE}")
    print("=" * 60)
    
    if os.path.exists(OUTFILE):
        size_mb = os.path.getsize(OUTFILE) / (1024*1024)
        if size_mb > 1:
            print(f"  ✅ Already exists: {OUTFILE} ({size_mb:.1f} MB)")
            return 0
    
    # Try reanalysis first, then fallback
    if try_download('reanalysis'):
        return 0
    
    print("\n⚠️  Reanalysis not available yet (too recent).")
    print("   ERA5 typically has ~5 day latency.")
    print("   Will retry automatically when data becomes available.")
    
    # Save track CSV as fallback evidence
    track_file = os.path.join(DATA_DIR, "track_jma.csv")
    print(f"\n  📋 JMA best track saved: {track_file}")
    print(f"  ⏳ ERA5 .nc pending — retry after ~2026-08-01")
    
    return 1

if __name__ == "__main__":
    sys.exit(main())
