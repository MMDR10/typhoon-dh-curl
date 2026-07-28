# 0-Lag Typhoon Genesis Detection via 850 hPa Helicity

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21650650.svg)](https://doi.org/10.5281/zenodo.21650650)

**A geometric caliper that detects tropical cyclone organization at zero temporal lag — built and validated on a 6 GB RAM consumer NAS.**

---

## What This Is

Traditional NWP models smooth away the vorticity structures that precede typhoon genesis. This project demonstrates an alternative: a **core-shell helicity decomposition** (`dH_curl`) that operates on raw ERA5 850 hPa vorticity and captures the characteristic **W-shaped pre-landfall trajectory** — pre-conditioning → relaxation → rebound — with zero lag relative to JTWC TD declaration.

## Key Results

| Typhoon | U-strength | Mode | JTWC Peak |
|---------|:----------:|------|:---------:|
| Meranti | **1.66** | Deep-U | 120 kt (Cat 4) |
| Mangkhut | **1.63** | Flat | 155 kt (Cat 5) |
| Saola | **1.08** | Flat | 140 kt (Cat 5) |
| Goni | 0.38 | Deep-U | 120 kt (Cat 4) |
| Hato | 0.34 | Asymmetric-Post | 100 kt (Cat 3) |
| Hagibis | **0.10** | Flat | 100 kt (Cat 3) |

- **54 grid-search experiments** across operator × level × sigma × core
- **Three fingerprint types**: Deep-U, Asymmetric-Post, Flat
- **Negative result**: No atmosphere–seismic coupling (|r| < 0.18) — physically meaningful null
- **All compute**: 6 GB RAM NAS. No GPU. No HPC.

## Structure

```
typhoon-dh-curl/
├── paper.md                         # Full paper (English)
├── scripts/
│   ├── analysis_6typhoon.py         # 6-typhoon cross-analysis
│   ├── coupling_analysis.py         # Atmosphere–flood–seismic v1
│   ├── coupling_significance.py     # Statistical tests
│   ├── coupling_v2.py               # Spatial filtering + precipitation
│   └── download_era5.py             # ERA5 CDS API pipeline
├── data/6typhoon_results/
│   ├── *_dh_curl.csv                # Per-typhoon dH_curl time series (6 files)
│   ├── 6typhoon_summary.json        # Aggregate mode + U-strength
│   ├── *_helicity_summary.json      # Grid search intermediates (23 files)
│   └── *_seismic_batch_summary.json # Earthquake coupling raw data (2 files)
└── output/
    ├── coupling_analysis_v1.json     # Full-domain coupling results
    └── coupling_analysis_v2.json     # Spatially-filtered coupling results
```

## Quick Start

```bash
pip install numpy scipy xarray netCDF4 cdsapi

# Run the 6-typhoon cross-analysis
python scripts/analysis_6typhoon.py

# Run the coupling survey
python scripts/coupling_analysis.py
python scripts/coupling_significance.py
```

**Note:** ERA5 .nc files are not included in this repo (35–105 MB each). Download them via `scripts/download_era5.py` with a valid CDS API key.

## Hardware Requirements

- **Minimum:** 4 GB RAM, 2 GB disk
- **Tested on:** Synology DS218 (6 GB RAM, ARM dual-core)
- Full pipeline runs in <30 minutes on Raspberry Pi 4 equivalent

## Citation

```bibtex
@misc{dr2026typhoon,
  title={0-Lag Phase Locking: Detecting Tropical Cyclone Genesis via 850 hPa Helicity on Edge Hardware},
  author={DR},
  year={2026},
  doi={10.5281/zenodo.21650650},
  note={Independent Research}
}
```

## Related Work

- [Noise Topography](https://github.com/MMDR10/noise-topography) — The Ô-HAT helicity framework
- [ENSO Watch](https://mmdr10.github.io/enso-watch) — Real-time ENSO phase-space monitoring

---

*Built by DR under MKP's direction. All errors are mine; all insights are from the data.*
