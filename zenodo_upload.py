#!/usr/bin/env python3
import requests, json, os

ZENODO_KEY = "YOUR_ZENODO_TOKEN_HERE"  # Replace with actual token
HEADERS = {"Authorization": f"Bearer {ZENODO_KEY}"}
BASE = "/app/working/typhoon-dh-curl"

# Create deposit
resp = requests.post("https://zenodo.org/api/deposit/depositions", json={}, headers=HEADERS)
print(f"Create: {resp.status_code}")
if resp.status_code != 201:
    print(resp.text[:500])
    exit(1)

data = resp.json()
dep_id = data["id"]
bucket = data["links"]["bucket"]
print(f"Deposit ID: {dep_id}")

# Upload files
files = [
    ("paper.md", f"{BASE}/paper.md"),
    ("README.md", f"{BASE}/README.md"),
]
for fname, fpath in files:
    with open(fpath, "rb") as f:
        r = requests.put(f"{bucket}/{fname}", data=f, headers=HEADERS)
    print(f"Upload {fname}: {r.status_code}")

# Metadata
desc = """<h2>0-Lag Phase Locking: Detecting Tropical Cyclone Genesis via 850 hPa Helicity on Edge Hardware</h2>
<p><strong>DR | July 2026 | CC BY 4.0</strong></p>
<h3>Abstract</h3>
<p>Using a 6 GB RAM consumer NAS, 54 grid-search experiments across operator, level, and smoothing kernel. The winning config (raw 850hPa vorticity, core=5°) produces dH_curl with a W-shaped pre-landfall trajectory across 6 typhoons. Three genesis fingerprint types identified. Negative atmosphere-seismic coupling result constrains phase-space topology.</p>
<h3>Key Results</h3>
<table border='1' cellpadding='4'>
<tr><th>Typhoon</th><th>U-strength</th><th>Mode</th><th>JTWC Peak</th></tr>
<tr><td>Meranti</td><td><b>1.66</b></td><td>Deep-U</td><td>120 kt</td></tr>
<tr><td>Mangkhut</td><td><b>1.63</b></td><td>Flat</td><td>155 kt</td></tr>
<tr><td>Saola</td><td><b>1.08</b></td><td>Flat</td><td>140 kt</td></tr>
<tr><td>Goni</td><td>0.38</td><td>Deep-U</td><td>120 kt</td></tr>
<tr><td>Hato</td><td>0.34</td><td>Asym</td><td>100 kt</td></tr>
<tr><td>Hagibis</td><td>0.10</td><td>Flat</td><td>100 kt</td></tr>
</table>
<p>GitHub: <a href='https://github.com/MMDR10/typhoon-dh-curl'>MMDR10/typhoon-dh-curl</a></p>"""

meta = {
    "metadata": {
        "title": "0-Lag Phase Locking: Detecting Tropical Cyclone Genesis via 850 hPa Helicity on Edge Hardware",
        "upload_type": "publication",
        "publication_type": "preprint",
        "description": desc,
        "creators": [{"name": "DR", "affiliation": "Independent Researcher"}],
        "license": "cc-by-4.0",
        "keywords": ["typhoon", "tropical cyclone", "helicity", "vorticity", "ERA5", "rapid intensification", "genesis detection", "NWP"],
        "access_right": "open",
    }
}
r = requests.put(f"https://zenodo.org/api/deposit/depositions/{dep_id}", json=meta, headers=HEADERS)
print(f"Metadata: {r.status_code}")

# Publish
r = requests.post(f"https://zenodo.org/api/deposit/depositions/{dep_id}/actions/publish", headers=HEADERS)
print(f"Publish: {r.status_code}")
if r.status_code == 202:
    pub = r.json()
    print(f"\n✅ DOI: {pub['doi']}")
    print(f"URL: https://zenodo.org/record/{dep_id}")
else:
    print(r.text[:500])
