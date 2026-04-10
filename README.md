# geemap_ndvi_download / Sentinel-2 NDVI Download

A Python tool for downloading cloud-masked Sentinel-2 NDVI images at 10 m spatial resolution for any area of interest, using Google Earth Engine (GEE) and `geemap`.

Each Sentinel-2 acquisition that passes the cloud/shadow filters is saved as an individual GeoTIFF — no temporal compositing, no averaging. The output is ready for downstream analysis such as NDVI time series construction or data fusion workflows.

---

## Repository structure

```
.
├── sentinel2_ndvi_download.py   # Core class — do not edit for day-to-day use
├── example.py                   # Edit this to define your AOI and run the download
├── requirements.txt             # All Python dependencies (pip)
└── README.md
```

Output files are written to the `output/` folder (git-ignored).

---

## Requirements

- A Google account registered for **Google Earth Engine** (see [GEE access](#google-earth-engine-access) below)
- Python 3.11
- Dependencies listed in `requirements.txt`

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/bravemaster3/geemap_ndvi_download.git
cd path_on_your_PC/to/geemap_ndvi_download
```

### 2. Create a virtual environment

**Option A — conda (recommended)**
```bash
conda create -n sentinel2_ndvi python=3.11
conda activate sentinel2_ndvi
```

**Option B — venv**
```bash
python -m venv sentinel2_ndvi
# Windows:
sentinel2_ndvi\Scripts\activate
# macOS/Linux:
source sentinel2_ndvi/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

---

## Google Earth Engine access

Since 2023, GEE requires a registered Google Cloud project even for non-commercial/research use.

1. Go to https://code.earthengine.google.com/register
2. Click **Register a Noncommercial or Commercial Cloud project** and follow the steps
3. Note your **project ID** (e.g. `ee-yourname`) — you will need it in `example.py`

To authenticate for the first time, run this once in Python:

```python
import ee
ee.Authenticate()   # opens a browser — log in with the account registered for GEE
```

This saves a token on your machine so you only need to do it once.

---

## Usage

Open [`example.py`](example.py) and edit the AOI bounding box, output directory, date range, and GEE project ID, then run it in Spyder or from the command line.

The file is fully commented and ready to use as a starting point.

Each image that passes cloud filtering is saved as:

```
output/kulbacksliden_2022_2022/sentinel2_ndvi_2022-05-14.tif
output/kulbacksliden_2022_2022/sentinel2_ndvi_2022-05-24.tif
...
```

---

## Cloud masking parameters

There are some default settings. If you need stricter or more permissive filtering, pass any of the following optional parameters:

```python
downloader = Sentinel2NDVIDownload(
    aoi=aoi,
    output_dir="./output/krycklan_2022",
    start_date="2022-05-01",
    end_date="2022-09-30",
    crs="EPSG:3006",

    # --- Optional cloud masking settings (defaults shown) ---
    cloud_filter=80,        # Max CLOUDY_PIXEL_PERCENTAGE in S2 metadata (%)
                            # Lower = fewer but cleaner images (e.g. 60 for stricter filtering)

    cld_prb_thresh=50,      # Cloud probability threshold from s2cloudless (%)
                            # Lower = more aggressive cloud masking

    nir_drk_thresh=0.15,    # NIR reflectance threshold for shadow detection
                            # Fraction of surface reflectance scale (10 000)

    cld_prj_dist=2,         # Maximum cloud projection distance for shadow detection (km)

    buffer=100,             # Dilation radius applied to cloud/shadow mask (metres)
                            # Larger = more conservative, masks more around cloud edges

    nodata_thresh=80,       # Max tolerated no-data percentage in the NDVI band (%)
                            # Images with more missing data than this are discarded

    scale=10,               # Output pixel resolution in metres (default: 10 m)
)
```

**Tip for Sweden:** cloud cover is frequent, especially early and late in the season. If you are getting too few images, try relaxing `cloud_filter` to 90 and `nodata_thresh` to 90. If images look noisy with cloud artefacts, lower `cld_prb_thresh` to 40.

---

## Output coordinate system

By default, images are saved in **EPSG:3006** (SWEREF99 TM), the Swedish national standard. This ensures the 10 m pixel size is meaningful (metres, not degrees).

To use UTM zone 33N instead:
```python
crs="EPSG:32633"
```

Avoid `EPSG:4326` (WGS84 geographic) for metric work — pixel size in degrees makes the 10 m scale parameter meaningless.

---

## Running from the command line

The script can also be run directly from a terminal without editing any file:

```bash
python sentinel2_ndvi_download.py \
    --output_dir ./output/kulbacksliden_2022 \
    --start_date 2022-05-01 \
    --end_date   2022-09-30 \
    --bbox 19.52 64.15 19.58 64.19 \
    --crs EPSG:3006
```

`--bbox` takes four values: `xmin ymin xmax ymax` (lon/lat, WGS84).

Run `python sentinel2_ndvi_download.py --help` for the full list of options.

> **Note:** when using the CLI, GEE initialisation is handled automatically inside the script. Make sure you have already authenticated once (`ee.Authenticate()`) so a valid token exists on your machine.

---

## Dependencies

| Package | Purpose |
|---|---|
| `earthengine-api` | Google Earth Engine Python API |
| `geemap` | Downloading GEE images as local GeoTIFFs |
| `rasterio` | Reading and writing GeoTIFF files |
| `numpy` | Array operations |
| `pyproj` | PROJ data path setup |

Install all at once:
```bash
pip install -r requirements.txt
```

---

## Troubleshooting

**`EEException: Earth Engine client library not initialized`**
Make sure `ee.Initialize(project='your-project-id')` is called before constructing any `ee.Geometry` or instantiating the class.

**`EEException: Not signed up for Earth Engine or project is not registered`**
Your Google account is not linked to a registered GEE project. Follow the steps at https://code.earthengine.google.com/register.

**`CondaSSLError` when creating the environment**
Run `conda create` from the `base` environment (deactivate any active env first with `conda deactivate`).

**No images found / empty output folder**
Try relaxing the cloud filters (`cloud_filter=90`, `nodata_thresh=90`) or extending the date range. Cloud cover in boreal Sweden is high, especially outside June–August.
