"""
sentinel2_ndvi_fetch.py
-----------------------
Fetches Sentinel-2 SR images for a given AOI and time range,
applies cloud/shadow masking, computes NDVI, and saves individual
GeoTIFF files at 10 m resolution.

Python API usage
----------------
from sentinel2_ndvi_fetch import Sentinel2NDVIDownload

downloader = Sentinel2NDVIDownload(
    aoi=my_ee_geometry,
    output_dir="./output",
    start_date="2022-05-01",
    end_date="2022-09-30",
    cloud_filter=80,
    cld_prb_thresh=50,
    nir_drk_thresh=0.15,
    cld_prj_dist=2,
    buffer=100,
    nodata_thresh=80,
)
downloader.run()

CLI usage
---------
python sentinel2_ndvi_fetch.py \\
    --output_dir ./output/krycklan_2022 \\
    --start_date 2022-05-01 \\
    --end_date   2022-09-30 \\
    --bbox 19.75 64.05 20.05 64.25   # xmin ymin xmax ymax (lon/lat WGS84)

All cloud-masking parameters are optional and have sensible defaults.
Run  python sentinel2_ndvi_fetch.py --help  for full usage.

Dependencies (install in a fresh venv)
---------------------------------------
    pip install earthengine-api geemap rasterio numpy pyproj
"""

import os
from datetime import datetime

import ee
import geemap
import numpy as np
import pyproj
import rasterio

# ── Make sure PROJ finds its data ─────────────────────────────────────────────
proj_data_path = pyproj.datadir.get_data_dir()
os.environ["PROJ_DATA"] = proj_data_path
os.environ["PROJ_LIB"] = proj_data_path

# ── Authenticate / initialise GEE ─────────────────────────────────────────────
# try:
#     ee.Initialize()
# except Exception:
#     ee.Authenticate()
#     ee.Initialize()


# ─────────────────────────────────────────────────────────────────────────────
class Sentinel2NDVIDownload:
    """
    Downloads individual cloud-masked Sentinel-2 NDVI GeoTIFFs for an AOI.

    Parameters
    ----------
    aoi : ee.Geometry
        Area of interest (pass a pre-built ee.Geometry).
    output_dir : str
        Directory where output TIFFs will be written.
    start_date, end_date : str  (YYYY-MM-DD)
        Temporal window for the image search.
    cloud_filter : int
        Maximum CLOUDY_PIXEL_PERCENTAGE accepted from the S2 metadata.
    cld_prb_thresh : int
        Cloud-probability threshold (%) for the s2cloudless band.
    nir_drk_thresh : float
        NIR reflectance threshold (fraction of SR_BAND_SCALE) for dark pixels
        used in shadow detection.
    cld_prj_dist : int
        Maximum cloud-projection distance for shadow detection (km).
    buffer : int
        Dilation radius (metres) applied to the combined cloud/shadow mask.
    nodata_thresh : float
        Maximum tolerated no-data percentage in the NDVI band; images above
        this threshold are discarded.
    scale : int
        Output pixel resolution in metres (default 10 m).
    crs : str
        Output coordinate reference system as an EPSG string (default
        'EPSG:32633', UTM zone 33N, suitable for most of Sweden).
        Use 'EPSG:3006' for SWEREF99 TM (Swedish national standard).
        Avoid 'EPSG:4326' for metric work: pixel size in degrees makes
        the 10 m scale parameter meaningless.
    """

    def __init__(
        self,
        aoi: ee.Geometry,
        output_dir: str,
        start_date: str,
        end_date: str,
        cloud_filter: int = 80,
        cld_prb_thresh: int = 50,
        nir_drk_thresh: float = 0.15,
        cld_prj_dist: int = 2,
        buffer: int = 100,
        nodata_thresh: float = 80,
        scale: int = 10,
        crs: str = "EPSG:32633",
    ):

        self.aoi = aoi
        self.output_dir = output_dir
        self.start_date = start_date
        self.end_date = end_date
        self.cloud_filter = cloud_filter
        self.cld_prb_thresh = cld_prb_thresh
        self.nir_drk_thresh = nir_drk_thresh
        self.cld_prj_dist = cld_prj_dist
        self.buffer = buffer
        self.nodata_thresh = nodata_thresh
        self.scale = scale
        self.crs = crs

        os.makedirs(self.output_dir, exist_ok=True)

    # ── Cloud / shadow masking ────────────────────────────────────────────────

    def _get_s2_sr_cld_col(self) -> ee.ImageCollection:
        """Join S2 SR collection with s2cloudless probability collection."""
        s2_sr = (
            ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(self.aoi)
            .filterDate(self.start_date, self.end_date)
            .filter(ee.Filter.lte("CLOUDY_PIXEL_PERCENTAGE", self.cloud_filter))
        )
        s2_cld = (
            ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
            .filterBounds(self.aoi)
            .filterDate(self.start_date, self.end_date)
        )
        return ee.ImageCollection(
            ee.Join.saveFirst("s2cloudless").apply(
                primary=s2_sr,
                secondary=s2_cld,
                condition=ee.Filter.equals(
                    leftField="system:index", rightField="system:index"
                ),
            )
        )

    def _add_cloud_bands(self, img: ee.Image) -> ee.Image:
        cld_prb = ee.Image(img.get("s2cloudless")).select("probability")
        is_cloud = cld_prb.gt(self.cld_prb_thresh).rename("clouds")
        return img.addBands(ee.Image([cld_prb, is_cloud]))

    def _add_shadow_bands(self, img: ee.Image) -> ee.Image:
        SR_BAND_SCALE = 1e4
        not_water = img.select("SCL").neq(6)
        dark_pixels = (
            img.select("B8")
            .lt(self.nir_drk_thresh * SR_BAND_SCALE)
            .multiply(not_water)
            .rename("dark_pixels")
        )
        shadow_azimuth = ee.Number(90).subtract(
            ee.Number(img.get("MEAN_SOLAR_AZIMUTH_ANGLE"))
        )
        cld_proj = (
            img.select("clouds")
            .directionalDistanceTransform(shadow_azimuth, self.cld_prj_dist * 10)
            .reproject(crs=img.select(0).projection(), scale=100)
            .select("distance")
            .mask()
            .rename("cloud_transform")
        )
        shadows = cld_proj.multiply(dark_pixels).rename("shadows")
        return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

    def _add_cld_shdw_mask(self, img: ee.Image) -> ee.Image:
        img = self._add_cloud_bands(img)
        img = self._add_shadow_bands(img)
        is_cld_shdw = img.select("clouds").add(img.select("shadows")).gt(0)
        is_cld_shdw = (
            is_cld_shdw.focalMin(2)
            .focalMax(self.buffer * 2 / 20)
            .reproject(crs=img.select([0]).projection(), scale=20)
            .rename("cloudmask")
        )
        return img.addBands(is_cld_shdw)

    def _apply_cld_shdw_mask(self, img: ee.Image) -> ee.Image:
        return img.select("B.*").updateMask(img.select("cloudmask").Not())

    # ── NDVI ─────────────────────────────────────────────────────────────────

    def _compute_ndvi(self, img: ee.Image) -> ee.Image:
        ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI")
        return img.addBands(ndvi)

    # ── No-data filter ────────────────────────────────────────────────────────

    def _filter_nodata(self, img: ee.Image) -> ee.Image:
        """Tag the image with its no-data percentage; mask images that exceed
        the threshold (GEE will exclude them via the notNull filter below)."""
        nodata_mask = img.select("NDVI").eq(0)
        masked_out = nodata_mask.mask().reduce(ee.Reducer.min()).Not()
        masked_area = (
            ee.Image.pixelArea()
            .updateMask(masked_out)
            .reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=img.geometry(),
                scale=10,
                maxPixels=1e13,
            )
            .getNumber("area")
        )
        total_area = img.geometry().area(10)
        nodata_pct = masked_area.divide(total_area).multiply(100)
        return img.set("nodata_percentage", nodata_pct)

    # ── Full pipeline ─────────────────────────────────────────────────────────

    def _build_collection(self) -> ee.ImageCollection:
        return (
            self._get_s2_sr_cld_col()
            .map(self._add_cld_shdw_mask)
            .map(self._apply_cld_shdw_mask)
            .map(self._compute_ndvi)
            .map(self._filter_nodata)
            .filter(ee.Filter.notNull(["nodata_percentage"]))
            .filter(ee.Filter.lt("nodata_percentage", self.nodata_thresh))
        )

    def run(self) -> list[str]:
        """
        Execute the fetch-mask-save pipeline.

        Returns
        -------
        list[str]
            Paths to all written NDVI GeoTIFFs.
        """
        col = self._build_collection()
        ndvi_list = col.select("NDVI").toList(col.size())
        n = ndvi_list.size().getInfo()
        print(f"Found {n} usable Sentinel-2 images.")

        saved_paths: list[str] = []
        for i in range(n):
            img = ee.Image(ndvi_list.get(i))
            ts_ms = img.get("system:time_start").getInfo()
            date_str = datetime.utcfromtimestamp(ts_ms / 1000).strftime("%Y-%m-%d")

            out_path = os.path.join(self.output_dir, f"sentinel2_ndvi_{date_str}.tif")
            print(f"  [{i+1}/{n}] Downloading {date_str} → {out_path}")
            geemap.ee_export_image(
                img,
                filename=out_path,
                scale=self.scale,
                region=self.aoi,
                crs=self.crs,
                file_per_band=False,
            )
            saved_paths.append(out_path)

        print(f"\nDone. {len(saved_paths)} file(s) written to: {self.output_dir}")
        return saved_paths


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Download cloud-masked Sentinel-2 NDVI GeoTIFFs for a bounding box.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required
    parser.add_argument("--output_dir",  required=True, help="Directory to write output TIFFs.")
    parser.add_argument("--start_date",  required=True, help="Start date (YYYY-MM-DD).")
    parser.add_argument("--end_date",    required=True, help="End date   (YYYY-MM-DD).")
    parser.add_argument(
        "--bbox",
        required=True,
        nargs=4,
        type=float,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        help="Bounding box in WGS84 lon/lat: xmin ymin xmax ymax.",
    )

    # Optional cloud-masking parameters
    parser.add_argument("--cloud_filter",   type=int,   default=80,   help="Max CLOUDY_PIXEL_PERCENTAGE.")
    parser.add_argument("--cld_prb_thresh", type=int,   default=50,   help="Cloud probability threshold (%%).")
    parser.add_argument("--nir_drk_thresh", type=float, default=0.15, help="NIR dark-pixel threshold (fraction).")
    parser.add_argument("--cld_prj_dist",   type=int,   default=2,    help="Cloud projection distance (km).")
    parser.add_argument("--buffer",         type=int,   default=100,  help="Cloud/shadow dilation buffer (m).")
    parser.add_argument("--nodata_thresh",  type=float, default=80,   help="Max tolerated no-data percentage.")
    parser.add_argument("--scale",          type=int,   default=10,   help="Output pixel resolution (m).")
    parser.add_argument("--crs",            type=str,   default="EPSG:32633", help="Output CRS (e.g. EPSG:32633, EPSG:3006).")

    args = parser.parse_args()

    xmin, ymin, xmax, ymax = args.bbox
    aoi = ee.Geometry.Polygon(
        [[xmin, ymax], [xmax, ymax], [xmax, ymin], [xmin, ymin], [xmin, ymax]]
    )

    downloader = Sentinel2NDVIDownload(
        aoi=aoi,
        output_dir=args.output_dir,
        start_date=args.start_date,
        end_date=args.end_date,
        cloud_filter=args.cloud_filter,
        cld_prb_thresh=args.cld_prb_thresh,
        nir_drk_thresh=args.nir_drk_thresh,
        cld_prj_dist=args.cld_prj_dist,
        buffer=args.buffer,
        nodata_thresh=args.nodata_thresh,
        scale=args.scale,
        crs=args.crs,
    )
    downloader.run()