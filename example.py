import ee
ee.Authenticate()  # only needed once, creates a token on your machine
ee.Initialize(project='ee-bravemaster102')  # this is your project ID... if you already have one on GEE, just use it. else... https://code.earthengine.google.com/register
from sentinel2_ndvi_download import Sentinel2NDVIDownload

aoi = ee.Geometry.Polygon([
    [19.52, 64.19],  # upper left
    [19.58, 64.19],  # upper right
    [19.58, 64.15],  # lower right
    [19.52, 64.15],  # lower left
    [19.52, 64.19],  # closing point
])

downloader = Sentinel2NDVIDownload(
    aoi=aoi,
    output_dir="./output/kulbacksliden_2022_single",
    start_date="2022-05-01",
    end_date="2022-05-31",
    crs="EPSG:3006",
)

# Option 1 — download individual images (one file per cloud-free acquisition)
downloader.run()

# Option 2 — download composites (one file per time window, median of all images in the window)
downloader = Sentinel2NDVIDownload(
    aoi=aoi,
    output_dir="./output/kulbacksliden_2022_composite",
    start_date="2022-05-01",
    end_date="2022-05-31",
    crs="EPSG:3006",
)
# Uncomment ONE of the lines below and comment out Option 1 above
# downloader.run_composites(period="week")
downloader.run_composites(period="dekad")
# downloader.run_composites(period="biweekly")
# downloader.run_composites(period="month")