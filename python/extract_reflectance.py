"""Sample satellite reflectance at field survey points.

For every survey point the script picks the image whose acquisition date is
closest to the survey date (within ``--max-days``, default 10, as in the
paper) and reads the pixel value of every band at the point location.

Inputs
------
points : CSV with columns ``sample_id, lon, lat, date`` (date as YYYY-MM-DD).
rasters: a folder of multiband GeoTIFFs whose file names contain the
         acquisition date as YYYYMMDD, e.g. ``S2_20230815.tif``.
bands  : band names in raster order, e.g. ``B2 B3 B4 B5 B6 B7 B8 B8A B11 B12``.

Example
-------
python python/extract_reflectance.py points.csv rasters/ S2_samples.csv \
    --bands B2 B3 B4 B5 B6 B7 B8 B8A B11 B12
"""

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio

DATE_RE = re.compile(r"(\d{8})")


def raster_dates(folder: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(folder.glob("*.tif*")):
        m = DATE_RE.search(path.name)
        if m:
            rows.append({"raster": path, "img_date": datetime.strptime(m.group(1), "%Y%m%d")})
    if not rows:
        raise SystemExit(f"No GeoTIFFs with a YYYYMMDD date in their name found in {folder}")
    return pd.DataFrame(rows)


def match_images(points: pd.DataFrame, images: pd.DataFrame, max_days: int) -> pd.DataFrame:
    points = points.copy()
    points["date"] = pd.to_datetime(points["date"])
    best = []
    for d in points["date"]:
        gap = (images["img_date"] - d).abs()
        i = gap.idxmin()
        best.append((images.at[i, "raster"], images.at[i, "img_date"], gap[i].days) if gap[i].days <= max_days
                    else (None, pd.NaT, gap[i].days))
    points[["raster", "img_date", "days_gap"]] = pd.DataFrame(best, index=points.index)
    return points


def sample(points: pd.DataFrame, bands: list[str]) -> pd.DataFrame:
    out = []
    for raster, group in points.dropna(subset=["raster"]).groupby("raster"):
        with rasterio.open(raster) as src:
            if src.count != len(bands):
                raise SystemExit(f"{raster.name}: {src.count} bands, but {len(bands)} band names were given")
            gdf = gpd.GeoDataFrame(group, geometry=gpd.points_from_xy(group.lon, group.lat), crs=4326).to_crs(src.crs)
            values = list(src.sample([(p.x, p.y) for p in gdf.geometry]))
        out.append(pd.concat([group.reset_index(drop=True),
                              pd.DataFrame(values, columns=bands)], axis=1))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("points", type=Path)
    ap.add_argument("rasters", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--bands", nargs="+", required=True)
    ap.add_argument("--max-days", type=int, default=10)
    args = ap.parse_args()

    points = pd.read_csv(args.points)
    matched = match_images(points, raster_dates(args.rasters), args.max_days)
    missing = matched["raster"].isna().sum()
    if missing:
        print(f"{missing} points have no image within {args.max_days} days and are skipped")

    result = sample(matched, args.bands).drop(columns=["raster"])
    result.to_csv(args.output, index=False)
    print(f"Wrote {len(result)} rows to {args.output}")


if __name__ == "__main__":
    main()
