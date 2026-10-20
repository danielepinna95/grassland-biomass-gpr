# Data

The field dataset used in the paper (950 AGB samples, three grasslands in north-east Italy, May 2023 – Nov 2024) is **not yet publicly available**. For access, contact the corresponding author.

To run the pipeline on your own data, save a CSV as `data/grassland_agb_s2_planetscope.csv` (or change `CONFIG$data_path` in `R/00_utils.R`) with these columns:

| Column | Description | Unit | Required |
|---|---|---|---|
| `sample_id` | Unique row identifier | – | yes |
| `point_id` | Survey point identifier within a site (points can be revisited across surveys) | – | yes |
| `lon`, `lat` | Point coordinates, WGS84 (EPSG:4326) | degrees | yes (Moran's I) |
| `site` | Site / field name | – | yes (field-level aggregation) |
| `survey` | Survey date or month, e.g. `2023-05` | – | yes (field-level aggregation) |
| `agb_kg_ha` | Above-ground dry matter, cut-and-dry (**model target**) | kg DM ha⁻¹ | yes |
| `B2`, `B3`, `B4`, `B5`, `B6`, `B7`, `B8`, `B8A`, `B11`, `B12` | Sentinel-2 L2A surface reflectance | reflectance × 10,000 | for `S2` models |
| `PLB2` … `PLB8` | PlanetScope SuperDove surface reflectance, bands 2–8 | reflectance × 10,000 | for `PS` models |

Missing values are empty cells. In the paper, each survey uses the cloud-free scene closest in time, within 10 days. `python/extract_reflectance.py` builds these columns from your own raster stacks.
