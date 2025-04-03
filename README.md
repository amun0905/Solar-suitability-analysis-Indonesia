# ☀️ Project Solaris: Solar Energy Site Suitability in Indonesia

Project Solaris is a spatial data science project focused on identifying and evaluating the most suitable locations for large-scale solar energy projects in Indonesia. Using a combination of land cover data, solar radiation, infrastructure proximity, environmental restrictions, and economic evaluation, the project performs a comprehensive multi-criteria assessment for solar potential.

---

## 📍 Objectives

- Identify land in Indonesia suitable for solar development
- Integrate environmental and infrastructure data (grid, roads, mountains, nature reserves)
- Use solar irradiance data to estimate energy production
- Apply AHP-style scoring and economic evaluation (NPV, LCOE) to prioritize development zones

---

## 🧰 Tools & Packages

- **R Language** with:
  - `sf`, `terra`, `raster` — spatial data processing
  - `tidyverse`, `ggplot2`, `tmap`, `ggspatial` — data wrangling and visualization
  - `ncdf4` — handling NetCDF climate data
  - `gstat` — spatial interpolation (IDW)
  - `readxl`, `writexl` — Excel I/O

---

## 🗂️ Data Inputs

- **Land Cover:** Copernicus LC, filtered for suitability
- **Solar Radiation:** NetCDF from Copernicus for 4 seasonal days in 2023
- **Infrastructure:** Electric grid and road buffers
- **Constraints:** Populated areas, mountainous terrain, nature reserves
- **Economic Variables:** CAPEX, grid connection costs, energy prices

---

## 📈 Workflow

1. **Preprocessing & Land Masking**
   - Extract LCCS types and reduce to suitable land
   - Clip to Indonesia boundary

2. **Solar Irradiance Analysis**
   - Combine multiple seasonal solar datasets
   - Estimate energy output using custom function

3. **Spatial Interpolation**
   - Use IDW to create smooth irradiance surface

4. **Buffering & Constraints**
   - Apply exclusion zones (mountains, protected areas, dense population)

5. **Scoring & Suitability Analysis**
   - Rank sites based on weighted criteria:
     - Solar power (65%)
     - Distance to grid (23%)
     - Distance to population (12%)

6. **Economic Modeling**
   - Estimate site capacity, CAPEX
   - Calculate Net Present Value (NPV) and Levelized Cost of Electricity (LCOE)

---

## 🗺️ Outputs

- Raster and vector outputs for suitability zones
- Ranked site matrix with energy potential and economic scores
- AHP-weighted scores for site selection
- Final plots of viable solar zones and economic viability

---

## 📊 Example Outputs

You can include some of these if available:
```r
ggplot() +
  geom_sf(data = boundaryU) +
  geom_sf(data = matrix2, aes(fill = "kwh")) +
  labs(title = "Solar Suitability Map")
