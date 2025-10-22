# Solar Power Plant Investment in Indonesia: A Decarbonized Electrical Supply Vision to 2040

The project applies Geographic Information Systems (GIS) and Analytical Hierarchy Process (AHP) techniques to identify suitable areas for solar power plant development in Indonesia by 2040.  
The analysis combines spatial, environmental, and socioeconomic datasets to assess land suitability, project electricity demand, and evaluate financial feasibility.

---

## Project Overview

Indonesia aims to achieve net-zero emissions by 2060 while maintaining strong economic growth.  
This study develops a spatial and financial model to identify suitable solar power plant locations that would meet projected 2040 electricity demand under a decarbonization pathway.

**Key objectives:**
- Estimate Indonesia’s 2040 electricity demand based on GDP and population growth.
- Identify spatially suitable areas for new solar power plants using AHP in GIS.
- Quantify the total renewable capacity required.
- Evaluate the financial viability of the proposed investment through NPV and LCOE analyses.

---

## Methodology

### 1. Demand Projection
Electricity demand for 2040 was projected using historical correlations between GDP per capita and power consumption.  
- Projected population (2040): 312 million  
- Projected GDP per capita (2040): \$8,895  
- Estimated total electricity demand: 678 TWh  

From this, total required generation capacity was estimated as 209 GW, with 79.4 GW expected from renewables.  
Given that 19 GW of renewable capacity already exists or is planned, an additional 60.4 GW is required — assumed to come from new solar plants.

---

## Spatial Analysis Framework

### Step 1 – Define Constraints (Unsuitable Areas)
Excluded regions include:
- Forests and tree cover  
- Water bodies and urban areas  
- Populated areas  
- Protected natural areas  
- Mountainous terrain (>1500 m elevation or within 500 m buffer)  

### Step 2 – Select Suitability Factors
| Factor | Description | Source |
|--------|--------------|---------|
| **Solar Radiation** | Average downward surface solar radiation (ERA5 data, 4 representative days & hours) | Muñoz Sabater (2019) |
| **Proximity to Electrical Grid** | Distance to Indonesia’s national grid | OpenStreetMap (2018) |
| **Proximity to Population Centers** | Distance to populated areas (proxy for transmission cost) | WorldPop (2020) |

### Step 3 – Apply Analytical Hierarchy Process (AHP)
Relative importance of factors:
- Solar radiation = 3× grid proximity  
- Solar radiation = 5× population proximity  
- Grid proximity = 2× population proximity  

**Final AHP weights:**
| Criterion | Weight |
|------------|---------|
| Solar radiation | 0.65 |
| Grid proximity | 0.23 |
| Population proximity | 0.12 |

Each raster cell (20×20 km) was scored from 1 to 5 and weighted to produce a **final suitability index**.

### Step 4 – Spatial Interpolation and Scoring
- Missing solar data were filled using **Inverse Distance Weighting (IDW)** interpolation.  
- Suitability map produced for all non-excluded areas.  
- Top-ranked cells representing **127 TWh** of production (≈60.4 GW capacity) were selected as optimal sites.

---

## Results

### Optimal Solar Sites
Nine 12 km² sites were identified as the most suitable for new solar power plants, offering:
- High solar irradiance
- Proximity to existing grid
- Close access to population centers  
These locations collectively satisfy Indonesia’s 2040 solar capacity target (≈60.4 GW).

### Suitability Overview
- 65 of 1,226 analyzed cells scored the **maximum suitability value (5)**.  
- Large areas across Indonesia show **high solar potential**, suggesting significant untapped renewable energy capacity.

---

## Financial Analysis

| Metric | Value | Description |
|--------|--------|-------------|
| Installed capacity | 60.4 GW | New solar required by 2040 |
| Annual generation | 127 TWh | Based on 0.24 capacity factor |
| Annual revenue | £9.7 billion | Assuming £0.0767/kWh |
| CAPEX | £56.2 billion | Installation + grid connection |
| NPV (5% discount) | £81 billion | Positive → profitable |
| LCOE | £0.03/kWh | Half of the assumed selling price |

Even under an 8% discount rate, LCOE remains below £0.04/kWh.  
The NPV turns negative only at 17%, confirming strong investment potential.

---

## Key Insights

- Indonesia has abundant high-quality solar sites, especially in central and eastern regions.  
- AHP-GIS integration effectively captures spatial and economic trade-offs.  
- Financial returns remain attractive under conservative assumptions.  
- IDW interpolation efficiently handled missing solar radiation data.  

---

## Future Work

- Incorporate additional suitability factors (land cost, road proximity, terrain slope, transmission infrastructure).  
- Use finer spatial resolution for higher precision in site identification.  
- Add operational expenditure (OPEX) and uncertainty analysis to financial modelling.  
- Explore Fuzzy AHP or **hybrid MCDM methods for sensitivity assessment.


