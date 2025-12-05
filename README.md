# Vegetation Cover Change Detection (1985-2025)

A Google Earth Engine script for detecting multi-decadal vegetation cover change using harmonized Landsat time series with scientifically validated thresholds.

## Methodology

### Approach

This script uses a **trend-first classification** approach: computing linear NDVI trends across the 40-year record before applying state-based thresholds. This enables detection of active succession at any vegetation stage.

### NDVI Thresholds

Classification thresholds are aligned with peer-reviewed literature:

| NDVI Range | Class | Scientific Basis |
|------------|-------|------------------|
| ≥ 0.60 | Dense Canopy | NDVI 0.6-1.0 = dense green vegetation (FAO, ResearchGate meta-analyses) |
| 0.40-0.60 | Transitional | Corresponds to FAO "open/fragmented forest" (10-40% canopy) |
| 0.20-0.40 | Sparse | Sparsely vegetated areas (Copernicus, ISU studies) |
| < 0.20 | Bare | Exposed soil, minimal vegetation |

**Sources:**
- FAO Forest Resources Assessment: 10% canopy = forest threshold, 40%+ = closed forest
- NDVI-canopy correlation: R² = 0.88-0.94 (UAV validation studies)
- Dense vegetation NDVI > 0.5; forest typically > 0.6 (ResearchGate synthesis)

### Trend Analysis

Linear regression is computed on summer NDVI (June-September) from 1985-2025:

```
NDVI(t) = slope × t + intercept
```

| Slope (per year) | Classification | Source |
|------------------|----------------|--------|
| > +0.005 | Gaining | MDPI Kunming study (2000-2020): moderate improvement |
| -0.005 to +0.005 | Stable | Same study: stable conditions range |
| < -0.005 | Losing | Indicates active decline |

**Validation:** The ±0.005/yr threshold represents ~0.2 NDVI change over 40 years, corresponding to a full vegetation class transition.

### Change Classification

Six change classes combine start/end state with trend direction:

| Class | Transition | Trend Requirement |
|-------|------------|-------------------|
| Canopy Loss | Dense → Sparse/Bare | — |
| Canopy Thinning | Dense → Transitional | Losing |
| Emerging Biomass | Sparse → Transitional | Gaining |
| Canopy Thickening | Transitional → Dense | Gaining |
| Canopy Densification | Dense → Dense | Gaining |
| Canopy Establishment | Sparse → Dense | 5-year epoch |

### Trajectory Projection

For gaining areas, linear extrapolation estimates years to reach the dense canopy threshold:

```
Years to threshold = (0.60 - current_NDVI) / slope
```

This provides a first-order approximation. Vegetation succession often follows sigmoid (logistic) growth, so actual trajectories may accelerate then plateau.

## Data Sources

- **Landsat 5 TM** (1984-2012): `LANDSAT/LT05/C02/T1_L2`
- **Landsat 7 ETM+** (1999-present): `LANDSAT/LE07/C02/T1_L2`
- **Landsat 8 OLI** (2013-present): `LANDSAT/LC08/C02/T1_L2`
- **Landsat 9 OLI-2** (2021-present): `LANDSAT/LC09/C02/T1_L2`

All imagery is Level-2 Surface Reflectance, Collection 2, Tier 1, with cloud/shadow masking via QA_PIXEL.

## Usage

1. Open in [Google Earth Engine Code Editor](https://code.earthengine.google.com/)
2. Define `roi` or use current map bounds
3. Run script
4. Click map to inspect individual pixels with NDVI time series chart

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DENSE_CANOPY` | 0.60 | NDVI threshold for closed canopy |
| `TRANSITIONAL` | 0.40 | NDVI threshold for woodland/shrub |
| `SPARSE` | 0.20 | NDVI threshold for open vegetation |
| `GAINING_SLOPE` | 0.005 | Slope threshold for active gain |
| `LOSING_SLOPE` | -0.005 | Slope threshold for active loss |

## Example Output

![Example output showing vegetation change detection](example_output_2025.jpg)

## Intended Applications

- Agricultural land abandonment and forest succession
- Post-disturbance vegetation recovery monitoring
- Landscape transition analysis in rural/peri-urban areas
- Temporal mapping of land-use change trajectories

## Citation

```
Pizzi, G. (2025). Vegetation Cover Change Detection Script.
GitHub repository. https://github.com/gbrlpzz/forest-cover-change
```

## References

1. FAO (2020). Global Forest Resources Assessment. Rome.
2. Kunming NDVI Study, MDPI Remote Sensing (2000-2020 analysis)
3. Copernicus Land Monitoring Service - CORINE Land Cover Technical Guide
4. Kennedy et al. (2010). Detecting trends in forest disturbance and recovery using yearly Landsat time series.

## License

Apache License 2.0. See LICENSE for details.

## Contact

Gabriele Pizzi | info@gabrielepizzi.com | gabrielepizzi.com
