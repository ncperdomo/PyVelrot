# PyVelrot

Python reimplementation of the GAMIT/GLOBK Fortran program **`velrot.f`**
(T. Herring et al.).

PyVelrot estimates a Helmert transformation (translations + rotations ± scale) between
two GNSS velocity fields using weighted least squares, then writes the transformed
System 1 velocities together with the System 2 sites.

---

## Features

- 3-parameter (T), 6-parameter (TR), or 7-parameter (TRS) Helmert transformation
- Weighted least squares with optional vertical down-weighting
- Flexible site-linking via a fundamental-site file (`eq_dist`, explicit pairs, exclusions)
- Output compatible with the GAMIT `.vel` file format
- Matches Fortran output to within display-format precision

---

## Repository layout

```
PyVelrot/
├── pyvelrot.py                     # Python module
├── requirements.txt
├── PyVelrot_validation.ipynb       # Validation notebook (no GAMIT needed)
├── data/
│   ├── reilinger_2006.vel          # System 1 — 429 published GNSS velocities
│   ├── serpelloni_2022.vel         # System 2 — 3252 European IGb14 velocities
│   └── velrot.lnk                  # Site-link file (eq_dist 1000)
└── reference/
    └── reilinger_2006_igb14.tmp    # Precomputed Fortran velrot output
```

---

## Installation

```bash
pip install numpy pandas matplotlib jupyter cartopy   # notebook dependencies
```

No compiled code; no GAMIT/GLOBK installation required.

---

## Usage

### Command line

```bash
python pyvelrot.py <sys1.vel> <sys1_frame> <sys2.vel> <sys2_frame> \
                   <output.tmp> <out_frame> <links.lnk> <height_weight> <param_opt>
```

**Example** (replicates the included test):

```bash
python pyvelrot.py data/reilinger_2006.vel NONE \
                   data/serpelloni_2022.vel NONE \
                   output/igb14.tmp NONE \
                   data/velrot.lnk 0 TR
```

### Python API

```python
from pyvelrot import PyVelrot

vr = PyVelrot()
vr.run(
    sys1_file="data/reilinger_2006.vel", sys1_frame="NONE",
    sys2_file="data/serpelloni_2022.vel", sys2_frame="NONE",
    out_file="output/igb14.tmp",          out_frame="NONE",
    fund_file="data/velrot.lnk",
    height_weight=0,    # horizontal only
    param_opt="TR",     # translations + rotations
)

# Access results
params = vr.get_parameters()
# Keys: Tx, Ty, Tz (mm/yr), Rx, Ry, Rz (mas/yr), Scale (ppb/yr),
#       sigma_*, rms_mm_yr, chi_fit, num_data, num_fund

vels = vr.get_sys1_velocities()
# List of dicts: name, lon, lat, ve, vn, vu (mm/yr), se, sn, su, rho_ne
```

### Arguments

| Argument | Description |
|----------|-------------|
| `sys1_file` | System 1 `.vel` file (to be transformed) |
| `sys1_frame` | Reference frame of System 1 (`NONE` if already in out_frame) |
| `sys2_file` | System 2 `.vel` file (the reference) |
| `sys2_frame` | Reference frame of System 2 |
| `out_file` | Output file path (or `"6"` to write to stdout) |
| `out_frame` | Output reference frame |
| `fund_file` | Fundamental-site file; links sites between sys1 and sys2 |
| `height_weight` | Weight for vertical component: 0 = horizontal only, 1 = equal |
| `param_opt` | `T` (3-param), `TR` (6-param, default), `TRS` (7-param), `L` (2-param local) |

### Link file format

```
eq_dist 1000          # auto-link all site pairs within 1000 m
SITE1_GPS SITE2_GPS   # explicit named pair
-BADSITE_GPS          # remove a site from the link list
```

---

## Validation

Run the included notebook from the `PyVelrot/` directory:

```bash
jupyter notebook PyVelrot_validation.ipynb
```

The notebook does **not** require GAMIT/GLOBK. It runs PyVelrot on the included test data
and compares the output against the precomputed Fortran reference file
`reference/reilinger_2006_igb14.tmp`.

### Test results

Test: 429 System 1 sites (Reilinger et al. 2006) + 3252 System 2 sites (Serpelloni et al. 2022),
178 linked pairs, 6-parameter TR transformation, horizontal-only weighting.

**Transformation parameters:**

| Parameter | Fortran | PyVelrot | Difference |
|-----------|---------|----------|------------|
| Tx (mm/yr) | +0.7551 ± 0.1766 | +0.7521 ± 0.1762 | −0.003 |
| Ty (mm/yr) | +0.0674 ± 0.2367 | +0.0667 ± 0.2365 | −0.001 |
| Tz (mm/yr) | +0.7002 ± 0.2621 | +0.7027 ± 0.2623 | +0.003 |
| Rx (mas/yr) | −0.0924 ± 0.0058 | −0.0924 ± 0.0058 | < 0.0001 |
| Ry (mas/yr) | −0.5290 ± 0.0077 | −0.5290 ± 0.0077 | 0.0000 |
| Rz (mas/yr) | +0.7490 ± 0.0085 | +0.7490 ± 0.0084 | < 0.0001 |
| RMS (mm/yr) | 0.72 | 0.72 | 0.00 |
| NRMS | 1.35 | 1.35 | 0.00 |

**Per-site velocities (429 sites):**

| Component | Max residual |
|-----------|-------------|
| East | ≤ 0.006 mm/yr |
| North | ≤ 0.006 mm/yr |
| Up | ≤ 0.058 mm/yr |

The horizontal residuals are dominated by display-precision rounding (the `.vel` format
stores velocities to 2 decimal places, ±0.005 mm/yr). The vertical residual (< 2% of
the 3 mm/yr vertical uncertainty) comes from floating-point differences between
`numpy.linalg.solve`/`inv` and Fortran's Gaussian-elimination solver `invert_vis`.

---

## Mapping to Fortran subroutines

| Python | Fortran | Purpose |
|--------|---------|---------|
| `read_vel_file()` | `read_sys` | Read `.vel` file → ECEF coords + NEU covariance |
| `_read_fund_sites()` | `read_fund` | Parse link file, build site-pair list |
| `_clear_norm()` | `clear_norm` | Initialise normal equations with parameter constraints |
| `_increment_norm()` | `increment_norm` | Accumulate Aᵀ W A and Aᵀ W d |
| `_transframe()` | `transframe` + `invert_vis` | Solve and invert normal equations |
| `get_parts()` | `get_parts` | Build 3×7 Helmert design matrix |
| `_update_tran()` | `update_tran` | Apply transformation + propagate uncertainty to sys1 |
| `_output_sum()` | `output_sum` | Write header, parameters, per-site residuals |
| `_output_frame()` | `output_frame` | Write transformed velocity field |

---

## `.vel` file format

```
* comment lines start with any character other than a space
 lon(°)    lat(°)    ve    vn    de    dn    se    sn    rho    vu    du    su   site
```

All velocities in mm/yr. `de/dn` = rate adjustment; `se/sn/su` = 1-sigma uncertainty
(mm/yr); `rho` = NE correlation coefficient; `vu/du/su` = vertical rate, adjustment, sigma.

---

## Reference

Reilinger, R., et al. (2006). GPS constraints on continental deformation in the
Africa‐Arabia‐Eurasia continental collision zone and implications for the dynamics of
plate interactions. *Journal of Geophysical Research*, 111(B5).
