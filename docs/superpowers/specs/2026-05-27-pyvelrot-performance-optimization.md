# PyVelrot Performance Optimization

**Date:** 2026-05-27  
**Scope:** `pyvelrot.py` in-place; no API changes; same numerical outputs.

---

## Diagnosis (profiled, 9.6 s total)

| Bottleneck | Time | Cause |
|---|---|---|
| `_output_frame` → `_check_cp` | 6.2 s (65%) | 3,681 Python calls each scanning all other sites; ~2.8 M `np.linalg.norm` on 3-element arrays |
| `_read_fund_sites` eq_dist loop | 1.2 s (13%) | 429×3,252 = 1.4 M Python-loop `np.linalg.norm` calls |
| `rotate_geod`/`xyz_to_geod` | 0.1 s (1%) | 8,325 individual scalar calls with 5-iter Bowring each |
| `_increment_norm` triple loop | negligible | 178× 2×7×7 = 17,248 Python ops |

Root cause: calling NumPy on O(1)-sized arrays millions of times from Python.
NumPy per-call overhead (~1 µs) dominates the nanosecond computation.

---

## Optimizations (Option B — Thorough)

### O1 — Add batch coordinate utilities
New module-level functions (`_geod_to_xyz_batch`, `_xyz_to_geod_batch`,
`_rotation_matrix_neu_batch`) process N sites in one NumPy pass.
Scalar originals are unchanged (public API).

### O2 — Vectorize `read_vel_file`
Parse all data rows into Python lists, then convert to ECEF and rotate to XYZ
in one batch pass. Assembles `coords` and `covs` with pre-allocated arrays.

### O3 — Vectorize `_frame_update`
Replace per-site Python loop with `coords[:, :, 1] += np.cross(rot_vec, coords[:, :, 0])`.

### O4 — Vectorize pairwise distance in `_read_fund_sites`
Replace O(N1×N2) Python loop with BLAS-based distance matrix:
`dist2 = ‖p1‖² + ‖p2‖² - 2 p1 p2ᵀ`, then `np.where(dist2 < eq_dist²)`.
Eliminates 1.4 M scalar `np.linalg.norm` calls. Same result as before.

### O5 — Vectorize `_increment_norm`
Replace triple loop `for i,j,k` with matrix products:
`norm_eq += (A*w).T @ A`,  `bvec += A.T @ (d*w)`.

### O6 — Vectorize `_update_tran`
Apply transformation to all N1 sites at once:
- Translation: `sys1_coord[:,:,1] += trans_parm[:3]`
- Rotation: `+= np.cross(omega, pos)` for all sites
- Scale (7-param): `+= scale * pos`

Build batched A_neu via the rotation-block formula (-Skew(pos)) and einsum
for the covariance propagation. Preserves the intentional Fortran N↔E swap.

### O7 — Vectorize `_output_frame` (check_cp + coordinate transforms)
Precompute sys1/sys2 proximity flags as boolean arrays using vectorized
distance matrix (same formula as O4). Replaces 3,681 Python `_check_cp` scans.
Batch-compute NEU velocities via `np.einsum('nij,nj->ni', R_batch, vel_batch)`.

---

## Invariants preserved
- All floating-point operations identical (same BLAS ops, same Bowring iterations).
- Fortran N↔E swap in `_update_tran` preserved exactly.
- Output format strings unchanged.
- Public API (`run`, `get_parameters`, `get_sys1_velocities`) unchanged.
- Scalar utility functions unchanged.

## Test definition
Compare `output/reilinger_2006_igb14_python.tmp` (optimized run) against
the baseline file (current Python output), ignoring the timestamp line.
Acceptance: zero differences outside the timestamp.
