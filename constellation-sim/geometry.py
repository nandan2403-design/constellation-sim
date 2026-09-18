"""Steps 2 to 7: build the constellation, move it, find when each satellite sees each target."""
import numpy as np
from skyfield.api import load
import config as C


def time_grid():
    """Times from 0 to 24 h every 10 s (8,641 values)."""
    return np.arange(0, C.DURATION_S + C.STEP_S, C.STEP_S, dtype=float)


def walker_elements():
    """Step 2. RAAN and starting position (argument of latitude) for all 24 satellites."""
    per_plane = C.TOTAL_SATS // C.PLANES
    raan, u0 = [], []
    for p in range(C.PLANES):
        for j in range(per_plane):
            raan.append(360.0 / C.PLANES * p)                         # 0, 120, 240
            u0.append(360.0 / per_plane * j                            # 45 deg apart
                      + 360.0 * C.PHASING / C.TOTAL_SATS * p)          # 15 deg per plane
    return np.radians(raan), np.radians(u0)


def orbit_constants():
    a = C.EARTH_RADIUS_KM + C.ALTITUDE_KM
    n = np.sqrt(C.MU_KM3_S2 / a**3)          # mean motion, rad/s
    return a, n, 2 * np.pi / n               # radius, mean motion, period


def propagate_eci(t):
    """Step 3. Exact circular two-body positions, no J2. Returns (24, len(t), 3) in km."""
    a, n, _ = orbit_constants()
    raan, u0 = walker_elements()
    inc = np.radians(C.INCLINATION_DEG)
    u = u0[:, None] + n * t[None, :]
    cO, sO = np.cos(raan)[:, None], np.sin(raan)[:, None]
    x = a * (cO * np.cos(u) - sO * np.sin(u) * np.cos(inc))
    y = a * (sO * np.cos(u) + cO * np.sin(u) * np.cos(inc))
    z = a * np.sin(u) * np.sin(inc)
    return np.stack([x, y, z], axis=-1)


def earth_angle(t):
    """Greenwich angle: value at the epoch (from Skyfield) plus sidereal rotation."""
    ts = load.timescale()
    theta0 = np.radians(ts.utc(*C.EPOCH_UTC).gmst * 15.0)   # gmst is in hours
    return theta0 + C.EARTH_ROTATION_RAD_S * t


def eci_to_ecef(r_eci, t):
    """Step 4. Rotate about the z-axis so positions are fixed to the spinning Earth."""
    th = earth_angle(t)
    c, s = np.cos(th), np.sin(th)
    x = c * r_eci[..., 0] + s * r_eci[..., 1]
    y = -s * r_eci[..., 0] + c * r_eci[..., 1]
    return np.stack([x, y, r_eci[..., 2]], axis=-1)


def targets():
    """Step 5. 25 target points on a spherical Earth. Returns (25, 3) km, lats, lons."""
    lat, lon = np.meshgrid(C.TARGET_LATS, C.TARGET_LONS, indexing="ij")
    lat, lon = lat.ravel(), lon.ravel()
    la, lo = np.radians(lat), np.radians(lon)
    r = C.EARTH_RADIUS_KM * np.stack(
        [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    return r, lat, lon


def elevation_deg(sat_ecef, target_ecef):
    """Step 6. Elevation of every satellite, at every time, seen from one target."""
    rho = sat_ecef - target_ecef                              # line of sight
    up = target_ecef / np.linalg.norm(target_ecef)            # local vertical
    sin_el = (rho @ up) / np.linalg.norm(rho, axis=-1)
    return np.degrees(np.arcsin(np.clip(sin_el, -1, 1)))


def _crossing(t0, t1, e0, e1, m):
    """Time where elevation crosses m, by straight-line interpolation."""
    return t0 + (m - e0) / (e1 - e0) * (t1 - t0)


def find_windows(t, el, m=C.MIN_ELEVATION_DEG):
    """Turn an elevation series into a list of (start, end, peak) visibility windows."""
    above = el > m
    if not above.any():
        return []
    d = np.diff(above.astype(int))
    starts = list(np.where(d == 1)[0] + 1)
    ends = list(np.where(d == -1)[0] + 1)
    if above[0]:
        starts.insert(0, 0)
    if above[-1]:
        ends.append(len(t))
    out = []
    for s, e in zip(starts, ends):
        ts_ = t[0] if s == 0 else _crossing(t[s - 1], t[s], el[s - 1], el[s], m)
        te_ = t[-1] if e == len(t) else _crossing(t[e - 1], t[e], el[e - 1], el[e], m)
        out.append((float(ts_), float(te_), float(el[s:e].max())))
    return out


def build_pass_table():
    """Everything in Phase 1-2. Returns time grid, satellite ECEF, targets, pass table."""
    t = time_grid()
    sat = eci_to_ecef(propagate_eci(t), t)
    tgt, lat, lon = targets()
    table = []                                   # one dict per (sat, target, window)
    for g in range(len(tgt)):
        el = elevation_deg(sat, tgt[g])          # (24, 8641)
        for s in range(C.TOTAL_SATS):
            for (a, b, pk) in find_windows(t, el[s]):
                table.append(dict(sat=s, tgt=g, start=a, end=b, peak=pk))
    _add_regional_pass_ids(table)
    return t, sat, (tgt, lat, lon), table


def _add_regional_pass_ids(table):
    """'One target per pass': merge each satellite's windows over ANY target into one pass."""
    next_id = 0
    for s in range(C.TOTAL_SATS):
        rows = sorted((r for r in table if r["sat"] == s), key=lambda r: r["start"])
        cur_end = -1.0
        for r in rows:
            if r["start"] > cur_end:           # a new regional pass begins
                next_id += 1
            cur_end = max(cur_end, r["end"])
            r["rp"] = next_id


def revisit_stats(table, n_targets):
    """Step 7. Merge all satellites' windows per target, measure the gaps between them."""
    per_target = []
    for g in range(n_targets):
        iv = sorted((r["start"], r["end"]) for r in table if r["tgt"] == g)
        merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        # gap = end of one covered period to start of the next (edges of the day excluded)
        gaps = [merged[i + 1][0] - merged[i][1] for i in range(len(merged) - 1)]
        covered = sum(b - a for a, b in merged) / C.DURATION_S
        per_target.append(dict(gaps=np.array(gaps), covered=covered))
    return per_target


def sub_satellite_points(sat_ecef):
    """Latitude and longitude directly under each satellite, for ground-track plots."""
    x, y, z = sat_ecef[..., 0], sat_ecef[..., 1], sat_ecef[..., 2]
    lat = np.degrees(np.arcsin(z / np.linalg.norm(sat_ecef, axis=-1)))
    lon = np.degrees(np.arctan2(y, x))
    return lat, lon
