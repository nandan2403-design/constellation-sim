"""Sanity checks. Run: python tests.py  (every line should say PASS)"""
import numpy as np
import config as C
import geometry as G

def check(name, ok):
    print(("PASS  " if ok else "FAIL  ") + name)

t = G.time_grid()
a, n, period = G.orbit_constants()
r = G.propagate_eci(t)
check("8,641 time steps", len(t) == 8641)
check("orbit radius constant at 6,928.137 km", np.allclose(np.linalg.norm(r, axis=-1), a))
check("period about 95.7 min", abs(period / 60 - 95.65) < 0.1)
rp = G.propagate_eci(np.array([0.0, period]))
check("satellite returns to start after one period", np.allclose(rp[:, 0], rp[:, 1], atol=1e-6))
lat = np.degrees(np.arcsin(r[..., 2] / a))
check("maximum latitude is 53 deg", abs(lat.max() - 53) < 0.05)
raan, u0 = G.walker_elements()
check("planes 120 deg apart", np.allclose(np.degrees(raan[[0, 8, 16]]), [0, 120, 240]))
check("45 deg in-plane spacing, 15 deg plane offset",
      np.isclose(np.degrees(u0[1] - u0[0]), 45) and np.isclose(np.degrees(u0[8] - u0[0]), 15))
tgt = np.array([C.EARTH_RADIUS_KM, 0, 0])
check("satellite straight overhead has 90 deg elevation",
      np.isclose(G.elevation_deg(np.array([[a, 0, 0]]), tgt)[0], 90))
th = G.earth_angle(np.array([0.0, 3600.0]))
check("Earth turns ~15.04 deg per hour", abs(np.degrees(th[1] - th[0]) - 15.041) < 0.001)
_, _, _, table = G.build_pass_table()
longest = max(w["end"] - w["start"] for w in table) / 60
check(f"longest pass about 8 min (got {longest:.2f})", 7.5 < longest < 8.6)

import tasking as T
C.MSG_DELAY_MIN_S = C.MSG_DELAY_MAX_S = 0.001
C.AUCTION_TIMEOUT_S = 0.01
w = T.index_windows(table)
same = True
for s in range(5):
    c = T.run_centralised(w, *T.make_scenario(25, 0.0, s)[:2])
    a = T.run_auction(w, *T.make_scenario(25, 0.0, s))
    same &= c[0] == a[0] and abs(c[1] - a[1]) < 0.1   # latency within 6 s

check("auction with zero delay gives exactly the centralised result", same)
