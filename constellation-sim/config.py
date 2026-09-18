"""All the settings for the simulator live here. Change a number, re-run run_all.py."""

# ---- Earth and physics (spherical Earth) ----
EARTH_RADIUS_KM = 6378.137
MU_KM3_S2 = 398600.4418            # Earth's gravitational parameter
EARTH_ROTATION_RAD_S = 7.2921159e-5  # sidereal rate, NOT 2*pi/86400

# ---- Constellation: Walker Delta 24/3/1 ----
TOTAL_SATS = 24
PLANES = 3
PHASING = 1
ALTITUDE_KM = 550.0
INCLINATION_DEG = 53.0

# ---- Simulation time ----
STEP_S = 10
DURATION_S = 24 * 3600
EPOCH_UTC = (2026, 1, 1, 0, 0, 0)  # start time, used for Earth's rotation angle

# ---- Targets: 5 x 5 grid over Southeast Asia ----
TARGET_LATS = [-10.0, -2.5, 5.0, 12.5, 20.0]
TARGET_LONS = [95.0, 102.5, 110.0, 117.5, 125.0]
MIN_ELEVATION_DEG = 10.0

# ---- Imaging requests ----
REQUESTS_PER_HOUR = 6.0       # tuned so the constellation is busy but not overloaded
REQUEST_WINDOW_H = 18         # requests arrive between 0 h and 18 h
DEADLINE_H = 3                # each request must be imaged within 3 h of arriving

# ---- Decentralised auction ----
MSG_DELAY_MIN_S = 5
MSG_DELAY_MAX_S = 30
AUCTION_TIMEOUT_S = 60        # >= worst case (25 s spread + 30 s bid delay)

# ---- Failures ----
FAILURE_FRACTIONS = [0.0, 0.1, 0.2, 0.3]
FAILURE_WINDOW_H = 18         # satellites fail at a random time in the first 18 h
DETECTION_DELAY_MIN = 180        # how long the ground takes to notice a dead satellite
DETECTION_SWEEP_MIN = [0, 90, 180, 360]   # bonus plot: vary the delay at 30% failure
SEEDS = 30                    # repeat each experiment 30 times with different randomness
