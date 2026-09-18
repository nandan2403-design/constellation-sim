"""Run everything: python run_all.py
Creates figures/coverage_map.png, figures/revisit_histogram.png,
figures/completion_vs_failure.png and results/results.md
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")                    # save images without opening windows
import matplotlib.pyplot as plt
import config as C
import geometry as G
import tasking as T

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
t0 = time.time()

# ---------------- Phase 1-2: geometry
print("1/4  Propagating orbits and finding access windows...")
t, sat_ecef, (tgt, lat, lon), table = G.build_pass_table()
a, n, period = G.orbit_constants()
print(f"     orbit radius {a:.3f} km, period {period/60:.2f} min, {len(table)} windows")

# ---------------- Step 7: revisit
print("2/4  Revisit analysis...")
rev = G.revisit_stats(table, len(tgt))
mean_gap = np.array([r["gaps"].mean() / 60 if len(r["gaps"]) else np.nan for r in rev])
max_gap = np.array([r["gaps"].max() / 60 if len(r["gaps"]) else np.nan for r in rev])
all_gaps = np.concatenate([r["gaps"] for r in rev]) / 60

# Plot 1: coverage map
fig, ax = plt.subplots(figsize=(8, 6))
slat, slon = G.sub_satellite_points(sat_ecef)
for s in range(C.TOTAL_SATS):                           # all 24 satellites, full day
    lo, la = slon[s].copy(), slat[s].copy()
    lo[np.abs(np.diff(lo, prepend=lo[0])) > 180] = np.nan   # break line at the date line
    ax.plot(lo, la, lw=0.5, color="grey", alpha=0.6,
            label="Ground tracks (all 24 satellites)" if s == 0 else None)
sc = ax.scatter(lon, lat, c=mean_gap, s=160, cmap="viridis_r", edgecolor="k", zorder=3)
fig.colorbar(sc, ax=ax, label="Mean revisit gap (min)")
ax.set(xlim=(85, 135), ylim=(-20, 30), xlabel="Longitude (deg E)", ylabel="Latitude (deg)",
       title="Walker 24/3/1 coverage over Southeast Asia (24 h)")
ax.legend(loc="lower left", fontsize=8)
ax.grid(alpha=0.3)
fig.savefig("figures/coverage_map.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# Plot 2: revisit gap histogram
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(all_gaps, bins=45, color="steelblue", edgecolor="white")
ax.set_yscale("log")   # most gaps are short; log scale keeps the long outages visible
ax.axvline(all_gaps.mean(), color="orange", ls="--", label=f"Mean {all_gaps.mean():.0f} min")
ax.axvline(all_gaps.max(), color="red", ls="--", label=f"Max {all_gaps.max():.0f} min")
ax.axvline(np.median(all_gaps), color="green", ls=":", label=f"Median {np.median(all_gaps):.0f} min")
ax.set(xlabel="Revisit gap (minutes)", ylabel="Count, log scale (all 25 targets)",
       title="Revisit gaps, Walker 24/3/1, 10 deg elevation mask")
ax.legend()
fig.savefig("figures/revisit_histogram.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------- Steps 8-12: tasking and failures
print(f"3/4  Tasking: {len(C.FAILURE_FRACTIONS)} failure levels x {C.SEEDS} seeds x 2 schemes...")
windows = T.index_windows(table)
res = {"central": {}, "auction": {}}
for f in C.FAILURE_FRACTIONS:
    rows = {"central": [], "auction": []}
    for seed in range(C.SEEDS):
        reqs, fail_time, rng = T.make_scenario(len(tgt), f, seed)
        rows["central"].append(T.run_centralised(windows, reqs, fail_time))
        rows["auction"].append(T.run_auction(windows, reqs, fail_time, rng))
    for k in rows:
        arr = np.array(rows[k])                                 # columns: rate, latency
        res[k][f] = dict(rate=arr[:, 0].mean(),
                         ci=1.96 * arr[:, 0].std(ddof=1) / np.sqrt(len(arr)),
                         lat=np.nanmean(arr[:, 1]))
    print(f"     {int(f*100):>2}% failed: centralised {res['central'][f]['rate']:.1%}, "
          f"auction {res['auction'][f]['rate']:.1%}")

# Plot 3: completion vs failure
fig, ax = plt.subplots(figsize=(8, 5))
x = [f * 100 for f in C.FAILURE_FRACTIONS]
for k, name in [("central", "Centralised greedy"), ("auction", "Decentralised auction")]:
    ax.errorbar(x, [res[k][f]["rate"] * 100 for f in C.FAILURE_FRACTIONS],
                yerr=[res[k][f]["ci"] * 100 for f in C.FAILURE_FRACTIONS],
                marker="o", capsize=4, label=name)
ax.set(xlabel="Satellites disabled (%)", ylabel="Task completion rate (%)",
       title=f"Completion vs failure ({C.SEEDS} seeds, 95% CI, "
             f"ground detection delay {C.DETECTION_DELAY_MIN} min)")
ax.legend()
ax.grid(alpha=0.3)
fig.savefig("figures/completion_vs_failure.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# Bonus plot 4: how slow detection hurts centralised (30% failure)
print("     Bonus: sweeping ground detection delay at 30% failure...")
worst = C.FAILURE_FRACTIONS[-1]
saved_delay = C.DETECTION_DELAY_MIN
sweep = {"central": [], "auction": []}
for d in C.DETECTION_SWEEP_MIN:
    C.DETECTION_DELAY_MIN = d
    c, a_ = [], []
    for seed in range(C.SEEDS):
        reqs, fail_time, rng = T.make_scenario(len(tgt), worst, seed)
        c.append(T.run_centralised(windows, reqs, fail_time)[0])
        a_.append(T.run_auction(windows, reqs, fail_time, rng)[0])
    sweep["central"].append(np.mean(c) * 100)
    sweep["auction"].append(np.mean(a_) * 100)
C.DETECTION_DELAY_MIN = saved_delay
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(C.DETECTION_SWEEP_MIN, sweep["central"], "o-", label="Centralised greedy")
ax.plot(C.DETECTION_SWEEP_MIN, sweep["auction"], "o-", label="Decentralised auction")
ax.set(xlabel="Time for ground to notice a failure (min)", ylabel="Task completion rate (%)",
       title=f"Effect of failure detection delay ({worst:.0%} of satellites disabled)")
ax.legend()
ax.grid(alpha=0.3)
fig.savefig("figures/bonus_detection_delay.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------- Results table
print("4/4  Writing results/results.md")
lo_f, hi_f = C.FAILURE_FRACTIONS[0], C.FAILURE_FRACTIONS[-1]
fmt = lambda k, f: f"{res[k][f]['rate']:.1%} ± {res[k][f]['ci']:.1%}"
lines = [
    "| Metric | Value |", "|---|---|",
    f"| Orbital period | {period/60:.2f} min |",
    f"| Mean revisit gap (all targets) | {all_gaps.mean():.1f} min |",
    f"| Median revisit gap | {np.median(all_gaps):.1f} min |",
    f"| Max revisit gap (worst target) | {all_gaps.max():.1f} min |",
    f"| Share of day each target is in view | {np.mean([r['covered'] for r in rev]):.0%} |",
    f"| Completion, centralised, {lo_f:.0%} failure | {fmt('central', lo_f)} |",
    f"| Completion, auction, {lo_f:.0%} failure | {fmt('auction', lo_f)} |",
    f"| Completion, centralised, {hi_f:.0%} failure | {fmt('central', hi_f)} |",
    f"| Completion, auction, {hi_f:.0%} failure | {fmt('auction', hi_f)} |",
    f"| Mean latency, centralised, {lo_f:.0%} / {hi_f:.0%} | "
    f"{res['central'][lo_f]['lat']:.0f} / {res['central'][hi_f]['lat']:.0f} min |",
    f"| Mean latency, auction, {lo_f:.0%} / {hi_f:.0%} | "
    f"{res['auction'][lo_f]['lat']:.0f} / {res['auction'][hi_f]['lat']:.0f} min |",
]
with open("results/results.md", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nDone in {time.time()-t0:.0f} s. Open the 'figures' folder to see the plots.")
