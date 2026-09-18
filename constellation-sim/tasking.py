"""Steps 8 to 12: imaging requests, two schedulers, satellite failures.

Model assumptions (also listed in the README):
- Capacity: a satellite can take ONE image per regional pass.
- Centralised: the ground sees everything instantly, but only learns a satellite has died
  DETECTION_DELAY_MIN after it happens. Tasks on a dead satellite are re-planned then.
- Decentralised: satellites talk over crosslinks with 5-30 s delays. A dead satellite
  stops bidding. If a winner dies, the losers notice no 'done' message and re-auction.
"""
import heapq
import itertools
from collections import defaultdict
import numpy as np
import config as C


# ---------------------------------------------------------------- scenario
def make_scenario(n_targets, fail_fraction, seed):
    """Same seed = same requests and same failures, so both schemes face identical conditions."""
    rng = np.random.default_rng(seed)
    horizon = C.REQUEST_WINDOW_H * 3600
    reqs, t = [], rng.exponential(3600 / C.REQUESTS_PER_HOUR)
    while t < horizon:
        reqs.append(dict(id=len(reqs), arrive=t, tgt=int(rng.integers(n_targets)),
                         deadline=t + C.DEADLINE_H * 3600))
        t += rng.exponential(3600 / C.REQUESTS_PER_HOUR)
    n_fail = int(round(fail_fraction * C.TOTAL_SATS))           # 10% -> 2, 20% -> 5, 30% -> 7
    fail_time = np.full(C.TOTAL_SATS, np.inf)
    dead = rng.choice(C.TOTAL_SATS, size=n_fail, replace=False)
    fail_time[dead] = rng.uniform(0, C.FAILURE_WINDOW_H * 3600, size=n_fail)
    delays = np.random.default_rng(seed + 10_000)                # message delays
    return reqs, fail_time, delays


def index_windows(table):
    """windows[(sat, tgt)] -> list of (start, end, regional_pass_id), sorted by start."""
    w = defaultdict(list)
    for r in table:
        w[(r["sat"], r["tgt"])].append((r["start"], r["end"], r["rp"]))
    for k in w:
        w[k].sort()
    return w


def earliest_free_window(windows, sat, tgt, ready, deadline, busy):
    """Earliest time `sat` can image `tgt`, no earlier than `ready`, using a free pass."""
    for start, end, rp in windows.get((sat, tgt), ()):
        img = max(start, ready)
        if img > deadline:
            break
        if img < end and rp not in busy:
            return img, rp
    return None


def score(reqs, successes):
    """Completion rate and mean latency (minutes) over completed requests."""
    lat = [successes[r["id"]] - r["arrive"] for r in reqs if r["id"] in successes]
    rate = len(lat) / len(reqs) if reqs else 0.0
    return rate, (np.mean(lat) / 60 if lat else np.nan)


# ---------------------------------------------------------------- centralised greedy
def run_centralised(windows, reqs, fail_time):
    detect = fail_time + C.DETECTION_DELAY_MIN * 60
    busy = set()                                # regional passes already used
    plan = defaultdict(list)                    # sat -> [(image_time, req)]
    successes = {}
    events = [(r["arrive"], 0, r["id"]) for r in reqs]           # (time, kind, data)
    events += [(detect[s], 1, s) for s in range(C.TOTAL_SATS) if np.isfinite(detect[s])]
    heapq.heapify(events)
    known_dead = set()

    def assign(req, now):
        best = None
        for s in range(C.TOTAL_SATS):
            if s in known_dead:
                continue
            opt = earliest_free_window(windows, s, req["tgt"], now, req["deadline"], busy)
            if opt and (best is None or opt[0] < best[0]):
                best = (opt[0], opt[1], s)
        if best:
            img, rp, s = best
            busy.add(rp)
            plan[s].append((img, req))

    while events:
        now, kind, data = heapq.heappop(events)
        if kind == 0:                                           # a request arrives
            assign(reqs[data], now)
        else:                                                   # ground notices sat died
            s = data
            known_dead.add(s)
            lost = [(img, rq) for img, rq in plan[s] if img >= fail_time[s]]
            plan[s] = [(img, rq) for img, rq in plan[s] if img < fail_time[s]]
            for _, rq in lost:                                  # re-plan the lost tasks
                assign(rq, now)

    for s, items in plan.items():                               # which images really happened
        for img, rq in items:
            if img < fail_time[s]:
                successes[rq["id"]] = min(img, successes.get(rq["id"], np.inf))
    return score(reqs, successes)


# ---------------------------------------------------------------- decentralised auction
def run_auction(windows, reqs, fail_time, rng):
    counter = itertools.count()                 # tie-breaker so the heap never compares dicts
    events = []

    def push(time, kind, **data):
        heapq.heappush(events, (time, next(counter), kind, data))

    def delay():
        return rng.uniform(C.MSG_DELAY_MIN_S, C.MSG_DELAY_MAX_S)

    alive = lambda s, now: fail_time[s] > now
    busy = [set() for _ in range(C.TOTAL_SATS)]         # each satellite's OWN used/locked passes
    my_bid = {}                                         # (sat, req, round) -> (bid, rp, img)
    heard = defaultdict(dict)                           # (sat, req, round) -> {other: bid}
    done_heard = defaultdict(set)                       # sat -> requests known finished
    successes = {}

    for r in reqs:                                      # ground broadcasts every request
        for s in range(C.TOTAL_SATS):
            push(r["arrive"] + delay(), "request", sat=s, req=r["id"], rnd=0)

    while events:
        now, _, kind, d = heapq.heappop(events)
        s = d["sat"]
        if not alive(s, now):
            continue                                    # dead satellites do nothing
        req = reqs[d["req"]]
        key = (s, d["req"], d.get("rnd", 0))

        if kind == "request":                           # decide whether to bid
            if d["req"] in done_heard[s]:
                continue
            opt = earliest_free_window(windows, s, req["tgt"], now + C.AUCTION_TIMEOUT_S,
                                       req["deadline"], busy[s])
            if opt is None:
                continue                                # cannot help, stays silent
            img, rp = opt
            busy[s].add(rp)                             # lock this pass while bidding
            my_bid[key] = (-img, rp, img)               # earlier image = higher bid
            for o in range(C.TOTAL_SATS):
                if o != s:
                    push(now + delay(), "bid", sat=o, req=d["req"], rnd=d["rnd"],
                         frm=s, bid=-img)
            push(now + C.AUCTION_TIMEOUT_S, "timeout", sat=s, req=d["req"], rnd=d["rnd"])

        elif kind == "bid":                             # remember other satellites' bids
            heard[key][d["frm"]] = d["bid"]

        elif kind == "timeout":                         # pick the winner on my own clock
            bid, rp, img = my_bid[key]
            others = heard[key]
            best_other = max(((b, -o) for o, b in others.items()), default=None)
            i_win = best_other is None or (bid, -s) > best_other
            if i_win:
                push(img, "image", sat=s, req=d["req"])
            else:
                busy[s].discard(rp)                     # release my locked pass
                w_bid, _ = best_other
                # if no 'done' arrives by then, assume the winner failed and re-auction
                push(-w_bid + C.MSG_DELAY_MAX_S + 1, "check", sat=s, req=d["req"],
                     rnd=d["rnd"])

        elif kind == "image":                           # take the picture
            successes[d["req"]] = min(now, successes.get(d["req"], np.inf))
            for o in range(C.TOTAL_SATS):
                if o != s:
                    push(now + delay(), "done", sat=o, req=d["req"])

        elif kind == "done":
            done_heard[s].add(d["req"])

        elif kind == "check":
            if d["req"] not in done_heard[s] and now < req["deadline"]:
                push(now, "request", sat=s, req=d["req"], rnd=d["rnd"] + 1)

    return score(reqs, successes)
