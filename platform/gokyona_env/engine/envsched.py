"""Environment schedules of Gokyona Environment (29 Oct 2026): w(t) in [0, 1] for a named environmental variable
(T, pH or salinity). Used by the trait (environment-dependent landscape: which of two landscapes is selected) and,
with couple_chemistry, by the chemistry (the variable itself moves between v0 at w = 0 and v1 at w = 1).
  constant    w = 0
  season      w = 0.5 - 0.5 cos(2 pi t / period), period default 1 simulated year (day/night is below the step
              resolution and averaged, like wet-dry cycling)
  trend       triangle 0 -> 1 -> 0, period default 200 simulated years
  catastrophe w = 0 (events are drawn by the engine; the trait may set w = 1 per site for shift_years after one)
Impacts (meteor / comet; option): rare GLOBAL events drawn by the engine; after one, w moves to pulse_w for
pulse_years, then relaxes back with time constant relax_years (0 = at once); an optional feed pulse multiplies the
feedstock by feed_mult for feed_years."""
import math

YEAR_H = 8766.0
VARIABLES = ("T", "pH", "salinity")


def impact_pulse(icfg, t_h, last_h):
    """Weight of the post-impact environment pulse in [0, 1]: 1 for pulse_years after the last impact, then (if
    relax_years > 0) exp(-(t - last - pulse) / relax), else 0."""
    if not icfg or last_h is None or last_h < -1e29:
        return 0.0
    age = (t_h - last_h) / YEAR_H
    if age < 0:
        return 0.0
    pul = float(icfg.get("pulse_years", 0) or 0)
    if age < pul:
        return 1.0
    rel = float(icfg.get("relax_years", 0) or 0)
    return math.exp(-(age - pul) / rel) if rel > 0 else 0.0


def w_total(cfg, t_h, icfg=None, last_h=None):
    """Schedule w(t), moved towards the impact pulse value pulse_w (default 1) by the pulse weight."""
    w = w_at(cfg, t_h)
    f = impact_pulse(icfg, t_h, last_h)
    return w + f * (float((icfg or {}).get("pulse_w", 1.0)) - w) if f > 0 else w


def w_at(cfg, t_h):
    sch = (cfg or {}).get("schedule", "constant")
    yr = t_h / YEAR_H
    if sch == "season":
        return 0.5 - 0.5 * math.cos(2 * math.pi * yr / float(cfg.get("period_years", 1.0)))
    if sch == "trend":
        ph = (yr / float(cfg.get("period_years", 200.0))) % 1.0
        return 2 * ph if ph < 0.5 else 2 - 2 * ph
    if sch in ("constant", "catastrophe"):
        return 0.0
    raise ValueError("unknown schedule " + str(sch))
