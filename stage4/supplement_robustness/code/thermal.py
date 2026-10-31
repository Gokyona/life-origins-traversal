"""Thermal guard for long stage-4 jobs on the Arc B70 (read-only: never changes fans or limits).

Source: the xe driver's hwmon in sysfs (/sys/class/drm/card*/device/hwmon/hwmon*/ with
name "xe"): temp*_input with labels pkg, vram, mctrl, pcie, vram_ch_*; each has its own
temp*_crit (pkg 100 C, vram 105 C on this card). Power is derived from energy*_input
(microjoules) between readings.

Rule: every >= 60 s (checked between batches) log time, max temperature, pkg temperature and
power to gpu_temp.log. If any sensor is at or above its own crit - 10 C (90 C if a sensor has
no crit), the caller's checkpoint function is run and the job pauses, re-reading every 30 s,
until the hottest sensor is below 75 C. Every pause is logged. The CPU package temperature
(coretemp "Package id 0") is checked the same way, pausing at 90 C.

GPU lock: the first Guard of a process takes an exclusive flock on lot-s4/.gpu.lock and holds
it until the process exits; a second GPU job exits at once (code 3). Several XPU jobs at once
oversubscribe the 32 GB of VRAM; the xe driver then moves buffers into system RAM (shmem, not
counted in any RSS or cgroup limit), which ran the machine out of memory on 2026-10-26.
"""
import fcntl, glob, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
RESUME_BELOW = 75.0
DEFAULT_LIMIT = 90.0
MARGIN = 10.0
CPU_LIMIT = 90.0
_LOCK = None


# Part III copy (P3-1): the lock is the Stage 4 lock file in lot-s4/, so Stage 4 and Part III jobs
# exclude each other (one GPU job at a time on the B70); the temperature log stays in this folder.
GPU_LOCK = os.environ.get("LOT_GPU_LOCK", os.path.join(HERE, "..", "..", "lot-s4", ".gpu.lock"))


def gpu_lock(path=GPU_LOCK):
    """Exclusive per-process GPU lock (see module docstring); exits with code 3 if taken."""
    global _LOCK
    if _LOCK is not None:
        return
    f = open(path, "a+")
    try:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        f.seek(0)
        sys.stderr.write(f"GPU lock {path} held by another job ({f.read().strip()}); exiting\n")
        sys.exit(3)
    f.seek(0); f.truncate(); f.write(f"pid {os.getpid()} {' '.join(sys.argv)}\n"); f.flush()
    _LOCK = f


def find_cpu_pkg():
    for f in sorted(glob.glob("/sys/class/hwmon/hwmon*/temp*_label")):
        if _read(f) == "Package id 0" and _read(os.path.join(os.path.dirname(f), "name")) == "coretemp":
            return f[:-len("_label")] + "_input"
    return None


def _read(path):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return None


def find_hwmon():
    for h in sorted(glob.glob("/sys/class/drm/card*/device/hwmon/hwmon*")):
        if _read(os.path.join(h, "name")) == "xe":
            return h
    return None


class Guard:
    def __init__(self, log_path=os.path.join(HERE, "gpu_temp.log"), every=60.0, job=""):
        gpu_lock()
        self.h = find_hwmon()
        self.log_path, self.every, self.job = log_path, every, job
        self.last = 0.0
        self.e_last = None
        self.temps = []
        self.pauses = 0
        self.sensors = []
        if self.h:
            for f in sorted(glob.glob(os.path.join(self.h, "temp*_input"))):
                base = f[:-len("_input")]
                crit = _read(base + "_crit")
                label = _read(base + "_label") or os.path.basename(base)
                limit = (float(crit) / 1000 - MARGIN) if crit else DEFAULT_LIMIT
                self.sensors.append((f, label, limit))
        cpu = find_cpu_pkg()
        if cpu:
            self.sensors.append((cpu, "cpu_pkg", CPU_LIMIT))
        src = f"hwmon {self.h}" if self.h else "no xe hwmon found (guard inactive)"
        limits = ", ".join(f"{l} pause>={lim:.0f}" for _, l, lim in self.sensors if l in ("pkg", "vram", "pcie", "mctrl", "cpu_pkg"))
        self._log(f"start job={job} source={src}; {limits}; vram_ch_* pause>=95 C; resume below {RESUME_BELOW:.0f} C")

    def _log(self, msg):
        line = time.strftime("%Y-%m-%d %H:%M:%S") + " " + msg
        with open(self.log_path, "a") as f:
            f.write(line + "\n")
            f.flush()
            os.fsync(f.fileno())

    def read(self):
        vals = []
        for f, label, limit in self.sensors:
            v = _read(f)
            if v is not None:
                vals.append((label, float(v) / 1000, limit))
        e = _read(os.path.join(self.h, "energy1_input")) if self.h else None
        now = time.time()
        power = None
        if e is not None:
            if self.e_last is not None and now > self.e_last[0]:
                power = (float(e) - self.e_last[1]) / 1e6 / (now - self.e_last[0])
            self.e_last = (now, float(e))
        return vals, power

    def check(self, checkpoint=None, force=False):
        """Call between batches. Logs every `every` seconds; pauses if too hot."""
        if not self.sensors or (not force and time.time() - self.last < self.every):
            return
        self.last = time.time()
        vals, power = self.read()
        tmax = max(v for _, v, _ in vals)
        pkg = next((v for l, v, _ in vals if l == "pkg"), float("nan"))
        cpu = next((v for l, v, _ in vals if l == "cpu_pkg"), float("nan"))
        self.temps.append(tmax)
        pw = f"{power:.0f} W" if power is not None else "n/a"
        self._log(f"job={self.job} max {tmax:.0f} C pkg {pkg:.0f} C cpu {cpu:.0f} C power {pw}")
        hot = [(l, v, lim) for l, v, lim in vals if v >= lim]
        if not hot:
            return
        self.pauses += 1
        self._log(f"PAUSE {self.pauses} job={self.job}: " + ", ".join(f"{l} {v:.0f}>={lim:.0f}" for l, v, lim in hot))
        if checkpoint:
            checkpoint()
        t0 = time.time()
        while True:
            time.sleep(30)
            vals, _ = self.read()
            tmax = max(v for _, v, _ in vals)
            if tmax < RESUME_BELOW:
                break
        self._log(f"RESUME job={self.job} after {time.time() - t0:.0f} s, max {tmax:.0f} C")

    def summary(self):
        t = self.temps
        return {"min_C": min(t) if t else None, "mean_C": sum(t) / len(t) if t else None,
                "max_C": max(t) if t else None, "readings": len(t), "pauses": self.pauses,
                "source": self.h}
