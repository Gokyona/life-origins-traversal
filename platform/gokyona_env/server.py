#!/usr/bin/env python3
"""Gokyona Environment - local web server (stdlib only, Python 3.10+).

Serves the browser UI (gokyona_env/web/) and a small JSON API over <data-dir>/runs/<name>/
(see CONTRACT.md). Binds 127.0.0.1:9300 by default. Never starts a simulation by itself:
the engine is launched only when the user presses Start in the UI (POST /api/start).

    gokyona-env [--data-dir DIR] [--host 127.0.0.1] [--port 9300] [--mem-max 24G] [--no-browser]
"""
import argparse
import datetime
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import threading
import time
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

PKG_DIR = pathlib.Path(__file__).resolve().parent
PKG_PARENT = PKG_DIR.parent          # cwd for the engine, so `-m gokyona_env.engine.run` also works un-installed
WEB_DIR = PKG_DIR / "web"
DATA_FILES = PKG_DIR / "data"
NAME_RE = re.compile(r"^[A-Za-z0-9_-]{1,40}$")
CONTROL_CMDS = ("pause", "resume", "stop")
ALIVE_STATES = ("running", "starting", "paused")
MAX_BODY = 1 << 20
MIME = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
        ".css": "text/css; charset=utf-8", ".json": "application/json", ".svg": "image/svg+xml",
        ".png": "image/png", ".ico": "image/x-icon", ".txt": "text/plain; charset=utf-8"}


class Cfg:
    """Runtime settings (filled from argparse in main())."""
    data_dir = pathlib.Path.home() / ".gokyona-env"
    host = "127.0.0.1"
    port = 9300
    mem_max = "24G"

    @classmethod
    def runs_dir(cls):
        return cls.data_dir / "runs"


_start_lock = threading.Lock()


class ApiError(Exception):
    def __init__(self, msg, code=400):
        super().__init__(msg)
        self.code = code


# ---------------------------------------------------------------- helpers
def now_iso():
    return datetime.datetime.now().replace(microsecond=0).isoformat()


def load_json(path, retries=3):
    """Read a JSON file (None if missing); retry briefly in case a writer is mid-replace."""
    for i in range(retries):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except FileNotFoundError:
            return None
        except (json.JSONDecodeError, OSError):
            if i == retries - 1:
                raise
            time.sleep(0.05)
    return None


def write_json_atomic(path, obj):
    path = pathlib.Path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")
    os.replace(tmp, path)


def run_dir(name, must_exist=True):
    if not isinstance(name, str) or not NAME_RE.match(name):
        raise ApiError("invalid run name (allowed: A-Z a-z 0-9 _ -, 1-40 characters)")
    runs = Cfg.runs_dir().resolve()
    d = (runs / name).resolve()
    if d.parent != runs:
        raise ApiError("invalid run name")
    if must_exist and not d.is_dir():
        raise ApiError("no such run: " + name, 404)
    return d


def presets():
    return load_json(DATA_FILES / "presets.json") or {}


def params():
    return load_json(DATA_FILES / "params.json") or []


def num(v, what, lo=None, hi=None, integer=False):
    if isinstance(v, bool):
        raise ApiError(f"{what}: not a number")
    if not isinstance(v, (int, float)):
        try:
            v = float(v)
        except (TypeError, ValueError):
            raise ApiError(f"{what}: not a number")
    if not math.isfinite(v):
        raise ApiError(f"{what}: not a finite number")
    if integer:
        if float(v) != int(v):
            raise ApiError(f"{what}: must be an integer")
        v = int(v)
    if lo is not None and v < lo:
        raise ApiError(f"{what}: must be >= {lo}")
    if hi is not None and v > hi:
        raise ApiError(f"{what}: must be <= {hi}")
    return v


# ---------------------------------------------------------------- one-GPU-job check
def systemd_units():
    """Active gk-*/lot-* user units, or None if systemctl is unavailable/failing."""
    if not shutil.which("systemctl"):
        return None
    cmd = ["systemctl", "--user", "list-units", "--no-pager", "--plain", "--no-legend", "gk-*", "lot-*"]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    except Exception:  # noqa: BLE001
        return None
    if out.returncode != 0:
        return None
    units = []
    for line in out.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[2] in ("active", "activating", "reloading", "deactivating"):
            units.append({"unit": parts[0], "active": parts[2], "sub": parts[3]})
    return units


def alive_runs():
    """Runs whose status.json says running|starting|paused."""
    out = []
    runs = Cfg.runs_dir()
    if not runs.is_dir():
        return out
    for d in sorted(runs.iterdir()):
        p = d / "status.json"
        if not p.is_file():
            continue
        try:
            st = load_json(p) or {}
        except Exception:  # noqa: BLE001
            continue
        if isinstance(st, dict) and st.get("state") in ALIVE_STATES:
            age = time.time() - p.stat().st_mtime
            out.append({"run": d.name, "state": st.get("state"), "status_age_s": round(age, 1)})
    return out


def start_check(name=None):
    blockers, notes = [], []
    units = systemd_units()
    if units is None:
        notes.append("systemctl not available: unit check skipped")
    else:
        for u in units:
            base = u["unit"].rsplit(".", 1)[0]
            if base.startswith("gk-"):
                blockers.append(f"{u['unit']} is active (another Gokyona run)")
            elif base.startswith("lot-") and base != "lot-tempwatch":
                blockers.append(f"{u['unit']} is active (another GPU job)")
    alive = alive_runs()
    for a in alive:
        blockers.append(f"run '{a['run']}' has status '{a['state']}' (status.json updated {a['status_age_s']:.0f} s ago;"
                        f" if that engine is dead, fix {Cfg.runs_dir() / a['run'] / 'status.json'})")
    if name:
        try:
            d = run_dir(name)
            if not (d / "config.json").is_file():
                blockers.append(f"run '{name}' has no config.json (press Create first)")
        except ApiError as e:
            blockers.append(str(e))
    ok = not blockers
    reason = ("allowed: no other GPU job active" if ok else "refused: " + "; ".join(blockers))
    if notes:
        reason += " (" + "; ".join(notes) + ")"
    return {"ok": ok, "reason": reason, "units": units or [], "alive_runs": alive,
            "systemd": units is not None, "systemd_run": bool(shutil.which("systemd-run"))}


def _reap(proc):
    try:
        proc.wait()
    except Exception:  # noqa: BLE001
        pass


def start_run(name):
    d = run_dir(name)
    with _start_lock:
        chk = start_check(name)
        if not chk["ok"]:
            raise ApiError(chk["reason"], 409)
        engine = [sys.executable, "-u", "-m", "gokyona_env.engine.run", "--run-dir", str(d)]
        if shutil.which("systemd-run") and chk["systemd"]:
            cmd = ["systemd-run", "--user", "--scope", f"--unit=gk-{name}",
                   "-p", f"MemoryMax={Cfg.mem_max}", "-p", "MemorySwapMax=2G"] + engine
        else:
            cmd = engine
        log_path = d / "engine.log"
        env = dict(os.environ)
        env["PYTHONPATH"] = str(PKG_PARENT) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        with open(log_path, "ab") as log:
            log.write(f"\n=== {now_iso()} server: {' '.join(cmd)}\n".encode())
            log.flush()
            proc = subprocess.Popen(cmd, cwd=str(PKG_PARENT), env=env, stdin=subprocess.DEVNULL, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True, close_fds=True)
        threading.Thread(target=_reap, args=(proc,), daemon=True).start()
        time.sleep(1.0)
        rc = proc.poll()
        if rc is not None and rc != 0:
            tail = ""
            try:
                tail = log_path.read_text(errors="replace")[-800:]
            except OSError:
                pass
            raise ApiError(f"engine exited immediately (rc={rc}). engine.log tail:\n{tail}", 500)
        return {"ok": True, "run": name, "pid": proc.pid, "cmd": " ".join(cmd), "log": str(log_path)}


# ---------------------------------------------------------------- create
def create_run(body):
    name = body.get("name")
    d = run_dir(name, must_exist=False)
    pr = presets()
    preset = body.get("preset")
    if preset not in pr:
        raise ApiError(f"unknown preset: {preset!r}")
    if (d / "status.json").exists() or (d / "engine.log").exists():
        raise ApiError(f"run '{name}' already exists and has been started; choose a new name", 409)

    pdefs = {p["name"]: p for p in params()}
    env = dict(pr[preset].get("env") or {})          # full preset env; user overrides below
    for k, v in (body.get("env") or {}).items():
        if k not in pdefs and k not in env:
            raise ApiError(f"env: unknown variable {k!r}")
        v = num(v, f"env.{k}")
        p = pdefs.get(k, {})
        if p.get("min") is not None:
            v = max(p["min"], v)
        if p.get("max") is not None:
            v = min(p["max"], v)
        if p.get("step") == 1:
            v = float(round(v))
        env[k] = v

    trait = dict(pr[preset].get("trait") or {})
    trait.setdefault("enabled", False)
    for k, v in (body.get("trait") or {}).items():
        if k in ("enabled", "fidelity_link", "epistatic", "pair_score"):
            if not isinstance(v, bool):
                raise ApiError(f"trait.{k}: must be true/false")
            trait[k] = v
        elif k == "mode":
            if v not in ("self", "shared", "membrane", "stability", "division", "retention"):
                raise ApiError("trait.mode: must be one of self, shared, membrane, stability, division, retention")
            trait[k] = v
        elif k == "s":
            trait[k] = float(num(v, "trait.s", 0.0, 4.0))
        elif k == "landscape":
            if v not in ("domain", "additive-mean"):
                raise ApiError("trait.landscape: must be 'domain' or 'additive-mean'")
            trait[k] = v
        elif k == "s_f":
            trait[k] = float(num(v, "trait.s_f", 0.0, 10.0))
        elif k == "L_rep":
            trait[k] = num(v, "trait.L_rep", 1, 8192, integer=True)
        elif k == "seed":
            trait[k] = num(v, "trait.seed", 0, 2**31 - 1, integer=True)
        elif k == "k":
            pass  # read-only, always 6
        else:
            raise ApiError(f"trait: unknown field {k!r}")
    trait["k"] = 6
    trait.setdefault("landscape", "domain")
    trait.setdefault("epistatic", False)
    trait.setdefault("pair_score", True)

    cfg = {
        "name": name,
        "preset": preset,
        "env": env,
        "trait": trait,
        "lmax": num(body.get("lmax", 1024), "lmax", 8, 16384, integer=True),
        "multi_start": bool(body.get("multi_start", False)),
        "ration": body.get("ration", "random") if body.get("ration", "random") in ("random", "slot", "vesicle") else "random",
        "census_years": float(num(body.get("census_years", 100), "census_years", 1, 1e6)),
        "vesicle_slots": num(body.get("vesicle_slots", 8), "vesicle_slots", 1, 127, integer=True),
        "migration_m": float(num(body.get("migration_m", 0), "migration_m", 0, 1)),
        "migration_mode": body.get("migration_mode", "random") if body.get("migration_mode", "random") in ("random", "ring") else "random",
        "half_pool_mix": bool(body.get("half_pool_mix", False)),
        "dt_max_h": float(num(body.get("dt_max_h", 0), "dt_max_h", 0, 1e6)),
        "sites": num(body.get("sites", 64), "sites", 2, 4096, integer=True),
        "slots": num(body.get("slots", 512), "slots", 1, 1 << 20, integer=True),
        "vram_budget_gb": 20,
        "report_every_s": 5,
        "checkpoint_every_s": 240,
        "created": now_iso(),
    }
    d.mkdir(parents=True, exist_ok=True)
    write_json_atomic(d / "config.json", cfg)
    return {"ok": True, "run": name, "config": cfg}


# ---------------------------------------------------------------- read endpoints
def list_runs():
    out = []
    runs = Cfg.runs_dir()
    if not runs.is_dir():
        return out
    for d in sorted(runs.iterdir()):
        if not d.is_dir() or not NAME_RE.match(d.name):
            continue
        state = None
        try:
            st = load_json(d / "status.json")
            if isinstance(st, dict):
                state = st.get("state")
            elif (d / "config.json").is_file():
                state = "created"
        except Exception:  # noqa: BLE001
            state = "unreadable"
        try:
            mtime = max([p.stat().st_mtime for p in d.iterdir()] + [d.stat().st_mtime])
        except OSError:
            mtime = 0
        out.append({"name": d.name, "state": state, "mtime": mtime})
    return out


def get_status(name):
    d = run_dir(name)
    st = load_json(d / "status.json")
    if st is None:
        if (d / "config.json").is_file():
            return {"name": name, "state": "created"}
        raise ApiError("run has no status yet", 404)
    return st


def get_history(name, since):
    d = run_dir(name)
    p = d / "history.jsonl"
    if not p.is_file():
        return []
    with open(p, "rb") as f:
        data = f.read()
    lines = data.split(b"\n")[:-1]  # drop "" or an incomplete last line still being written
    out = []
    for raw in lines[since:]:
        try:
            out.append(json.loads(raw))
        except (json.JSONDecodeError, UnicodeDecodeError):
            out.append(None)  # keep indices aligned with line numbers
    return out


def get_hardware():
    try:
        hw = load_json(Cfg.data_dir / "hardware.json")
    except Exception:  # noqa: BLE001
        hw = None
    return hw if isinstance(hw, dict) else {}


# ---------------------------------------------------------------- HTTP
def _sanitize(o):
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, dict):
        return {k: _sanitize(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_sanitize(v) for v in o]
    return o


class Handler(BaseHTTPRequestHandler):
    server_version = "Gokyona/0.1"

    def log_message(self, fmt, *args):
        if os.environ.get("GK_VERBOSE"):
            super().log_message(fmt, *args)

    def _send(self, code, body, ctype="application/json"):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code, obj):
        try:
            body = json.dumps(obj, allow_nan=False, default=str).encode()
        except ValueError:  # NaN/Infinity read from a file -> null
            body = json.dumps(_sanitize(obj), default=str).encode()
        self._send(code, body)

    def _body(self):
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_BODY:
            raise ApiError("request too large", 413)
        raw = self.rfile.read(n) if n else b"{}"
        try:
            obj = json.loads(raw or b"{}")
        except json.JSONDecodeError:
            raise ApiError("body is not valid JSON")
        if not isinstance(obj, dict):
            raise ApiError("body must be a JSON object")
        return obj

    def _host_ok(self):
        # DNS-rebinding guard when bound to loopback.
        if Cfg.host not in ("127.0.0.1", "localhost", "::1"):
            return True
        host = (self.headers.get("Host") or "").rsplit(":", 1)[0].strip("[]").lower()
        return host in ("127.0.0.1", "localhost", "::1", "")

    def do_GET(self):
        self._dispatch("GET")

    def do_POST(self):
        self._dispatch("POST")

    def _dispatch(self, method):
        try:
            if not self._host_ok():
                raise ApiError("bad Host header", 403)
            path = urlparse(self.path).path
            if path.startswith("/api/"):
                return self._json(200, self._api(method, path))
            if method != "GET":
                raise ApiError("method not allowed", 405)
            return self._static(path)
        except ApiError as e:
            self._json(e.code, {"ok": False, "error": str(e)})
        except Exception as e:  # noqa: BLE001
            self._json(500, {"ok": False, "error": f"server error: {type(e).__name__}: {e}"})

    def _api(self, method, path):
        q = {k: v[0] for k, v in parse_qs(urlparse(self.path).query).items()}
        if method == "GET":
            if path == "/api/runs":
                return list_runs()
            if path == "/api/status":
                return get_status(q.get("run"))
            if path == "/api/history":
                try:
                    since = max(0, int(q.get("since", 0)))
                except ValueError:
                    raise ApiError("since must be an integer")
                return get_history(q.get("run"), since)
            if path == "/api/presets":
                return presets()
            if path == "/api/params":
                return params()
            if path == "/api/check":
                return start_check(q.get("run") or None)
            if path == "/api/hardware":
                return get_hardware()
        elif method == "POST":
            origin = self.headers.get("Origin")
            if origin and urlparse(origin).hostname not in ("127.0.0.1", "localhost", "::1", Cfg.host):
                raise ApiError("cross-origin request refused", 403)
            body = self._body()
            if path == "/api/create":
                return create_run(body)
            if path == "/api/start":
                return start_run(body.get("run"))
            if path == "/api/control":
                d = run_dir(body.get("run"))
                cmd = body.get("cmd")
                if cmd not in CONTROL_CMDS:
                    raise ApiError("cmd must be one of " + ", ".join(CONTROL_CMDS))
                ctl = {"cmd": cmd, "ts": now_iso()}
                write_json_atomic(d / "control.json", ctl)
                return {"ok": True, "run": d.name, "control": ctl}
        raise ApiError(f"unknown endpoint {method} {path}", 404)

    def _static(self, path):
        if path in ("", "/"):
            path = "/index.html"
        web = WEB_DIR.resolve()
        target = (web / path.lstrip("/")).resolve()
        if web not in target.parents or not target.is_file():
            raise ApiError("not found", 404)
        self._send(200, target.read_bytes(), MIME.get(target.suffix, "application/octet-stream"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gokyona-env", description="Gokyona Environment local web UI")
    ap.add_argument("--data-dir", default=str(pathlib.Path.home() / ".gokyona-env"),
                    help="where runs/<name>/ and hardware.json live (default ~/.gokyona-env)")
    ap.add_argument("--host", default="127.0.0.1", help="bind address (default 127.0.0.1)")
    ap.add_argument("--port", type=int, default=9300, help="port (default 9300)")
    ap.add_argument("--mem-max", default="24G", help="MemoryMax for the systemd scope, if systemd is used (default 24G)")
    ap.add_argument("--no-browser", action="store_true", help="do not open a web browser")
    a = ap.parse_args(argv)
    if not re.match(r"^[0-9]+[KMGT]?$", a.mem_max):
        ap.error("--mem-max must look like 24G")

    Cfg.data_dir = pathlib.Path(a.data_dir).expanduser().resolve()
    Cfg.host, Cfg.port, Cfg.mem_max = a.host, a.port, a.mem_max
    Cfg.runs_dir().mkdir(parents=True, exist_ok=True)

    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    srv.daemon_threads = True
    shown = "localhost" if a.host in ("127.0.0.1", "::1") else a.host
    url = f"http://{shown}:{a.port}/"
    print(f"Gokyona Environment on {url}  (bound to {a.host}; data in {Cfg.data_dir}; Ctrl+C to stop)", flush=True)
    if not a.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()


if __name__ == "__main__":
    main()
