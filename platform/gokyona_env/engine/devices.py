"""Device selection and memory/temperature probes (portable; the Linux Intel-xe hwmon parts are optional)."""
import glob, os
import torch


def pick(name="auto"):
    if name and name != "auto":
        return torch.device(name)
    if hasattr(torch, "xpu") and torch.xpu.is_available():
        return torch.device("xpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def describe(dev):
    try:
        if dev.type == "xpu":
            return torch.xpu.get_device_name(0)
        if dev.type == "cuda":
            return torch.cuda.get_device_name(0)
    except Exception:  # noqa: BLE001
        pass
    return dev.type


def _host_mem():
    """(free, total) bytes of system RAM."""
    try:
        info = {}
        for line in open("/proc/meminfo"):
            k, v = line.split(":", 1)
            info[k] = int(v.split()[0]) * 1024
        return info.get("MemAvailable", info.get("MemFree", 0)), info.get("MemTotal", 0)
    except OSError:
        try:
            tot = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
            return tot // 2, tot
        except (ValueError, OSError, AttributeError):
            return 8 << 30, 16 << 30


def mem_info(dev):
    """(free, total) bytes on the device (host RAM for cpu / mps)."""
    try:
        if dev.type == "xpu":
            return torch.xpu.mem_get_info()
        if dev.type == "cuda":
            return torch.cuda.mem_get_info()
    except Exception:  # noqa: BLE001
        pass
    return _host_mem()


def allocated(dev):
    try:
        if dev.type == "xpu":
            return torch.xpu.memory_allocated(), torch.xpu.max_memory_allocated()
        if dev.type == "cuda":
            return torch.cuda.memory_allocated(), torch.cuda.max_memory_allocated()
    except Exception:  # noqa: BLE001
        pass
    return None, None


def reset_peak(dev):
    try:
        if dev.type == "xpu":
            torch.xpu.empty_cache(); torch.xpu.reset_peak_memory_stats()
        elif dev.type == "cuda":
            torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    except Exception:  # noqa: BLE001
        pass


def sync(dev):
    try:
        if dev.type == "xpu":
            torch.xpu.synchronize()
        elif dev.type == "cuda":
            torch.cuda.synchronize()
    except Exception:  # noqa: BLE001
        pass


def gpu_temp():
    """GPU temperature in C from Linux hwmon (Intel xe, amdgpu) if available, else None."""
    for h in glob.glob("/sys/class/drm/card*/device/hwmon/hwmon*"):
        try:
            if open(os.path.join(h, "name")).read().strip() in ("xe", "i915", "amdgpu"):
                vals = []
                for f in glob.glob(os.path.join(h, "temp*_input")):
                    vals.append(int(open(f).read()) / 1000.0)
                if vals:
                    return max(vals)
        except OSError:
            continue
    try:
        import subprocess
        out = subprocess.run(["nvidia-smi", "--query-gpu=temperature.gpu", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=3).stdout.strip()
        return float(out.splitlines()[0]) if out else None
    except Exception:  # noqa: BLE001
        return None
