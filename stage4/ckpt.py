"""Crash-safe saving for long stage-4 runs (power loss, reboot, kill).

  atomic_json(obj, path)   write to path.tmp, fsync, rename over path, fsync the folder:
                           after a power cut the file is either the old or the new
                           version, never half-written
  atomic_npz(path, **arr)  the same for numpy arrays (e.g. cached ground truth)
  load_json(path)          the saved result, or None if there is none

Pattern for a script: save after every finished unit of work (one rep, one map), and on
start with --resume load the saved result and skip the units already in it.
"""
import json, os
import numpy as np


def _fsync_dir(path):
    fd = os.open(os.path.dirname(os.path.abspath(path)) or ".", os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    _fsync_dir(path)


def atomic_npz(path, **arrays):
    tmp = path + ".tmp.npz"
    with open(tmp, "wb") as f:
        np.savez(f, **arrays)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    _fsync_dir(path)


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)
