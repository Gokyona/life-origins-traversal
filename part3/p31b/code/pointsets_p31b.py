"""Point sets for P3-1b (read-only sources in lot-s4): the 38 points of the 4c-5 point study (map4c5_5a/config.json;
10 coexistence groups: engine12538, p1558 x8, p3709 x8, C_engine10653 x8, rep1_1926 x8, 9c x5), 54-D 4c-4 coordinates
whose first 42 are the 4c-3 dimensions; and the 4c-3 box2 scan points (validate_p31._box2_points)."""
import json, os
import numpy as np

S4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lot-s4")


def coex_points():
    c = json.load(open(os.path.join(S4, "map4c5_5a", "config.json")))
    names = list(c["points"])
    X = np.array([c["points"][k] for k in names], dtype=np.float32)
    return names, X


def group_of(name):
    for g in ("p1558", "p3709", "C_engine10653", "rep1_1926"):
        if name.startswith(g):
            return g
    return name
