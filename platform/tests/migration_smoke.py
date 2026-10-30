"""Migration option checks (CPU): every move conserves chains and letters of the run, keeps duplexes symmetric,
keeps compartment labels consistent with alive vesicles, and never crosses between the halves.
    python tests/migration_smoke.py [m] [mode] [half_pool_mix 0/1]"""
import json, os, sys
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gokyona_env", "engine"))
import gk_chem
from trait import Trait
m = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
mode = sys.argv[2] if len(sys.argv) > 2 else "random"
hm = bool(int(sys.argv[3])) if len(sys.argv) > 3 else True
gk_chem.configure(4, 64)
N, S = 8, 128
pre = json.load(open(os.path.join(HERE, "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
orig = gk_chem._migrate
bad = {"letters": 0, "chains": 0, "pair": 0, "label": 0, "cross_half": 0, "checks": 0}


def checked(VS, L, seq, chem, pair, dlen, lab, busy, B, N_, S_, m_, mode_, mgen, mig, t_h):
    n0, l0, k0 = int((L > 0).sum()), int(L.sum()), len(mig["log"])
    orig(VS, L, seq, chem, pair, dlen, lab, busy, B, N_, S_, m_, mode_, mgen, mig, t_h)
    bad["checks"] += 1
    bad["chains"] += int((L > 0).sum()) != n0
    bad["letters"] += int(L.sum()) != l0
    R = L.shape[0]
    rr, cc = (pair >= 0).nonzero(as_tuple=True)
    bad["pair"] += int((pair[rr, pair[rr, cc]] != cc).sum())
    inside = lab > 0
    rr, cc = inside.nonzero(as_tuple=True)
    bad["label"] += int((~VS.alive[rr, lab[rr, cc].long() - 1]).sum())
    for (_, r, _, d, _, _) in mig["log"][k0:]:
        bad["cross_half"] += ((r % N_) >= N_ // 2) != ((d % N_) >= N_ // 2)


gk_chem._migrate = checked


class Stop:
    k = 0
    def due(self): return True
    def __call__(self, loc):
        Stop.k += 1
        if Stop.k >= 120:
            print("migration", {k: v for k, v in loc["mig"].items() if k != "log"}, "log entries", len(loc["mig"]["log"]))
            return "stop"


u = torch.tensor([pre["u"]], dtype=torch.float32)
ov = {k: float(v) for k, v in pre["env"].items() if k != "avail_time"}
ov.update({"endless": True, "migration": m, "migration_mode": mode, "half_pool_mix": hm})
g = torch.Generator(); g.manual_seed(1)
tr = Trait(dict(pre["trait"]), N, torch.device("cpu"))
gk_chem.simulate(u, g, (N, S, 300), override=ov, trait=tr, hook=Stop())
print("violations", bad, "PASS" if all(v == 0 for k, v in bad.items() if k != "checks") and bad["checks"] > 0 else "FAIL")
