"""W = 1 identity: gk_chem (multi-word engine, 1 word) against chem4c4 (md5 b79360f6), the Stage-4 engine of paper 2, on CPU.
chem4c4 is in the paper repository (https://github.com/Gokyona/life-origins-traversal, code/lot-s4/).
    python tests/identity_w1.py <folder containing chem4c4.py>      (or set GK_CHEM4C4_DIR)"""
import json, os, sys, time
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
LOT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("GK_CHEM4C4_DIR")
if not LOT:
    sys.exit("give the folder that contains chem4c4.py (paper repository, code/lot-s4/) as argument or GK_CHEM4C4_DIR")
sys.path.insert(0, os.path.join(HERE, "..", "gokyona_env", "engine"))
import gk_chem
sys.path.insert(0, LOT)
import chem4c4
gk_chem.configure(1)
pre = json.load(open(os.path.join(HERE, "..", "gokyona_env", "data", "presets.json")))
torch.set_num_threads(4)
res = {}
for name in ("engine12538", "9c_engine15233"):
    u = torch.tensor([pre[name]["u"]] * 2, dtype=torch.float32)
    for seed in (1, 2):
        out = []
        for mod in (chem4c4, gk_chem):
            g = torch.Generator(); g.manual_seed(seed)
            t = time.time()
            ov = {"return_state": True}
            if mod is gk_chem:
                ov["ration"] = "slot"          # chem4c4 order explicitly (the platform default is random)
            o = mod.simulate(u, g, (4, 128, 120), override=ov)
            out.append((o, time.time() - t))
        (a, ta), (b, tb) = out
        bad = [k for k, v in a.items() if isinstance(v, torch.Tensor) and k in b and isinstance(b[k], torch.Tensor)
               and v.shape == b[k].shape and not torch.equal(torch.nan_to_num(v, nan=-7.0), torch.nan_to_num(b[k], nan=-7.0))]
        st = torch.equal(a["state"]["L"], b["state"]["L"]) and torch.equal(a["state"]["seq"], b["state"]["seq"][..., 0]) \
            and torch.equal(a["state"]["chem"], b["state"]["chem"])
        res[f"{name}_s{seed}"] = {"mismatch": bad, "state_equal": bool(st), "copies": float(a["copies"].sum()),
                                  "vesicles": float(a.get("vesicles_formed", torch.zeros(1)).sum()), "t_old": round(ta, 1), "t_new": round(tb, 1)}
        print(name, seed, res[f"{name}_s{seed}"], flush=True)
print("ALL IDENTICAL" if all(not r["mismatch"] and r["state_equal"] for r in res.values()) else "NOT IDENTICAL")
