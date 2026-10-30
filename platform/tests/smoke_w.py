"""Multi-word smoke: W words (L_MAX = 64 W) on CPU; reports lengths and checks sequence invariants
(no letters set beyond each chain's length)."""
import json, os, sys, time
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gokyona_env", "engine"))
import gk_chem, words as wd
Wn = int(sys.argv[1]) if len(sys.argv) > 1 else 4
gk_chem.configure(Wn)
pre = json.load(open(os.path.join(HERE, "..", "gokyona_env", "data", "presets.json")))
torch.set_num_threads(4)
u = torch.tensor([pre["engine12538"]["u"]] * 2, dtype=torch.float32)
g = torch.Generator(); g.manual_seed(1)
t = time.time()
o = gk_chem.simulate(u, g, (4, 128, 120), override={"return_state": True})
L, seq = o["state"]["L"], o["state"]["seq"]
beyond = (seq & ~wd.mask(L)).ne(0).any(-1)
print("W", Wn, "L_MAX", gk_chem.L_MAX, "sec", round(time.time() - t, 1), "max L", int(L.max()), "mean L", float(L[L > 0].float().mean()),
      "chains>64", int((L > 64).sum()), "copies", float(o["copies"].sum()), "copy_len", o["copy_mean_len"].tolist(),
      "bits beyond length", int(beyond.sum()), "vesicles", float(o["vesicles_formed"].sum()))
