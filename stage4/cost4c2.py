"""Cost per run of chem4c2 at screening fidelity on 512 Sobol points (with the paired p0 = 0
control for L1_loose runs, as in the map). python cost4c2.py [out.json]"""
import json, os, sys, time
import torch
from scipy.stats import qmc
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c2 as chem  # noqa: E402
import ckpt  # noqa: E402
dev = torch.device("xpu")
B = 512
u = torch.tensor(qmc.Sobol(chem.DIM, scramble=True, seed=77).random(B), dtype=torch.float32, device=dev)
g = torch.Generator(device=dev); g.manual_seed(61)
torch.xpu.synchronize(); t = time.time()
o = chem.simulate(u, g, "screen", record=True)
torch.xpu.synchronize(); t_main = time.time() - t
t = time.time()
strict, a, c = chem.paired_control(u, 61, "screen", o)
torch.xpu.synchronize(); t_ctl = time.time() - t
res = {"batch": B, "seconds_main": t_main, "seconds_control": t_ctl, "ms_per_run": 1000 * (t_main + t_ctl) / B,
       "E1": float(o["E1"].float().mean()), "E2": float(o["E2"].float().mean()),
       "L1_loose": float(o["L1_loose"].float().mean()), "L1_strict": float(strict.float().mean()),
       "reached_end": float(o["reached_end"].float().mean()), "steps": int(o["steps"][0])}
print(res)
ckpt.atomic_json(res, sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "cost4c2.json"))
