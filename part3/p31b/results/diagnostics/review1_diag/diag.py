import sys, math, json, importlib
sys.path.insert(1, "/run/media/bern/Data2T/myClaude/lot-p3/p31b")
import torch
dev = torch.device("cpu")
BASE = dict(T=22.85, pH=7.5, Mg=0.1, salinity=0.25, s_exp=0.0, P=100.0, feed=10.0, tau=1e12, k_c=1e-12, p0=0.0,
            depth=0.0, thermo=1.0, mineral_type=0.0, x_nuc=1.0, k_UV=0.0, avail_time=1e8, L_t=12.0, act_a=100.0,
            T_p=22.85, tau_p=1 / 60, f_p=1e-3, delta=1.1, k_on=3.16e6, m_Mg=1.0, E_ext=50.0, eps_0=1e-3, eps_1=1e-3,
            t_half_act=3.0, E_act_h=75.0, k_h_fixed=0.0, init_factor=0.0)
BASE.update(E_ext=18.4, act_a=0.0, K_I=4.0)
modname = sys.argv[1]
chem = importlib.import_module(modname)
def gen(s):
    g = torch.Generator(device=dev); g.manual_seed(s); return g
km = chem.PAIRS["GC"]["km_kmax"]
Ms = 780.0; c = Ms / 2000; r = km[1][1] * c / (km[1][0] + c)
fid = (64, 64, 1); nt = 32; B = 128
cases = {}
kf = 100 * r; q = r / (r + kf); lam = r * q; dt = 0.02 / lam
ov = dict(BASE); ov.update({"clamp_mstar": Ms, "koff_fixed": kf, "t_end": dt, "init_templates": {"n": nt, "L": 12, "seq": 0}})
cases["2a"] = ("GC", chem.DIM_P31, ov)
ov = {"T": 25.0, "T_m8": 25.0, "k_c": 1.0, "k_t": 1e5, "thermo": 1.0, "depth": 0.0, "mineral_type": 0.0,
      "x_nuc": 1.0, "init_factor": 0.0, "k_h_fixed": 0.0, "k_UV": 0.0, "tau": 1e12, "eps": 1e-3, "p0": 0.0,
      "avail_time": 1e8, "P": 100.0, "feed": 10.0, "L_t": 12.0, "pH": 7.5, "clamp_mstar": 1000.0,
      "t_end": 0.02 / 1e5, "init_templates": {"n": nt, "L": 12, "seq": 0}}
cases["2c"] = (None, chem.DIM, ov)
res = {}
for name, (pair, dim, ov) in cases.items():
    for multi in (False, True):
        chem.PAIR = pair; chem.MULTI_START = multi; chem.DBG.clear()
        u = torch.full((B, dim), 0.5, dtype=torch.float32, device=dev)
        o = chem.simulate(u, gen(9000 + multi), fid, override=dict(ov))
        d = dict(chem.DBG); d["copies"] = float(o["copies"].sum()); d["dropped"] = float(o["copies_dropped"].sum())
        e = d["eligible"]
        d = {k: round(v / e, 6) if k not in ("eligible", "room_sum") else v for k, v in d.items()}
        res[f"{name}_{'on' if multi else 'off'}"] = d
        print(name, "on" if multi else "off", d, flush=True)
json.dump(res, open(f"diag_{modname}.json", "w"), indent=1)
