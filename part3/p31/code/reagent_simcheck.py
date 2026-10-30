"""Cold-gap change: the simulator's reagent + monomer network (chem_p31, CPU) against the calibration ODE
(reagent_fit.py, fit B) for the Fig. 5 protocol: 20 mM monomer, 200 mM reagent once, no supply, no chains,
22 C and -3 C (liquid via the validation solute term), 400 min. Prints the mean activated blocks per site (mM)
on the 32-point record grid and the ODE value; output reagent_simcheck.json. CPU only.
    python reagent_simcheck.py [module]
"""
import json, math, os, sys
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
mod = sys.argv[1] if len(sys.argv) > 1 else "chem_p31"
chem = __import__(mod)
import reagent_fit as RF
FIT = json.load(open(os.path.join(HERE, "reagent_fit.json")))["B_with_Ea"]
chem.PAIR = "AS2U"
F = chem.REAGENT_FIT
out = {}
for T in (22.0, -3.0):
    ov = dict(T=T, depth=0.0, salinity=0.25, solute_osm=1.2, P=100.0, pH=8.0, Mg=0.03, x_nuc=1.0, feed=1e-9, tau=1e12,
              k_c=1e-12, init_factor=0.0, k_UV=0.0, UV=1e-4, H2=1e-6, Y_UV=1e-4, Y_H2=1e-4, mix=1e-4,
              act_a=0.0, k_act=F["k_act"], E_actv=F["E_actv"], k_r=F["k_r"], E_r=F["E_r"], y_ref=F["y_ref"], dE_y=F["dE_y"],
              t_half_act=math.log(2) / F["k_h"], E_act_h=F["E_h"], T_p=T, f_p=1e-3, tau_p=1 / 3600,
              reagent_R0=200.0, init_monomer_mM=20.0, t_end=400 / 60, k_t=0.0)
    u = torch.full((2, chem.DIM_P31), 0.5)
    g = torch.Generator(); g.manual_seed(5)
    o = chem.simulate(u, g, (16, 64, 400), record=True, override=ov)
    ts = np.arange(33) / 32 * 400 / 60
    sim = (o["mstar_series"][0] / chem.V).numpy()
    per = {"RT": FIT["RT"], "m3C": FIT["m3C"]}["RT" if T > 0 else "m3C"]
    ode = RF.sim((per["k_act"], per["k_r"], per["y"], per["k_h"]), ts)[2]
    out[str(T)] = {"t_h": ts.tolist(), "sim_mM": sim.tolist(), "ode_mM": ode.tolist(),
                   "max_abs_diff_mM": float(np.nanmax(np.abs(sim - ode))), "liquid": bool(o["liquid"][0]),
                   "reagent_end_mM": float(o["reagent_end_mM"][0])}
    print(T, "liquid", bool(o["liquid"][0]), "max|sim-ode| mM", round(out[str(T)]["max_abs_diff_mM"], 3),
          "sim@2h", round(float(np.interp(2, ts, sim)), 3), "ode@2h", round(float(np.interp(2, ts, ode)), 3),
          "sim@400min", round(float(sim[-1]), 3), "ode", round(float(ode[-1]), 3))
json.dump(out, open(os.path.join(HERE, "reagent_simcheck.json"), "w"), indent=1)
