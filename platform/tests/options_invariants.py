"""Invariant (conservation) checks of the Gokyona options "ration" and "protected_domain" (29 Oct 2026). CPU only.
Darwin demo, 8 sites x 256 slots, 64 vesicle slots, L_MAX 256; checked after every step:
  monomer pools F, M, Ms >= 0; 0 <= L <= L_MAX; no bit set beyond a chain's length; empty slots are all-zero with no
  label, pair or tag; duplex pairs symmetric and within the row; a label points to a live vesicle; protected-domain mark
  true on empty slots.
    python tests/options_invariants.py out.json [seconds per config] [name filter, comma-separated substrings]"""
import json, os, sys, time
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gokyona_env", "engine"))
import gk_chem, words as wd  # noqa: E402
from trait import Trait  # noqa: E402

torch.set_num_threads(4)
OUT = sys.argv[1]
SEC = float(sys.argv[2]) if len(sys.argv) > 2 else 45
PRE = json.load(open(os.path.join(HERE, "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
N, S = 8, 256
gk_chem.configure(4, 64)
NV = gk_chem.NV


class Check:
    def __init__(self):
        self.k, self.t0, self.fail, self.prev = 0, time.time(), {}, None

    def due(self):
        return True

    def bad(self, name, cond):
        if bool(cond):
            self.fail[name] = self.fail.get(name, 0) + 1

    def __call__(self, loc):
        self.k += 1
        L, seq, chem, pair, lab, VS = loc["L"], loc["seq"], loc["chem"], loc["pair"], loc["lab"], loc["VS"]
        F, M, Ms = loc["F"], loc["M"], loc["Ms"]
        R = L.shape[0]
        self.bad("negative_monomer_pool", (F < -1e-6).any() or (M < -1e-6).any() or (Ms < -1e-6).any())
        self.bad("length_out_of_range", (L < 0).any() or (L > gk_chem.L_MAX).any())
        self.bad("bits_beyond_length", ((seq & ~wd.mask(L)) != 0).any())
        emp = L == 0
        self.bad("empty_slot_not_clear", (seq[emp] != 0).any() or (lab[emp] != 0).any() or (pair[emp] >= 0).any() or (loc["ftag"][emp] != 0).any())
        pr, pc = (pair >= 0).nonzero(as_tuple=True)
        if pr.numel():
            self.bad("pair_out_of_row", (pair[pr, pc] >= S).any())
            self.bad("pair_not_symmetric", (pair[pr, pair[pr, pc]] != pc).any())
        lr, lc = (lab > 0).nonzero(as_tuple=True)
        if lr.numel():
            self.bad("label_to_dead_vesicle", (~VS.alive[lr, lab[lr, lc].long() - 1]).any())
        if loc.get("prom") is not None:
            self.bad("prom_false_on_empty", (~loc["prom"][emp]).any())
        if time.time() - self.t0 > SEC:
            return "stop"


res = {"size": f"{N} sites x {S} slots, 64 vesicle slots, L_MAX 256, Darwin demo, seed 1", "configs": {}}
for name, extra in (("off", {}), ("ration_random", {"ration": "random"}), ("ration_vesicle", {"ration": "vesicle"}),
                    ("protected", {"protected_domain": True}), ("random+protected", {"ration": "random", "protected_domain": True}),
                    ("random+guard", {"ration": "random", "domain_break_guard": 24}),
                    ("random+P1", {"ration": "random", "promoter_required": True}),
                    ("random+env_season", {"ration": "random", "_trait": {"env_landscape": {"schedule": "season"}}}),
                    ("random+env_trend", {"ration": "random", "_trait": {"env_landscape": {"schedule": "trend", "period_years": 20}}}),
                    ("random+catastrophe+shift", {"ration": "random", "catastrophe_rate": 1 / (2 * 8766.0), "catastrophe_kill": 0.95,
                                                  "_trait": {"env_landscape": {"schedule": "catastrophe", "shift_years": 2}}}),
                    ("random+env_chem_T_trend", {"ration": "random", "env_chem": {"variable": "T", "v0": 4.47, "v1": 60.0, "schedule": "trend", "period_years": 20},
                                                 "_trait": {"env_landscape": {"schedule": "trend", "period_years": 20, "variable": "T"}}}),
                    ("random+env_chem_pH_season", {"ration": "random", "env_chem": {"variable": "pH", "v0": 7.58, "v1": 6.0, "schedule": "season"}}),
                    ("random+env_chem_salinity_trend", {"ration": "random", "env_chem": {"variable": "salinity", "v0": 0.74, "v1": 0.05, "schedule": "trend", "period_years": 20}}),
                    ("random+gradient_moving_line_local", {"ration": "random", "migration": 1e-2, "migration_mode": "line", "local_pool_mix": True,
                                       "_trait": {"env_landscape": {"schedule": "constant", "gradient": True, "band_amp": 0.5, "band_period_years": 10}}}),
                    ("random+gradient_season_uv", {"ration": "random", "migration": 1e-3, "migration_mode": "line", "local_pool_mix": True,
                                       "uv_gradient": {"uv_pole": 0.017, "uv_eq": 0.068},
                                       "_trait": {"env_landscape": {"schedule": "constant", "gradient": True, "band_amp": 0.5, "band_period_years": 10, "season_amp": 0.3}}}),
                    ("sg_B_kmax2_noreencap", {"ration": "random", "k_max": 2, "no_reencapsulation": True, "vesicle_stats": True}),
                    ("sg_C_content", {"ration": "random", "k_max": 2, "no_reencapsulation": True, "division_mode": "content", "vesicle_stats": True}),
                    ("sg_D_segregate", {"ration": "random", "k_max": 2, "no_reencapsulation": True, "division_mode": "content", "segregate": True, "vesicle_stats": True}),
                    ("sg_Dp_divtrait2", {"ration": "random", "k_max": 2, "no_reencapsulation": True, "division_mode": "content", "segregate": True,
                                         "division_n0": "kmax", "division_trait_s": 2.0, "vesicle_stats": True}),
                    ("motif_local", {"ration": "random", "p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64},
                                     "division_mode": "content", "local_conc": True, "vesicle_stats": True}),
                    ("motif_local_minsize", {"ration": "random", "p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64},
                                             "division_mode": "content", "local_conc": True, "vesicle_min_size": 0.5, "vesicle_stats": True}),
                    ("growth_uncapped", {"ration": "random", "growth_uncapped": True}),
                    ("copy_kinetic_uncapped", {"ration": "random", "growth_uncapped": True, "copy_kinetic": True}),
                    ("copy_kinetic_uncapped_dt10", {"ration": "random", "growth_uncapped": True, "copy_kinetic": True, "dt_max": 10.0}),
                    ("kinetic_competition_dt10", {"ration": "random", "growth_uncapped": True, "copy_kinetic": True, "monomer_competition": True, "dt_max": 10.0}),
                    ("copy_physical_full", {"ration": "random", "growth_uncapped": True, "copy_kinetic": True, "monomer_competition": True, "copy_physical": True}),
                    ("growthphys_mineral", {"ration": "random", "growth_physical": True, "mineral_oligo": {"yield": 0.01}, "growth_uncapped": True,
                                            "copy_kinetic": True, "monomer_competition": True, "copy_physical": True, "long_primers": {"min": 8, "max": 16}, "dt_max": 30.0}),
                    ("random+impact", {"ration": "random", "impact": {"rate_per_hour": 1 / (3 * 8766.0), "kill": 0.99, "pulse_years": 1, "relax_years": 1, "feed_mult": 3, "feed_years": 1},
                                       "env_chem": {"variable": "T", "v0": 4.47, "v1": 40.0, "schedule": "constant"},
                                       "_trait": {"env_landscape": {"schedule": "constant", "impact": {"pulse_years": 1, "relax_years": 1}}}})):
    if len(sys.argv) > 3 and not any(k in name for k in sys.argv[3].split(",")):
        continue
    extra = dict(extra)
    tx = extra.pop("_trait", {})
    u = torch.tensor([PRE["u"]], dtype=torch.float32)
    ov = {k: float(v) for k, v in PRE["env"].items() if k != "avail_time"}
    ov["mineral_type"] = float(round(ov.get("mineral_type", 0)))
    ov["endless"] = True
    ov["census_h"] = 20 * 8766.0
    ov.update(extra)
    tr = Trait(dict(PRE["trait"], enabled=True, **tx), N, torch.device("cpu"))
    ck = Check()
    g = torch.Generator(); g.manual_seed(1)
    gk_chem.simulate(u, g, (N, S, 300), override=ov, trait=tr, hook=ck)
    res["configs"][name] = {"steps_checked": ck.k, "failures": ck.fail, "pass": not ck.fail}
    print(name, res["configs"][name], flush=True)
res["all_pass"] = all(c["pass"] for c in res["configs"].values())
json.dump(res, open(OUT, "w"), indent=1)
print("ALL PASS" if res["all_pass"] else "FAILURES", flush=True)
