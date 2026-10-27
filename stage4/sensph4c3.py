"""Stage 4c-3: acid-branch sensitivity at the reproducible L2 points of the box2 scan (#3709,
#813; coordinator request 27 Oct 2026). The hydrolysis law (Li & Breaker) has no acid branch,
so at pH < 5 lifetimes may be overestimated; pH enters chem4c3 only through hydrolysis, so a
pH floor of 5 for hydrolysis is the override pH = max(pH, 5). Per point and arm (16 reps,
validation fidelity, both copy_random controls): L2 strict/loose, mean margin, mean copy length,
share of chains at L_MAX (copies bounded by the length cap?), and the L_H-mer turnover time.
    python sensph4c3.py [out.json]
"""
import os, sys
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402
import preflight_mem as pf  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c3 pH-floor sensitivity")
FID, REPS = (16, 256, 600), 16


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sensph4c3.json")
    allx = pf.box2(4096, 0)
    res = ckpt.load_json(path) or {}
    for i in (3709, 813):
        x = allx[i]
        ph = float(chem.decode(torch.tensor(x[None, :]))["pH"][0])
        for arm, ov in (("as_is", {}), ("pH_floor_5", {"pH": max(ph, 5.0)})):
            key = f"p{i}_{arm}"
            if key in res:
                continue
            u = torch.tensor(np.repeat(x[None, :], REPS, 0), dtype=torch.float32, device=dev)
            g = torch.Generator(device=dev)
            g.manual_seed(80000 + i)
            o = chem.simulate(u, g, FID, record=True, override=dict(ov))
            strict, s_m, s_c = chem.paired_control_l2(u, 80000 + i, FID, o, override=ov)
            guard.check()
            h = o["length_hist"]
            mk = o["L2_loose"]
            res[key] = {"pH_model": ph, "pH_used": ov.get("pH", ph), "L2_loose": int(mk.sum()), "L2_strict": int(strict.sum()),
                        "margin_mean": float(torch.nan_to_num(s_m - s_c, nan=0.0)[mk].mean()) if bool(mk.any()) else None,
                        "copy_len_mean": float((o["copy_mean_len"] * o["copies"]).sum() / o["copies"].sum().clamp(min=1)),
                        "share_chains_at_L_MAX": float(h[:, chem.L_MAX].sum() / h[:, 1:].sum().clamp(min=1)),
                        "tau_h_hours": float(o["L2_tau_h"].mean()),
                        "R_long_per_template_per_lifetime_median": float(o["long_copies_per_template_per_lifetime"].median())}
            ckpt.atomic_json(res, path)
            print(key, res[key], flush=True)
    res["thermal"] = guard.summary()
    ckpt.atomic_json(res, path)


if __name__ == "__main__":
    main()
