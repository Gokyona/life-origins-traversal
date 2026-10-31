"""Review 6 gates for the Cornish-Fisher Poisson sampler (candidates chem_p31b_r6.py, chem4c4b_r6.py), on the XPU:
  agreement with the current sampler (chem_p31b._poisson_inv, same seeds) on >= 1e8 draws, >= 99.9999%;
  2e(ii) at lam = 0.02, 0.5, 5, 50, 1e3, 1e4 (validate_p31b criteria); search iterations logged (cap 200 raises).
    python gate_r6.py out.json
"""
import json, math, os, sys, time
import importlib.util
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ckpt  # noqa: E402
import validate_p31b as VB  # noqa: E402

dev = torch.device("xpu")


def load(n, p):
    s = importlib.util.spec_from_file_location(n, os.path.join(HERE, p))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def g(sd):
    x = torch.Generator(device=dev)
    x.manual_seed(sd)
    return x


def main():
    out = sys.argv[1]
    old = load("old_p31b", "chem_p31b.py")
    res = {"candidates": {f: os.popen(f"md5sum {os.path.join(HERE, f)}").read().split()[0] for f in ("chem_p31b_r6.py", "chem4c4b_r6.py")}}
    ok = True
    for cand in ("chem_p31b_r6.py", "chem4c4b_r6.py"):
        new = load("new_" + cand[:-3], cand)
        agree, total, per = 0, 0, {}
        n = 1 << 22
        for lam in (0.02, 0.5, 5.0, 50.0, 1e3, 1e4, "mixed"):
            a_l = 0
            for sd in range(4):
                if lam == "mixed":
                    gg = torch.Generator(device=dev); gg.manual_seed(500 + sd)
                    L = torch.pow(10.0, -3 + 7 * torch.rand(n, device=dev, generator=gg))
                else:
                    L = torch.full((n,), lam, device=dev)
                a = old._poisson_inv(L, g(7100 + sd))
                b = new._poisson_inv(L, g(7100 + sd))
                a_l += int((a == b).sum())
            per[str(lam)] = a_l / (4 * n)
            agree += a_l
            total += 4 * n
        frac = agree / total
        d = {"draws": total, "agreement": frac, "per_lambda": per, "max_iter": new.INV_STATS["max_iter"],
             "pass_agreement": frac >= 0.999999}
        if cand == "chem_p31b_r6.py":
            crit = lambda r: abs(r["z_mean"]) < 3 and abs(r["z_var"]) < 3 and r["p_chi2_mean"] > 0.01 and r["p_chi2_var"] > 0.01 \
                and abs(r["rel_mean"]) <= 0.005
            e = {}
            for lam in (0.02, 0.5, 5.0, 50.0, 1e3, 1e4):
                r = VB._seedcheck(lambda sd, nn: new._poisson_inv(torch.full((nn,), lam, device=dev), g(sd)), lam, lam, lam * (1 + 3 * lam))
                r["pass"] = crit(r)
                e[str(lam)] = r
                print("2e(ii)", lam, {k: v for k, v in r.items() if k != "seed_means_rel"}, flush=True)
            d["2e_ii"] = e
            d["pass_2e"] = all(v["pass"] for v in e.values())
            d["max_iter_after_2e"] = new.INV_STATS["max_iter"]
        d["pass"] = d["pass_agreement"] and d.get("pass_2e", True)
        ok &= d["pass"]
        res[cand] = d
        print(cand, {k: v for k, v in d.items() if k != "2e_ii"}, flush=True)
        ckpt.atomic_json(res, out)
    res["pass"] = bool(ok)
    ckpt.atomic_json(res, out)
    print("GATE R6", "PASS" if ok else "FAIL", flush=True)


if __name__ == "__main__":
    main()
