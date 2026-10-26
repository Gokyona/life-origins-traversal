"""Compare the 4c-1 ground truth of map4c1 (GT batches 0-65 run with chem 0c9d1264, the rest with
19da9687) and map4c1_gt2 (all with 19da9687, same seeds): per-point E1 / E2 majority changes,
split into the early points (batches 0-65 = points 0..1055) and the rest.

  python compare_gt4c1.py map4c1 map4c1_gt2
"""
import json, sys
import numpy as np

sys.path.insert(0, ".")
import run4c1_map as rm  # noqa: E402


def probs(d, n_test=2000, seeds=32):
    recs = rm.gather(rm.load_lines(f"{d}/gt.jsonl"))
    e1, e2 = rm.outcomes(recs)
    return e1.reshape(n_test, seeds).mean(1), e2.reshape(n_test, seeds).mean(1)


def main(a, b):
    a1, a2 = probs(a)
    b1, b2 = probs(b)
    early = np.arange(2000) < 1056
    out = {}
    for name, x, y in (("E1", a1, b1), ("E2", a2, b2)):
        ch = (x > 0.5) != (y > 0.5)
        out[name] = {"changed_early_1056": int(ch[early].sum()), "changed_late_944": int(ch[~early].sum()),
                     "mean_abs_dp_early": float(np.abs(x - y)[early].mean()),
                     "mean_abs_dp_late": float(np.abs(x - y)[~early].mean()),
                     "life_old": float((x > 0.5).mean()), "life_new": float((y > 0.5).mean())}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
