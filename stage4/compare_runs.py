"""Compare two run4c1_map.py output folders (kill-and-resume test).

Simulations (every batch's per-run records), COMPASS selections and accuracy curves must be
bit-identical. The final-map analysis (importance, partial dependence, slices) is refitted
on the GPU, whose float reductions are not bit-reproducible even between two uninterrupted
runs, so it is compared with a tolerance of 0.002.
"""
import json, sys


def recs(d, f):
    ls = [json.loads(l) for l in open(f"{d}/{f}")]
    ls.sort(key=lambda x: (x["round"], x["i0"]))
    return [(x["round"], x["i0"], json.dumps(x["runs"])) for x in ls]


def maxdiff(a, b):
    if isinstance(a, dict):
        return max([maxdiff(a[k], b[k]) for k in a] or [0.0])
    if isinstance(a, list):
        return max([maxdiff(x, y) for x, y in zip(a, b)] or [0.0])
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b)
    return 0.0 if a == b else float("inf")


def main(da, db):
    ok = True
    for f in ("gt.jsonl", "sobol_rep0.jsonl", "engine_rep0.jsonl"):
        a, b = recs(da, f), recs(db, f)
        same = a == b
        print(f"{f}: batches {len(a)} / {len(b)}, runs identical: {same}")
        ok &= same
    same = json.load(open(f"{da}/engine_rep0_rounds.json")) == json.load(open(f"{db}/engine_rep0_rounds.json"))
    print("engine selections identical:", same)
    ok &= same
    A, B = json.load(open(f"{da}/result.json")), json.load(open(f"{db}/result.json"))
    for k in ("life_fraction", "composition_E2_positive"):
        print(f"result {k} identical:", A[k] == B[k])
        ok &= A[k] == B[k]
    same = all(A["methods"][m][0]["curves"] == B["methods"][m][0]["curves"] for m in A["methods"])
    print("accuracy curves identical:", same)
    ok &= same
    md = maxdiff(A["maps"], B["maps"])
    print(f"final-map analysis max abs difference: {md:.4g} (tolerance 0.002)")
    ok &= md <= 0.002
    print("KILL-RESUME TEST", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2]) else 1)
