"""Review 9: the copy concentration record must not change the dynamics: base_conc / A2b_conc / A2c_conc must equal
base / A2b / A2c in every field stored in both (the record fields exist only in the *_conc lines)."""
import json, sys
L = [json.loads(l) for l in open(sys.argv[1])]
D = {d["variant"]: d for d in L}
ok = True
for a, b in (("base", "base_conc"), ("A2b", "A2b_conc"), ("A2c", "A2c_conc")):
    ra, rb = D[a]["runs"], D[b]["runs"]
    diff = [k for k in ra if k in rb and ra[k] != rb[k]]
    print(a, b, "shared fields", sum(k in rb for k in ra), "differing", diff)
    ok &= not diff
print("CONC IDENTITY", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
