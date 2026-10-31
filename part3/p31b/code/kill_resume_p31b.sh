#!/bin/bash
# P3-1b validation test 5: kill and resume (MULTI_START on, PAIR=None, the 38 coexistence points tiled, 5 batches of 32,
# FID (16, 256, 300)). An uninterrupted run (ka) and a run killed with kill -9 after its second batch and resumed (kb)
# must give identical results (all fields except seconds).
cd "$(dirname "$0")"
S=${1:-kill_test_p31b}
PY=${PY:-/home/bern/ComfyUI/venv/bin/python3}
rm -rf $S/ka $S/kb; mkdir -p $S
A="--set coex --tile 5x32 --fid 16,256,300 --multi 1 --tag kill"
$PY -u run_p31b.py --out-dir $S/ka $A > $S/ka.log 2>&1
$PY -u run_p31b.py --out-dir $S/kb $A > $S/kb1.log 2>&1 & P=$!
until [ -f $S/kb/runs.jsonl ] && [ $(wc -l < $S/kb/runs.jsonl) -ge 2 ]; do sleep 1; done
kill -9 $P; wait $P 2>/dev/null
echo "killed after $(wc -l < $S/kb/runs.jsonl) batches"
$PY -u run_p31b.py --out-dir $S/kb $A --resume > $S/kb2.log 2>&1
$PY - $S/ka/runs.jsonl $S/kb/runs.jsonl <<'PYEOF'
import json, sys
load = lambda p: sorted((json.loads(l) for l in open(p)), key=lambda d: d["batch"])
a, b = load(sys.argv[1]), load(sys.argv[2])
strip = lambda d: {k: v for k, v in d.items() if k != "seconds"}
same = len(a) == len(b) and all(strip(x) == strip(y) for x, y in zip(a, b))
print("batches", len(a), len(b), "runs", sum(d["n"] for d in a), "IDENTICAL" if same else "DIFFERENT")
PYEOF
