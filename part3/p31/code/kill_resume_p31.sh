#!/bin/bash
# P3-1 validation test 9: kill and resume. An uninterrupted run (ka) and a run killed with kill -9 after
# its second batch and resumed (kb) must give identical results (all fields except seconds).
cd "$(dirname "$0")"
S=${1:-kill_test_p31}
PAIR=${PAIR:-AS2U}
PY=${PY:-/home/bern/ComfyUI/venv/bin/python3}
rm -rf $S/ka $S/kb; mkdir -p $S
A="--pair $PAIR --n 640 --batch 128"
$PY -u run_p31.py --out-dir $S/ka $A > $S/ka.log 2>&1
$PY -u run_p31.py --out-dir $S/kb $A > $S/kb1.log 2>&1 & P=$!
until [ -f $S/kb/runs_$PAIR.jsonl ] && [ $(wc -l < $S/kb/runs_$PAIR.jsonl) -ge 2 ]; do sleep 1; done
kill -9 $P; wait $P 2>/dev/null
echo "killed after $(wc -l < $S/kb/runs_$PAIR.jsonl) batches"
$PY -u run_p31.py --out-dir $S/kb $A --resume > $S/kb2.log 2>&1
$PY - $S/ka/runs_$PAIR.jsonl $S/kb/runs_$PAIR.jsonl <<'EOF'
import json, sys
load = lambda p: sorted((json.loads(l) for l in open(p)), key=lambda d: d["batch"])
a, b = load(sys.argv[1]), load(sys.argv[2])
strip = lambda d: {k: v for k, v in d.items() if k != "seconds"}
same = len(a) == len(b) and all(strip(x) == strip(y) for x, y in zip(a, b))
print("batches", len(a), len(b), "runs", sum(d["n"] for d in a), "IDENTICAL" if same else "DIFFERENT")
EOF
