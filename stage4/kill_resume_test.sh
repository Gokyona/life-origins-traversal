#!/bin/bash
# Kill-and-resume test for run4c1_map.py: an uninterrupted run (ka) and a run killed with
# kill -9 in the middle of a COMPASS round and resumed (kb) must give identical results.
cd "$(dirname "$0")"
S=${1:-kill_test}
PY=${PY:-/home/bern/ComfyUI/venv/bin/python3}
rm -rf $S/ka $S/kb; mkdir -p $S
A="--n-test 40 --gt-seeds 4 --budget 4096 --batch 512"
$PY -u run4c1_map.py --out-dir $S/ka $A > $S/ka.log 2>&1
$PY -u run4c1_map.py --out-dir $S/kb $A > $S/kb1.log 2>&1 & P=$!
until [ -f $S/kb/engine_rep0.jsonl ] && [ $(wc -l < $S/kb/engine_rep0.jsonl) -ge 5 ]; do sleep 1; done
kill -9 $P; wait $P 2>/dev/null
echo "killed after $(wc -l < $S/kb/engine_rep0.jsonl) engine batches"
$PY -u run4c1_map.py --out-dir $S/kb $A --resume > $S/kb2.log 2>&1
$PY compare_runs.py $S/ka $S/kb
