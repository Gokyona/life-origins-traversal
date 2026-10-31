#!/bin/bash
# P3-1b: after queue_s3.sh (PID $1) succeeded (chem4c5b validation PASS): the coex variants, then the L2 detector replays.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
while kill -0 $1 2>/dev/null; do sleep 15; done
grep -q "chem4c5b validation PASS" logs/queue_s3.log || { echo "STOP: queue_s3 did not finish with the chem4c5b validation passing"; exit 1; }
B=$($PY -c "import json; d=json.load(open('validate_4c5b.json')); print(d['test_c']['B'])")
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_coex.py s3_coex --variants base,A0,A2,A2a,A3,A4 --batch $B > logs/s3_coex_run.log 2>&1 || { echo "STOP: coex runner failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants det_LH5,det_LH7,det_SPAN4,det_SPAN6 --resume > logs/s3_l2_det_run.log 2>&1 || { echo "STOP: L2 detector replays failed"; exit 1; }
echo "QUEUE S3B DONE"
