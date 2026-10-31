#!/bin/bash
# P3-1b queue after Review 10: chem4c5b validation (a) and (c) on chem4c4b with EMPTY_CACHE_EVERY, then coex (8 variants),
# L2 detector replays, A5 part 2. Stops at the first failure.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_4c5b.py validate_4c5b.json ac > logs/validate_4c5b_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_4c5b.json')); sys.exit(0 if all(d.get('test_'+t,{}).get('pass') for t in 'ac') else 1)" || { echo "STOP: chem4c5b validation (a or c) failed"; exit 1; }
echo "chem4c5b validation PASS"
B=$($PY -c "import json; d=json.load(open('validate_4c5b.json')); print(d['test_c']['B'])")
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_coex.py s3_coex --variants base,A2b,A2c,A0,A2,A2a,A3,A4 --batch $B > logs/s3_coex_run.log 2>&1 || { echo "STOP: coex runner failed"; exit 1; }
echo "COEX DONE"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants det_LH5,det_LH7,det_SPAN4,det_SPAN6 --resume > logs/s3_l2_det_run.log 2>&1 || { echo "STOP: L2 detector replays failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants A5_part2 --points engine_rep1_1926,p1558 --resume > logs/s3_l2_a5p2_run.log 2>&1 || { echo "STOP: A5 part 2 failed"; exit 1; }
echo "QUEUE S3G DONE"
