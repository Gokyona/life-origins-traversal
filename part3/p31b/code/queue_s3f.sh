#!/bin/bash
# P3-1b queue after Review 9 (replaces queue_s3e.sh from the end of A5 part 1 on). Waits for PID $1.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
while kill -0 $1 2>/dev/null; do sleep 15; done
grep -q "^done" logs/s3_l2_a5p1_run.log || { echo "STOP: A5 part 1 did not finish"; exit 1; }
echo "A5 PART 1 DONE"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_p31b.py validate_p31b_t1_r9.json 1 > logs/validate_p31b_t1_r9_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_p31b_t1_r9.json')); sys.exit(0 if d['test_1']['pass'] else 1)" || { echo "STOP: test 1 failed on the Review 9 code"; exit 1; }
echo "test 1 (Review 9 code) PASS"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants A2e,base_conc,A2b_conc,A2c_conc --resume > logs/s3_l2_r9_run.log 2>&1 || { echo "STOP: Review 9 L2 runs failed"; exit 1; }
echo "A2E DONE"
$PY check_conc_identity.py s3_l2/runs.jsonl > logs/check_conc_identity.log 2>&1 || { echo "STOP: copy concentration record changed the dynamics"; exit 1; }
echo "CONC IDENTITY PASS"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_4c5b.py validate_4c5b.json abc > logs/validate_4c5b_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_4c5b.json')); sys.exit(0 if all(d.get('test_'+t,{}).get('pass') for t in 'abc') else 1)" || { echo "STOP: chem4c5b validation failed"; exit 1; }
echo "chem4c5b validation PASS"
B=$($PY -c "import json; d=json.load(open('validate_4c5b.json')); print(d['test_c']['B'])")
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_coex.py s3_coex --variants base,A2b,A2c,A0,A2,A2a,A3,A4 --batch $B > logs/s3_coex_run.log 2>&1 || { echo "STOP: coex runner failed"; exit 1; }
echo "COEX DONE"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants det_LH5,det_LH7,det_SPAN4,det_SPAN6 --resume > logs/s3_l2_det_run.log 2>&1 || { echo "STOP: L2 detector replays failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants A5_part2 --points engine_rep1_1926,p1558 --resume > logs/s3_l2_a5p2_run.log 2>&1 || { echo "STOP: A5 part 2 failed"; exit 1; }
echo "QUEUE S3F DONE"
