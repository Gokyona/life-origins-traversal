#!/bin/bash
# P3-1b queue after Reviews 7 and 8 (replaces queue_s3c.sh from the end of the L2 r6 reruns on). Waits for PID $1.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
while kill -0 $1 2>/dev/null; do sleep 15; done
grep -q "^done" logs/s3_l2_r6_run.log || { echo "STOP: L2 r6 reruns did not finish"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_p31b.py validate_p31b_t1_r7.json 1 > logs/validate_p31b_t1_r7_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_p31b_t1_r7.json')); sys.exit(0 if d['test_1']['pass'] else 1)" || { echo "STOP: test 1 failed on the Review 7 switches"; exit 1; }
echo "test 1 (Review 7 switches off) PASS"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants base_rec,A2b,A2c,A2d --resume > logs/s3_l2_r7_run.log 2>&1 || { echo "STOP: Review 7 diagnostics failed"; exit 1; }
echo "REVIEW 7 DIAGNOSTICS DONE"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants A5_part1 --points p3709,engine_rep1_9823,engine10653 --resume > logs/s3_l2_a5p1_run.log 2>&1 || { echo "STOP: A5 part 1 failed"; exit 1; }
echo "A5 PART 1 DONE"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_4c5b.py validate_4c5b.json abc > logs/validate_4c5b_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_4c5b.json')); sys.exit(0 if all(d.get('test_'+t,{}).get('pass') for t in 'abc') else 1)" || { echo "STOP: chem4c5b validation failed"; exit 1; }
echo "chem4c5b validation PASS"
B=$($PY -c "import json; d=json.load(open('validate_4c5b.json')); print(d['test_c']['B'])")
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_coex.py s3_coex --variants base,A0,A2,A2a,A3,A4 --batch $B > logs/s3_coex_run.log 2>&1 || { echo "STOP: coex runner failed"; exit 1; }
echo "COEX DONE"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants det_LH5,det_LH7,det_SPAN4,det_SPAN6 --resume > logs/s3_l2_det_run.log 2>&1 || { echo "STOP: L2 detector replays failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants A5_part2 --points engine_rep1_1926,p1558 --resume > logs/s3_l2_a5p2_run.log 2>&1 || { echo "STOP: A5 part 2 failed"; exit 1; }
echo "QUEUE S3D DONE"
