#!/bin/bash
# P3-1b after section 1: test 1 on chem_p31b with the section 3 switches (off path), then section 3 on the 7 L2 points,
# then the chem4c5b validation (Review 4 R4.3 a, b, c). One GPU job at a time, user scope MemoryMax 24G; stops at the first
# failure. Waits for queue_s1.sh (PID $1) and requires its TEST 6 PASS.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
while kill -0 $1 2>/dev/null; do sleep 15; done
grep -q "TEST 6 PASS" logs/queue_s1.log || { echo "STOP: test 6 did not pass"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_p31b.py validate_p31b_t1_s3switches.json 1 > logs/validate_p31b_t1_s3switches_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_p31b_t1_s3switches.json')); sys.exit(0 if d['test_1']['pass'] else 1)" || { echo "STOP: test 1 failed on the s3 switches"; exit 1; }
echo "test 1 (s3 switches off) PASS"
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants base,A0,A2,A2a,A4,A6,A1 > logs/s3_l2_run.log 2>&1 || { echo "STOP: s3 L2 runner failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_4c5b.py validate_4c5b.json abc > logs/validate_4c5b_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_4c5b.json')); sys.exit(0 if all(d.get('test_'+t,{}).get('pass') for t in 'abc') else 1)" || { echo "STOP: chem4c5b validation failed"; exit 1; }
echo "chem4c5b validation PASS"
