#!/bin/bash
# P3-1b queue after Review 6 (replaces queue_s3.sh / queue_s3b.sh from the end of the running L2 variants on).
# Waits for the running run_s3_l2.py (PID $1); then the R6 sampler gates; swap if they pass (else keep the current sampler
# and continue); L2 A0_r6/A2_r6 if swapped; chem4c5b validation a, b, c; coex variants; L2 detector replays.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
while kill -0 $1 2>/dev/null; do sleep 15; done
grep -q "^done" logs/s3_l2_run.log || { echo "STOP: s3 L2 variants did not finish"; exit 1; }
SWAP=0
run $PY -u gate_r6.py gate_r6.json > logs/gate_r6_run.log 2>&1
if grep -q "GATE R6 PASS" logs/gate_r6_run.log; then
  cp chem_p31b.py logs/chem_p31b_60c68e27_pre_r6.py; cp chem4c4b.py logs/chem4c4b_fa9e35b3_pre_r6.py
  cp chem_p31b_r6.py chem_p31b.py; cp chem4c4b_r6.py chem4c4b.py
  LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_p31b.py validate_p31b_t1_r6.json 1 > logs/validate_p31b_t1_r6_run.log 2>&1
  if $PY -c "import json,sys; d=json.load(open('validate_p31b_t1_r6.json')); sys.exit(0 if d['test_1']['pass'] else 1)"; then
    SWAP=1; echo "R6 sampler swapped in (gates and test 1 PASS)"
  else
    cp logs/chem_p31b_60c68e27_pre_r6.py chem_p31b.py; cp logs/chem4c4b_fa9e35b3_pre_r6.py chem4c4b.py
    echo "R6: test 1 FAILED after swap, reverted to the current sampler"
  fi
else
  echo "R6: gate FAILED, keeping the current sampler"
fi
if [ $SWAP = 1 ]; then
  LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants A0_r6,A2_r6 --resume > logs/s3_l2_r6_run.log 2>&1 || { echo "STOP: L2 r6 reruns failed"; exit 1; }
fi
LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_4c5b.py validate_4c5b.json abc > logs/validate_4c5b_run.log 2>&1
$PY -c "import json,sys; d=json.load(open('validate_4c5b.json')); sys.exit(0 if all(d.get('test_'+t,{}).get('pass') for t in 'abc') else 1)" || { echo "STOP: chem4c5b validation failed"; exit 1; }
echo "chem4c5b validation PASS"
B=$($PY -c "import json; d=json.load(open('validate_4c5b.json')); print(d['test_c']['B'])")
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_coex.py s3_coex --variants base,A0,A2,A2a,A3,A4 --batch $B > logs/s3_coex_run.log 2>&1 || { echo "STOP: coex runner failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run $PY -u run_s3_l2.py --out-dir s3_l2 --variants det_LH5,det_LH7,det_SPAN4,det_SPAN6 --resume > logs/s3_l2_det_run.log 2>&1 || { echo "STOP: L2 detector replays failed"; exit 1; }
echo "QUEUE S3C DONE"
