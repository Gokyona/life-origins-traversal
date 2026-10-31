#!/bin/bash
# P3-1b section 1 queue after tests e,2,3,4,1: test 5 (kill/resume), then test 6 (pre-flight, batch 256, 128 if 256
# fails). One GPU job at a time, each in a user scope MemoryMax 24G; stops at the first failure.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; rc=$?; echo "$(date '+%F %T') END rc=$rc $*"; return $rc; }
while pgrep -f "python3 -u validate_p31b.py validate_p31b.json" > /dev/null; do sleep 10; done
$PY -c "import json,sys; d=json.load(open('validate_p31b.json')); ok=all(d[f'test_{t}']['pass'] for t in 'e2341'); print('tests e2341', ok); sys.exit(0 if ok else 1)" || { echo "STOP: validation failed"; exit 1; }
LOT_EMPTY_CACHE_EVERY=1 run bash kill_resume_p31b.sh kill_test_p31b > logs/kill_resume_p31b_run.log 2>&1
grep -q IDENTICAL logs/kill_resume_p31b_run.log || { echo "STOP: test 5 failed"; exit 1; }
for B in 256 128; do
  LOT_EMPTY_CACHE_EVERY=1 run $PY -u preflight_p31b.py --pair none --multi 1 --sizes $B --fid screen \
      --sets coex,box2,coex,coex,coex,coex --out preflight_p31b_$B.json --job p31b-t6-$B > logs/preflight_p31b_${B}_run.log 2>&1
  if $PY -c "import json,sys; d=json.load(open('preflight_p31b_$B.json')); sys.exit(0 if d['sizes']['$B'].get('pass') else 1)"; then
    echo "TEST 6 PASS at batch $B"; exit 0
  fi
  echo "test 6 failed at batch $B"
done
echo "STOP: test 6 failed at 256 and 128"; exit 1
