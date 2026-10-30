#!/bin/bash
# Post-cold-gap GPU queue (chem_p31 md5 8b6da3c1...): waits for the running old-code tests 7/8, then one job at a
# time, each in its own user scope (MemoryMax 24G).
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
while pgrep -f "validate_p31.py validate_p31.json 78" >/dev/null; do sleep 20; done
echo "$(date '+%F %T') old tests 7/8 finished"
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; echo "$(date '+%F %T') END rc=$? $*"; }
run $PY -u validate_p31.py validate_p31_coldgap.json 24 > validate_p31_cg_t24_run.log 2>&1
run $PY -u validate_p31.py validate_p31_coldgap.json 78 > validate_p31_cg_t78_run.log 2>&1
run bash kill_resume_p31.sh kill_test_p31 > kill_resume_p31_run.log 2>&1
LOT_EMPTY_CACHE_EVERY=1 run $PY -u preflight_p31.py --pair AS2U --sizes 256,512 --out preflight_p31_AS2U_cg.json --job as2u-cg > preflight_p31_AS2U_cg_run.log 2>&1
echo "$(date '+%F %T') QUEUE DONE"
