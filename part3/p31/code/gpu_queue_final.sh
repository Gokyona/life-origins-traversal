#!/bin/bash
# Final P3-1 GPU queue (chem_p31 after cold-gap + inhibition): one job at a time, each in a user scope MemoryMax 24G.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; echo "$(date '+%F %T') END rc=$? $*"; }
run $PY -u validate_p31.py validate_p31_final.json 24 > validate_p31_final_t24_run.log 2>&1
run $PY -u validate_p31.py validate_p31_final.json 78 > validate_p31_final_t78_run.log 2>&1
run bash kill_resume_p31.sh kill_test_p31 > kill_resume_p31_run.log 2>&1
LOT_EMPTY_CACHE_EVERY=1 run $PY -u preflight_p31.py --pair AS2U --sizes 256,512 --out preflight_p31_AS2U_final.json --job as2u-final > preflight_p31_AS2U_final_run.log 2>&1
echo "$(date '+%F %T') QUEUE DONE"
