#!/bin/bash
# P3-1 GPU queue: one job at a time, each in its own user scope with MemoryMax 24G.
cd "$(dirname "$0")"
PY=/home/bern/ComfyUI/venv/bin/python3
run() { echo "$(date '+%F %T') START $*"; systemd-run --user --scope -q -p MemoryMax=24G -p MemorySwapMax=0 "$@"; echo "$(date '+%F %T') END rc=$? $*"; }
for step in "$@"; do
  case $step in
    pf_as2u) run $PY -u preflight_p31.py --pair AS2U --sizes 128,256,512 --out preflight_p31_AS2U.json --job as2u > preflight_p31_AS2U_run.log 2>&1 ;;
    pf_gc)   run $PY -u preflight_p31.py --pair GC --sizes 512 --out preflight_p31_GC.json --job gc > preflight_p31_GC_run.log 2>&1 ;;
    pf_au)   run $PY -u preflight_p31.py --pair AU --sizes 512 --out preflight_p31_AU.json --job au > preflight_p31_AU_run.log 2>&1 ;;
    v24)     run $PY -u validate_p31.py validate_p31.json 24 > validate_p31_t24_run.log 2>&1 ;;
    v78)     run $PY -u validate_p31.py validate_p31.json 78 > validate_p31_t78_run.log 2>&1 ;;
    v4ec)    LOT_EMPTY_CACHE_EVERY=1 run $PY -u validate_p31.py validate_p31_t4_emptycache.json 4 > validate_p31_t4_emptycache_run.log 2>&1 ;;
    pf2)     LOT_EMPTY_CACHE_EVERY=1 run $PY -u preflight_p31.py --pair AS2U --sizes 256,512 --out preflight_p31_AS2U_ec.json --job as2u-ec > preflight_p31_AS2U_ec_run.log 2>&1 ;;
    pf2gc)   LOT_EMPTY_CACHE_EVERY=1 run $PY -u preflight_p31.py --pair GC --sizes 512 --out preflight_p31_GC_ec.json --job gc-ec > preflight_p31_GC_ec_run.log 2>&1 ;;
    pf2au)   LOT_EMPTY_CACHE_EVERY=1 run $PY -u preflight_p31.py --pair AU --sizes 512 --out preflight_p31_AU_ec.json --job au-ec > preflight_p31_AU_ec_run.log 2>&1 ;;
    kr)      run bash kill_resume_p31.sh kill_test_p31 > kill_resume_p31_run.log 2>&1 ;;
  esac
done
echo "$(date '+%F %T') QUEUE DONE"
