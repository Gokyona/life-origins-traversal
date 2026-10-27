#!/bin/bash
# Watchdog for the queued stage-4c jobs (in order: map4c1, sens4c1, map4c1_gt2, map4c2b, sens4c2,
# map4c2_focus, map4c3; the focus map uses the detector settings adopted by sens4c2; map4c3 is the
# full-space 4c-3 map, coordinator decision 26 Oct 2026). map4c2b is the fresh
# restart of map4c2 (killed by the 2026-10-26 07:07 out-of-memory event; map4c2 is kept as is).
# Each job runs in its own systemd user scope (lot-<dir>, MemoryMax 24G, swap 2G), outside the
# Claude app's cgroup. Start the watchdog itself the same way:
#   systemd-run --user --scope --unit=lot-watchdog -p MemoryMax=1G setsid nohup ./watchdog.sh &
# map4c1_gt2 = full ground-truth rerun with the fixed chem4c1, reusing map4c1's method runs
# (its folder is pre-filled with map4c1's config and sobol/engine files, so it runs with --resume).
# Every 60 s: if the current job has no result.json and is not running, or has written
# nothing for STALL_MIN minutes, it is (killed and) restarted with --resume. Logs to
# resume_all.log. Exits when every job has its result.json.
cd "$(dirname "$0")"
PY=${PY:-/home/bern/ComfyUI/venv/bin/python3}
STALL_MIN=${STALL_MIN:-25}
JOBS=("map4c1|run4c1_map.py|--n-test 2000 --gt-seeds 32 --budget 64000 --reps 1"
      "sens4c1|sens4c1_hydrolysis.py|--n-test 2000 --seeds 8"
      "map4c1_gt2|run4c1_map.py|--n-test 2000 --gt-seeds 32 --budget 64000 --reps 1"
      "map4c2b|run4c2_map.py|--n-test 1000 --gt-seeds 16 --budget 32000 --reps 1 --methods engine"
      "sens4c2|sens4c2_detector.py|"
      "map4c2_focus|run4c2_map.py|--n-test 1000 --gt-seeds 16 --budget 32000 --reps 1 --methods engine --focus --detector sens4c2/result.json"
      "map4c3|run4c3_map.py|--n-test 1000 --gt-seeds 16 --budget 32000 --reps 1 --methods engine --batch 1024")
log() { echo "$(date '+%F %T') watchdog: $*" >> resume_all.log; }
log "started"
while true; do
  current=""
  for j in "${JOBS[@]}"; do
    IFS='|' read -r dir script args <<< "$j"
    [ -f "$dir/result.json" ] && continue
    current=$dir; break
  done
  [ -z "$current" ] && { log "all jobs finished"; exit 0; }
  pid=$(pgrep -f -- "-u $script --out-dir $dir( |\$)" | head -1)
  if [ -n "$pid" ] && [ -d "$dir" ]; then
    newest=$(find "$dir" -type f -printf '%T@\n' | sort -n | tail -1)
    age=$(( $(date +%s) - ${newest%.*} ))
    if [ "$age" -gt $((STALL_MIN * 60)) ]; then
      log "$dir: no output for $((age / 60)) min, killing pid $pid"
      kill "$pid"; sleep 10; kill -9 "$pid" 2>/dev/null; pid=""
    fi
  fi
  if [ -z "$pid" ]; then
    resume=""; [ -f "$dir/config.json" ] && resume="--resume"
    log "starting $script --out-dir $dir $args $resume"
    systemctl --user reset-failed "lot-$dir.scope" 2>/dev/null
    setsid nohup systemd-run --user --scope --unit="lot-$dir" -p MemoryMax=24G -p MemorySwapMax=2G \
      $PY -u $script --out-dir $dir $args $resume >> "${dir}_run.log" 2>&1 < /dev/null &
    sleep 30
  fi
  sleep 60
done
