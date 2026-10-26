#!/bin/bash
# Night watchdog for the queued stage-4c-1 jobs (in order: map4c1, sens4c1, map4c1_gt2).
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
      "map4c1_gt2|run4c1_map.py|--n-test 2000 --gt-seeds 32 --budget 64000 --reps 1")
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
  pid=$(pgrep -f -- "-u $script --out-dir $dir " | head -1)
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
    setsid nohup $PY -u $script --out-dir $dir $args $resume >> "${dir}_run.log" 2>&1 < /dev/null &
    sleep 30
  fi
  sleep 60
done
