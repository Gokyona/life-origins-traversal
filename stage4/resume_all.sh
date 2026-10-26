#!/bin/bash
# Restart every unfinished stage-4 job with --resume (after a power cut or reboot).
# A job folder is unfinished when it has config.json but no result.json.
# Run by hand:   bash /run/media/bern/Data2T/myClaude/lot-s4/resume_all.sh
# (Starting this automatically at boot changes the machine and needs the owner's approval.)
cd "$(dirname "$0")"
PY=${PY:-/home/bern/ComfyUI/venv/bin/python3}
for cfg in */config.json; do
  d=$(dirname "$cfg")
  [ -f "$d/result.json" ] && continue
  case "$d" in sens4c1*) script=sens4c1_hydrolysis.py ;; map4c2*) script=run4c2_map.py ;; sens4c2*) script=sens4c2_detector.py ;; *) script=run4c1_map.py ;; esac
  if pgrep -f "$script --out-dir $d " >/dev/null; then echo "$d: already running"; continue; fi
  args=$($PY -c "
import json; c = json.load(open('$cfg'))
keys = ('n_test', 'gt_seeds', 'budget', 'reps', 'batch', 'fidelity', 'methods', 'seeds', 'detector')
a = [f'--{k.replace(\"_\", \"-\")} {c[k]}' for k in keys if k in c]
print(' '.join(a + (['--focus'] if c.get('focus') else [])))")
  echo "$(date '+%F %T') resuming $d with $script: $args" | tee -a resume_all.log
  nohup $PY -u $script --out-dir "$d" $args --resume >> "$d/run.log" 2>&1 &
done
