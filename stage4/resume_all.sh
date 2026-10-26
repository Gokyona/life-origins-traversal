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
  if pgrep -f "run4c1_map.py --out-dir $d " >/dev/null; then echo "$d: already running"; continue; fi
  args=$($PY -c "import json;c=json.load(open('$cfg'));print(' '.join(f'--{k.replace(\"_\",\"-\")} {v}' for k,v in c.items()))")
  echo "$(date '+%F %T') resuming $d: $args" | tee -a resume_all.log
  nohup $PY -u run4c1_map.py --out-dir "$d" $args --resume >> "$d/run.log" 2>&1 &
done
