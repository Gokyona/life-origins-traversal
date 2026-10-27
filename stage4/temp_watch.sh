#!/bin/bash
# Overnight temperature watch for the stage-4 jobs (read-only on hardware: never touches fans,
# limits or power). Every 60 s: reads the Arc B70 (xe hwmon) and CPU package temperatures and
# logs them; if the GPU's hottest sensor >= 90 C or the CPU package >= 90 C, it pauses (SIGSTOP)
# the running stage-4 python jobs and resumes them (SIGCONT) once both are below 75 C.
# Run in its own scope: systemd-run --user --scope --unit=lot-tempwatch -p MemoryMax=256M ./temp_watch.sh
cd "$(dirname "$0")"
LOG=temp_watch.log
paused=0
gpu_h=$(for h in /sys/class/drm/card*/device/hwmon/hwmon*; do [ "$(cat $h/name 2>/dev/null)" = xe ] && echo $h && break; done)
cpu_f=$(for f in /sys/class/hwmon/hwmon*/temp*_label; do [ "$(cat $f)" = "Package id 0" ] && echo ${f%_label}_input && break; done)
while true; do
  g=0; for f in $gpu_h/temp*_input; do v=$(( $(cat $f) / 1000 )); [ $v -gt $g ] && g=$v; done
  c=$(( $(cat $cpu_f) / 1000 ))
  pids=$(pgrep -f "lot-s4|validate4c3.py|scan4c3.py|sens4c3.py|run4c3_map.py|run4c2_map.py" | grep -v $$)
  if [ $paused -eq 0 ] && { [ $g -ge 90 ] || [ $c -ge 90 ]; }; then
    for p in $(pgrep -f "python3 -u (validate4c3|scan4c3|sens4c3|sensph4c3|s2stats4c3|run4c3_map|run4c2_map)"); do kill -STOP $p; done
    paused=1; echo "$(date '+%F %T') PAUSE gpu ${g}C cpu ${c}C" >> $LOG
  elif [ $paused -eq 1 ] && [ $g -lt 75 ] && [ $c -lt 75 ]; then
    for p in $(pgrep -f "python3 -u (validate4c3|scan4c3|sens4c3|sensph4c3|s2stats4c3|run4c3_map|run4c2_map)"); do kill -CONT $p; done
    paused=0; echo "$(date '+%F %T') RESUME gpu ${g}C cpu ${c}C" >> $LOG
  fi
  echo "$(date '+%F %T') gpu ${g}C cpu ${c}C paused=$paused" >> $LOG
  sleep 60
done
