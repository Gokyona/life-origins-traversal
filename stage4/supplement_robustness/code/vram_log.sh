#!/bin/bash
# Logs VRAM used, power (xpu-smi), hottest xe sensor and RAM every 10 s while a process pattern runs; ALERT at >= 85 C.
PAT="$1"; OUT="$2"
while pgrep -f "$PAT" > /dev/null; do
  mp=$(xpu-smi stats -d 0 -j 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(round(d['memory']['used_mib']['tile_0']['current']), round(d['power']['gpu_power_w']['tile_0']['current']))" 2>/dev/null)
  t=$(cat /sys/class/drm/card0/device/hwmon/hwmon6/temp*_input 2>/dev/null | sort -n | tail -1); t=$((t/1000))
  r=$(free -m | awk '/Mem:/{print $3}')
  a=""; [ "$t" -ge 85 ] && a=" ALERT"
  echo "$(date '+%F %T') vram_MiB,power_W=$mp tmax_C=$t ram_used_MiB=$r$a" >> "$OUT"
  sleep 10
done
echo "$(date '+%F %T') process gone" >> "$OUT"
