#!/bin/bash
# GPU/RAM monitor while PID $1 lives; writes ALERT at >= 85 C.
P=$1; OUT=$2
while kill -0 $P 2>/dev/null; do
  t=$(cat /sys/class/drm/card0/device/hwmon/hwmon6/temp*_input 2>/dev/null | sort -n | tail -1); t=$((t/1000))
  m=$(xpu-smi stats -d 0 -j 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(round(d['memory']['used_mib']['tile_0']['current']))" 2>/dev/null)
  a=""; [ "$t" -ge 85 ] && a=" ALERT"
  echo "$(date '+%F %T') tmax_C=$t vram_MiB=$m ram_used_MiB=$(free -m | awk '/Mem:/{print $3}')$a" >> "$OUT"
  sleep 20
done
