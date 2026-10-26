#!/bin/bash
# Resume the stage-4 queue after a power cut or reboot: (re)start the watchdog in its own
# systemd user scope if it is not already running. The watchdog's JOBS list (watchdog.sh) is
# the only source of jobs; it runs one job at a time, each in its own scope (lot-<dir>,
# MemoryMax 24G), with --resume, so this script never starts a GPU job itself.
# Run by hand:   bash /run/media/bern/Data2T/myClaude/lot-s4/resume_all.sh
# (Starting this automatically at boot changes the machine and needs the owner's approval.)
cd "$(dirname "$0")"
if systemctl --user is-active -q lot-watchdog.scope || pgrep -f -- "watchdog.sh" >/dev/null; then
  echo "watchdog already running"
  exit 0
fi
echo "$(date '+%F %T') resume_all: starting watchdog" | tee -a resume_all.log
systemctl --user reset-failed lot-watchdog.scope 2>/dev/null
setsid nohup systemd-run --user --scope --unit=lot-watchdog -p MemoryMax=1G ./watchdog.sh > /dev/null 2>&1 < /dev/null &
