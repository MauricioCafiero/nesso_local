#!/bin/zsh
# Queue a multi-ligand screen: launches detached (survives closing the terminal)
# with a memory watchdog that aborts the run if the machine runs out of memory.
#
# usage: zsh code/run_screen.sh <protein.fasta> <ligands.smi> <name> [extra nesso flags]
#
# Watch it with:  tail -f outputs/screens/<name>/run.log
# Results:        outputs/screens/<name>/results.csv
#
# Safe to re-run: finished ligands are skipped, so an aborted screen resumes.
# Note this does NOT survive closing the lid -- macOS suspends the job on sleep.
# Start it once the machine is awake and staying awake.

set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
protein="$1"; ligands="$2"; name="$3"; shift 3
[ -z "$name" ] && { echo "usage: zsh code/run_screen.sh <protein.fasta> <ligands.smi> <name> [flags]"; exit 2; }

dir="$root/outputs/screens/$name"
mkdir -p "$dir"
# One log per attempt, with run.log pointing at the latest. Appending every
# attempt to a single file leaves stale START/ABORT/END markers behind, which
# makes a resumed run look like it failed when reading the log.
log="$dir/run-$(date +%Y%m%d-%H%M%S).log"
ln -sf "$(basename "$log")" "$dir/run.log"

cat > "$dir/.driver.sh" <<DRIVER
#!/bin/zsh
caffeinate -w \$\$ &            # hold off idle sleep for the whole run
export HF_HUB_DISABLE_XET=1     # plain HTTPS LFS; xet stalls on weak wifi
echo "START \$(date)"
"$root/.venv/bin/python" "$root/code/screen_nesso.py" "$protein" "$ligands" --name "$name" $@ &
SPID=\$!
ustreak=0
while kill -0 \$SPID 2>/dev/null; do
  free=\$(memory_pressure 2>/dev/null | awk '/free percentage/{gsub(/%/,"",\$NF); print \$NF}')
  npid=\$(pgrep -f "nesso predict" | head -1)
  st=\$(ps -o stat= -p "\$npid" 2>/dev/null | tr -d ' ')
  if [ "\${free:-100}" -lt 12 ] 2>/dev/null; then
    echo "!! ABORT \$(date +%H:%M:%S): free memory \${free}% -- target too large for this machine"
    kill -9 \$SPID 2>/dev/null; pkill -9 -f "nesso predict" 2>/dev/null; break
  fi
  # state U alone is just swapping; only abort if it never leaves U
  if [ "\$st" = "U" ]; then ustreak=\$((ustreak+1)); else ustreak=0; fi
  if [ \$ustreak -ge 60 ]; then
    echo "!! ABORT \$(date +%H:%M:%S): no progress for 10 min -- wedged"
    kill -9 \$SPID 2>/dev/null; pkill -9 -f "nesso predict" 2>/dev/null; break
  fi
  sleep 10
done
wait \$SPID 2>/dev/null
echo "END \$(date)"
"$root/.venv/bin/python" "$root/code/screen_nesso.py" "$protein" "$ligands" --name "$name" --collect-only || true
DRIVER

chmod +x "$dir/.driver.sh"
"$root/.venv/bin/python" -c "import subprocess,sys; subprocess.Popen(['zsh', sys.argv[1]], stdout=open(sys.argv[2],'ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True, cwd=sys.argv[3])" "$dir/.driver.sh" "$log" "$root"
echo "screen '$name' launched detached"
echo "  log:     $log"
echo "  results: $dir/results.csv"
