#!/bin/bash
# Sequential, resumable batch runner for the 100 complete hyperonic EoSs.

set -uo pipefail

EXE="./logtov_seq_geom_tidal.out"
EOS_DIR="eos/hyperonic_dataset/hyperonic_100/complete_tables"
OUT_DIR="results/hyperonic_100/coarse"
LOG_DIR="results/hyperonic_100/logs"
STATUS="results/hyperonic_100/batch_status.csv"
UNRESOLVED="results/hyperonic_100/unresolved.txt"

RHO0="1.e14"
DRHO="0.05e14"
N_FIRST=400
N_EXTENDED=600

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [ ! -x "$EXE" ]; then
    echo "ERROR: executable not found: $EXE"
    exit 1
fi

COUNT=$(find "$EOS_DIR" -maxdepth 1 -name '*_complete.dat' | wc -l | tr -d ' ')
if [ "$COUNT" -ne 100 ]; then
    echo "ERROR: expected 100 complete EoSs, found $COUNT"
    exit 1
fi

# Preserve the user's current infile.
HAD_INFILE=0
if [ -f infile ]; then
    cp infile /tmp/infile.hyperonic100.backup
    HAD_INFILE=1
fi

restore_infile() {
    if [ "$HAD_INFILE" -eq 1 ] && [ -f /tmp/infile.hyperonic100.backup ]; then
        cp /tmp/infile.hyperonic100.backup infile
    fi
}
trap restore_infile EXIT INT TERM

echo "eos,status,nmodels,rows,mmax,index_max" > "$STATUS"
: > "$UNRESOLVED"

resolved_check () {
python3 - "$1" <<'PY'
import sys, numpy as np
p = sys.argv[1]
d = np.loadtxt(p)
if d.ndim == 1:
    d = d.reshape(1, -1)
M = d[:,2]
i = int(np.nanargmax(M))
resolved = i < len(M) - 2
print(f"{len(M)},{M[i]:.10g},{i}")
sys.exit(0 if resolved else 1)
PY
}

run_one () {
    eos="$1"
    nmodels="$2"
    outfile="$3"
    logfile="$4"

    cat > infile <<EOF
$RHO0 $DRHO $nmodels
$eos
EOF

    rm -f logtov_seq_geom_tidal.dat

    {
        echo "EoS: $eos"
        echo "nmodels requested: $nmodels"
        echo "start: $(date)"
        /usr/bin/time -p "$EXE"
        rc=$?
        echo "solver exit code: $rc"
        echo "end: $(date)"
    } > "$logfile" 2>&1

    if [ "$rc" -ne 0 ]; then
        return "$rc"
    fi

    if [ ! -s logtov_seq_geom_tidal.dat ]; then
        return 1
    fi

    cp logtov_seq_geom_tidal.dat "$outfile"
    return 0
}

i=0
for eos in "$EOS_DIR"/*_complete.dat; do
    i=$((i+1))
    stem=$(basename "$eos" .dat)
    outfile="$OUT_DIR/${stem}_TOV.dat"
    log1="$LOG_DIR/${stem}_400.log"
    log2="$LOG_DIR/${stem}_600.log"

    echo "============================================================"
    echo "[$i/100] $stem"

    # Resume safely.
    if [ -s "$outfile" ]; then
        info=$(resolved_check "$outfile" 2>/dev/null)
        if [ $? -eq 0 ]; then
            echo "Existing resolved output found; skipping."
            echo "$stem,skipped_existing,unknown,$info" >> "$STATUS"
            continue
        fi
    fi

    echo "Running initial $N_FIRST-model scan..."
    if ! run_one "$eos" "$N_FIRST" "$outfile" "$log1"; then
        echo "FAILED initial run."
        echo "$stem,failed,$N_FIRST,0,nan,nan" >> "$STATUS"
        continue
    fi

    info=$(resolved_check "$outfile")
    if [ $? -eq 0 ]; then
        echo "Turning point resolved."
        echo "$stem,resolved,$N_FIRST,$info" >> "$STATUS"
        continue
    fi

    echo "Maximum not securely resolved; extending to $N_EXTENDED models..."
    if ! run_one "$eos" "$N_EXTENDED" "$outfile" "$log2"; then
        echo "FAILED extended run."
        echo "$stem,failed_extended,$N_EXTENDED,0,nan,nan" >> "$STATUS"
        continue
    fi

    info=$(resolved_check "$outfile")
    if [ $? -eq 0 ]; then
        echo "Turning point resolved in extended run."
        echo "$stem,resolved_extended,$N_EXTENDED,$info" >> "$STATUS"
    else
        echo "WARNING: maximum still unresolved."
        echo "$stem" >> "$UNRESOLVED"
        echo "$stem,unresolved,$N_EXTENDED,$info" >> "$STATUS"
    fi
done

echo
echo "Batch complete."
echo "Outputs    : $OUT_DIR"
echo "Logs       : $LOG_DIR"
echo "Status CSV : $STATUS"
echo "Unresolved : $UNRESOLVED"
