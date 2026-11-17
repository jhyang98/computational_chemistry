#!/usr/bin/env bash
# Additional ED peak runs for the hypothesis analysis (T = 1).
#   step 1: finer eps grid at L = 8..14, same s range as results/peak_T1.dat
#           -> results/peak_T1_fine.dat
#   step 2: L = 16 for every (eps, g)          -> results/peak_T1_L16.dat
#   step 3: L = 18 for selected (eps, g)       -> results/peak_T1_L18.dat
#   step 4: strict line at L = 16, 18          -> results/strict_T1_L16_18.dat
# For L >= 16 the full s range is too slow, so the coarse s grid is a window
# around the L = 14 s* of the same (eps, g); points that are nan at L = 14
# (chi_KK maximum at the edge of the full s range) are skipped.
#
# Run from the project root:
#   bash scripts/run_ed_fine.sh            # all steps
#   bash scripts/run_ed_fine.sh 2          # one step (1, 2, 3 or 4)
#   bash scripts/run_ed_fine.sh merge      # rebuild the three tables from the job files
# Independent jobs run in parallel: NPROC=4 bash scripts/run_ed_fine.sh
set -euo pipefail

# ---- parameters ------------------------------------------------------------
NPROC=${NPROC:-$(nproc)}
EPS_FINE="0.0005 0.0015 0.002 0.004 0.005 0.007 0.015 0.02"
G_LIST="-0.15:0.05:0.15"
S_FULL="-0.5:0.005:0.8"
EPS_L18="0.0003 0.001 0.003 0.01 0.03"
G_L18="-0.15 0 0.15"
WIN_LO=0.02          # window: [s*(14) - WIN_LO, s*(14) + WIN_HI]
WIN_HI=0.03
WIN_STEP=0.001
# ------------------------------------------------------------------------------

TMP=results/tmp_ed_fine
mkdir -p "$TMP"

merge() { # merge job files $3.. into $1 (header of the first job, $2 = note), rows sorted
    local out=$1 note=$2
    shift 2
    {
        awk -v note="$note" '/^#/ {
            if ($0 ~ /^# command:/) print "# merged by scripts/run_ed_fine.sh (one east_ed run per job)"
            else if ($0 ~ /^#   L /) hdr = $0
            else print }
            END { if (note != "") print "# " note; print hdr }' "$1"
        cat "$@" | grep -v '^#' | sort -k1,1n -k2,2g -k3,3g
    } > "$out"
}

s14() { # s* at L = 14 for eps $1, g $2 from the stored tables (nan if absent)
    awk -v e="$1" -v g="$2" '!/^#/ && $1 == 14 && $2 + 0 == e + 0 && $3 + 0 == g + 0 { print $4; f = 1; exit }
        END { if (!f) print "nan" }' results/peak_T1.dat results/peak_T1_fine.dat 2>/dev/null
}

window_jobs() { # print one east_ed command per (eps, g) for size $1
    local L=$1 e g s
    shift 1
    for e in $EPS_LIST_W; do
        for g in $G_LIST_W; do
            s=$(s14 "$e" "$g")
            if [ "$s" = "nan" ]; then
                echo "skip L=$L eps=$e g=$g (nan at L=14)" >&2
                continue
            fi
            echo "./bin/east_ed peak -L $L --eps $e --g $g --s $(awk -v s="$s" -v a="$WIN_LO" \
                -v b="$WIN_HI" -v d="$WIN_STEP" 'BEGIN { printf "%.6f:%g:%.6f", s - a, d, s + b }') \
                --out $TMP/L${L}_e${e}_g${g}.dat"
        done
    done
}

step1() {
    echo "[$(date '+%H:%M:%S')] step 1: L = 8..14, eps = $EPS_FINE"
    for e in $EPS_FINE; do
        echo "./bin/east_ed peak -L 8,10,12,14 --eps $e --g $G_LIST --s $S_FULL --out $TMP/fine_e$e.dat"
    done | xargs -P "$NPROC" -I{} bash -c {} 2> "$TMP/step1.err"
    merge results/peak_T1_fine.dat "s grid $S_FULL" $TMP/fine_e*.dat
}

step2() {
    echo "[$(date '+%H:%M:%S')] step 2: L = 16"
    rm -f $TMP/L16_*.dat
    EPS_LIST_W=$(awk '!/^#/ { print $2 }' results/peak_T1.dat results/peak_T1_fine.dat | sort -gu | tr '\n' ' ')
    G_LIST_W=$(awk '!/^#/ { print $3 }' results/peak_T1.dat | sort -gu | tr '\n' ' ')
    EPS_LIST_W=$EPS_LIST_W G_LIST_W=$G_LIST_W window_jobs 16 \
        | xargs -P "$NPROC" -I{} bash -c {} 2> "$TMP/step2.err"
    merge results/peak_T1_L16.dat "s grid: window [s*(L=14) - $WIN_LO, s*(L=14) + $WIN_HI], step $WIN_STEP" $TMP/L16_*.dat
}

step3() {
    echo "[$(date '+%H:%M:%S')] step 3: L = 18"
    rm -f $TMP/L18_*.dat
    EPS_LIST_W=$EPS_L18 G_LIST_W=$G_L18 window_jobs 18 \
        | xargs -P "$NPROC" -I{} bash -c {} 2> "$TMP/step3.err"
    merge results/peak_T1_L18.dat "s grid: window [s*(L=14) - $WIN_LO, s*(L=14) + $WIN_HI], step $WIN_STEP" $TMP/L18_*.dat
}

step4() {
    echo "[$(date '+%H:%M:%S')] step 4: strict line, L = 16, 18"
    for L in 16 18; do
        echo "./bin/east_ed strict -L $L --g $G_LIST --out $TMP/strict_L$L.dat"
    done | xargs -P "$NPROC" -I{} bash -c {} 2> "$TMP/step4.err"
    merge results/strict_T1_L16_18.dat "" $TMP/strict_L*.dat
}

merge_all() {
    local w="s grid: window [s*(L=14) - $WIN_LO, s*(L=14) + $WIN_HI], step $WIN_STEP"
    merge results/peak_T1_fine.dat "s grid $S_FULL" $TMP/fine_e*.dat
    merge results/peak_T1_L16.dat "$w" $TMP/L16_*.dat
    merge results/peak_T1_L18.dat "$w" $TMP/L18_*.dat
    merge results/strict_T1_L16_18.dat "" $TMP/strict_L*.dat
}

case "${1:-all}" in
    all) step1; step2; step3; step4 ;;
    merge) merge_all ;;
    1) step1 ;;
    2) step2 ;;
    3) step3 ;;
    4) step4 ;;
    *) echo "usage: bash scripts/run_ed_fine.sh [all | 1 | 2 | 3 | 4 | merge]" >&2; exit 1 ;;
esac
echo "[$(date '+%H:%M:%S')] done"
