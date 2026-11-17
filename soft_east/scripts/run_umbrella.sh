#!/usr/bin/env bash
# Umbrella sampling in the activity K at fixed (eps, g, L, tobs, s), combined
# with MBAR. The windows cover k = K/(L tobs) from KMIN_ALL to KMAX_ALL; window
# i spans [KMIN_ALL + i*D, KMIN_ALL + (i+2)*D] with D = (KMAX_ALL-KMIN_ALL)/(NWIN+1),
# so neighbouring windows overlap by half a window. Every window starts from
# INIT and first moves into its window (approach phase), so windows are
# independent and can run in parallel.
#
# Run from the project root:
#   bash scripts/run_umbrella.sh              # all windows in turn, then MBAR
#   bash scripts/run_umbrella.sh window <i>   # window i only (0 .. NWIN-1)
#   bash scripts/run_umbrella.sh mbar         # MBAR over the finished windows
#   bash scripts/run_umbrella.sh list         # print the windows and exit
set -euo pipefail

# ---- parameters to set before running -------------------------------------
EPS=0.01
G=0.0
L=100
TOBS=600
S=0.025                 # s of the umbrella runs; near the expected s* (not critical)
KMIN_ALL=0.0            # k range covered by the windows: must contain both peaks
KMAX_ALL=0.08
NWIN=22
INIT=empty              # starting trajectory: empty, equilibrium or a file
T=1
N_APPROACH=1000000      # maximum moves to enter a window
N_EQUIL=20000           # moves inside the window before sampling
N_PROD=1000000          # sampled moves per window
EVERY=10
# MBAR / analysis
S_TARGET="0.015:0.0005:0.040"   # reweighting grid in s (also the s* search range)
KCUT=auto               # auto, a value of k, or none
HIST_S="0.020,0.025,0.030"      # s values for reweighted P(k) output
HIST_BINS=60
SMOOTH_RES=30           # smoothing resolution m: window ~ 1/m of the k range
N_BOOT=100
BLOCK=100
OUT_ROOT=umb
# ---------------------------------------------------------------------------

DIR=$(printf "%s/eps%s_g%+.3f/L%d_t%d_s%+.4f" "$OUT_ROOT" "$EPS" "$G" "$L" "$TOBS" "$S")

window_bounds() { # prints "kmin kmax" of window $1
    awk -v a="$KMIN_ALL" -v b="$KMAX_ALL" -v n="$NWIN" -v i="$1" \
        'BEGIN { d = (b - a) / (n + 1); printf "%.10g %.10g\n", a + i * d, a + (i + 2) * d }'
}

run_window() {
    local i=$1 TAG KLO KHI
    TAG=$(printf "w%02d" "$i")
    read -r KLO KHI <<< "$(window_bounds "$i")"
    mkdir -p "$DIR"
    cat > "$DIR/$TAG.inp" <<EOT
L            $L
T            $T
eps          $EPS
s            $S
g            $G
tobs         $TOBS
seed         0
init         $INIT
umb_kmin     $KLO
umb_kmax     $KHI
n_approach   $N_APPROACH
n_equil      $N_EQUIL
n_prod       $N_PROD
sample_every $EVERY
save         $DIR/$TAG.traj
output       $DIR/$TAG.dat
EOT
    echo "[$(date '+%H:%M:%S')] window $i: k in [$KLO, $KHI]"
    ./bin/east_tps "$DIR/$TAG.inp" > "$DIR/$TAG.log"
}

run_mbar() {
    local i
    {
        for ((i = 0; i < NWIN; i++)); do echo "window $DIR/$(printf "w%02d" "$i").dat"; done
        echo "target_s    $S_TARGET"
        echo "target_g    $G"
        if [ "$KCUT" != "none" ]; then echo "kcut        $KCUT"; fi
        echo "hist_s      $HIST_S"
        echo "hist_bins   $HIST_BINS"
        echo "smooth_res  $SMOOTH_RES"
        echo "hist_output $DIR/mbar_hist.dat"
        echo "n_boot      $N_BOOT"
        echo "block       $BLOCK"
        echo "seed        1"
        echo "output      $DIR/mbar.dat"
    } > "$DIR/mbar.inp"
    ./bin/east_mbar "$DIR/mbar.inp"
    echo "[$(date '+%H:%M:%S')] MBAR done: $DIR/mbar.dat, $DIR/mbar_hist.dat"
    grep -E "s_equal|^# peaks " "$DIR/mbar.dat" || true
}

case "${1:-all}" in
    all)
        for ((i = 0; i < NWIN; i++)); do run_window "$i"; done
        run_mbar
        ;;
    window)
        if [ $# -ne 2 ] || [ "$2" -lt 0 ] || [ "$2" -ge "$NWIN" ]; then
            echo "usage: bash scripts/run_umbrella.sh window <0..$((NWIN - 1))>" >&2
            exit 1
        fi
        run_window "$2"
        ;;
    mbar)
        run_mbar
        ;;
    list)
        for ((i = 0; i < NWIN; i++)); do echo "window $i: k in [$(window_bounds "$i")]"; done
        ;;
    *)
        echo "usage: bash scripts/run_umbrella.sh [all | window <i> | mbar | list]" >&2
        exit 1
        ;;
esac
