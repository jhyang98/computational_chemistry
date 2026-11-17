#!/usr/bin/env bash
# TPS for the finite-size / observation-time scaling of the interim report.
#
# For every combination of eps, g, L and tobs below, the s windows in S_LIST
# are run from low to high s, each window starting from the last trajectory of
# the previous one, and the windows are then combined with MBAR.
#
# Run from the project root:
#   bash scripts/run_tps_scaling.sh                      # all combinations
#   bash scripts/run_tps_scaling.sh <eps> <g> <L> <tobs> # one combination
#
# Each combination writes only to its own folder
#   <OUT_ROOT>/eps<eps>_g<g>/L<L>_t<tobs>/
# so single combinations can run in parallel (separate shells or batch jobs).
#
# Hysteresis check: list S_LIST from high to low s, set INIT_FIRST=empty (the
# first window then starts from an inactive trajectory) and use another
# OUT_ROOT so the results of the ascending run are kept.
set -euo pipefail

# ---- parameters to set before running -------------------------------------
EPS_LIST="0.01"
G_LIST="0.0"
L_LIST="25 50 100"
TOBS_LIST="200 400 600"
# EXAMPLE windows only: choose them around s* from the ED results and a pilot
# run. Neighbouring windows must overlap (see README, section 4).
S_LIST="0.010 0.015 0.020 0.025 0.030 0.035 0.040"
S_TARGET="0.010:0.0005:0.040"   # MBAR evaluation grid
T=1
N_EQUIL=20000                   # shooting moves before sampling, per window
N_PROD=1000000                  # shooting moves while sampling, per window
EVERY=10                        # record K, E every EVERY moves
N_BOOT=100                      # MBAR bootstrap replicas (0: no error bars)
BLOCK=100                       # bootstrap block length in samples
INIT_FIRST=equilibrium          # first window starts from: equilibrium, empty,
                                # or a saved trajectory file
OUT_ROOT=tps                    # results go to OUT_ROOT/eps.._g../L.._t../
# ---------------------------------------------------------------------------

run_combo() {
    local EPS=$1 G=$2 L=$3 TOBS=$4
    local DIR INIT TAG S
    DIR=$(printf "%s/eps%s_g%+.3f/L%d_t%d" "$OUT_ROOT" "$EPS" "$G" "$L" "$TOBS")
    mkdir -p "$DIR"
    INIT=$INIT_FIRST
    for S in $S_LIST; do
        TAG=$(printf "s%+.4f" "$S")
        cat > "$DIR/$TAG.inp" <<EOT
L            $L
T            $T
eps          $EPS
s            $S
g            $G
tobs         $TOBS
seed         0
n_equil      $N_EQUIL
n_prod       $N_PROD
sample_every $EVERY
init         $INIT
save         $DIR/$TAG.traj
output       $DIR/$TAG.dat
EOT
        echo "[$(date '+%H:%M:%S')] eps=$EPS g=$G L=$L tobs=$TOBS s=$S"
        ./bin/east_tps "$DIR/$TAG.inp" > "$DIR/$TAG.log"
        INIT="$DIR/$TAG.traj"
    done
    {
        for S in $S_LIST; do echo "window $DIR/$(printf "s%+.4f" "$S").dat"; done
        echo "target_s $S_TARGET"
        echo "target_g $G"
        echo "n_boot   $N_BOOT"
        echo "block    $BLOCK"
        echo "seed     1"
        echo "output   $DIR/mbar.dat"
    } > "$DIR/mbar.inp"
    ./bin/east_mbar "$DIR/mbar.inp"
    echo "[$(date '+%H:%M:%S')] done: $DIR"
}

if [ $# -eq 4 ]; then
    run_combo "$1" "$2" "$3" "$4"
    exit 0
elif [ $# -ne 0 ]; then
    echo "usage: bash scripts/run_tps_scaling.sh [<eps> <g> <L> <tobs>]" >&2
    exit 1
fi

for EPS in $EPS_LIST; do
    for G in $G_LIST; do
        for L in $L_LIST; do
            for TOBS in $TOBS_LIST; do
                run_combo "$EPS" "$G" "$L" "$TOBS"
            done
        done
    done
done
