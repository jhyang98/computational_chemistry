#!/usr/bin/env bash
# Replica-exchange umbrella sampling (east_reus) in K at fixed (eps, g, L,
# tobs, s), followed by MBAR and the equal-weight point.
# Windows: [KMIN_ALL + i*D, KMIN_ALL + (i+2)*D], D = (KMAX_ALL-KMIN_ALL)/(NWIN+1).
#
# Run from the project root:
#   bash scripts/run_reus.sh          # REUS run (or resume), then MBAR
#   bash scripts/run_reus.sh mbar     # MBAR only (after a finished run)
#   bash scripts/run_reus.sh params   # print the parameters as KEY=VALUE lines
#
# Every parameter below can be overridden from the environment, e.g.
#   EPS=0.003 G=0.05 bash scripts/run_reus.sh
# On a SLURM cluster submit with scripts/submit_reus.sh instead.
#
# Checkpoints: the REUS state is saved in <dir>/reus.ckpt every
# CHECKPOINT_EVERY seconds and when the run stops early (STOP_AFTER seconds,
# SIGTERM or SIGUSR1). Running the same command again resumes from it. Exit
# status 3 means "stopped with a checkpoint, run again". A finished run leaves
# <dir>/reus.done (its parameters); a later run with the same parameters skips
# REUS and only redoes MBAR, one with different parameters stops with an error.
#
# Parallel windows: build with "make OMP=1" and set OMP_NUM_THREADS, e.g.
#   OMP_NUM_THREADS=4 bash scripts/run_reus.sh
set -euo pipefail

# ---- parameters (defaults; environment variables of the same name win) -----
EPS=${EPS:-0.01}
G=${G:-0.0}
L=${L:-100}
TOBS=${TOBS:-600}
S=${S:-0.018}                       # s of the runs; near the expected s* (not critical)
KMIN_ALL=${KMIN_ALL:-0.0}           # k range covered by the windows: must contain both peaks
KMAX_ALL=${KMAX_ALL:-0.11}
NWIN=${NWIN:-30}
# starting trajectories: empty, equilibrium, or a file pattern with %02d for the
# window index (missing files start empty), e.g. the .traj files of an earlier
# umbrella run: INIT="umb/eps0.01_g+0.000/L100_t600_s+0.0250/w%02d.traj"
INIT=${INIT:-empty}
T=${T:-1}
N_APPROACH=${N_APPROACH:-1000000}   # maximum moves for a replica to enter its window, per attempt
N_EQUIL=${N_EQUIL:-20000}           # moves per window before sampling (with swaps)
N_PROD=${N_PROD:-1000000}           # sampled moves per window
SWAP_EVERY=${SWAP_EVERY:-10}        # moves between swap attempts
EVERY=${EVERY:-10}                  # moves between samples (multiple of SWAP_EVERY)
PROGRESS=${PROGRESS:-10000}         # progress report every this many cycles (0: off)
CHECKPOINT_EVERY=${CHECKPOINT_EVERY:-1800}  # seconds of wall time between checkpoints
STOP_AFTER=${STOP_AFTER:-0}         # stop with a checkpoint after this many seconds (0: never)
# MBAR / analysis
S_TARGET=${S_TARGET:-0.010:0.00025:0.035}   # reweighting grid in s (also the s* search range)
KCUT=${KCUT:-auto}                  # auto, a value of k, or none
HIST_S=${HIST_S:-0.015,0.0175,0.020}        # s values for reweighted P(k) output
HIST_BINS=${HIST_BINS:-60}
SMOOTH_RES=${SMOOTH_RES:-30}        # smoothing resolution m: window ~ 1/m of the k range
N_BOOT=${N_BOOT:-30}
BLOCK=${BLOCK:-100}
OUT_ROOT=${OUT_ROOT:-reus}
# ---------------------------------------------------------------------------
PARAM_NAMES="EPS G L TOBS S KMIN_ALL KMAX_ALL NWIN INIT T N_APPROACH N_EQUIL N_PROD SWAP_EVERY
EVERY PROGRESS CHECKPOINT_EVERY STOP_AFTER S_TARGET KCUT HIST_S HIST_BINS SMOOTH_RES N_BOOT BLOCK
OUT_ROOT"

DIR=$(printf "%s/eps%s_g%+.3f/L%d_t%d_s%+.4f" "$OUT_ROOT" "$EPS" "$G" "$L" "$TOBS" "$S")

reus_params() { # the parameters that define the REUS samples
    cat <<EOT
L             $L
T             $T
eps           $EPS
s             $S
g             $G
tobs          $TOBS
nwin          $NWIN
k_lo          $KMIN_ALL
k_hi          $KMAX_ALL
init          $INIT
n_approach    $N_APPROACH
n_equil       $N_EQUIL
n_prod        $N_PROD
swap_every    $SWAP_EVERY
sample_every  $EVERY
EOT
}

run_reus() {
    local rc=0
    mkdir -p "$DIR"
    if [ -f "$DIR/reus.done" ]; then
        if reus_params | cmp -s - "$DIR/reus.done"; then
            echo "[$(date '+%F %T')] REUS already finished in $DIR; skipping to MBAR"
            return 0
        fi
        echo "error: $DIR holds a finished REUS run with other parameters (see $DIR/reus.done);" >&2
        echo "       set OUT_ROOT to another folder, or remove $DIR/reus.done to run again" >&2
        exit 1
    fi
    {
        reus_params
        echo "seed          0"
        echo "progress      $PROGRESS"
        echo "checkpoint    $DIR/reus.ckpt"
        echo "checkpoint_every $CHECKPOINT_EVERY"
        echo "stop_after    $STOP_AFTER"
        echo "output_prefix $DIR/"
    } > "$DIR/reus.inp"
    if [ -f "$DIR/reus.ckpt" ]; then
        echo "[$(date '+%F %T')] REUS resume from checkpoint: $DIR"
    else
        echo "[$(date '+%F %T')] REUS start: $DIR"
    fi
    ./bin/east_reus "$DIR/reus.inp" >> "$DIR/reus.log" || rc=$?
    if [ "$rc" -eq 3 ]; then
        echo "[$(date '+%F %T')] REUS stopped with a checkpoint; run the same command again to resume"
        exit 3
    elif [ "$rc" -ne 0 ]; then
        echo "[$(date '+%F %T')] REUS failed (exit status $rc)" >&2
        exit "$rc"
    fi
    reus_params > "$DIR/reus.done"
    echo "[$(date '+%F %T')] REUS done: summary in $DIR/summary.dat"
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
    echo "[$(date '+%F %T')] MBAR done: $DIR/mbar.dat, $DIR/mbar_hist.dat"
    grep -E "s_equal|^# peaks " "$DIR/mbar.dat" || true
}

case "${1:-all}" in
    all)
        run_reus
        run_mbar
        ;;
    mbar)
        run_mbar
        ;;
    params)
        for k in $PARAM_NAMES; do printf '%s=%q\n' "$k" "${!k}"; done
        echo "# output folder: $DIR"
        ;;
    *)
        echo "usage: bash scripts/run_reus.sh [all | mbar | params]" >&2
        exit 1
        ;;
esac
