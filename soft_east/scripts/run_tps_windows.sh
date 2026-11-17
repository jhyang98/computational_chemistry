#!/usr/bin/env bash
# TPS over a set of s windows at fixed g, chained as in the earlier East-model
# code: each window starts from the last trajectory of the previous window.
# The windows are then combined with MBAR.
# Run from the project root:  bash scripts/run_tps_windows.sh
set -euo pipefail

L=100
T=1
EPS=0.01
G=0.0
TOBS=400
S_LIST="0.000 0.010 0.020 0.030 0.040"
N_EQUIL=20000
N_PROD=1000000
EVERY=10

mkdir -p tps
INIT=equilibrium
for S in $S_LIST; do
    TAG=$(printf "s%+.3f_g%+.3f" "$S" "$G")
    cat > system.inp <<EOT
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
save         tps/$TAG.traj
output       tps/$TAG.dat
EOT
    echo "== window s=$S g=$G (init: $INIT)"
    ./bin/east_tps
    INIT=tps/$TAG.traj
done

{
    for S in $S_LIST; do echo "window tps/$(printf "s%+.3f_g%+.3f" "$S" "$G").dat"; done
    echo "target_s 0.0:0.001:0.04"
    echo "target_g $G"
    echo "n_boot   100"
    echo "block    100"
    echo "seed     1"
    echo "output   tps/mbar_g$(printf "%+.3f" "$G").dat"
} > tps/mbar.inp
./bin/east_mbar tps/mbar.inp
echo "done: tps/*.dat, MBAR result in tps/mbar_g$(printf "%+.3f" "$G").dat"
