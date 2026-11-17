#!/usr/bin/env bash
# SLURM batch script that runs scripts/run_reus.sh. Submit it with
# scripts/submit_reus.sh, which sets the node, cores and time limit and
# freezes the parameters in the file named by JOB_ENV.
#
# Time limit: REUS stops by itself STOP_AFTER seconds into the run and leaves
# a checkpoint (run_reus.sh exits with status 3); this script then requeues
# the job and the next run resumes from the checkpoint. If the run is still
# busy when SLURM sends SIGUSR1 shortly before the limit (for example in
# MBAR), east_reus is asked to write its checkpoint and the job is requeued.
# A requeued job keeps its job ID and appends to the same log file.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
if [ -z "${JOB_ENV:-}" ] || [ ! -f "$JOB_ENV" ]; then
    echo "error: JOB_ENV is unset or missing; submit with scripts/submit_reus.sh" >&2
    exit 1
fi
set -a
# shellcheck disable=SC1090
. "$JOB_ENV"
set +a
MODE=${MODE:-all}
MAX_REQUEUE=${MAX_REQUEUE:-10}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

stamp() { date '+%F %T'; }

requeue() {
    local n=${SLURM_RESTART_COUNT:-0}
    if [ -z "${SLURM_JOB_ID:-}" ]; then
        echo "[$(stamp)] not under SLURM: run again to resume"
        exit 3
    fi
    if [ "$n" -ge "$MAX_REQUEUE" ]; then
        echo "[$(stamp)] not requeued: already requeued $n times (MAX_REQUEUE = $MAX_REQUEUE)"
        exit 1
    fi
    echo "[$(stamp)] requeue job $SLURM_JOB_ID (restart $((n + 1)) of at most $MAX_REQUEUE)"
    if ! scontrol requeue "$SLURM_JOB_ID"; then
        echo "[$(stamp)] error: scontrol requeue failed" >&2
        exit 1
    fi
    exit 0
}

on_usr1() {
    echo "[$(stamp)] time limit near"
    # east_reus writes a checkpoint and stops at the next cycle on SIGTERM
    if pkill -TERM -P "$RUN_PID" -x east_reus 2> /dev/null; then
        for _ in $(seq 1 180); do
            pgrep -P "$RUN_PID" -x east_reus > /dev/null || break
            sleep 1
        done
    fi
    requeue
}
trap on_usr1 USR1

echo "[$(stamp)] job ${SLURM_JOB_ID:-local} on $(hostname), restart ${SLURM_RESTART_COUNT:-0}," \
    "$OMP_NUM_THREADS threads, mode $MODE, parameters $JOB_ENV"
bash scripts/run_reus.sh "$MODE" &
RUN_PID=$!
wait "$RUN_PID"
rc=$?
if [ "$rc" -eq 3 ]; then requeue; fi
echo "[$(stamp)] finished with exit status $rc"
exit "$rc"
