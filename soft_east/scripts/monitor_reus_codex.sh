#!/usr/bin/env bash
# Monitor one REUS SLURM job and invoke Codex only on a health event or completion.
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "usage: bash scripts/monitor_reus_codex.sh JOB_ID JOB_ENV" >&2
    exit 1
fi

JOB_ID=$1
JOB_ENV=$2
POLL_SECONDS=${POLL_SECONDS:-300}
STALE_SECONDS=${STALE_SECONDS:-5400}
CODEX_RETRIES=${CODEX_RETRIES:-3}
SYNC_RETRIES=${SYNC_RETRIES:-5}
LOG_DIR=monitor_logs

cd "$(dirname "$0")/.."
mkdir -p "$LOG_DIR"
MONITOR_LOG="$LOG_DIR/reus_${JOB_ID}.log"
CODEX_LOG="$LOG_DIR/reus_${JOB_ID}_codex.log"

stamp() { date '+%F %T'; }

record() {
    printf '[%s] %s\n' "$(stamp)" "$*" | tee -a "$MONITOR_LOG"
}

job_line() {
    squeue -h -j "$JOB_ID" -o '%T|%M|%R' 2>/dev/null || return 2
}

output_dir() {
    awk -v id="$JOB_ID" '$1 == id { print $5; exit }' slurm_jobs/jobs.txt
}

run_codex() {
    local purpose=$1 prompt=$2 attempt rc=1
    for ((attempt = 1; attempt <= CODEX_RETRIES; attempt++)); do
        record "starting Codex $purpose attempt $attempt of $CODEX_RETRIES"
        if codex exec -C . --skip-git-repo-check --approve-for-me \
            --add-dir /tmp \
            -o "$LOG_DIR/reus_${JOB_ID}_${purpose}_last.txt" "$prompt" \
            >> "$CODEX_LOG" 2>&1; then
            record "Codex $purpose completed"
            return 0
        else
            rc=$?
        fi
        record "Codex $purpose failed with status $rc"
        sleep 60
    done
    return "$rc"
}

continue_after_requeue_limit() {
    local submit_output new_job new_env
    record "automatic requeue limit reached; submitting from the saved parameters"
    submit_output=$(bash scripts/submit_reus.sh --from "$JOB_ENV" 2>&1) || {
        printf '%s\n' "$submit_output" >> "$MONITOR_LOG"
        return 1
    }
    printf '%s\n' "$submit_output" >> "$MONITOR_LOG"
    new_job=$(printf '%s\n' "$submit_output" | sed -n 's/^submitted  : job \([0-9][0-9]*\).*/\1/p')
    new_env=$(printf '%s\n' "$submit_output" | sed -n 's/^parameters : //p')
    if [ -z "$new_job" ] || [ -z "$new_env" ]; then
        record "could not parse the continuation submission"
        return 1
    fi
    JOB_ID=$new_job
    JOB_ENV=$new_env
    MONITOR_LOG="$LOG_DIR/reus_${JOB_ID}.log"
    CODEX_LOG="$LOG_DIR/reus_${JOB_ID}_codex.log"
    record "monitoring continuation job $JOB_ID"
}

find_replacement_job() {
    local out=$1 line new_job new_env
    line=$(awk -v out="$out" '$5 == out { line=$0 } END { print line }' slurm_jobs/jobs.txt)
    new_job=$(printf '%s\n' "$line" | awk '{ print $1 }')
    new_env=$(printf '%s\n' "$line" | awk '{ print $6 }')
    if [ -n "$new_job" ] && [ "$new_job" != "$JOB_ID" ] && [ -n "$new_env" ]; then
        JOB_ID=$new_job
        JOB_ENV=$new_env
        MONITOR_LOG="$LOG_DIR/reus_${JOB_ID}.log"
        CODEX_LOG="$LOG_DIR/reus_${JOB_ID}_codex.log"
        record "Codex submitted replacement job $JOB_ID"
        return 0
    fi
    return 1
}

finalize_job() {
    local out=$1 prompt attempt
    prompt=$(cat <<EOF
Work autonomously on completed soft_east SLURM job $JOB_ID.
Read AGENTS.md, STATUS.md, and README.md first. The frozen parameter file is $JOB_ENV,
the scheduler log is slurm_logs/soft_east_${JOB_ID}.out, and the result directory is $out.
Follow the completed-run checklist in STATUS.md exactly. Validate summary.dat, mbar.dat,
and mbar_hist.dat; rerun only the same-parameter MBAR analysis if needed. Read the frozen
parameters to identify the campaign condition and compare it with the scientifically relevant
completed runs. Interpret whether P(k) has one or two peaks. At the chi_KK maximum, report
chi_KK, its s location, the width of P(k), and Binder U when the available samples or histogram
support it. Also report s_equal, delta k, delta e, round trips, worst swap acceptance, worst KS,
minimum N_eff, elapsed time, requeue behavior, and any caveat. Update STATUS.md and README.md
in Korean with stable results and next-step implications. Do not submit a new scientific
condition. Do not perform any Git operation; the monitor synchronizes GitHub after you finish.
Keep code, comments, and command output in English.
EOF
)
    run_codex finalize "$prompt" || return 1
    for ((attempt = 1; attempt <= SYNC_RETRIES; attempt++)); do
        record "starting GitHub sync attempt $attempt of $SYNC_RETRIES"
        if bash scripts/sync_github_soft_east.sh "Complete REUS job $JOB_ID analysis" \
            >> "$MONITOR_LOG" 2>&1; then
            record "GitHub sync completed"
            touch "$LOG_DIR/reus_${JOB_ID}.finalized"
            return 0
        fi
        record "GitHub sync attempt $attempt failed"
        sleep 120
    done
    return 1
}

diagnose_job() {
    local out=$1 reason=$2 prompt
    prompt=$(cat <<EOF
Diagnose soft_east SLURM job $JOB_ID autonomously because the monitor observed: $reason.
Read AGENTS.md, STATUS.md, and README.md first. Inspect $JOB_ENV,
slurm_logs/soft_east_${JOB_ID}.out, squeue/scontrol, and $out. Follow the documented recovery
rules. You may resume the same frozen parameters, rerun same-parameter MBAR, fix a clear
implementation or operational error, and submit a documented recoverable continuation without
asking. Never change the scientific parameters or delete results. If a new scientific decision
is required, stop without submitting it and record the blocking reason in STATUS.md. If you
submit a replacement job, use bash scripts/submit_reus.sh and keep the same output directory.
Do not perform any Git operation; the monitor handles synchronization after final completion.
Keep code, comments, and command output in English.
EOF
)
    run_codex diagnose "$prompt"
}

if [ ! -f "$JOB_ENV" ]; then
    echo "error: parameter file not found: $JOB_ENV" >&2
    exit 1
fi

OUT_DIR=$(output_dir)
if [ -z "$OUT_DIR" ]; then
    echo "error: job $JOB_ID is not listed in slurm_jobs/jobs.txt" >&2
    exit 1
fi

record "monitor started for job $JOB_ID, output $OUT_DIR"
missing_checks=0
stale_marker=0

while :; do
    state=$(job_line) || {
        record "squeue query failed; retrying"
        sleep "$POLL_SECONDS"
        continue
    }
    log_file="slurm_logs/soft_east_${JOB_ID}.out"
    if [ -n "$state" ]; then
        missing_checks=0
        if [ -f "$log_file" ]; then
            now=$(date +%s)
            modified=$(stat -c %Y "$log_file")
            if [ $((now - modified)) -ge "$STALE_SECONDS" ] && [ "$modified" -ne "$stale_marker" ]; then
                stale_marker=$modified
                record "job is active but its log is stale: $state"
                diagnose_job "$OUT_DIR" "the active job log has not changed for at least $STALE_SECONDS seconds" || true
                find_replacement_job "$OUT_DIR" || true
            fi
        fi
        sleep "$POLL_SECONDS"
        continue
    fi

    missing_checks=$((missing_checks + 1))
    if [ "$missing_checks" -lt 4 ]; then
        sleep 60
        continue
    fi

    if [ -f "$log_file" ] && grep -q 'not requeued: already requeued' "$log_file"; then
        continue_after_requeue_limit || {
            diagnose_job "$OUT_DIR" "the automatic requeue limit was reached and continuation submission failed" || true
            find_replacement_job "$OUT_DIR" || exit 1
        }
        missing_checks=0
        continue
    fi

    if [ -f "$log_file" ] && grep -q 'finished with exit status 0' "$log_file" && \
       [ -f "$OUT_DIR/reus.done" ] && [ -f "$OUT_DIR/mbar.dat" ]; then
        record "job completed successfully"
        if ! finalize_job "$OUT_DIR"; then
            record "finalization failed; manual intervention required"
            exit 1
        fi
        record "monitor finished"
        exit 0
    fi

    diagnose_job "$OUT_DIR" "the job left the queue without a verified successful result" || true
    if find_replacement_job "$OUT_DIR"; then
        missing_checks=0
        continue
    fi
    record "no recoverable replacement job was found"
    exit 1
done
