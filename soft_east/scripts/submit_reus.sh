#!/usr/bin/env bash
# Submit scripts/run_reus.sh as a SLURM job (batch script: scripts/slurm_reus.sh).
#
# Run from the project root:
#   bash scripts/submit_reus.sh                            # defaults of run_reus.sh
#   bash scripts/submit_reus.sh EPS=0.003 G=0.05 TOBS=400  # override parameters for this job
#   bash scripts/submit_reus.sh MODE=mbar EPS=0.003        # MBAR only (REUS already finished)
#   bash scripts/submit_reus.sh --dry-run EPS=0.003        # show the job, do not submit
#   bash scripts/submit_reus.sh --from slurm_jobs/<file>.env  # same parameters as that job
#                                                          # (e.g. to continue it after
#                                                          # MAX_REQUEUE was reached)
#
# KEY=VALUE pairs set parameters of scripts/run_reus.sh for this job only. All
# parameters are frozen at submission in slurm_jobs/<file>.env, so editing
# run_reus.sh later does not change jobs that are already queued. Every
# submission is listed in slurm_jobs/jobs.txt.
set -euo pipefail

# ---- cluster settings --------------------------------------------------------
NODE="node3"        # SLURM name of compute node 3
PARTITION=""        # partition (empty: the cluster default)
CPUS=4              # cores per job = OpenMP threads
TIME=24:00:00       # requested limit; capped at the SLURM partition maximum
JOB_NAME=soft_east  # jobs with this name run one at a time (--dependency=singleton),
                    # so at most CPUS cores are in use at once
MAX_REQUEUE=10      # automatic requeues allowed per job
MARGIN=1800         # REUS checkpoints and stops this many seconds before the limit
SIGNAL_LEAD=600     # SIGUSR1 to the batch script this many seconds before the limit
LOG_DIR=slurm_logs
JOB_DIR=slurm_jobs
# ------------------------------------------------------------------------------

to_seconds() { # D-HH:MM:SS, HH:MM:SS, or MM:SS -> seconds
    local t=$1 d=0 h=0 m=0 s=0
    local -a fields
    if [[ $t == *-* ]]; then
        d=${t%%-*}
        t=${t#*-}
    fi
    IFS=: read -r -a fields <<< "$t"
    case ${#fields[@]} in
        3) h=${fields[0]}; m=${fields[1]}; s=${fields[2]} ;;
        2) m=${fields[0]}; s=${fields[1]} ;;
        *) return 1 ;;
    esac
    [[ $d =~ ^[0-9]+$ && ${h:-x} =~ ^[0-9]+$ && ${m:-x} =~ ^[0-9]+$ && ${s:-x} =~ ^[0-9]+$ ]] ||
        return 1
    echo $(((((10#$d * 24) + 10#$h) * 60 + 10#$m) * 60 + 10#$s))
}

to_slurm_time() { # seconds -> D-HH:MM:SS or HH:MM:SS
    local total=$1 d h m s
    d=$((total / 86400))
    total=$((total % 86400))
    h=$((total / 3600))
    total=$((total % 3600))
    m=$((total / 60))
    s=$((total % 60))
    if [ "$d" -gt 0 ]; then
        printf '%d-%02d:%02d:%02d\n' "$d" "$h" "$m" "$s"
    else
        printf '%02d:%02d:%02d\n' "$h" "$m" "$s"
    fi
}

partition_limit() {
    local row
    if [ -n "$PARTITION" ]; then
        row=$(sinfo -h -p "$PARTITION" -o '%P|%l' |
            awk -F'|' 'NF >= 2 {gsub(/[[:space:]]/, "", $1); gsub(/[[:space:]]/, "", $2); print $1 "|" $2; exit}')
    else
        row=$(sinfo -h -o '%P|%l' |
            awk -F'|' '$1 ~ /\*/ {gsub(/[[:space:]*]/, "", $1); gsub(/[[:space:]]/, "", $2); print $1 "|" $2; exit}')
    fi
    [ -n "$row" ] || return 1
    printf '%s\n' "$row"
}

dry=0
mode=""
overrides=()
known=" $(bash scripts/run_reus.sh params | grep -v '^#' | cut -d= -f1 | tr '\n' ' ') "
if [ "${1:-}" = "--from" ]; then
    if [ $# -lt 2 ] || [ ! -f "$2" ]; then
        echo "error: --from needs an existing parameter file (slurm_jobs/*.env)" >&2
        exit 1
    fi
    keep_max=$MAX_REQUEUE
    set -a
    # shellcheck disable=SC1090
    . "$2"
    set +a
    mode=${MODE:-}
    unset STOP_AFTER MODE
    MAX_REQUEUE=$keep_max # the cluster setting above, not the one in the file
    shift 2
fi
for a in "$@"; do
    case "$a" in
        --dry-run) dry=1 ;;
        MODE=*) mode=${a#MODE=} ;;
        *=*)
            k=${a%%=*}
            if [[ $known != *" $k "* || $k == STOP_AFTER ]]; then
                echo "error: unknown parameter '$k' (see: bash scripts/run_reus.sh params)" >&2
                exit 1
            fi
            overrides+=("$a")
            ;;
        *)
            echo "usage: bash scripts/submit_reus.sh [--from FILE.env] [--dry-run] [MODE=all|mbar]" \
                "[KEY=VALUE ...]" >&2
            exit 1
            ;;
    esac
done
mode=${mode:-all}
if [[ $mode != all && $mode != mbar ]]; then
    echo "error: MODE must be all or mbar" >&2
    exit 1
fi
if [ -z "$NODE" ]; then
    echo "error: set NODE at the top of scripts/submit_reus.sh to the SLURM name of compute" >&2
    echo "       node 3 (list the nodes with: sinfo -N -l)" >&2
    exit 1
fi
if ! requested_limit=$(to_seconds "$TIME"); then
    echo "error: TIME must look like 24:00:00 or 1-00:00:00" >&2
    exit 1
fi
if ! partition_info=$(partition_limit); then
    if [ -n "$PARTITION" ]; then
        echo "error: cannot read the SLURM time limit for partition '$PARTITION'" >&2
    else
        echo "error: cannot find the default SLURM partition and its time limit" >&2
        echo "       set PARTITION at the top of scripts/submit_reus.sh" >&2
    fi
    exit 1
fi
job_partition=${partition_info%%|*}
server_time=${partition_info#*|}
limit=$requested_limit
effective_time=$TIME
case ${server_time^^} in
    INFINITE|UNLIMITED)
        ;;
    *)
        if ! server_limit=$(to_seconds "$server_time"); then
            echo "error: unsupported SLURM time limit '$server_time' for partition '$job_partition'" >&2
            exit 1
        fi
        if [ "$server_limit" -lt "$limit" ]; then
            limit=$server_limit
            effective_time=$(to_slurm_time "$limit")
        fi
        ;;
esac
if [ "$limit" -lt 120 ]; then
    echo "error: effective SLURM time limit must be at least 120 seconds" >&2
    exit 1
fi
effective_margin=$MARGIN
max_margin=$((limit / 4))
if [ "$effective_margin" -gt "$max_margin" ]; then effective_margin=$max_margin; fi
effective_signal_lead=$SIGNAL_LEAD
max_signal_lead=$((effective_margin / 2))
if [ "$effective_signal_lead" -gt "$max_signal_lead" ]; then
    effective_signal_lead=$max_signal_lead
fi
stop_after=$((limit - effective_margin))
if [ "$stop_after" -le 0 ] || [ "$effective_signal_lead" -le 0 ] ||
   [ "$effective_signal_lead" -ge "$effective_margin" ]; then
    echo "error: cannot fit checkpoint and signal margins inside effective time $effective_time" >&2
    exit 1
fi
[ -x bin/east_reus ] && [ -x bin/east_mbar ] || {
    echo "error: build first (make OMP=1)" >&2
    exit 1
}

mkdir -p "$LOG_DIR" "$JOB_DIR"
envfile="$JOB_DIR/${JOB_NAME}_$(date +%Y%m%d_%H%M%S)_$$.env"
{
    env ${overrides[@]+"${overrides[@]}"} STOP_AFTER="$stop_after" bash scripts/run_reus.sh params |
        grep -v '^#'
    printf 'MODE=%q\n' "$mode"
    printf 'MAX_REQUEUE=%q\n' "$MAX_REQUEUE"
} > "$envfile"
outdir=$(env ${overrides[@]+"${overrides[@]}"} bash scripts/run_reus.sh params | sed -n 's/^# output folder: //p')

args=(--job-name="$JOB_NAME" --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --time="$effective_time"
      --nodelist="$NODE" --requeue --dependency=singleton --open-mode=append
      --signal=B:USR1@"$effective_signal_lead" --output="$LOG_DIR/%x_%j.out"
      --export=ALL,JOB_ENV="$envfile")
if [ -n "$PARTITION" ]; then args+=(--partition="$PARTITION"); fi

echo "parameters : $envfile"
echo "output     : $outdir"
echo "partition  : $job_partition (maximum $server_time)"
echo "mode       : $mode, $CPUS cores on $NODE, time limit $effective_time, REUS stops after ${stop_after} s"
echo "command    : sbatch ${args[*]} scripts/slurm_reus.sh"
if [ "$dry" -eq 1 ]; then
    echo "(dry run: not submitted)"
    rm -f "$envfile"
    exit 0
fi
jobid=$(sbatch --parsable "${args[@]}" scripts/slurm_reus.sh)
jobid=${jobid%%;*}
printf '%s  %s  %s  %s  %s\n' "$jobid" "$(date '+%F %T')" "$mode" "$outdir" "$envfile" \
    >> "$JOB_DIR/jobs.txt"
echo "submitted  : job $jobid (log: $LOG_DIR/${JOB_NAME}_$jobid.out)"
