import argparse
import fcntl
import getpass
import json
import os
import re
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANALYSIS = ROOT / "analysis" / HERE.name
CONFIG = json.loads((HERE / "config.json").read_text())
CONFIG["scratch"] = CONFIG["scratch"].format(user=getpass.getuser())
STATE_PATH = HERE / "campaign_state.json"
EVENTS = HERE / "events.jsonl"
LOCK_PATH = HERE / "controller.lock"
CHILDREN = {}


def event(kind, **fields):
    record = {"time": datetime.now().isoformat(timespec="seconds"), "event": kind, **fields}
    with EVENTS.open("a") as handle:
        handle.write(json.dumps(record) + "\n")
    print(json.dumps(record), flush=True)


def save_state(state):
    temporary = STATE_PATH.with_suffix(".tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    temporary.replace(STATE_PATH)


def segment_name(segment):
    return f"{segment}_segment_{segment * CONFIG['segment_fs']}fs"


def trajectory_name(trajectory):
    return f"traj_{trajectory:03d}"


def segment_dir(trajectory, segment):
    return HERE / trajectory_name(trajectory) / segment_name(segment)


def trajectory_count(state):
    return state.get("n_trajectories", CONFIG["n_trajectories"])


def result_state(work):
    result = work / "qm/results/result.json"
    if not result.exists():
        return None
    try:
        return json.loads(result.read_text()).get("state")
    except json.JSONDecodeError:
        return None


def submission_job_id(work):
    submission = work / ".maestro/submission.json"
    if not submission.exists():
        return None
    try:
        tasks = json.loads(submission.read_text()).get("tasks", [])
    except json.JSONDecodeError:
        return None
    return str(tasks[0]["job_id"]) if tasks else None


def job_active(job_id):
    if not job_id:
        return False
    try:
        query = subprocess.run(
            ["squeue", "-h", "-j", str(job_id), "-o", "%T"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return True
    return bool(query.stdout.strip())


def recover_scratch(work, restart_dir):
    scratch = work / "qm/run/scratch.json"
    if not scratch.exists():
        return
    userscr = json.loads(scratch.read_text())["userscr"]
    restart_dir.mkdir(parents=True, exist_ok=True)
    try:
        listing = subprocess.run(
            ["ssh", CONFIG["node"], "find", userscr, "-maxdepth", "1", "-type", "f", "-print"],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return
    for remote in listing.stdout.splitlines():
        if Path(remote).suffix not in {".namdrst", ".rst", ".trj", ".dat"}:
            continue
        try:
            subprocess.run(
                ["scp", f"{CONFIG['node']}:{remote}", str(restart_dir / Path(remote).name)],
                timeout=20,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            pass


def latest_valid_packet(restart_dir):
    for packet in sorted(restart_dir.glob("*.namdrst"), key=lambda path: path.stat().st_mtime, reverse=True):
        text = packet.read_text(errors="replace")
        if "QM PARTICLE COORDINATES FOR $DATA GROUP" in text and re.search(r"(?m)^ \$MD READ=", text):
            return packet
    return None


def packet_time(packet):
    match = re.search(r"TTOTAL=\s*([0-9.]+) FS", packet.read_text(errors="replace"))
    return float(match.group(1)) if match else 0.0


def soscf_input(source, output):
    text = source.read_text()
    text = re.sub(
        r"^ \$SCF.*?\$END$",
        " $SCF DIRSCF=.TRUE. DIIS=.FALSE. SOSCF=.TRUE. $END",
        text,
        count=1,
        flags=re.M,
    )
    output.write_text(text)


def full_restart_input(trajectory, attempt):
    name = trajectory_name(trajectory)
    source = segment_dir(trajectory, 1) / f"{name}_segment_1.inp"
    output = segment_dir(trajectory, 2) / f"{name}_full_20fs_attempt_{attempt:02d}.inp"
    text = source.read_text()
    text = re.sub(r"NSTEPS=\s*\d+", f"NSTEPS={2 * CONFIG['steps_per_segment']}", text, count=1)
    if attempt % 3 == 1:
        scf = " $SCF DIRSCF=.TRUE. DIIS=.FALSE. SOSCF=.TRUE. $END"
        maxit = 400
    elif attempt % 3 == 2:
        scf = " $SCF DIRSCF=.TRUE. DIIS=.TRUE. SOSCF=.FALSE. DAMP=.TRUE. SHIFT=.TRUE. $END"
        maxit = 400
    else:
        scf = " $SCF DIRSCF=.TRUE. DIIS=.TRUE. SOSCF=.FALSE. $END"
        maxit = 200
    text = re.sub(r"MAXIT=\d+", f"MAXIT={maxit}", text, count=1)
    text = re.sub(r"^ \$SCF.*?\$END$", scf, text, count=1, flags=re.M)
    if f"NSTEPS={2 * CONFIG['steps_per_segment']}" not in text:
        raise RuntimeError(f"failed to prepare full 20 fs input from {source}")
    output.write_text(text)
    return output


def restart_input(packet, output, target_time, use_soscf, new_random, robust_scf=False):
    text = packet.read_text(errors="replace")
    coordinate_match = re.search(
        r"QM PARTICLE COORDINATES FOR \$DATA GROUP -----\n(.*?)(?=\n \$EFRAG)", text, re.S
    )
    md_match = re.search(r"(?ms)^ \$MD READ=.*?^ \$END\s*$", text)
    electronic_match = re.search(r"(?ms)^ \$TDDFT\s*$.*\Z", text)
    if coordinate_match is None or md_match is None or electronic_match is None:
        raise ValueError("restart packet lacks coordinates, $MD, or electronic restart groups")
    remaining = max(1, round((target_time - packet_time(packet)) / CONFIG["timestep_fs"]))
    md = re.sub(r"NSTEPS=\s*\d+", f"NSTEPS={remaining:8d}", md_match.group(0), count=1)
    md = re.sub(r"NEWRND=\.[TFtf]\.\s*", f"NEWRND={'.t.' if new_random else '.f.'}\n", md, count=1)
    md = re.sub(r" DELR=\s*0\.020", " DELR=   0.020 TRJFMT=F10o5", md, count=1)
    if robust_scf:
        scf = "DIIS=.TRUE. SOSCF=.FALSE. DAMP=.TRUE. SHIFT=.TRUE."
        maxit = 400
    else:
        scf = "DIIS=.FALSE. SOSCF=.TRUE." if use_soscf else "DIIS=.TRUE. SOSCF=.FALSE."
        maxit = 200
    header = (
        " $CONTRL SCFTYP=ROHF RUNTYP=MD DFTTYP=B3LYP TDDFT=MRSF\n"
        f"  MULT=3 COORD=UNIQUE MAXIT={maxit} $END\n"
        " $SYSTEM MWORDS=500 $END\n"
        " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 DIFFSP=.TRUE. $END\n"
        f" $SCF DIRSCF=.TRUE. {scf} $END\n"
        " $DATA\nO3 MRSF-TDDFT FSSH automatic restart\nC1\n"
        f"{coordinate_match.group(1).strip()}\n $END\n"
    )
    output.write_text(header + md + "\n" + electronic_match.group(0).rstrip() + "\n")


def prepare_input(trajectory, segment, attempt, run_mode="segment"):
    if run_mode == "full_restart":
        return full_restart_input(trajectory, attempt)
    directory = segment_dir(trajectory, segment)
    base = directory / f"{trajectory_name(trajectory)}_segment_{segment}.inp"
    if attempt == 0:
        return base
    packet = latest_valid_packet(directory / "restart")
    output = directory / f"{trajectory_name(trajectory)}_segment_{segment}_attempt_{attempt:02d}.inp"
    if packet is None:
        soscf_input(base, output)
        return output
    restart_input(
        packet,
        output,
        segment * CONFIG["segment_fs"],
        use_soscf=attempt % 2 == 1,
        new_random=False,
        robust_scf=attempt >= 8,
    )
    return output


def initialize_state():
    if STATE_PATH.exists():
        return json.loads(STATE_PATH.read_text())
    segment = 1
    trajectories = {}
    for trajectory in range(1, CONFIG["n_trajectories"] + 1):
        complete = (segment_dir(trajectory, segment) / "COMPLETE").exists()
        trajectories[str(trajectory)] = {
            "status": "complete" if complete else "pending",
            "attempt": 0,
            "work": None,
            "job_id": None,
            "submitted_at": None,
        }
    state = {
        "version": 2,
        "campaign_id": datetime.now().strftime("%Y%m%d_%H%M%S"),
        "segment": segment,
        "n_trajectories": CONFIG["n_trajectories"],
        "target_segments": CONFIG["target_segments"],
        "trajectories": trajectories,
    }
    save_state(state)
    return state


def launch(state, trajectory):
    record = state["trajectories"][str(trajectory)]
    segment = state["segment"]
    attempt = record["attempt"]
    run_mode = record.get("run_mode", "segment")
    input_file = prepare_input(trajectory, segment, attempt, run_mode)
    token = datetime.now().strftime("%Y%m%d%H%M%S%f")
    mode_tag = "full_20fs_" if run_mode == "full_restart" else ""
    work = segment_dir(trajectory, segment) / f"campaign_{state['campaign_id']}_{mode_tag}attempt_{attempt:02d}_{token}"
    log_path = segment_dir(trajectory, segment) / f"campaign_{state['campaign_id']}_{mode_tag}attempt_{attempt:02d}.log"
    log_handle = log_path.open("ab")
    process = subprocess.Popen(
        ["maestro-python", str(HERE / "run_segment.py"), str(input_file), str(work)],
        cwd=HERE,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    CHILDREN[str(work)] = (process, log_handle)
    record.update({
        "status": "running",
        "work": str(work),
        "job_id": None,
        "submitted_at": time.time(),
    })
    save_state(state)
    event("submit", trajectory=trajectory, segment=segment, attempt=attempt, run_mode=run_mode, work=str(work))


def finish_attempt(state, trajectory, success):
    record = state["trajectories"][str(trajectory)]
    work = Path(record["work"])
    attempt_restart = work.parent / f"{work.name}_restart"
    recover_scratch(work, attempt_restart)
    canonical_restart = segment_dir(trajectory, state["segment"]) / "restart"
    canonical_restart.mkdir(parents=True, exist_ok=True)
    for packet in attempt_restart.glob("*"):
        shutil.copy2(packet, canonical_restart / packet.name)
    child = CHILDREN.pop(str(work), None)
    if child:
        child[1].close()
    if success:
        record["status"] = "complete"
        record["completed_work"] = str(work)
        (segment_dir(trajectory, state["segment"]) / "COMPLETE").write_text(datetime.now().isoformat() + "\n")
        event("trajectory_complete", trajectory=trajectory, segment=state["segment"], attempt=record["attempt"])
    else:
        record["attempt"] += 1
        if record["attempt"] >= CONFIG["max_attempts"]:
            record["status"] = "failed"
            event("trajectory_failed", trajectory=trajectory, segment=state["segment"], attempts=record["attempt"])
        else:
            record["status"] = "pending"
            event("retry_queued", trajectory=trajectory, segment=state["segment"], attempt=record["attempt"])
    record["job_id"] = None
    record["work"] = None
    save_state(state)


def refresh_running(state):
    for key, record in state["trajectories"].items():
        if record["status"] != "running":
            continue
        trajectory = int(key)
        work = Path(record["work"])
        current_result = result_state(work)
        if not record["job_id"]:
            record["job_id"] = submission_job_id(work)
            if record["job_id"]:
                save_state(state)
        if current_result in {"success", "failed"} and record["job_id"] and job_active(record["job_id"]):
            continue
        if current_result == "success":
            finish_attempt(state, trajectory, True)
            continue
        if current_result == "failed":
            finish_attempt(state, trajectory, False)
            continue
        child = CHILDREN.get(str(work))
        if child and child[0].poll() is not None:
            finish_attempt(state, trajectory, False)
            continue
        if not child and record["job_id"] and not job_active(record["job_id"]):
            finish_attempt(state, trajectory, False)
            continue
        if not child and not record["job_id"] and time.time() - record["submitted_at"] > 120:
            finish_attempt(state, trajectory, False)


def fill_slots(state):
    running = sum(record["status"] == "running" for record in state["trajectories"].values())
    for key, record in state["trajectories"].items():
        if running >= CONFIG["concurrent_trajectories"]:
            break
        if record["status"] == "pending":
            launch(state, int(key))
            running += 1


def segment_logs(trajectory, segment, campaign_id):
    directory = segment_dir(trajectory, segment)
    candidates = []
    for work in directory.glob(f"campaign_{campaign_id}_attempt_*"):
        log = work / "qm/run/job.out"
        if log.exists() and b"*** AT T=" in log.read_bytes():
            candidates.append(log)
    if not candidates:
        for work in directory.glob("work_attempt_*"):
            log = work / "qm/run/job.out"
            if log.exists() and b"*** AT T=" in log.read_bytes():
                candidates.append(log)
    if not candidates:
        return []
    return sorted(candidates, key=lambda path: path.parents[2].stat().st_mtime)


def authoritative_logs(state, trajectory):
    record = state["trajectories"][str(trajectory)]
    if record.get("run_mode") == "full_restart" and record.get("completed_work"):
        log = Path(record["completed_work"]) / "qm/run/job.out"
        if not log.exists():
            raise RuntimeError(f"missing successful full-restart log for trajectory {trajectory}")
        return [log]
    logs = []
    for completed_segment in range(1, state["segment"] + 1):
        logs.extend(segment_logs(trajectory, completed_segment, state["campaign_id"]))
    return logs


def analyze_segment(state):
    segment = state["segment"]
    name = segment_name(segment)
    csv_files = []
    for trajectory in range(1, trajectory_count(state) + 1):
        logs = authoritative_logs(state, trajectory)
        if not logs:
            raise RuntimeError(f"missing all trajectory logs for trajectory {trajectory}")
        output = ANALYSIS / trajectory_name(trajectory) / name
        subprocess.run(
            ["python3", str(ANALYSIS / "analyze_segments.py"), str(output), *map(str, logs)],
            check=True,
        )
        csv_files.append(output / "observables.csv")
    total = ANALYSIS / "total" / state.get("analysis_output", name)
    analyzer = "analyze_censored_ensemble.py" if any(
        record["status"] == "failed" for record in state["trajectories"].values()
    ) else "analyze_ensemble.py"
    command = ["python3", str(ANALYSIS / analyzer)]
    if state.get("analysis_output"):
        command.extend([
            "--target-time",
            str(segment * CONFIG["segment_fs"]),
            "--highlight-trajectory",
            trajectory_name(trajectory_count(state)),
        ])
    subprocess.run([*command, str(total), *map(str, csv_files)], check=True)
    event("analysis_complete", segment=segment, output=str(total))


def advance_segment(state):
    previous = state["segment"]
    following = previous + 1
    for trajectory in range(1, trajectory_count(state) + 1):
        packet = latest_valid_packet(segment_dir(trajectory, previous) / "restart")
        if packet is None:
            raise RuntimeError(f"trajectory {trajectory} has no valid restart packet for segment {following}")
        directory = segment_dir(trajectory, following)
        directory.mkdir(parents=True, exist_ok=True)
        restart_input(
            packet,
            directory / f"{trajectory_name(trajectory)}_segment_{following}.inp",
            following * CONFIG["segment_fs"],
            use_soscf=False,
            new_random=False,
        )
        state["trajectories"][str(trajectory)] = {
            "status": "pending",
            "attempt": 0,
            "work": None,
            "job_id": None,
            "submitted_at": None,
        }
    state["segment"] = following
    save_state(state)
    event("segment_start", segment=following)


def configure_full_restart(state, targets):
    if state["segment"] != 2:
        raise RuntimeError("full 0-20 fs restart is only valid while segment 2 is active")
    for trajectory in targets:
        if str(trajectory) not in state["trajectories"]:
            raise ValueError(f"unknown trajectory {trajectory}")
        state["trajectories"][str(trajectory)] = {
            "status": "pending",
            "attempt": 0,
            "run_mode": "full_restart",
            "work": None,
            "job_id": None,
            "submitted_at": None,
            "completed_work": None,
        }
        complete_marker = segment_dir(trajectory, 2) / "COMPLETE"
        if complete_marker.exists():
            complete_marker.unlink()
    state["phase"] = "full_restart"
    state["full_restart_targets"] = targets
    state["full_restart_started_at"] = datetime.now().isoformat(timespec="seconds")
    save_state(state)
    event("full_restart_configured", trajectories=targets)


def extend_segments(state, additional_segments):
    if additional_segments < 1:
        raise ValueError("additional segment count must be positive")
    current_target = state.get("target_segments", CONFIG["target_segments"])
    state["target_segments"] = current_target + additional_segments
    save_state(state)
    event("segment_target_extended", old_target=current_target, new_target=state["target_segments"])


def add_trajectories(state, additional_trajectories, analysis_output=None):
    if additional_trajectories < 1:
        raise ValueError("additional trajectory count must be positive")
    if state["segment"] != 1:
        raise RuntimeError("trajectory extension is supported while segment 1 is active or complete")
    first = trajectory_count(state) + 1
    subprocess.run(
        ["python3", str(HERE / "setup_wigner.py"), "--start", str(first), "--count", str(additional_trajectories)],
        check=True,
    )
    for trajectory in range(first, first + additional_trajectories):
        state["trajectories"][str(trajectory)] = {
            "status": "pending",
            "attempt": 0,
            "work": None,
            "job_id": None,
            "submitted_at": None,
        }
    state["n_trajectories"] = first + additional_trajectories - 1
    if analysis_output:
        state["analysis_output"] = analysis_output
    save_state(state)
    event("trajectories_added", first=first, count=additional_trajectories, total=state["n_trajectories"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-restart", type=int, nargs="+")
    parser.add_argument("--extend-segments", type=int)
    parser.add_argument("--add-trajectories", type=int)
    parser.add_argument("--analysis-output")
    args = parser.parse_args()
    lock_handle = LOCK_PATH.open("w")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        event("controller_already_running")
        return
    lock_handle.write(str(os.getpid()) + "\n")
    lock_handle.flush()
    state = initialize_state()
    if args.full_restart:
        configure_full_restart(state, args.full_restart)
    if args.extend_segments:
        extend_segments(state, args.extend_segments)
    if args.add_trajectories:
        add_trajectories(state, args.add_trajectories, args.analysis_output)
    event("controller_start", campaign_id=state["campaign_id"], segment=state["segment"], phase=state.get("phase", "production"))
    while True:
        refresh_running(state)
        fill_slots(state)
        statuses = [record["status"] for record in state["trajectories"].values()]
        if all(status == "complete" for status in statuses):
            event("segment_barrier_complete", segment=state["segment"])
            analyze_segment(state)
            if state.get("phase") == "full_restart":
                state["phase"] = "production"
                state["full_restart_completed_at"] = datetime.now().isoformat(timespec="seconds")
                save_state(state)
                event("full_restart_complete", trajectories=state.get("full_restart_targets", []))
            if state["segment"] >= state.get("target_segments", CONFIG["target_segments"]):
                event("campaign_complete", segment=state["segment"])
                return
            advance_segment(state)
            continue
        if "running" not in statuses and "pending" not in statuses:
            analyze_segment(state)
            event("campaign_partial_complete", segment=state["segment"])
            event("campaign_blocked", segment=state["segment"])
            return
        time.sleep(CONFIG["poll_seconds"])


if __name__ == "__main__":
    main()
