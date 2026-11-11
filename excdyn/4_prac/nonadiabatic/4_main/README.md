# 4_main: O3 MRSF-TDDFT/FSSH at 0.5 fs

Campaign `20261106_161236` is complete with censored attrition: 51 of 52 Wigner trajectories reached 30 fs, while `traj_004` remains valid through 16 fs and failed after 10 SCF attempts.

`traj_031`-`traj_052` use consecutive deterministic Wigner seeds and the same 0.5 fs, 60-step protocol. The latest ensemble is summarized in `analysis/4_main/total/3_add2/`, retaining `traj_004` as censored data and highlighting `traj_052` in red.

- Timestep: 0.5 fs (`DT=5.0D-16 s`)
- Segment: 30 fs, 60 steps
- Initial state: S3 (`IROOT=3`)
- Parallel layout used: 30 jobs at one core each
- Sampling: deterministic zero-temperature Wigner distribution, seeds 20261106-20261157

The retained attempt directories contain `qm/run/job.out` and `qm/results/result.json`; duplicated MAESTRO metadata, duplicate output copies, scratch intermediates, and attempt-specific restart copies were removed after analysis. Canonical restart packets remain under each trajectory's `1_segment_30fs/restart/` directory.

Run `python3 monitor.py` for status. The controller supports `--extend-segments N` and `--add-trajectories N`, but extension requires an explicit decision about how to handle the censored `traj_004`.
