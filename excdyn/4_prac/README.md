# O3 MRSF-TDDFT Surface-Hopping Study

## Current production result

`nonadiabatic/4_main` is the retained 52-trajectory production campaign. All trajectories use zero-temperature Wigner sampling and the same MRSF-TDDFT/FSSH settings.

- Method: MRSF-TDDFT, B3LYP/6-31+G(d), five singlet states
- Dynamics: GAMESS FSSH/NAMD with classical nuclei
- Timestep: 0.5 fs
- Segment: one 30 fs segment, 60 steps
- Parallel layout: 30 trajectories at one core each on node3
- Outcome: 51 trajectories reached 30 fs; `traj_004` is censored at 16 fs after 10 SCF attempts
- Ensemble analysis: `analysis/4_main/total/3_add2/`

For the 52-member ensemble, the 51 survivors at 30 fs have active-state fractions S1=0.059, S2=0.235, S3=0.686, and S4=0.020. The new `traj_052` finishes on S2 without crossing the 2.0 A dissociation threshold.

## Layout

- `geometry/`: optimized O3 structure and normal modes
- `initial/`: electronic-state selection calculations
- `nonadiabatic/4_main/`: retained inputs, Wigner samples, attempt logs, restart packets, state, and controller
- `analysis/4_main/`: per-trajectory and censored ensemble analysis
- `nonadiabatic/backup/`: clean recovery copy of essential scripts and a configuration template

## Recovery and extension

The controller recovers GAMESS restart packets into each trajectory's `restart/` directory and retries from the latest valid packet. If an initial SCF fails before a restart packet exists, that attempt necessarily restarts from 0 fs. The retry limit is 10.

`python3 controller.py --extend-segments 1` increases the target by one additional 30 fs segment. `python3 controller.py --add-trajectories 5` creates five additional deterministic Wigner trajectories while segment 1 is current. Review the failed/censored trajectory policy before extending this completed campaign.

Use `cd nonadiabatic/4_main && python3 monitor.py` to write and display the current campaign state.
