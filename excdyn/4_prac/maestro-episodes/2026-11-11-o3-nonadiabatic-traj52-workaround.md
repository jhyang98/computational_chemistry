---
task: unsupported
engine: gamess
error_class: scf_not_converged
failure_class: scientific
outcome: workaround
wall_time_s: 1162
n_atoms: 3
method: MRSF-TDDFT/B3LYP
basis: 6-31+G(d)
cores: 1
scheduler: slurm
attempts: 2
---
Symptom: The initial DIIS calculation stopped before dynamics because the ROHF SCF did not converge, and no valid NAMD restart packet existed.
Attempts: The project controller restarted trajectory 52 from its deterministic Wigner sample with DIIS disabled and SOSCF enabled. The second attempt completed all 60 steps.
Result: `traj_052` reached 30.0 fs on S2 without dissociation. The censored 52-trajectory analysis completed in `analysis/4_main/total/3_add2`, with trajectory 52 coordinates highlighted in red.
Context: Three-atom ozone MRSF-TDDFT/FSSH from S3, 0.5 fs timestep, one core on node3. MAESTRO lacks a registered nonadiabatic surface-hopping task, so the existing raw QMTask/GAMESS campaign controller was used.
