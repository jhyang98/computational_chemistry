---
task: unsupported
engine: none
error_class: support_gap
failure_class: unclassified
outcome: workaround
wall_time_s: 42309
cores: 4
scheduler: slurm
attempts: 1
---
Symptom: MAESTRO has no task for replica-exchange umbrella sampling and MBAR analysis of the softened East kinetic model.
Attempts: The already documented support-gap decision was followed by one project-specific C99 REUS run on SLURM. The run completed without requeue, and its MBAR analysis converged.
Result: The standalone workflow completed at epsilon 0.001, g 0, L 100, and observation time 1200 MCS. The activity distribution remained unimodal, while the susceptibility maximum moved to s 0.0042 and increased to 33.64.
Context: The measured wall time includes REUS and MBAR. This result is a workaround using the existing project implementation, not a MAESTRO calculation.
