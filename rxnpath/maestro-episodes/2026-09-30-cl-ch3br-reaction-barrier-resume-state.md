---
task: ReactionBarrierTask
engine: orca
error_class: submission_error
failure_class: infra_transient
outcome: workaround
wall_time_s: 2
n_atoms: 6
method: PBE0
basis: aug-cc-pVDZ
cores: 4
mode: slurm
---
Symptom: The first submission stopped before creating a scheduler job with `resources.cores is unanswered`; no graph node ran.
Attempts: The script explicitly supplied `Resources(cores=4)`. Inspection showed the failed run's serialized resume binding omitted runtime resource blocks, so a retry from that directory would load an incomplete binding.
Result: Archived the zero-result state and retried cleanly with the same confirmed four-core SLURM configuration.
Context: Gas-phase chloride plus methyl bromide SN2 reaction barrier; PBE0/aug-cc-pVDZ; six atoms at each endpoint.
