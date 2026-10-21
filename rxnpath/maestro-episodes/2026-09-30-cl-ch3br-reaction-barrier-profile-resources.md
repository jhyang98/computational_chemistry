---
task: ReactionBarrierTask
engine: none
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
Symptom: A clean retry failed before scheduler submission because the built-in Python reaction-profile node had an unanswered resource core count.
Attempts: Traced MAESTRO's planner. All seven ORCA nodes retained the explicit four-core allocation; only the Python profile node inherited an unanswered default.
Result: Bound the built-in profile node to the same four-core SLURM allocation and resubmitted.
Context: Gas-phase chloride plus methyl bromide SN2 reaction barrier; PBE0/aug-cc-pVDZ; six atoms at each endpoint.
