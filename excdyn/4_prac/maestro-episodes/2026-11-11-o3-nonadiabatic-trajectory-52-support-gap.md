---
task: unsupported
engine: none
error_class: support_gap
outcome: support_gap
---
Symptom: The user requested trajectory 52 for the existing MRSF-TDDFT/FSSH nonadiabatic dynamics campaign.
Attempts: Searched the live MAESTRO catalog for nonadiabatic dynamics, trajectory, molecular dynamics, excited-state tasks, and nonadiabatic recipes. Only the classical MDTask was relevant, and it does not perform excited-state surface hopping.
Result: The existing project-specific raw QMTask/GAMESS campaign controller is used as the workaround.
Context: O3, 30 fs trajectory, 0.5 fs timestep, one core on node3, followed by censored ensemble analysis with the new trajectory highlighted.
