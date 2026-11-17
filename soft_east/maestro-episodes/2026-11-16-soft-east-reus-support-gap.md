---
task: unsupported
engine: none
error_class: support_gap
outcome: support_gap
---
Symptom: The requested calculation is replica-exchange umbrella sampling with MBAR for a custom softened East kinetically constrained model.
Attempts: The live MAESTRO catalog was searched for replica exchange umbrella sampling, kinetically constrained model, and Monte Carlo tasks. No matching task or recipe was found. The molecular dynamics search returned MDTask, which applies to classical molecular dynamics of prepared molecular systems and does not represent this discrete-time trajectory-ensemble calculation.
Result: Use the existing project-specific C99 executables and SLURM submission workflow.
Context: The calculation uses hard activity windows, replica exchange, trajectory Monte Carlo, and MBAR analysis.
