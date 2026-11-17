---
task: unsupported
engine: none
error_class: support_gap
outcome: support_gap
---
Symptom: The user requested replica-exchange umbrella sampling for the softened East model at epsilon = 0.001, g = 0, L = 100, and t_obs = 600 to test for phase coexistence.
Attempts: The live MAESTRO catalog was searched for East model, replica-exchange umbrella sampling, umbrella sampling, replica exchange, trajectory sampling, phase coexistence, and Monte Carlo tasks, quantities, and recipes. No matching capability was found.
Result: MAESTRO does not support this model or sampling workflow, so the calculation must use the project's existing standalone C REUS implementation and SLURM submission scripts.
Context: This is a nonequilibrium statistical-mechanics trajectory calculation rather than a registered quantum-chemistry or molecular-dynamics task.
