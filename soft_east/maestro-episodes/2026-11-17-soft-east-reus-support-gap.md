---
task: unsupported
engine: none
error_class: support_gap
failure_class: unclassified
outcome: support_gap
---
Symptom: The user requested a replica-exchange umbrella-sampling calculation for the softened East kinetic model with epsilon 0.001 and an observation time of 1200 MCS.
Attempts: Searched the live MAESTRO capability catalog for kinetic Monte Carlo, Markov-chain trajectory sampling, Monte Carlo, replica-exchange umbrella sampling, umbrella sampling, the East model, and related tasks and recipes. No matching task or recipe was available. MDTask was the nearest result, but it is for classical molecular dynamics of atomic coordinates and does not implement this stochastic lattice-model calculation.
Result: Use the project's existing C99 east_reus executable and SLURM submission workflow outside MAESTRO.
Context: The calculation samples biased trajectory ensembles of a softened East model and combines overlapping activity windows with MBAR. It is not an atomistic molecular-dynamics run.
