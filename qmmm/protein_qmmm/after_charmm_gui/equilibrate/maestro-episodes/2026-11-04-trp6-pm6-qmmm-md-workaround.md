---
task: unsupported
engine: none
error_class: support_gap
outcome: workaround
wall_time_s: 17
n_atoms: 5975
method: PM6
cores: 1
scheduler: slurm
attempts: 1
---
Symptom: MAESTRO's QM/MM MD loop requires PySCF, whose live capability matrix does not include PM6.

Attempts: One native Amber sander job was submitted to node3 with the equilibrated restart, CHARMM-GUI Amber topology, and the requested Trp6 QM mask.

Result: The workaround completed 200 steps at 0.5 fs, reaching 100 fs total. The final reported temperature was 301.10 K, and Amber produced the restart and NetCDF trajectory.

Context: The run used NVT at 300 K, a neutral 18-atom QM region, and the CA-CB link boundary. Amber reported 16.96 seconds elapsed on one CPU.
