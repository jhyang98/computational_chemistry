---
task: unsupported
engine: none
error_class: support_gap
outcome: support_gap
---
Symptom: The requested QM/MM NVT molecular dynamics uses PM6 for the Trp6 side chain.

Attempts: The live catalog showed MDTask supports QMMMEngine, but its dynamics loop requires PySCFEngine as the QM sub-engine. The PySCF capability matrix contains HF, DFT, MP2, and CCSD but not PM6.

Result: MAESTRO cannot represent this method/backend combination. The calculation was prepared through Amber's native sander QM/MM implementation instead.

Context: The requested QM region is neutral and has a CA-CB boundary, with CB and the Trp6 side-chain atoms in the QM region.
