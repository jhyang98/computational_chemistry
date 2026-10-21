"""Gas-phase Cl- + CH3Br SN2 reaction barrier at PBE0/aug-cc-pVDZ."""

from maestro import Maestro, SystemQM, DFT, OrcaEngine, PythonEngine, ReactionBarrierTask
from maestro.engines.jobspec import Resources


mae = Maestro(mode="slurm", workdir=".")

# Atom order is identical at both endpoints: Cl, C, H, H, H, Br.
reactant = SystemQM(
    geometry="sn2_reactant_cl_ch3br_start.xyz", charge=-1, spin=0
)
product = SystemQM(
    geometry="sn2_reactant_cl_ch3br_end.xyz", charge=-1, spin=0
)

# No solvent is passed: this is a gas-phase calculation.
theory = DFT(functional="PBE0", basis="aug-cc-pVDZ")

# The site runinfo config supplies the validated ORCA module, executable,
# 100-hour walltime, and compute partition.  Four cores is this run's limit.
qm_engine = OrcaEngine(resources=Resources(cores=4))
# The final reaction-profile aggregation is a built-in Python node but still
# runs through SLURM, so it needs an explicit scheduler allocation too.
profile_engine = PythonEngine(resources=Resources(cores=4))

task = ReactionBarrierTask(
    system_react=reactant,
    system_prod=product,
    theory=theory,
    nimages=6,
)

res = mae.run(
    rundir="cl_ch3br_reaction_barrier_work",
    task=task,
    engines={"*": qm_engine, "python": profile_engine},
)

print("activation free energy in", res.unit("activation_free_energy"), ":",
      res.load("activation_free_energy"))
print("reaction free energy in", res.unit("reaction_free_energy"), ":",
      res.load("reaction_free_energy"))
print("profile:", res.path("figure"))
