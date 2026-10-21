"""Gas-phase Gibbs free energy for 2 H2 + O2 -> 2 H2O at 298.15 K, 1 atm."""

from dataclasses import dataclass

from maestro import CompositeTask, Maestro, MP2, OrcaEngine, SystemQM, ThermoTask
from maestro.engines.jobspec import Resources


@dataclass(kw_only=True)
class WaterFormationThermoTask(CompositeTask):
    """Three independent molecular thermochemistry calculations."""

    h2: SystemQM
    o2: SystemQM
    h2o: SystemQM
    theory: MP2
    temperature: float = 298.15
    pressure: float = 1.0

    def produces(self):
        return ThermoTask.can_produce()

    @classmethod
    def can_produce(cls):
        return ThermoTask.can_produce()

    def expand(self, g):
        for name, system in (("h2", self.h2), ("o2", self.o2), ("h2o", self.h2o)):
            g.add(
                name,
                ThermoTask(
                    system=system,
                    theory=self.theory,
                    temperature=self.temperature,
                    pressure=self.pressure,
                ),
                outputs=True,
            )


task = WaterFormationThermoTask(
    h2=SystemQM(geometry="h2.xyz", charge=0, spin=0),
    o2=SystemQM(geometry="o2.xyz", charge=0, spin=2),
    h2o=SystemQM(geometry="h2o.xyz", charge=0, spin=0),
    theory=MP2(basis="6-31G(d)"),
)

# The site runinfo config loads orca/v6.0.1; Resources controls the Slurm request.
engine = OrcaEngine(resources=Resources(cores=1, time="24:00:00"))
mae = Maestro(mode="slurm", workdir=".", runinfo_path="runinfo_orca.toml")
result = mae.run(
    rundir="water_formation_gibbs_work",
    task=task,
    engines=engine,
)

g_h2 = result.load("h2/gibbs_free_energy")
g_o2 = result.load("o2/gibbs_free_energy")
g_h2o = result.load("h2o/gibbs_free_energy")
delta_g = 2.0 * g_h2o - (2.0 * g_h2 + g_o2)

print("Reaction: 2 H2(g) + O2(g) -> 2 H2O(g)")
print(f"Delta G degrees (Hartree): {delta_g:.12f}")
print(f"Delta G degrees (kcal mol^-1): {delta_g * 627.509474:.6f}")
