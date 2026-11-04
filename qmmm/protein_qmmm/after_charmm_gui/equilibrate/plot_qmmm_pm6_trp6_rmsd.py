#!/usr/bin/env python3
"""Plot Trp6 heavy-atom RMSD for the 100 fs PM6 QM/MM trajectory."""

import csv
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
TOPOLOGY = ROOT.parent / "step3_input.parm7"
REFERENCE = ROOT / "equilibrated.pdb"
TRAJECTORY = ROOT / "qmmm_pm6_100fs.nc"
CPPTRAJ_INPUT = ROOT / "qmmm_pm6_trp6_rmsd_cpptraj.in"
CPPTRAJ_LOG = ROOT / "qmmm_pm6_trp6_rmsd_cpptraj.log"
RMSD_DATA = ROOT / "qmmm_pm6_trp6_heavy_atom_rmsd.dat"
CSV_OUTPUT = ROOT / "qmmm_pm6_trp6_rmsd_vs_time.csv"
PLOT_OUTPUT = ROOT / "qmmm_pm6_trp6_rmsd_vs_time.png"

ALIGNMENT_MASK = ":1-20@CA,C,N,O"
TRP6_HEAVY_MASK = ":6&!@H="
FRAME_INTERVAL_FS = 5.0


def calculate_rmsd() -> np.ndarray:
    commands = [
        f"parm {TOPOLOGY}",
        f"trajin {TRAJECTORY}",
        f"reference {REFERENCE} [qmmm_start]",
        "autoimage",
        f"rms reference {ALIGNMENT_MASK}",
        f"rms reference nofit {TRP6_HEAVY_MASK} out {RMSD_DATA}",
        "run",
    ]
    CPPTRAJ_INPUT.write_text("\n".join(commands) + "\n", encoding="utf-8")
    with CPPTRAJ_INPUT.open("r", encoding="utf-8") as input_file, CPPTRAJ_LOG.open(
        "w", encoding="utf-8"
    ) as log_file:
        subprocess.run(
            ["/appl/share/amber/bin/cpptraj"],
            stdin=input_file,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return np.loadtxt(RMSD_DATA, comments="#", usecols=1)


def write_outputs(rmsd_angstrom: np.ndarray) -> None:
    time_fs = np.arange(1, len(rmsd_angstrom) + 1) * FRAME_INTERVAL_FS
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["simulation_time_fs", "trp6_heavy_atom_rmsd_A"])
        writer.writerows(zip(time_fs, rmsd_angstrom))

    figure, axis = plt.subplots(figsize=(7.5, 5.0))
    axis.plot(time_fs, rmsd_angstrom, color="tab:purple", marker="o", linewidth=1.6)
    axis.set_xlabel("Simulation time (fs)")
    axis.set_ylabel("Trp6 heavy-atom RMSD (Å)")
    axis.set_title("Trp6 Heavy-Atom RMSD During PM6 QM/MM NVT MD")
    axis.set_xlim(0.0, 100.0)
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(PLOT_OUTPUT, dpi=200)
    plt.close(figure)

    print(f"Frames analyzed: {len(rmsd_angstrom)}")
    print(f"Time range: {time_fs[0]:.1f}-{time_fs[-1]:.1f} fs")
    print(f"RMSD range: {rmsd_angstrom.min():.4f}-{rmsd_angstrom.max():.4f} A")
    print(f"CSV: {CSV_OUTPUT}")
    print(f"Plot: {PLOT_OUTPUT}")


if __name__ == "__main__":
    write_outputs(calculate_rmsd())
