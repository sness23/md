# OpenMM "forever-ish" runner: writes topology.pdb + traj.dcd and keeps appending.
# Requirements: conda-forge::openmm
#
# Create & activate env:
#   mamba create -n live-md -c conda-forge openmm mdtraj -y
#   conda activate live-md
#
# Run:
#   python python/openmm_run.py
#
# Output files end up under ./data relative to project root.

import os, sys, time
from sys import stdout

# Configure
DATA_DIR = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)

TOPOLOGY_PATH = os.path.join(DATA_DIR, "topology.pdb")
TRAJ_PATH = os.path.join(DATA_DIR, "traj.dcd")

# --- Minimal alanine dipeptide in implicit solvent ---
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

# Build alanine dipeptide from PDB (built-in example via modeller)
# We'll create a simple peptide (Ace-Ala-Nme) by loading a known PDB-like string
ala2 = \"\"\"\
ATOM      1  N   ALA A   1      -1.207   1.207   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       1.207   0.000   1.207  1.00  0.00           C
ATOM      4  O   ALA A   1       2.121  -0.828   1.207  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.000   0.000  -1.540  1.00  0.00           C
ATOM      6  H   ALA A   1      -1.540   2.070   0.000  1.00  0.00           H
ATOM      7  HA  ALA A   1       0.000  -0.945   0.540  1.00  0.00           H
ATOM      8  HB1 ALA A   1       0.945   0.540  -2.070  1.00  0.00           H
ATOM      9  HB2 ALA A   1      -0.945   0.540  -2.070  1.00  0.00           H
ATOM     10  HB3 ALA A   1       0.000  -0.945  -2.080  1.00  0.00           H
TER
END
\"\"\"

from io import StringIO
pdb = app.PDBFile(StringIO(ala2))

# Force field & system
ff = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# We'll do an implicit-like setup by not adding solvent; keep it simple
system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

# Integrator
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
timestep = 2.0 * unit.femtoseconds
integrator = mm.LangevinIntegrator(temperature, friction, timestep)

# Context
platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# Minimize & write initial PDB
simulation.minimizeEnergy()
with open(TOPOLOGY_PATH, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

# Reporters: write a DCD every N steps
simulation.reporters.append(app.DCDReporter(TRAJ_PATH, 100))   # one frame each 100 steps
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, speed=True, separator=','))

print(f'[openmm_run] Writing to {TOPOLOGY_PATH} and {TRAJ_PATH}', file=sys.stderr)

# Run "forever-ish": batches of steps
while True:
    simulation.step(5000)  # ~10 ps per 5000 steps at 2 fs/frame (reporters decimate to every 100 steps)
    # small sleep to free CPU/GPU a bit if desired
    time.sleep(0.2)