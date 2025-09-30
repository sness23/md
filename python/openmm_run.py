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
from openmm import app
import openmm as mm
from openmm import unit

# Use OpenMM's built-in alanine dipeptide test system
from openmmtools import testsystems
test_system = testsystems.AlanineDipeptideImplicit()
system = test_system.system
pdb = test_system.topology
positions = test_system.positions

# Integrator
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
timestep = 2.0 * unit.femtoseconds
integrator = mm.LangevinIntegrator(temperature, friction, timestep)

# Context
platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb, system, integrator, platform)
simulation.context.setPositions(positions)

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