# OpenMM "forever-ish" runner: writes
# topology.pdb + traj.dcd and keeps appending.
# Requirements: conda-forge::openmm
#
# Create & activate env:
#   mamba create -n live-md -c conda-forge openmm mdtraj -y
#   conda activate live-md
#
# Run:
#   python python/openmm_run.py
#
# Output files end up under ./data relative to
# project root.

import os, sys, time, argparse
from sys import stdout

# Parse command line arguments
parser = argparse.ArgumentParser(description='Run OpenMM molecular dynamics simulation with live trajectory output')
parser.add_argument('--timestep', type=float, default=2.0, help='Integration timestep in femtoseconds (default: 2.0)')
parser.add_argument('--temperature', type=float, default=300.0, help='Temperature in Kelvin (default: 300.0)')
parser.add_argument('--friction', type=float, default=1.0, help='Friction coefficient in 1/ps (default: 1.0)')
parser.add_argument('--steps-per-batch', type=int, default=5000, help='MD steps per batch (default: 5000)')
parser.add_argument('--report-interval', type=int, default=100, help='DCD frame output interval in steps (default: 100)')
parser.add_argument('--log-interval', type=int, default=1000, help='Console log interval in steps (default: 1000)')
parser.add_argument('--batch-sleep', type=float, default=0.2, help='Sleep time between batches in seconds (default: 0.2)')
args = parser.parse_args()

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
temperature = args.temperature * unit.kelvin
friction = args.friction / unit.picosecond
timestep = args.timestep * unit.femtoseconds
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
simulation.reporters.append(app.DCDReporter(TRAJ_PATH, args.report_interval))
simulation.reporters.append(app.StateDataReporter(stdout, args.log_interval, step=True, potentialEnergy=True, temperature=True, speed=True, separator=','))

print(f'[openmm_run] Configuration:', file=sys.stderr)
print(f'  Timestep: {args.timestep} fs', file=sys.stderr)
print(f'  Temperature: {args.temperature} K', file=sys.stderr)
print(f'  Friction: {args.friction} 1/ps', file=sys.stderr)
print(f'  Steps per batch: {args.steps_per_batch}', file=sys.stderr)
print(f'  DCD report interval: {args.report_interval} steps', file=sys.stderr)
print(f'  Log interval: {args.log_interval} steps', file=sys.stderr)
print(f'[openmm_run] Writing to {TOPOLOGY_PATH} and {TRAJ_PATH}', file=sys.stderr)

# Run "forever-ish": batches of steps
while True:
    simulation.step(args.steps_per_batch)
    time.sleep(args.batch_sleep)
