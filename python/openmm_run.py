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

import os, sys, time, argparse, platform
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
parser.add_argument('--platform', type=str, default='auto', help='OpenMM platform: auto, CPU, OpenCL, Metal, CUDA (default: auto)')
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

# Platform selection with GPU preference on macOS
def select_platform(preferred='auto'):
    """Select best available OpenMM platform, preferring GPU on macOS."""
    if preferred.lower() != 'auto':
        try:
            plat = mm.Platform.getPlatformByName(preferred)
            print(f'[openmm_run] Using requested platform: {preferred}', file=sys.stderr)
            return plat, {}
        except Exception as e:
            print(f'[openmm_run] ⚠️  Platform {preferred} not available: {e}', file=sys.stderr)
            print(f'[openmm_run] Falling back to auto-selection', file=sys.stderr)

    # Auto-select: prefer GPU platforms on macOS
    available = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    print(f'[openmm_run] Available platforms: {available}', file=sys.stderr)

    is_macos = platform.system() == 'Darwin'

    # Priority order: GPU platforms first (Metal > OpenCL > CUDA), then CPU
    priority = ['Metal', 'OpenCL', 'CUDA', 'CPU', 'Reference'] if is_macos else ['CUDA', 'OpenCL', 'Metal', 'CPU', 'Reference']

    for plat_name in priority:
        if plat_name in available:
            plat = mm.Platform.getPlatformByName(plat_name)

            # Try different property configurations for each platform
            prop_attempts = []
            if plat_name == 'OpenCL':
                # Try multiple OpenCL configurations
                prop_attempts = [
                    {'OpenCLPrecision': 'mixed'},
                    {'OpenCLPrecision': 'single'},
                    {},  # default
                ]
            elif plat_name == 'CUDA':
                prop_attempts = [
                    {'CudaPrecision': 'mixed'},
                    {},
                ]
            else:
                prop_attempts = [{}]

            for props in prop_attempts:
                try:
                    # Test if platform works by creating a minimal context
                    test_sys = mm.System()
                    test_sys.addParticle(1.0)
                    test_int = mm.VerletIntegrator(1.0 * unit.femtoseconds)
                    test_ctx = mm.Context(test_sys, test_int, plat, props)
                    del test_ctx, test_int

                    props_str = f' with {props}' if props else ''
                    print(f'[openmm_run] Selected platform: {plat_name}{props_str}', file=sys.stderr)
                    return plat, props
                except Exception as e:
                    if props == prop_attempts[-1]:  # Last attempt for this platform
                        print(f'[openmm_run] Platform {plat_name} failed: {str(e).split(chr(10))[0]}', file=sys.stderr)
                    continue

    # Fallback to CPU if nothing else works
    print(f'[openmm_run] ⚠️  No GPU platform available, using CPU', file=sys.stderr)
    return mm.Platform.getPlatformByName('CPU'), {}

selected_platform, platform_props = select_platform(args.platform)

# Context
simulation = app.Simulation(pdb, system, integrator, selected_platform, platform_props)
simulation.context.setPositions(positions)

# Minimize & write initial PDB
simulation.minimizeEnergy()
with open(TOPOLOGY_PATH, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

# Reporters: write a DCD every N steps
simulation.reporters.append(app.DCDReporter(TRAJ_PATH, args.report_interval))
simulation.reporters.append(app.StateDataReporter(stdout, args.log_interval, step=True, potentialEnergy=True, temperature=True, speed=True, separator=','))

print(f'[openmm_run] Configuration:', file=sys.stderr)
print(f'  Platform: {simulation.context.getPlatform().getName()}', file=sys.stderr)
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
