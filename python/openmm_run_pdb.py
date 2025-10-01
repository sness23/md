#!/usr/bin/env python
# OpenMM runner for arbitrary PDB/CIF structures
# Reads structure files and runs MD with implicit solvent

import os, sys, time, argparse, platform
from sys import stdout

# Parse command line arguments
parser = argparse.ArgumentParser(description='Run OpenMM MD simulation from PDB/CIF file')
parser.add_argument('--input', type=str, required=True, help='Input structure file (PDB or CIF)')
parser.add_argument('--timestep', type=float, default=2.0, help='Integration timestep in femtoseconds (default: 2.0)')
parser.add_argument('--temperature', type=float, default=300.0, help='Temperature in Kelvin (default: 300.0)')
parser.add_argument('--friction', type=float, default=1.0, help='Friction coefficient in 1/ps (default: 1.0)')
parser.add_argument('--steps-per-batch', type=int, default=5000, help='MD steps per batch (default: 5000)')
parser.add_argument('--report-interval', type=int, default=100, help='DCD frame output interval in steps (default: 100)')
parser.add_argument('--log-interval', type=int, default=1000, help='Console log interval in steps (default: 1000)')
parser.add_argument('--batch-sleep', type=float, default=0.2, help='Sleep time between batches in seconds (default: 0.2)')
parser.add_argument('--platform', type=str, default='auto', help='OpenMM platform: auto, CPU, OpenCL, Metal, CUDA (default: auto)')
parser.add_argument('--solvent', type=str, default='implicit', choices=['implicit', 'explicit'], help='Solvent model (default: implicit)')
args = parser.parse_args()

# Configure output paths
DATA_DIR = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)

# Use standard filenames for web viewer compatibility
TOPOLOGY_PATH = os.path.join(DATA_DIR, "topology.pdb")
TRAJ_PATH = os.path.join(DATA_DIR, "traj.dcd")

# --- Load structure and prepare system ---
from openmm import app
import openmm as mm
from openmm import unit

try:
    from pdbfixer import PDBFixer
    use_pdbfixer = True
except ImportError:
    print(f'[openmm_run_pdb] ⚠️  PDBFixer not available, limited structure cleanup', file=sys.stderr)
    use_pdbfixer = False

print(f'[openmm_run_pdb] Loading structure from {args.input}', file=sys.stderr)

# Load and fix structure using PDBFixer if available
if use_pdbfixer:
    try:
        fixer = PDBFixer(filename=args.input)
        print(f'[openmm_run_pdb] ✓ Structure loaded: {fixer.topology.getNumAtoms()} atoms', file=sys.stderr)

        # Find and fix issues
        print(f'[openmm_run_pdb] Finding and fixing structure issues...', file=sys.stderr)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)  # Remove ligands and waters
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)  # pH 7.0

        print(f'[openmm_run_pdb] ✓ Structure fixed: {fixer.topology.getNumAtoms()} atoms', file=sys.stderr)

        topology = fixer.topology
        positions = fixer.positions

    except Exception as e:
        print(f'[openmm_run_pdb] ✗ PDBFixer failed: {e}', file=sys.stderr)
        sys.exit(1)
else:
    # Fallback: simple loading without fixing
    try:
        if args.input.lower().endswith('.cif'):
            pdb = app.PDBxFile(args.input)
        elif args.input.lower().endswith('.pdb'):
            pdb = app.PDBFile(args.input)
        else:
            raise ValueError(f"Unsupported file format: {args.input}. Use .pdb or .cif")
        print(f'[openmm_run_pdb] ✓ Structure loaded: {pdb.topology.getNumAtoms()} atoms', file=sys.stderr)

        # Basic cleanup
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.deleteWater()

        topology = modeller.topology
        positions = modeller.positions

    except Exception as e:
        print(f'[openmm_run_pdb] ✗ Failed to load structure: {e}', file=sys.stderr)
        sys.exit(1)

# Choose forcefield and create system
print(f'[openmm_run_pdb] Creating system with {args.solvent} solvent...', file=sys.stderr)
forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')

try:
    if args.solvent == 'implicit':
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds
        )
    else:
        # Explicit solvent - need water forcefield and modeller
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = app.Modeller(topology, positions)
        modeller.addSolvent(forcefield, padding=1.0*unit.nanometers)
        print(f'[openmm_run_pdb] ✓ Added solvent: {modeller.topology.getNumAtoms()} total atoms', file=sys.stderr)

        topology = modeller.topology
        positions = modeller.positions

        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds
        )
    print(f'[openmm_run_pdb] ✓ System created successfully', file=sys.stderr)
except Exception as e:
    print(f'[openmm_run_pdb] ✗ Failed to create system: {e}', file=sys.stderr)
    print(f'[openmm_run_pdb] This may be due to:', file=sys.stderr)
    print(f'  - Non-standard residues or ligands', file=sys.stderr)
    print(f'  - Missing atoms in the structure', file=sys.stderr)
    print(f'  - Incompatible forcefield parameters', file=sys.stderr)
    sys.exit(1)

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
            print(f'[openmm_run_pdb] Using requested platform: {preferred}', file=sys.stderr)
            return plat, {}
        except Exception as e:
            print(f'[openmm_run_pdb] ⚠️  Platform {preferred} not available: {e}', file=sys.stderr)
            print(f'[openmm_run_pdb] Falling back to auto-selection', file=sys.stderr)

    # Auto-select: prefer GPU platforms on macOS
    available = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    print(f'[openmm_run_pdb] Available platforms: {available}', file=sys.stderr)

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
                    print(f'[openmm_run_pdb] Selected platform: {plat_name}{props_str}', file=sys.stderr)
                    return plat, props
                except Exception as e:
                    if props == prop_attempts[-1]:  # Last attempt for this platform
                        print(f'[openmm_run_pdb] Platform {plat_name} failed: {str(e).split(chr(10))[0]}', file=sys.stderr)
                    continue

    # Fallback to CPU if nothing else works
    print(f'[openmm_run_pdb] ⚠️  No GPU platform available, using CPU', file=sys.stderr)
    return mm.Platform.getPlatformByName('CPU'), {}

selected_platform, platform_props = select_platform(args.platform)

# Create simulation
simulation = app.Simulation(topology, system, integrator, selected_platform, platform_props)
simulation.context.setPositions(positions)

# Minimize energy
print(f'[openmm_run_pdb] Minimizing energy...', file=sys.stderr)
simulation.minimizeEnergy()

# Write initial topology
with open(TOPOLOGY_PATH, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

# Use atomic write wrapper for DCD to prevent partial reads
class AtomicDCDReporter(app.DCDReporter):
    """DCD reporter that writes to temp file and atomically renames."""
    def __init__(self, file, reportInterval, append=False, enforcePeriodicBox=None):
        self._temp_file = file + '.tmp'
        self._final_file = file
        super().__init__(self._temp_file, reportInterval, append, enforcePeriodicBox)

    def report(self, simulation, state):
        super().report(simulation, state)
        # Flush and atomically rename after each frame
        self._out.flush()
        os.fsync(self._out.fileno())
        # Copy temp to final (atomic on most filesystems)
        import shutil
        shutil.copy2(self._temp_file, self._final_file)

# Add reporters
simulation.reporters.append(AtomicDCDReporter(TRAJ_PATH, args.report_interval))
simulation.reporters.append(app.StateDataReporter(stdout, args.log_interval, step=True, potentialEnergy=True, temperature=True, speed=True, separator=','))

print(f'[openmm_run_pdb] Configuration:', file=sys.stderr)
print(f'  Platform: {simulation.context.getPlatform().getName()}', file=sys.stderr)
print(f'  Solvent: {args.solvent}', file=sys.stderr)
print(f'  Atoms: {topology.getNumAtoms()}', file=sys.stderr)
print(f'  Timestep: {args.timestep} fs', file=sys.stderr)
print(f'  Temperature: {args.temperature} K', file=sys.stderr)
print(f'  Friction: {args.friction} 1/ps', file=sys.stderr)
print(f'  Steps per batch: {args.steps_per_batch}', file=sys.stderr)
print(f'  DCD report interval: {args.report_interval} steps', file=sys.stderr)
print(f'  Log interval: {args.log_interval} steps', file=sys.stderr)
print(f'[openmm_run_pdb] Writing to {TOPOLOGY_PATH} and {TRAJ_PATH}', file=sys.stderr)
print(f'[openmm_run_pdb] Running simulation... Press Ctrl+C to stop', file=sys.stderr)

# Run forever
try:
    while True:
        simulation.step(args.steps_per_batch)
        time.sleep(args.batch_sleep)
except KeyboardInterrupt:
    print(f'\n[openmm_run_pdb] Simulation stopped by user', file=sys.stderr)
