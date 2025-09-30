# OpenMM Runner Documentation

## Overview

`python/openmm_run.py` is a continuous molecular dynamics simulation runner that uses OpenMM's alanine dipeptide test system. It writes trajectory data to DCD format while the simulation runs, enabling live visualization in the browser.

## Current Implementation

The script currently uses OpenMM's built-in `testsystems.AlanineDipeptideImplicit()` which provides:
- Pre-configured alanine dipeptide topology
- Implicit solvent (GBSA)
- Force field parameters already applied
- Initial coordinates

### Output Files

- `data/topology.pdb` - Static structure file (written once after energy minimization)
- `data/traj.dcd` - Growing trajectory file (appended every 100 simulation steps)

### Simulation Parameters

```python
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
timestep = 2.0 * unit.femtoseconds
integrator = mm.LangevinIntegrator(temperature, friction, timestep)
```

**Frame Output:**
- Saves every 100 steps (200 fs per frame)
- Continuous batches of 5000 steps (~10 ps per batch)
- 0.2s sleep between batches

## Using Custom PDB Files

To run simulations on your own structures, you need to replace the test system with PDB loading code:

### Basic Example (Implicit Solvent)

```python
from openmm.app import PDBFile, ForceField
import openmm as mm
from openmm import unit

# Load your structure
pdb = PDBFile('path/to/your/structure.pdb')

# Choose force field
forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')

# Create system
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
    implicitSolvent=app.GBn2
)

# Setup integrator and simulation
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
```

### PDB File Requirements

Your PDB file must:
1. **Have all atoms** - Missing atoms will cause force field matching errors
2. **Use standard residue names** - Force fields expect standard nomenclature
3. **Be properly capped** - Protein termini need capping groups (ACE/NME for AMBER)
4. **Match the force field** - Use AMBER-compatible naming for AMBER force fields

### Common Force Field Options

**AMBER force fields:**
```python
# All-atom protein with implicit solvent
forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')

# All-atom protein with explicit water
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
```

**CHARMM force fields:**
```python
# CHARMM36 protein
forcefield = ForceField('charmm36.xml')

# CHARMM36 with explicit water
forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
```

### Explicit Solvent Example

```python
from openmm.app import PDBFile, ForceField, Modeller
import openmm as mm
from openmm import unit

# Load structure
pdb = PDBFile('protein.pdb')

# Add solvent
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=1.0*unit.nanometers)

# Create system with PME for long-range electrostatics
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=app.HBonds
)

# Integrator and simulation
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
```

### Handling Force Field Errors

**Error: "No template found for residue X"**

Solutions:
- Ensure residue names match force field conventions (e.g., HIS vs HSD/HSE/HSP)
- Add missing atoms using a tool like `pdbfixer`
- Use a different force field that supports your structure
- Manually add custom residue templates

**Using PDBFixer:**
```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='input.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

PDBFile.writeFile(fixer.topology, fixer.positions, open('fixed.pdb', 'w'))
```

## Platform Selection

The script currently uses CPU platform for compatibility:
```python
platform = mm.Platform.getPlatformByName('CPU')
```

**For faster simulations, use CUDA or OpenCL:**
```python
# CUDA (NVIDIA GPUs)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(topology, system, integrator, platform, properties)

# OpenCL (AMD/NVIDIA/Intel)
platform = mm.Platform.getPlatformByName('OpenCL')
properties = {'OpenCLPrecision': 'mixed'}
simulation = app.Simulation(topology, system, integrator, platform, properties)
```

## Modifying Simulation Parameters

### Change Temperature
```python
temperature = 310 * unit.kelvin  # 37°C for body temperature
integrator = mm.LangevinIntegrator(temperature, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
```

### Change Timestep
```python
timestep = 1.0 * unit.femtoseconds  # Smaller timestep for stability
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, timestep)
```

### Change Frame Output Frequency
```python
# Save every 50 steps instead of 100
simulation.reporters.append(app.DCDReporter(TRAJ_PATH, 50))
```

### Change Simulation Batch Size
```python
# Run 10,000 steps per batch instead of 5,000
while True:
    simulation.step(10000)
    time.sleep(0.2)
```

## Environment Variables

- `DATA_DIR` - Override output directory (default: `./data`)

Example:
```bash
DATA_DIR=/path/to/output python python/openmm_run.py
```

## Performance Tips

1. **Use GPU platforms** - 10-100x faster than CPU
2. **Adjust precision** - Use `mixed` precision for speed/accuracy balance
3. **Increase batch size** - Reduce Python overhead (but increases latency for live viewing)
4. **Remove state reporter** - The CSV output slows things down slightly
5. **Use constraints** - `constraints=app.HBonds` allows 2fs timestep

## Troubleshooting

### Simulation is unstable / crashes
- Reduce timestep to 1.0 fs
- Run longer energy minimization
- Check for steric clashes in input structure
- Add constraints: `constraints=app.AllBonds`

### Very slow performance
- Switch from CPU to CUDA/OpenCL platform
- Reduce system size (fewer atoms/smaller solvent box)
- Increase nonbonded cutoff (if using PME)

### Memory issues
- Reduce DCD output frequency
- Restart with new trajectory file periodically
- Use smaller solvent box for explicit solvent

## References

- [OpenMM User Guide](http://docs.openmm.org/latest/userguide/)
- [OpenMM Force Fields](http://docs.openmm.org/latest/userguide/application/02_running_sims.html#force-fields)
- [PDBFixer Documentation](https://github.com/openmm/pdbfixer)