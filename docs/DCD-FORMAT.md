# DCD File Format Specification

## Overview

DCD (CHARMM/NAMD Binary Trajectory Format) is a binary file format for storing molecular dynamics trajectory coordinates. It's widely used by CHARMM, NAMD, X-PLOR, and OpenMM.

**Key Features:**
- Binary format (compact, fast I/O)
- Little-endian or big-endian (system-dependent)
- Single or double precision coordinates
- Fortran-style record markers (block sizes)
- Appendable (new frames added to end)

## File Structure

```
┌────────────────────────────────────┐
│  Header Block (84 bytes + markers) │
├────────────────────────────────────┤
│  Title Block (variable size)       │
├────────────────────────────────────┤
│  N-Atoms Block (4 bytes + markers) │
├────────────────────────────────────┤
│  Frame 1 Coordinates               │
│    - X coordinates                 │
│    - Y coordinates                 │
│    - Z coordinates                 │
├────────────────────────────────────┤
│  Frame 2 Coordinates               │
│    - X coordinates                 │
│    - Y coordinates                 │
│    - Z coordinates                 │
├────────────────────────────────────┤
│  ...                               │
├────────────────────────────────────┤
│  Frame N Coordinates               │
└────────────────────────────────────┘
```

## Fortran Record Markers

DCD uses Fortran unformatted I/O conventions where each "record" is bracketed by 4-byte integers indicating the record size:

```
[4 bytes: size] [data...] [4 bytes: size]
```

Both size markers should match. This allows reading/writing records as atomic units.

## Header Block

### Structure (92 bytes total)

```
Offset | Size | Type    | Description
-------|------|---------|-----------------------------------
0      | 4    | int32   | Block size marker (84)
4      | 4    | char[4] | Magic: "CORD" (CHARMM DCD)
8      | 4    | int32   | NSET: Number of frames
12     | 4    | int32   | ISTART: Starting timestep
16     | 4    | int32   | NSAVC: Timesteps between frames
20     | 4    | int32   | NSTEP: Total simulation timesteps
24     | 4    | int32   | NFREE: Degrees of freedom (unused)
28     | 20   | int32×5 | Reserved/unused
48     | 8    | float64 | DELTA: Timestep in AKMA units
56     | 36   | int32×9 | Reserved/version info
92     | 4    | int32   | Block end marker (84)
```

### Field Details

**Magic "CORD":**
- CHARMM format: "CORD"
- X-PLOR format: "VELD" (velocities)
- If magic != "CORD", may need byte-swapping

**NSET (Number of Frames):**
- Total frames in file
- Updated when appending frames
- Some writers leave as 0, requiring parsing entire file

**ISTART (Starting Step):**
- Simulation timestep of first frame
- Often 0 for trajectories starting from t=0

**NSAVC (Save Frequency):**
- Number of timesteps between saved frames
- Example: NSAVC=100 means save every 100 steps

**DELTA (Timestep):**
- Time between simulation steps in AKMA units
- AKMA: 20.45482706 fs
- **OpenMM writes this as double (8 bytes)**
- **Some implementations use float (4 bytes)** + 4 bytes padding

**Byte Order:**
- Little-endian (Intel/AMD): Most common
- Big-endian (older systems): Rare
- Detect by checking if magic reads correctly

### Example Header (Hex Dump)

```
00000000  54 00 00 00 43 4f 52 44  a5 30 00 00 64 00 00 00
          └─ 84      └─ CORD      └─ 12453   └─ 100
          (size)     (magic)       (nFrames)  (istart)

00000010  64 00 00 00 74 00 13 00  00 00 00 00 00 00 00 00
          └─ 100      └─ 1262708   (unused)
          (nsavc)     (nstep)

00000020  00 00 00 00 00 00 00 00  00 00 00 00 e3 90 27 3d
          (unused)                  └─ 2.87e-8 (delta)

00000030  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00
          (unused/version info)

00000050  00 00 00 00 18 00 00 00  54 00 00 00
          (unused)    └─ 24        └─ 84
                      (version?)   (end marker)
```

## Title Block

### Structure (Variable Size)

```
Offset | Size      | Type    | Description
-------|-----------|---------|---------------------------
0      | 4         | int32   | Block size
4      | 4         | int32   | NTITLE: Number of title lines
8      | 80×NTITLE | char    | Title strings (80 chars each)
8+80N  | 4         | int32   | Block end marker
```

### Title Format

- Each title line is exactly 80 ASCII characters
- Null-padded if shorter
- Common content:
  - Line 1: Software info ("Created by OpenMM")
  - Line 2: Creation timestamp
  - Additional lines: User-defined metadata

### Example

```
┌───────────────────────────────────────────────────────────────┐
│ Size: 164 bytes                                               │
├───────────────────────────────────────────────────────────────┤
│ NTITLE: 2                                                     │
├───────────────────────────────────────────────────────────────┤
│ "Created by OpenMM" + 66 null bytes                           │
│ "Created Mon Sep 29 21:36:35 2025" + 48 null bytes           │
└───────────────────────────────────────────────────────────────┘
```

## N-Atoms Block

### Structure (12 bytes total)

```
Offset | Size | Type  | Description
-------|------|-------|---------------------------
0      | 4    | int32 | Block size (4)
4      | 4    | int32 | NATOM: Number of atoms
8      | 4    | int32 | Block end marker (4)
```

**Important:** This is the **only** place the atom count is stored. Must read this before parsing coordinate frames.

## Coordinate Frames

Each frame contains 3 separate blocks (X, Y, Z coordinates):

### X-Coordinates Block

```
Offset | Size      | Type    | Description
-------|-----------|---------|---------------------------
0      | 4         | int32   | Block size (NATOM × 4)
4      | NATOM × 4 | float32 | X coordinates
4+N×4  | 4         | int32   | Block end marker
```

### Y-Coordinates Block

Same structure as X, immediately following.

### Z-Coordinates Block

Same structure as X and Y, immediately following.

### Coordinate Units

- **Standard DCD:** Angstroms (Å)
- **OpenMM:** Angstroms (Å)
- **CHARMM:** Angstroms (Å)
- **GROMACS (via conversion):** May be nm, check metadata

### Frame Layout

```
Frame N:
┌─────────────────────┐
│ [4] Block size      │
│ [4×NATOM] X coords  │ ← float32 array
│ [4] End marker      │
├─────────────────────┤
│ [4] Block size      │
│ [4×NATOM] Y coords  │ ← float32 array
│ [4] End marker      │
├─────────────────────┤
│ [4] Block size      │
│ [4×NATOM] Z coords  │ ← float32 array
│ [4] End marker      │
└─────────────────────┘

Total frame size: 3 × (8 + NATOM × 4) bytes
For 22 atoms: 3 × 96 = 288 bytes per frame
```

## Extended DCD Format

Some DCD variants include additional data:

### Unit Cell Information

If present, each frame starts with a unit cell block:

```
[4] Block size (48)
[8] double: a (cell length)
[8] double: gamma (angle)
[8] double: b (cell length)
[8] double: beta (angle)
[8] double: alpha (angle)
[8] double: c (cell length)
[4] End marker (48)
```

**Detection:** Check header flags or assume present if first block != expected coordinate size.

### Fixed Atoms

Some formats support "fixed" atoms that aren't written every frame:

```
NAMNF = NATOM - number_of_free_atoms
```

Fixed atom positions taken from first frame only.

## Parsing Algorithm

### Pseudocode

```python
def parse_dcd(file_path):
    with open(file_path, 'rb') as f:
        # Read header
        block_size = read_int32(f)
        magic = read_string(f, 4)
        nframes = read_int32(f)
        # ... read rest of header ...
        skip_to_offset(f, 4 + block_size)
        read_int32(f)  # End marker

        # Read title
        title_size = read_int32(f)
        skip(f, title_size)
        read_int32(f)  # End marker

        # Read natoms
        read_int32(f)  # Block size
        natoms = read_int32(f)
        read_int32(f)  # End marker

        # Read frames
        frames = []
        for frame_idx in range(nframes):
            x_coords = read_coord_block(f, natoms)
            y_coords = read_coord_block(f, natoms)
            z_coords = read_coord_block(f, natoms)
            frames.append({'x': x_coords, 'y': y_coords, 'z': z_coords})

    return frames, natoms

def read_coord_block(f, natoms):
    block_size = read_int32(f)
    coords = read_float32_array(f, natoms)
    end_marker = read_int32(f)
    assert block_size == end_marker == natoms * 4
    return coords
```

### JavaScript Implementation

```javascript
function parseDCD(arrayBuffer) {
    const view = new DataView(arrayBuffer)
    let offset = 0

    // Header
    const hdrSize = view.getInt32(offset, true); offset += 4
    const magic = String.fromCharCode(
        view.getUint8(offset),
        view.getUint8(offset+1),
        view.getUint8(offset+2),
        view.getUint8(offset+3)
    ); offset += 4

    if (magic !== 'CORD') {
        throw new Error('Not a CHARMM DCD file')
    }

    const nFrames = view.getInt32(offset, true); offset += 4

    // Skip to end of header
    offset = 4 + hdrSize
    offset += 4  // End marker

    // Title block
    const titleSize = view.getInt32(offset, true); offset += 4
    offset += titleSize  // Skip title content
    offset += 4  // End marker

    // N-atoms
    offset += 4  // Block size
    const nAtoms = view.getInt32(offset, true); offset += 4
    offset += 4  // End marker

    // Frames
    const frames = []
    for (let f = 0; f < nFrames; f++) {
        const coords = {x: [], y: [], z: []}

        // X
        offset += 4  // Block size
        for (let i = 0; i < nAtoms; i++) {
            coords.x.push(view.getFloat32(offset, true))
            offset += 4
        }
        offset += 4  // End marker

        // Y
        offset += 4
        for (let i = 0; i < nAtoms; i++) {
            coords.y.push(view.getFloat32(offset, true))
            offset += 4
        }
        offset += 4

        // Z
        offset += 4
        for (let i = 0; i < nAtoms; i++) {
            coords.z.push(view.getFloat32(offset, true))
            offset += 4
        }
        offset += 4

        frames.push(coords)
    }

    return frames
}
```

## Byte Order Detection

```javascript
function detectByteOrder(arrayBuffer) {
    const view = new DataView(arrayBuffer, 0, 8)

    // Try little-endian
    const magicLE = String.fromCharCode(
        view.getUint8(4),
        view.getUint8(5),
        view.getUint8(6),
        view.getUint8(7)
    )

    if (magicLE === 'CORD') return 'little'

    // Try big-endian
    const magicBE = String.fromCharCode(
        view.getUint8(7),
        view.getUint8(6),
        view.getUint8(5),
        view.getUint8(4)
    )

    if (magicBE === 'CORD') return 'big'

    throw new Error('Cannot determine byte order')
}
```

## Writing DCD Files

### Python (MDTraj)

```python
import mdtraj as md

traj = md.load('input.pdb')
# ... run simulation ...
traj.save_dcd('output.dcd')
```

### Python (OpenMM)

```python
from openmm.app import DCDReporter

simulation.reporters.append(
    DCDReporter('trajectory.dcd', 100)  # Save every 100 steps
)
```

### Appending Frames

DCD files can be opened in append mode:

```python
# Initial write
reporter = DCDReporter('traj.dcd', 100)

# Later: append to existing file
reporter = DCDReporter('traj.dcd', 100, append=True)
```

**Note:** When appending, NSET in header should be updated, but some writers don't. Parsers should count frames instead of trusting NSET.

## Common Pitfalls

### 1. Timestep Field Size

**Problem:** OpenMM uses `double` (8 bytes), but original spec uses `float` (4 bytes) + padding.

**Solution:** Always read 8 bytes at offset 48.

### 2. NSET = 0

**Problem:** Some writers set NSET=0 and never update it.

**Solution:** Parse until EOF, don't rely on NSET.

### 3. Block Size Mismatches

**Problem:** Start and end markers don't match (corrupted file).

**Solution:** Validate markers, throw error if mismatch.

### 4. Partial Frames

**Problem:** File truncated mid-frame (simulation crashed).

**Solution:** Catch EOF exceptions, return complete frames only.

### 5. Unit Cell Data

**Problem:** Extended format with unit cell causes offset errors.

**Solution:** Check header flags or detect unexpectedly large first block.

## File Size Estimation

For N frames with M atoms:

```
Header: ~300 bytes
Per frame: 3 × (8 + M × 4) bytes

Total ≈ 300 + N × (24 + 12M) bytes
```

**Example:** 10,000 frames, 22 atoms
```
300 + 10000 × (24 + 264) = 300 + 2,880,000 = ~2.88 MB
```

## Conversion Tools

**MDTraj:**
```python
import mdtraj as md
traj = md.load('traj.dcd', top='topology.pdb')
traj.save('traj.xtc')  # Convert to XTC
traj.save('traj.h5')   # Convert to HDF5
```

**MDAnalysis:**
```python
import MDAnalysis as mda
u = mda.Universe('topology.pdb', 'traj.dcd')
with mda.Writer('output.xtc', u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u.atoms)
```

**VMD:**
```tcl
mol load pdb topology.pdb dcd trajectory.dcd
animate write xtc output.xtc
```

## References

- [CHARMM DCD Format](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html)
- [MDAnalysis DCD Reader](https://docs.mdanalysis.org/stable/documentation_pages/coordinates/DCD.html)
- [VMD Molfile Plugin](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/)
- [OpenMM DCDReporter](http://docs.openmm.org/latest/api-python/generated/openmm.app.dcdreporter.DCDReporter.html)

## Version History

- **Original CHARMM:** Single-precision, fixed format
- **X-PLOR Extension:** Added velocity support ("VELD")
- **NAMD Extension:** Added unit cell information
- **OpenMM:** Double-precision timestep, standard otherwise