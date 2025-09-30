# Frontend Architecture

## Overview

The frontend is a single-page application (`public/index.html`) that combines Babylon.js for 3D rendering with custom parsers for PDB and DCD file formats. All processing happens client-side in the browser.

## Technology Stack

- **Babylon.js 8.29+** - 3D rendering engine
- **Vanilla JavaScript** - No framework dependencies
- **ES6 Modules** - Modern JavaScript features
- **WebGL 2.0** - GPU-accelerated graphics

## File Structure

```javascript
public/index.html
├── HTML structure (canvas, sidebar, controls)
├── CSS (inline styles)
└── JavaScript (ES6 module)
    ├── Configuration & State
    ├── Babylon.js Scene Setup
    ├── PDB Parser
    ├── DCD Parser
    ├── Atom Rendering (CPK)
    ├── Animation Loop
    ├── Polling Mechanism
    └── Event Handlers
```

## Core Components

### 1. State Management

```javascript
// Global state variables
let polling = false        // Auto-detection enabled?
let playing = true         // Animation playing?
let currentFrame = 0       // Current frame index
let frames = []            // Array of frame coordinate objects
let atoms = []             // Array of atom objects {element, x, y, z}
let atomSpheres = []       // Babylon.js mesh references
let currentBytes = 0       // Last known DCD file size
```

### 2. Babylon.js Scene

#### Scene Setup
```javascript
const engine = new BABYLON.Engine(canvas, true)
const scene = new BABYLON.Scene(engine)
scene.clearColor = new BABYLON.Color3(0.043, 0.071, 0.125)  // Dark blue
```

#### Camera Configuration
```javascript
const camera = new BABYLON.ArcRotateCamera(
  'camera',
  0,                    // Alpha (horizontal rotation)
  Math.PI / 3,          // Beta (vertical rotation)
  30,                   // Radius (distance from target)
  BABYLON.Vector3.Zero(), // Target point
  scene
)
camera.wheelPrecision = 50      // Zoom sensitivity
camera.lowerRadiusLimit = 5     // Min zoom distance
camera.upperRadiusLimit = 100   // Max zoom distance
```

#### Lighting
```javascript
const light = new BABYLON.HemisphericLight(
  'light',
  new BABYLON.Vector3(0, 1, 0),  // Light direction (from above)
  scene
)
light.intensity = 0.8
```

#### Render Loop
```javascript
engine.runRenderLoop(() => scene.render())
window.addEventListener('resize', () => engine.resize())
```

### 3. PDB Parser

Parses ATOM and HETATM lines from PDB text format:

```javascript
function parsePDB(text) {
  const lines = text.split('\n')
  const parsedAtoms = []

  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      // Column-based parsing per PDB spec
      const element = line.substring(76, 78).trim() ||
                      line.substring(12, 14).trim().charAt(0)
      const x = parseFloat(line.substring(30, 38))
      const y = parseFloat(line.substring(38, 46))
      const z = parseFloat(line.substring(46, 54))

      parsedAtoms.push({ element, x, y, z })
    }
  }

  return parsedAtoms
}
```

**PDB Column Format:**
- Columns 1-6: Record type (ATOM/HETATM)
- Columns 13-14: Atom name
- Columns 31-38: X coordinate (Å)
- Columns 39-46: Y coordinate (Å)
- Columns 47-54: Z coordinate (Å)
- Columns 77-78: Element symbol

### 4. DCD Parser

Parses binary CHARMM/NAMD DCD trajectory files:

```javascript
function parseDCD(arrayBuffer) {
  const view = new DataView(arrayBuffer)
  let offset = 0

  // Block 1: Main header (84 bytes)
  const hdrSize = view.getInt32(offset, true); offset += 4
  const magic = new TextDecoder('ascii').decode(
    new Uint8Array(arrayBuffer, offset, 4)
  ); offset += 4
  const nFrames = view.getInt32(offset, true); offset += 4
  // ... skip to end of block
  offset = 4 + hdrSize
  offset += 4  // Skip end marker

  // Block 2: Title block (variable size)
  const titleSize = view.getInt32(offset, true); offset += 4
  offset += titleSize  // Skip entire title content
  offset += 4  // Skip end marker

  // Block 3: Number of atoms (4 bytes)
  offset += 4  // Skip block size marker
  const nAtoms = view.getInt32(offset, true); offset += 4
  offset += 4  // Skip end marker

  // Read coordinate frames
  const parsedFrames = []
  for (let f = 0; f < nFrames; f++) {
    const coords = { x: [], y: [], z: [] }

    // X coordinates
    offset += 4  // Block size
    for (let i = 0; i < nAtoms; i++) {
      coords.x.push(view.getFloat32(offset, true))
      offset += 4
    }
    offset += 4  // End marker

    // Y coordinates (same pattern)
    // Z coordinates (same pattern)

    parsedFrames.push(coords)
  }

  return parsedFrames
}
```

**DCD File Structure:**
```
[4 bytes] Block size
[84 bytes] Header data (magic, nFrames, timestep info)
[4 bytes] End marker

[4 bytes] Block size
[variable] Title strings
[4 bytes] End marker

[4 bytes] Block size (4)
[4 bytes] nAtoms
[4 bytes] End marker

For each frame:
  [4 bytes] Block size (nAtoms * 4)
  [nAtoms * 4 bytes] X coordinates (float32)
  [4 bytes] End marker

  [4 bytes] Block size (nAtoms * 4)
  [nAtoms * 4 bytes] Y coordinates (float32)
  [4 bytes] End marker

  [4 bytes] Block size (nAtoms * 4)
  [nAtoms * 4 bytes] Z coordinates (float32)
  [4 bytes] End marker
```

**Important Notes:**
- All integers are little-endian (true flag in DataView)
- Coordinates are single-precision floats (Float32)
- Block size markers before/after each block
- OpenMM's DCD uses double-precision timestep (8 bytes)

### 5. CPK Rendering

#### Color Scheme
```javascript
const CPK_COLORS = {
  'H': [1.0, 1.0, 1.0],    // White
  'C': [0.5, 0.5, 0.5],    // Gray
  'N': [0.2, 0.2, 1.0],    // Blue
  'O': [1.0, 0.2, 0.2],    // Red
  'S': [1.0, 1.0, 0.0],    // Yellow
  'P': [1.0, 0.5, 0.0],    // Orange
}
```

#### Atomic Radii (Angstroms)
```javascript
const CPK_RADII = {
  'H': 0.25,
  'C': 0.7,
  'N': 0.65,
  'O': 0.6,
  'S': 1.0,
  'P': 1.0
}
```

#### Sphere Creation
```javascript
function createAtoms() {
  // Dispose old spheres
  atomSpheres.forEach(s => s.dispose())
  atomSpheres = []

  for (let i = 0; i < atoms.length; i++) {
    const atom = atoms[i]
    const radius = CPK_RADII[atom.element] || 0.7
    const color = CPK_COLORS[atom.element] || [0.8, 0.8, 0.8]

    // Create sphere mesh
    const sphere = BABYLON.MeshBuilder.CreateSphere(
      `atom${i}`,
      { diameter: radius * 2, segments: 12 },
      scene
    )

    // Create material
    const material = new BABYLON.StandardMaterial(`mat${i}`, scene)
    material.diffuseColor = new BABYLON.Color3(...color)
    material.specularColor = new BABYLON.Color3(0.3, 0.3, 0.3)
    sphere.material = material

    // Set initial position
    sphere.position.x = atom.x
    sphere.position.y = atom.y
    sphere.position.z = atom.z

    atomSpheres.push(sphere)
  }
}
```

**Performance Note:** Using 12 segments keeps polygon count low (~144 triangles per sphere). For 22 atoms that's ~3,168 triangles total, easily handled by modern GPUs.

### 6. Animation Loop

Uses `requestAnimationFrame` for smooth rendering with time-based frame advancement:

```javascript
let lastFrameTime = 0
const frameInterval = 50  // 20 fps (50ms per frame)

function animate(timestamp) {
  // Update at 20fps regardless of render rate
  if (timestamp - lastFrameTime >= frameInterval) {
    if (playing && frames.length > 0) {
      currentFrame = (currentFrame + 1) % frames.length
      updateFrame()
    }
    lastFrameTime = timestamp
  }

  requestAnimationFrame(animate)
}
```

**Why 20 fps?**
- Matches typical MD visualization framerates
- Reduces CPU load vs 60fps
- Smooth enough for molecular motion perception
- Configurable via `frameInterval`

#### Frame Update
```javascript
function updateFrame() {
  if (frames.length === 0) return

  const frame = frames[currentFrame]

  // Update all atom positions
  for (let i = 0; i < atomSpheres.length && i < frame.x.length; i++) {
    atomSpheres[i].position.x = frame.x[i]
    atomSpheres[i].position.y = frame.y[i]
    atomSpheres[i].position.z = frame.z[i]
  }

  // Update UI
  frameInfoEl.textContent = `${currentFrame + 1}/${frames.length}`
}
```

### 7. Polling Mechanism

Detects file growth and reloads trajectory:

```javascript
async function checkForNewFrames() {
  const cfg = await getConfig()
  const dcdUrl = `${cfg.mdsrvUrl}/data/${cfg.trajectory}`
  const size = await headSize(dcdUrl)  // HEAD request for Content-Length

  if (size > currentBytes) {
    const oldFrameCount = frames.length
    currentBytes = size

    // Reload entire DCD file
    const dcdRes = await fetch(dcdUrl + '?t=' + Date.now())  // Cache bust
    const dcdBuffer = await dcdRes.arrayBuffer()
    frames = parseDCD(dcdBuffer)

    const newFrameCount = frames.length - oldFrameCount
    if (newFrameCount > 0) {
      statusEl.textContent = `loaded ${newFrameCount} new frames`
    }
  }
}
```

**Polling Strategy:**
- Uses HEAD request to check file size (lightweight)
- Only fetches full file when size increases
- Cache-busting with timestamp query parameter
- Default interval: 5 seconds (user-configurable)

**Known Limitation:** Reloads entire DCD file each time. For better performance, could:
- Parse only new bytes (requires byte-offset tracking)
- Use server-side frame extraction API
- Implement WebSocket streaming

### 8. Event Handlers

```javascript
// Play/Pause toggle
playPauseBtn.onclick = () => {
  playing = !playing
  playPauseBtn.textContent = playing ? 'Pause' : 'Play'
}

// Polling toggle
togglePollBtn.onclick = () => {
  polling = !polling
  togglePollBtn.textContent = polling ? 'Stop Polling' : 'Start Polling'

  if (polling) {
    const secs = Math.max(1, Number(pollSecsInput.value) || 5)
    pollTimer = setInterval(checkForNewFrames, secs * 1000)
  } else if (pollTimer) {
    clearInterval(pollTimer)
    pollTimer = null
  }
}
```

## Data Flow

```
┌──────────────┐
│  Load Page   │
└──────┬───────┘
       │
       ↓
┌──────────────────────┐
│ Fetch /api/config    │ → Get server URLs
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ Fetch topology.pdb   │ → Parse atoms
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ createAtoms()        │ → Babylon.js spheres
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ Fetch traj.dcd       │ → Parse frames
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ updateFrame()        │ → Set initial positions
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ animate() loop       │ ←─────┐
└──────┬───────────────┘        │
       │                        │
       ↓                        │
  [playing?] ──yes──→ advance frame
       │                        │
       └────────────────────────┘


┌──────────────────────┐
│ Polling timer        │ (every N seconds)
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ checkForNewFrames()  │ → HEAD request
└──────┬───────────────┘
       │
   [size grew?]
       │
       yes
       ↓
┌──────────────────────┐
│ Reload DCD           │ → Parse new frames
└──────┬───────────────┘
       │
       ↓
┌──────────────────────┐
│ Continue animation   │ → Show new frames
└──────────────────────┘
```

## Performance Characteristics

### Memory Usage

**Per atom:**
- Atom data: ~64 bytes (object + properties)
- Sphere mesh: ~8 KB (geometry + material)
- Per frame coordinate: 12 bytes (3x float32)

**For alanine dipeptide (22 atoms, 10k frames):**
- Atoms: ~1.4 KB
- Spheres: ~176 KB
- Frames: ~2.6 MB
- **Total: ~2.8 MB**

### CPU Usage

- **Idle (not playing):** Near zero
- **Playing animation:** ~2-5% CPU (20fps updates)
- **Parsing DCD:** Spikes to 20-50% during load
- **Rendering:** GPU-accelerated, minimal CPU

### Network Usage

- **Initial load:** ~2-3 MB (topology + trajectory)
- **HEAD polls:** <1 KB per request
- **Trajectory reload:** Full file size (grows with simulation)

## Browser Limitations

### Memory
- Modern browsers: ~2 GB per tab
- Practical limit: ~1M frames (depends on atom count)
- Exceeding causes slowdown or crash

### WebGL
- Max texture size: 8192x8192 (GPU dependent)
- Max vertex count: ~16M triangles
- Our usage: ~3K triangles << limits

### File Size
- Fetch API: No hard limits
- Parse time increases linearly with file size
- Large files (>100 MB) may cause UI freeze during parse

## Extending the Viewer

### Adding New Atom Types

```javascript
CPK_COLORS['Fe'] = [0.88, 0.40, 0.20]  // Iron (orange-brown)
CPK_RADII['Fe'] = 1.4
```

### Changing Representation

**Stick model (bonds):**
```javascript
function createBond(atom1, atom2) {
  const distance = BABYLON.Vector3.Distance(
    new BABYLON.Vector3(atom1.x, atom1.y, atom1.z),
    new BABYLON.Vector3(atom2.x, atom2.y, atom2.z)
  )

  const cylinder = BABYLON.MeshBuilder.CreateCylinder('bond', {
    diameter: 0.2,
    height: distance
  }, scene)

  // Position and orient cylinder between atoms
  // ...
}
```

**Ball-and-stick:**
```javascript
// Reduce sphere sizes
CPK_RADII['C'] = 0.35  // Half size
// Add bonds between connected atoms
```

### Adding Measurements

```javascript
function measureDistance(atom1Idx, atom2Idx) {
  const frame = frames[currentFrame]
  const dx = frame.x[atom1Idx] - frame.x[atom2Idx]
  const dy = frame.y[atom1Idx] - frame.y[atom2Idx]
  const dz = frame.z[atom1Idx] - frame.z[atom2Idx]
  return Math.sqrt(dx*dx + dy*dy + dz*dz)
}
```

### Custom Coloring

**By residue:**
```javascript
function getResidueColor(atomIndex) {
  const residueId = atoms[atomIndex].residue
  return colorPalette[residueId % colorPalette.length]
}
```

**By temperature factor:**
```javascript
function getBFactorColor(value) {
  // Gradient from blue (low) to red (high)
  const t = (value - minB) / (maxB - minB)
  return [t, 0, 1-t]
}
```

## Debugging

### Browser Console

Useful logging is included:
```javascript
console.log('Loaded atoms:', atoms.length)
console.log('DCD buffer size:', arrayBuffer.byteLength)
console.log('DCD header:', { magic, nFrames, nAtoms })
console.log('Parsed frames:', frames.length)
```

### Common Issues

**Atoms not visible:**
- Check camera position/target
- Verify atom coordinates are reasonable
- Check scene.clearColor isn't same as sphere colors

**Animation not smooth:**
- Reduce frameInterval
- Check CPU usage (may be maxed)
- Reduce atom count or sphere segments

**Parsing errors:**
- Check DCD file integrity: `python verify_traj.py`
- Verify byte order (little-endian)
- Check for partial file reads

## Future Enhancements

1. **Incremental DCD loading** - Parse only new frames
2. **Web Workers** - Move parsing off main thread
3. **Instanced rendering** - Single mesh for all atoms
4. **LOD (Level of Detail)** - Reduce segments when zoomed out
5. **Selection system** - Click to select atoms
6. **Measurement tools** - Distances, angles, dihedrals
7. **Trajectory controls** - Seek bar, frame stepping
8. **Export** - Screenshot, video recording
9. **Multiple trajectories** - Compare simulations
10. **Protein cartoon** - Secondary structure rendering