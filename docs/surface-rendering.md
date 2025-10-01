# Molecular Surface Rendering Implementation Plan

## Goal
Add transparent molecular surface overlay that shows the protein envelope while allowing the ligand (HETATM) to be visible inside.

## Approach

### 1. Surface Calculation (Backend Option - Recommended)

Use Python to pre-calculate surfaces and serve them:

**Libraries:**
- **MSMS** (Molecular Surface package) - fast, high quality
- **FreeSASA** - simple, good for SASA
- **PyMOL API** - if available, can generate surfaces
- **MDTraj** - has shrake_rupley for SASA

**Python service endpoint:**
```python
# New endpoint in server.js → Python microservice
@app.route('/api/surface/<frame>')
def get_surface(frame):
    # Load frame from trajectory
    # Calculate surface using MSMS/FreeSASA
    # Return vertices + faces as JSON/binary
    return {
        'vertices': [...],  # Nx3 array
        'faces': [...],     # Mx3 triangle indices
        'normals': [...]    # Nx3 normals
    }
```

### 2. Surface Calculation (Frontend Option)

Calculate in JavaScript using marching cubes:

**Libraries:**
- **three-bmfont-text** - Has marching cubes implementation
- **isosurface** - Fast marching cubes library
- **Custom implementation** - Roll our own

**Algorithm:**
1. Create 3D grid around protein
2. For each grid point, calculate distance to nearest atom (using van der Waals radii)
3. Use marching cubes to extract isosurface at probe radius (1.4Å for water)
4. Generate mesh with vertices and faces

**Pros/Cons:**
- ✅ No backend dependency
- ✅ Real-time updates as simulation runs
- ❌ Computationally expensive (slow for large proteins)
- ❌ Lower quality than MSMS

### 3. Mesh Rendering in Babylon.js

```javascript
// Create custom mesh from surface data
const surfaceMesh = new BABYLON.Mesh('surface', scene)

const vertexData = new BABYLON.VertexData()
vertexData.positions = surfaceVertices  // Flattened [x,y,z,x,y,z,...]
vertexData.indices = surfaceFaces       // Flattened [i1,i2,i3,i1,i2,i3,...]
vertexData.normals = surfaceNormals     // For smooth shading

vertexData.applyToMesh(surfaceMesh)

// Transparent material
const surfaceMaterial = new BABYLON.StandardMaterial('surfaceMat', scene)
surfaceMaterial.alpha = 0.3  // 30% opacity
surfaceMaterial.diffuseColor = new BABYLON.Color3(0.7, 0.7, 0.9)  // Light blue
surfaceMaterial.backFaceCulling = false  // Render both sides
surfaceMesh.material = surfaceMaterial
```

### 4. Filtering ATOM vs HETATM

In PDB parser, track atom types:

```javascript
function parsePDB(text) {
  const lines = text.split('\n')
  const parsedAtoms = []
  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const recordType = line.substring(0, 6).trim()
      const element = line.substring(76, 78).trim() || line.substring(12, 14).trim().charAt(0)
      const x = parseFloat(line.substring(30, 38))
      const y = parseFloat(line.substring(38, 46))
      const z = parseFloat(line.substring(46, 54))

      parsedAtoms.push({
        element, x, y, z,
        isProtein: recordType === 'ATOM'  // Track if it's protein
      })
    }
  }
  return parsedAtoms
}

// Filter for surface calculation
const proteinAtoms = atoms.filter(a => a.isProtein)
```

## Implementation Steps

### Phase 1: Basic Static Surface (Prototype)
1. Add surface calculation library (recommend: try frontend first for simplicity)
2. Calculate surface once on load for frame 0
3. Render as transparent mesh
4. Add keyboard toggle ('S' key to show/hide surface)

### Phase 2: Dynamic Updates
1. Recalculate surface every N frames (e.g., every 10 frames)
2. Smooth interpolation between surface updates
3. Show loading indicator during calculation

### Phase 3: Optimization
1. Move calculation to Web Worker (non-blocking)
2. Cache surfaces for frames
3. Use lower resolution for real-time, high-res on demand
4. Implement LOD (level of detail) based on camera distance

### Phase 4: Advanced Features
1. Multiple surface types (SAS, SES, VDW)
2. Color by hydrophobicity or electrostatics
3. Cavity detection and highlighting
4. Surface area calculations displayed in HUD

## Recommended Libraries

### For Surface Calculation:
1. **3Dmol.js** - Has built-in surface generation we could adapt
2. **NGL Viewer** - Open source, high-quality surfaces
3. **EDTSurf** - Fast Euclidean distance transform method
4. **msms.js** - Port of MSMS to JavaScript (if it exists)

### Alternative: Use Existing Viewer Libraries
Instead of rolling our own, we could:
- Integrate **3Dmol.js** surface code into our Babylon.js viewer
- Use **NGL Viewer** for surface, Babylon for trajectory playback
- Hybrid approach: Generate surface mesh externally, import as OBJ/STL

## Performance Considerations

For a 4000-atom protein like 1ERM:
- **Surface calculation**: 100-500ms (frontend) or 10-50ms (Python/MSMS)
- **Mesh generation**: 10-50ms
- **Rendering**: Real-time (60 fps) even with transparent overlay

**Optimization strategies:**
1. Calculate surface at lower temporal resolution (every 10th frame)
2. Use coarse grid for real-time, fine grid for snapshots
3. Adaptive quality based on FPS
4. Pre-calculate surfaces for entire trajectory (background process)

## UI Controls

Add to keyboard controls:
- **S** - Toggle surface visibility
- **[** / **]** - Decrease/increase surface transparency
- **Shift+S** - Switch surface type (SAS/SES/VDW)
- **T** - Toggle surface color scheme

Add to HUD (when visible):
- Surface area (Å²)
- Volume (Å³)
- Surface quality setting

## Example Code Snippet

```javascript
// Marching cubes surface calculation
function calculateProteinSurface(proteinAtoms, probeRadius = 1.4, gridSpacing = 0.5) {
  // 1. Create 3D grid
  const bounds = getBoundingBox(proteinAtoms)
  const grid = create3DGrid(bounds, gridSpacing)

  // 2. Calculate signed distance field
  for (let i = 0; i < grid.length; i++) {
    const point = grid[i]
    let minDist = Infinity

    for (const atom of proteinAtoms) {
      const dist = distance(point, atom) - VDW_RADII[atom.element] - probeRadius
      minDist = Math.min(minDist, dist)
    }

    grid[i].value = minDist
  }

  // 3. Extract isosurface at zero-crossing
  const surface = marchingCubes(grid, 0.0)

  return surface
}

// Toggle with keyboard
window.addEventListener('keydown', e => {
  if (e.key.toLowerCase() === 's') {
    surfaceVisible = !surfaceVisible
    surfaceMesh.setEnabled(surfaceVisible)
  }
})
```

## Testing Plan

1. **Small system** (alanine dipeptide, ~20 atoms) - verify algorithm correctness
2. **Medium system** (1ERM0, 1500 atoms) - check performance
3. **Large system** (1ERM, 4000 atoms) - stress test
4. **With ligand** (find PDB with HETATM) - verify filtering works

## References

- [MSMS algorithm paper](https://doi.org/10.1145/237170.237229)
- [EDTSurf paper](https://doi.org/10.1002/jcc.21297)
- [Marching cubes algorithm](https://en.wikipedia.org/wiki/Marching_cubes)
- [3Dmol.js surface code](https://github.com/3dmol/3Dmol.js/blob/master/3Dmol/surfaces.js)

## Next Steps

To implement this feature:
1. Start with frontend marching cubes (simpler, no backend changes)
2. Get basic working version with static surface
3. Add transparency and keyboard toggle
4. Optimize performance with Web Workers
5. If too slow, move to Python backend with MSMS

Would you like me to start implementing the basic version?
