# Live Molecular Dynamics Visualization

Real-time streaming visualization of molecular dynamics simulations using OpenMM and Babylon.js.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![Node](https://img.shields.io/badge/node-18+-green.svg)

## Overview

This project enables live visualization of running MD simulations in a web browser. As OpenMM generates trajectory frames, they're automatically detected and rendered in real-time using a custom Babylon.js-based CPK sphere viewer.

**Key Features:**
- ✨ Real-time trajectory streaming
- 🎨 CPK color-coded atomic spheres
- 🔄 Automatic frame detection and loading
- 🖱️ Interactive 3D camera controls
- 📊 Simple architecture with minimal dependencies
- 🚀 No external visualization servers required

## Quick Start

### Prerequisites

- **Python 3.9+** with conda/mamba
- **Node.js 18+**
- **OpenMM** (via conda-forge)
- **MDTraj** (for trajectory verification)

### Installation

```bash
# 1. Clone repository
git clone <your-repo-url>
cd md

# 2. Create Python environment
mamba create -n live-md -c conda-forge openmm mdtraj openmmtools -y
conda activate live-md

# 3. Install Node dependencies
npm install

# 4. Configure environment (optional)
cp .env.example .env
```

### Running

**Terminal 1 - Start the simulation:**
```bash
conda activate live-md
python python/openmm_run.py
```

**Terminal 2 - Start the web server:**
```bash
npm run dev
```

**Open browser:**
- Navigate to http://127.0.0.1:5173
- Click "Pause" to stop animation or "Start Polling" to auto-detect new frames
- Use mouse to rotate (left), pan (right), zoom (wheel)

## Architecture

### Components

```
┌─────────────────┐
│  OpenMM Runner  │  Generates trajectory frames
│ (Python/OpenMM) │  → data/topology.pdb (once)
└────────┬────────┘  → data/traj.dcd (continuous)
         │
         │ writes to filesystem
         ↓
┌─────────────────┐
│   data/ dir     │  Shared filesystem
│  topology.pdb   │
│  traj.dcd       │
└────────┬────────┘
         │
         │ served via HTTP
         ↓
┌─────────────────┐
│ Express Server  │  Serves files + frontend
│   (Node.js)     │  → /data/* (static files)
└────────┬────────┘  → /api/config (metadata)
         │
         │ HTTP
         ↓
┌─────────────────┐
│  Web Frontend   │  Babylon.js CPK viewer
│ (Babylon.js)    │  - Parses PDB/DCD in browser
└─────────────────┘  - Polls for new frames
```

### Data Flow

1. **OpenMM** runs molecular dynamics simulation
2. **DCDReporter** appends frames to `data/traj.dcd` every 100 steps
3. **Express server** serves files from `data/` directory with CORS
4. **Frontend** polls trajectory file size via HEAD requests
5. When file grows, frontend fetches and parses new DCD data
6. **Babylon.js** updates atom positions and renders frames

### File Formats

- **PDB** (Protein Data Bank): Text format containing topology, atom names, and initial coordinates
- **DCD** (CHARMM/NAMD trajectory): Binary format with 3D coordinates for each frame

## Project Structure

```
md/
├── data/                    # Generated at runtime
│   ├── topology.pdb        # Structure file (22 atoms)
│   └── traj.dcd            # Trajectory frames (growing)
├── public/
│   └── index.html          # Babylon.js viewer (PDB/DCD parsers, CPK renderer)
├── python/
│   └── openmm_run.py       # MD simulation runner
├── server.js               # Express server (file serving + CORS)
├── verify_traj.py          # Trajectory validation tool
├── package.json            # Node dependencies
├── .env.example            # Configuration template
├── CLAUDE.md               # AI agent instructions
├── README.md               # This file
└── README-openmm_run.md    # OpenMM runner documentation
```

## Configuration

### Environment Variables

Create `.env` file (or use defaults):

```bash
# Server configuration
PORT=5173                           # Web server port
DATA_DIR=./data                     # Output directory
```

### OpenMM Parameters

Edit `python/openmm_run.py`:

```python
temperature = 300 * unit.kelvin     # Simulation temperature
friction = 1.0 / unit.picosecond    # Langevin friction
timestep = 2.0 * unit.femtoseconds  # Integration timestep
```

DCD output frequency:
```python
simulation.reporters.append(app.DCDReporter(TRAJ_PATH, 100))  # Every 100 steps
```

Simulation batch size:
```python
simulation.step(5000)  # 5000 steps per batch (~10 ps)
```

## Usage

### Basic Workflow

1. **Start simulation** - Begins generating frames
2. **Start web server** - Serves viewer and data files
3. **Open browser** - Loads initial topology and trajectory
4. **Auto-polling** - Detects new frames every 5 seconds (configurable)

### Controls

- **Play/Pause** - Toggle trajectory animation
- **Poll interval** - Adjust polling frequency (seconds)
- **Start/Stop Polling** - Enable/disable auto-detection
- **Mouse controls:**
  - Left drag: Rotate camera
  - Right drag: Pan camera
  - Scroll: Zoom in/out

### Using Custom Structures

See [README-openmm_run.md](./README-openmm_run.md) for detailed instructions on:
- Loading custom PDB files
- Choosing force fields
- Setting up explicit solvent
- Troubleshooting structure issues

## Verification

Check trajectory integrity:

```bash
python verify_traj.py
```

Output:
```
Loading trajectory...

=== Trajectory Info ===
Frames: 6165
Atoms: 22
Time span: 0.00 - 6164.00 ps
Timestep: 1.00 ps

=== First frame coordinates (first 3 atoms) ===
Atom 0: (0.178, 0.135, 0.064) nm
...

✓ Trajectory is valid!
```

## Performance

### Current Performance (CPU)

- **System:** Alanine dipeptide (22 atoms, implicit solvent)
- **Rate:** ~100-200 steps/second on modern CPU
- **Frame generation:** ~0.5-1.0 fps (100 steps per frame)
- **File size:** ~30 KB per 100 frames

### Optimization Tips

1. **Use GPU acceleration:**
   ```python
   platform = mm.Platform.getPlatformByName('CUDA')
   properties = {'CudaPrecision': 'mixed'}
   simulation = app.Simulation(topology, system, integrator, platform, properties)
   ```
   Expected speedup: 10-100x

2. **Reduce frame output frequency:**
   ```python
   app.DCDReporter(TRAJ_PATH, 500)  # Every 500 steps
   ```

3. **Increase batch size:**
   ```python
   simulation.step(50000)  # Larger batches
   ```

## Browser Compatibility

- ✅ Chrome/Edge 90+
- ✅ Firefox 88+
- ✅ Safari 14+
- ⚠️ Mobile browsers (limited performance)

Requires:
- WebGL 2.0
- ES6 modules
- Fetch API
- DataView/ArrayBuffer

## Known Limitations

1. **Full trajectory reload** - Currently reloads entire DCD file when detecting growth (simple but inefficient)
2. **No frame interpolation** - Displays discrete frames without smoothing
3. **Limited to ~10k frames** - Browser memory constraints
4. **CPU-only by default** - Simulation uses CPU platform for compatibility
5. **Single trajectory** - No support for switching between trajectories

## Troubleshooting

### Simulation won't start

```bash
# Check OpenMM installation
python -c "import openmm; print(openmm.version.version)"

# Check for openmmtools
pip install openmmtools
```

### No frames appearing

```bash
# Check if DCD file exists and is growing
ls -lh data/traj.dcd

# Verify with MDTraj
python verify_traj.py
```

### Viewer shows blank screen

- Open browser console (F12) and check for errors
- Verify files are accessible: http://127.0.0.1:5173/data/topology.pdb
- Check CORS headers in Network tab

### Animation not playing

- Click "Pause" button (should change to "Play")
- Check console for frame parsing errors
- Verify frames have different coordinates: `python verify_traj.py`

## Development

### Adding Features

**New visualization modes:**
- Edit `public/index.html` sphere creation in `createAtoms()`
- Modify CPK_COLORS and CPK_RADII dictionaries

**Different parsers:**
- Add parser functions alongside `parsePDB()` and `parseDCD()`
- Update `loadAll()` to use new parsers

**Server-side processing:**
- Add routes in `server.js`
- Create API endpoints for frame data

### Testing

```bash
# Test trajectory integrity
python verify_traj.py

# Test web server
curl http://127.0.0.1:5173/api/config
curl -I http://127.0.0.1:5173/data/topology.pdb
```

## References

- [OpenMM Documentation](http://docs.openmm.org/)
- [Babylon.js Documentation](https://doc.babylonjs.com/)
- [CHARMM DCD Specification](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html)

## License

MIT License - see LICENSE file for details

## Acknowledgments

- Built with [OpenMM](https://openmm.org/)
- Visualization powered by [Babylon.js](https://www.babylonjs.com/)
- Test system from [OpenMMTools](https://openmmtools.readthedocs.io/)