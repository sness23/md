# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a live molecular dynamics (MD) visualization demo that streams running MD simulations to a web browser in near-real-time. The architecture consists of three independent components that work together:

1. **OpenMM Producer** (`python/openmm_run.py`) - Continuously runs an MD simulation and appends frames to `data/traj.dcd`
2. **MDsrv** - Serves the `data/` directory files over HTTP at port 8080
3. **Node/Express Web Server** (`server.js`) - Serves the Mol* viewer frontend and provides configuration via `/api/config`

The frontend (`public/index.html`) uses Mol* viewer to visualize the trajectory and polls for file size changes to detect and reload new frames.

## Development Setup

### Prerequisites
- Python 3.9+ with OpenMM and MDTraj (via conda/mamba)
- MDsrv (`pip install mdsrv` or `conda install -c ngl mdsrv`)
- Node 18+

### Installation
```bash
# 1. Create Python environment
mamba create -n live-md -c conda-forge openmm mdtraj -y
conda activate live-md

# 2. Install MDsrv
pip install mdsrv

# 3. Install Node dependencies
npm install

# 4. Configure environment
cp .env.example .env
```

### Running the System

The system requires three processes running simultaneously:

**Terminal A - OpenMM simulation:**
```bash
conda activate live-md
python python/openmm_run.py
```

**Terminal B - MDsrv file server:**
```bash
mdsrv --cfg scripts/app.cfg
```

**Terminal C - Web application:**
```bash
npm run dev
# Or: npm start
```

Then open http://127.0.0.1:5173 in a browser and click "Start Polling" to watch the live trajectory.

## Architecture Details

### Data Flow
1. OpenMM writes topology to `data/topology.pdb` (once) and continuously appends frames to `data/traj.dcd` (every 100 simulation steps)
2. MDsrv serves files from `./data` at `http://127.0.0.1:8080/data/`
3. Express server serves the frontend and provides MDsrv URL via `/api/config`
4. Frontend polls trajectory file size via HEAD requests and reloads when file grows

### Configuration
- `.env` - Controls ports and URLs (see `.env.example`)
- `scripts/app.cfg` - MDsrv configuration (host, port, data directories)
- OpenMM simulation parameters are hardcoded in `python/openmm_run.py`

### Key Implementation Notes
- The demo uses simple file growth detection (HEAD request for Content-Length) rather than true streaming
- DCD trajectory file is continuously appended by OpenMM's DCDReporter
- Frontend reloads entire trajectory when growth is detected (naive but simple approach)
- For production use, consider MDsrv's REST endpoints or WebSocket streaming for frame-accurate updates
- Current simulation is alanine dipeptide in implicit solvent (CPU-based)

## File Structure
```
data/                    # Output directory (created at runtime)
  topology.pdb           # Static structure file
  traj.dcd               # Growing trajectory file
public/
  index.html             # Mol* viewer + polling logic
python/
  openmm_run.py          # MD simulation runner
scripts/
  app.cfg                # MDsrv configuration
server.js                # Express server
.env                     # Local configuration (not tracked)
```