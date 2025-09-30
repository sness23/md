# Live MD Demo — OpenMM → MDsrv → Mol* (Node backend)

This is a minimal, pragmatic starter kit for your idea: keep an MD simulation running and keep showing the growing trajectory in a web app.

## What’s inside

- **OpenMM runner** (`python/openmm_run.py`) — continuously appends to `./data/traj.dcd` and writes `./data/topology.pdb`.
- **MDsrv** — serves the files under `/data` so the browser can fetch structure/trajectory.
- **Node/Express** (`server.js`) — serves the web UI and exposes `/api/config` with the MDsrv URL.
- **Mol*** frontend (`public/index.html`) — loads PDB + DCD and plays the trajectory. A simple polling loop HEADs the DCD and reloads if the file grows.

> This demo takes the simplest approach (reload on file growth). For sub-second true live streaming and on-the-fly alignment, use MDsrv’s REST endpoints in a smarter polling loop, or build a small WebSocket bridge that feeds frames into Mol*’s custom trajectory API.

---

## Prereqs

- **Python 3.9+** with **OpenMM** and (optionally) **MDTraj**
- **MDsrv** (installable via `pip install mdsrv` or `conda install -c ngl mdsrv`)
- **Node 18+**

### Quick install

```bash
# 1) Python env (mamba recommended)
mamba create -n live-md -c conda-forge openmm mdtraj -y
conda activate live-md

# 2) MDsrv
pip install mdsrv  # or: conda install -c ngl mdsrv

# 3) Node deps
npm install
cp .env.example .env
```

Edit `.env` if needed:
```
PORT=5173
MDSRV_URL=http://127.0.0.1:8080
DATA_DIR=./data
```

---

## Run everything

Open **three terminals** in the project root:

**A) Start the MD producer (OpenMM):**
```bash
conda activate live-md
python python/openmm_run.py
```
This writes `data/topology.pdb` once and appends frames to `data/traj.dcd` forever(ish).

**B) Start MDsrv (serves /data):**
```bash
mdsrv --cfg scripts/app.cfg
```
This serves your `./data` directory at `http://127.0.0.1:8080/data/`.

**C) Start the Node web app:**
```bash
npm run dev
```
Open http://127.0.0.1:5173 — you should see the viewer. Press **Start Polling** to auto-reload when new frames arrive.

---

## Notes & tips

- **File rotation**: For very long runs, periodically rotate `traj.dcd` (e.g., `traj_00.dcd`, `traj_01.dcd`) and change the frontend to follow the latest. A trivial approach is to write a symlink `traj.dcd` pointing at the newest file.
- **Decimation**: We write one DCD frame every 100 steps. Tweak the reporter interval for your GPU/CPU and network.
- **Alignment/PBC**: This minimal demo does no on-the-fly alignment. MDsrv/NGL can handle centering/superposition; integrate MDsrv’s processing in production.
- **Security**: If you expose this beyond localhost, run MDsrv behind your web server (e.g., Apache mod_wsgi), and lock `DATA_DIRS` to a dedicated directory. Add auth if desired.
- **Switch to XTC**: DCD works out of the box with OpenMM. For XTC (smaller), use MDTraj or MDAnalysis reporters.
- **Custom streaming**: For near-real-time, send frames via WebSocket from Python and push into Mol*’s custom trajectory provider — this avoids reloading the entire DCD.

---

## Project layout

```
live-md-demo/
  data/                 # output: topology.pdb + traj.dcd (grows as sim runs)
  public/
    index.html          # Mol* viewer + simple polling logic
  python/
    openmm_run.py       # OpenMM producer (alanine dipeptide)
  scripts/
    app.cfg             # MDsrv config (serves ./data at http://127.0.0.1:8080/data/)
  .env.example
  package.json
  server.js             # Node/Express server and /api/config
  README.md
```

Have fun — and ping me if you want the WebSocket streaming version next.