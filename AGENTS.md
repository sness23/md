# Repository Guidelines

## Project Structure & Module Organization
- `server.js` runs the Express service that serves `public/` and proxies trajectory assets; treat it as the single entry point for Node changes.
- `public/` hosts the viewer shell (`index.html`) while MD artifacts live in `data/` (`topology.pdb`, `traj.dcd`); keep generated files out of commits.
- Simulation helpers sit in `python/` (notably `openmm_run.py` for generating live data) and validation tools in `bin/` such as `verify_traj.py`.
- Reference material, deployment notes, and troubleshooting guides reside in `docs/`; sync any workflow changes with these files.

## Build, Test, and Development Commands
- Install deps once with `npm install`; Python tooling expects a conda env described in `python/openmm_run.py`.
- Use `npm run dev` for hot reload via nodemon, and `npm start` for production-like runs.
- Generate trajectories locally with `python python/openmm_run.py` (ensure `DATA_DIR` is writable); sanity-check outputs via `python bin/verify_traj.py`.
- Confirm the server wiring with `curl http://localhost:5173/api/config` before shipping changes.

## Coding Style & Naming Conventions
- Follow the existing ESM pattern: 2-space indentation, single quotes, and no semicolons in JavaScript.
- Prefer descriptive camelCase for variables and uppercase snake case for environment variables (`MDSRV_URL`, `DATA_DIR`).
- Keep static assets lowercase-with-dashes and colocate feature-specific files under `public/` to avoid path rewrites.

## Testing Guidelines
- No automated test suite yet; rely on manual validation by streaming OpenMM output and loading the viewer end-to-end.
- When touching simulation code, regenerate `data/` assets and compare trajectory stats with `verify_traj.py` to spot regressions.
- Document any bespoke test steps in PR descriptions so others can reproduce them quickly.

## Commit & Pull Request Guidelines
- Craft concise, imperative commit subjects (`Add nodemon hot reload`) and include focused changes per commit.
- PRs should link to relevant issues or docs, describe validation steps (`npm run dev`, trajectory replay), and include screenshots or GIFs if the UI changes.
- Flag configuration changes (env vars, ports) prominently and update `docs/DEPLOYMENT.md` when deployment expectations shift.

## Configuration Tips
- Load secrets from `.env`; default values ship in `server.js`, but avoid hardcoding credentials.
- If you point `DATA_DIR` elsewhere, update both the server and simulation scripts to keep paths aligned.
