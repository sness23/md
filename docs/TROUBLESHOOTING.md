# Troubleshooting Guide

## Quick Diagnostics

Run through these checks first:

```bash
# 1. Check if simulation is running
ps aux | grep openmm_run.py

# 2. Check if files are being created
ls -lh data/

# 3. Check if web server is running
curl http://localhost:5173/api/config

# 4. Check if files are accessible
curl -I http://localhost:5173/data/topology.pdb
curl -I http://localhost:5173/data/traj.dcd

# 5. Verify trajectory integrity
python verify_traj.py
```

## Installation Issues

### Cannot import openmm

**Symptoms:**
```
ModuleNotFoundError: No module named 'openmm'
```

**Solutions:**

1. **Verify conda environment:**
   ```bash
   conda activate live-md
   python -c "import openmm; print(openmm.version.version)"
   ```

2. **Reinstall OpenMM:**
   ```bash
   conda install -c conda-forge openmm -y
   ```

3. **Check Python version:**
   ```bash
   python --version  # Should be 3.9+
   ```

### Cannot import openmmtools

**Symptoms:**
```
ModuleNotFoundError: No module named 'openmmtools'
```

**Solution:**
```bash
conda activate live-md
conda install -c conda-forge openmmtools -y
# OR
pip install openmmtools
```

### npm install fails

**Symptoms:**
```
npm ERR! code EACCES
npm ERR! syscall access
```

**Solutions:**

1. **Fix permissions:**
   ```bash
   sudo chown -R $USER:$USER ~/.npm
   sudo chown -R $USER:$USER node_modules
   ```

2. **Clear cache:**
   ```bash
   npm cache clean --force
   rm -rf node_modules package-lock.json
   npm install
   ```

3. **Use different Node version:**
   ```bash
   nvm install 18
   nvm use 18
   npm install
   ```

## Simulation Issues

### Simulation crashes immediately

**Symptoms:**
```python
ValueError: No template found for residue...
```

**Solution:** The PDB file has incompatible residue names. See [README-openmm_run.md](../README-openmm_run.md) for force field troubleshooting.

**Quick fix:**
The default test system should always work:
```python
# Use OpenMM's built-in alanine dipeptide
from openmmtools import testsystems
test_system = testsystems.AlanineDipeptideImplicit()
```

### Simulation runs but no DCD file

**Symptoms:**
```bash
$ ls data/
topology.pdb
# traj.dcd is missing
```

**Causes:**
1. Simulation hasn't run 100 steps yet (first frame not written)
2. DATA_DIR environment variable is wrong
3. Permission issues

**Solutions:**

1. **Wait longer:**
   ```bash
   # Watch for file creation
   watch -n 1 ls -lh data/
   ```

2. **Check DATA_DIR:**
   ```bash
   # In openmm_run.py logs, look for:
   [openmm_run] Writing to /path/to/data/topology.pdb and /path/to/data/traj.dcd
   ```

3. **Fix permissions:**
   ```bash
   chmod 755 data
   ```

### Simulation very slow

**Symptoms:**
```
#"Step","Potential Energy (kJ/mole)","Temperature (K)","Speed (ns/day)"
1000,-187.2,295.3,0.05
```
Speed < 1 ns/day is very slow.

**Solutions:**

1. **Use GPU platform:**
   ```python
   platform = mm.Platform.getPlatformByName('CUDA')
   properties = {'CudaPrecision': 'mixed'}
   simulation = app.Simulation(topology, system, integrator, platform, properties)
   ```

   Expected speeds:
   - CPU: 0.1-1 ns/day (alanine dipeptide)
   - GPU (V100): 50-200 ns/day
   - GPU (A100): 100-500 ns/day

2. **Check CUDA availability:**
   ```python
   import openmm as mm
   for i in range(mm.Platform.getNumPlatforms()):
       platform = mm.Platform.getPlatform(i)
       print(f"{i}: {platform.getName()}")
   ```

   Should show: `CUDA`, `OpenCL`, `CPU`, `Reference`

3. **Verify GPU usage:**
   ```bash
   nvidia-smi -l 1
   # Should show Python process using GPU
   ```

### Out of memory

**Symptoms:**
```
RuntimeError: CUDA out of memory
```

**Solutions:**

1. **Use mixed precision:**
   ```python
   properties = {'CudaPrecision': 'mixed'}  # Instead of 'double'
   ```

2. **Reduce system size:**
   - Smaller molecule
   - Smaller solvent box
   - Fewer explicit water molecules

3. **Use CPU fallback:**
   ```python
   try:
       platform = mm.Platform.getPlatformByName('CUDA')
   except:
       platform = mm.Platform.getPlatformByName('CPU')
   ```

## Web Server Issues

### Server won't start

**Symptoms:**
```
Error: listen EADDRINUSE: address already in use :::5173
```

**Solution:**

1. **Find process using port:**
   ```bash
   # macOS/Linux
   lsof -i :5173

   # Or
   netstat -an | grep 5173
   ```

2. **Kill process:**
   ```bash
   kill -9 <PID>
   ```

3. **Use different port:**
   ```bash
   PORT=5174 npm start
   ```

### Cannot access server from other machines

**Symptoms:**
- Works on `http://localhost:5173`
- Fails on `http://192.168.1.5:5173`

**Solutions:**

1. **Check server is listening on 0.0.0.0:**
   ```javascript
   // server.js - should be:
   app.listen(PORT, '0.0.0.0', ...)

   // NOT:
   app.listen(PORT, 'localhost', ...)
   ```

2. **Check firewall:**
   ```bash
   # macOS
   sudo /usr/libexec/ApplicationFirewall/socketfilterfw --getglobalstate

   # Linux (ufw)
   sudo ufw status

   # Linux (firewalld)
   sudo firewall-cmd --list-all
   ```

3. **Allow port through firewall:**
   ```bash
   # macOS: System Preferences → Security & Privacy → Firewall

   # Linux (ufw)
   sudo ufw allow 5173/tcp

   # Linux (firewalld)
   sudo firewall-cmd --add-port=5173/tcp --permanent
   sudo firewall-cmd --reload
   ```

### 404 Not Found for data files

**Symptoms:**
```
GET http://localhost:5173/data/topology.pdb 404 (Not Found)
```

**Solutions:**

1. **Verify files exist:**
   ```bash
   ls -lh data/topology.pdb data/traj.dcd
   ```

2. **Check DATA_DIR in server.js:**
   ```javascript
   const DATA_DIR = process.env.DATA_DIR || './data'
   console.log('Serving data from:', DATA_DIR)
   ```

3. **Check file permissions:**
   ```bash
   chmod 644 data/topology.pdb data/traj.dcd
   chmod 755 data
   ```

4. **Restart server:**
   ```bash
   # Stop server (Ctrl+C)
   npm run dev
   ```

## Viewer Issues

### Blank screen

**Symptoms:**
- Page loads but canvas is blank
- No errors in console

**Solutions:**

1. **Check WebGL support:**
   ```
   Visit: https://get.webgl.org/
   Should show spinning cube
   ```

2. **Check browser console (F12):**
   ```
   Look for JavaScript errors
   ```

3. **Check Babylon.js loaded:**
   ```javascript
   // In console:
   typeof BABYLON
   // Should return "object", not "undefined"
   ```

4. **Try different browser:**
   - Chrome/Edge (best compatibility)
   - Firefox
   - Safari (may have issues)

### Error: Cannot read properties of undefined (reading 'data')

**Symptoms:**
```
TypeError: Cannot read properties of undefined (reading 'data')
```

**This was the old Mol* viewer error. Current Babylon.js viewer should not have this issue.**

If you see this:
1. Clear browser cache (Cmd/Ctrl + Shift + R)
2. Verify you're using the latest `index.html`
3. Check that Babylon.js is loading from CDN

### Atoms not visible

**Symptoms:**
- Scene loads
- No errors
- But no spheres visible

**Solutions:**

1. **Check atom count:**
   ```
   Open console, should see:
   "Loaded atoms: 22"
   "Atoms: 22" in sidebar
   ```

2. **Check coordinates are reasonable:**
   ```bash
   python verify_traj.py
   # Coordinates should be <100 Angstroms
   ```

3. **Reset camera:**
   ```javascript
   // In browser console:
   camera.setPosition(new BABYLON.Vector3(0, 0, 30))
   camera.setTarget(BABYLON.Vector3.Zero())
   ```

4. **Check sphere creation:**
   ```
   Console should show no errors in createAtoms()
   ```

### Animation not playing

**Symptoms:**
- Atoms are visible
- Frame counter shows "0/0" or "1/X" (stuck)
- Pause button doesn't change

**Solutions:**

1. **Check if frames loaded:**
   ```
   Console should show:
   "Loaded frames: XXXX"

   Sidebar should show:
   "Frame: 1/XXXX"
   ```

2. **Click Pause button:**
   - Should toggle between "Pause" and "Play"
   - Should change `playing` variable

3. **Check animation loop:**
   ```
   Console every few seconds should show:
   "Animation frame: XXX counter: XXX"
   ```

4. **Verify frames have different coordinates:**
   ```bash
   python verify_traj.py
   # Should show first and last frames differ
   ```

5. **Check requestAnimationFrame:**
   ```javascript
   // In console:
   let frames_checked = 0
   setInterval(() => {
       console.log('Current frame:', currentFrame)
       frames_checked++
       if (frames_checked > 5) console.log('Animation running!')
   }, 1000)
   ```

### Polling not detecting new frames

**Symptoms:**
- "Start Polling" clicked
- Button says "Stop Polling"
- Status stays "ready"
- New frames not loading

**Solutions:**

1. **Check poll interval:**
   - Should be >= 1 second
   - Try increasing to 10 seconds

2. **Verify simulation is running:**
   ```bash
   # Watch file size grow
   watch -n 1 ls -lh data/traj.dcd
   ```

3. **Check network tab (F12):**
   - Should see HEAD requests to `/data/traj.dcd` every N seconds
   - Should return 200 with Content-Length header

4. **Manual test:**
   ```javascript
   // In console:
   checkForNewFrames()
   // Should see "loading new frames..." if file grew
   ```

5. **Check for CORS errors:**
   ```
   Console should NOT show:
   "blocked by CORS policy"
   ```

### DCD parsing errors

**Symptoms:**
```
DCD parse error: RangeError: Offset is outside the bounds of the DataView
```

**Solutions:**

1. **Verify file integrity:**
   ```bash
   python verify_traj.py
   # Should complete without errors
   ```

2. **Check for partial writes:**
   ```bash
   # Simulation still running?
   ps aux | grep openmm_run

   # Wait a few seconds and try again
   ```

3. **Check file size:**
   ```bash
   ls -lh data/traj.dcd
   # Should be > 1 KB
   # Should be growing if simulation running
   ```

4. **Clear browser cache:**
   ```
   Cmd/Ctrl + Shift + R (hard reload)
   ```

5. **Test with known-good file:**
   ```bash
   # Make backup
   cp data/traj.dcd data/traj.dcd.bak

   # Stop simulation
   # Delete current file
   rm data/traj.dcd

   # Restart simulation
   # Wait for new file to be created
   ```

## Performance Issues

### High CPU usage

**Symptoms:**
```bash
top
# Shows >50% CPU for browser or Node
```

**Solutions:**

1. **Reduce animation framerate:**
   ```javascript
   // In index.html, change:
   const frameInterval = 50  // 20 fps
   // To:
   const frameInterval = 100  // 10 fps
   ```

2. **Reduce sphere segments:**
   ```javascript
   // Change:
   { diameter: radius * 2, segments: 12 }
   // To:
   { diameter: radius * 2, segments: 8 }
   ```

3. **Limit frames in memory:**
   ```javascript
   // Keep only last 1000 frames
   if (frames.length > 1000) {
       frames = frames.slice(-1000)
   }
   ```

### Browser tab freezes

**Symptoms:**
- Page unresponsive
- Browser shows "Page Unresponsive" dialog

**Causes:**
- Parsing very large DCD file (>100 MB)
- Too many frames loaded (>10k)

**Solutions:**

1. **Reduce output frequency:**
   ```python
   # In openmm_run.py, change:
   app.DCDReporter(TRAJ_PATH, 100)
   # To:
   app.DCDReporter(TRAJ_PATH, 500)  # Fewer frames
   ```

2. **Rotate trajectory files:**
   ```bash
   # Stop simulation
   # Archive old trajectory
   mv data/traj.dcd data/archive/traj_$(date +%Y%m%d).dcd

   # Restart simulation (creates new traj.dcd)
   ```

3. **Use Web Workers (future enhancement):**
   ```javascript
   // Move DCD parsing to worker thread
   // (not currently implemented)
   ```

### Slow file serving

**Symptoms:**
- Long wait for trajectory to load
- Network tab shows slow transfer

**Solutions:**

1. **Enable compression:**
   ```javascript
   // Add to server.js
   const compression = require('compression')
   app.use(compression())
   ```

2. **Use faster disk:**
   - Move `data/` to SSD
   - Use RAM disk for active trajectory

3. **Increase Node memory:**
   ```bash
   node --max-old-space-size=4096 server.js
   ```

## Network Issues

### CORS errors

**Symptoms:**
```
Access to fetch at 'http://...' has been blocked by CORS policy
```

**Solution:**
This shouldn't happen with current setup (Express adds CORS headers). If it does:

```javascript
// Verify in server.js:
app.use((req, res, next) => {
  res.header('Access-Control-Allow-Origin', '*')
  res.header('Access-Control-Allow-Methods', 'GET, HEAD, OPTIONS')
  res.header('Access-Control-Allow-Headers', 'Content-Type')
  next()
})
```

### Timeout errors

**Symptoms:**
```
TypeError: Failed to fetch
net::ERR_TIMED_OUT
```

**Solutions:**

1. **Check server is running:**
   ```bash
   curl http://localhost:5173/api/config
   ```

2. **Check network connectivity:**
   ```bash
   ping 192.168.1.5  # Your server IP
   ```

3. **Increase timeout:**
   ```javascript
   // Add to fetch calls:
   fetch(url, { signal: AbortSignal.timeout(30000) })  // 30 second timeout
   ```

## Data Issues

### Trajectory validation fails

**Symptoms:**
```bash
$ python verify_traj.py
✗ Error loading trajectory: ...
```

**Solutions:**

1. **Check file exists:**
   ```bash
   ls -lh data/traj.dcd
   ```

2. **Check file not corrupted:**
   ```bash
   file data/traj.dcd
   # Should say: "data"
   ```

3. **Try with MDTraj:**
   ```python
   import mdtraj as md
   try:
       traj = md.load('data/traj.dcd', top='data/topology.pdb')
       print(f"Loaded {traj.n_frames} frames")
   except Exception as e:
       print(f"Error: {e}")
   ```

4. **Regenerate files:**
   ```bash
   # Backup
   mv data data.old

   # Create fresh directory
   mkdir data

   # Restart simulation
   python python/openmm_run.py
   ```

### Coordinates look wrong

**Symptoms:**
- Atoms very far apart (>1000 Å)
- Atoms overlapping
- Strange geometries

**Cause:** Unit mismatch or parsing error

**Solutions:**

1. **Check units in PDB:**
   ```bash
   head -5 data/topology.pdb
   # Coordinates should be <100
   ```

2. **Check DCD units:**
   ```python
   # In verify_traj.py output:
   # Coordinates should be 0.1-10 nm (1-100 Å)
   ```

3. **Verify parser:**
   ```
   Console logs should show reasonable values:
   "Frame 0, atom 0: 1.78 1.35 0.64"
   NOT: "Frame 0, atom 0: 1780000 1350000 640000"
   ```

## Getting Help

### Collect diagnostic information

```bash
# System info
uname -a
python --version
node --version

# Check processes
ps aux | grep -E "(openmm|node)"

# Check files
ls -lh data/

# Check logs
tail -50 data/simulation.log  # If you have logging enabled

# Test trajectory
python verify_traj.py

# Test web server
curl -v http://localhost:5173/api/config
```

### Enable debug logging

**Server:**
```javascript
// Add to server.js
app.use((req, res, next) => {
    console.log(`[${new Date().toISOString()}] ${req.method} ${req.url}`)
    next()
})
```

**Browser:**
```javascript
// In index.html, add more console.log statements
console.log('State:', { playing, frames: frames.length, currentFrame })
```

### Report issues

When reporting issues, include:
1. Operating system and version
2. Python version and OpenMM version
3. Node.js version
4. Browser and version
5. Error messages (full text)
6. Console logs
7. Steps to reproduce

## Common Error Messages

| Error | Likely Cause | Solution |
|-------|-------------|----------|
| `No template found for residue` | Wrong force field | See README-openmm_run.md |
| `CUDA out of memory` | GPU memory full | Use mixed precision or CPU |
| `404 Not Found` | Files not served | Check DATA_DIR and permissions |
| `Failed to fetch` | Network/CORS issue | Check CORS headers |
| `Offset outside bounds` | DCD parse error | Verify file with verify_traj.py |
| `Address in use` | Port already used | Kill process or use different port |
| `Cannot read 'data'` | Old Mol* error | Update to Babylon.js version |

## Still Stuck?

1. **Check existing issues:** https://github.com/yourusername/md/issues
2. **Search documentation:** All docs in `docs/` folder
3. **Ask for help:** Create new issue with diagnostic info
4. **Email:** your.email@example.com