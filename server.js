import express from 'express'
import path from 'path'
import { fileURLToPath } from 'url'
import dotenv from 'dotenv'

dotenv.config()

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const app = express()
const PORT = process.env.PORT || 5173
const MDSRV_URL = process.env.MDSRV_URL || 'http://127.0.0.1:8080'
const DATA_DIR = process.env.DATA_DIR || './data'

// Static frontend
app.use(express.static(path.join(__dirname, 'public')))

// Simple config endpoint the frontend can read
app.get('/api/config', (req, res) => {
  res.json({
    mdsrvUrl: MDSRV_URL,
    // We expect two files served by MDsrv under /data alias (see app.cfg): topology.pdb + traj.dcd
    topology: 'topology.pdb',
    trajectory: 'traj.dcd'
  })
})

app.listen(PORT, () => {
  console.log(`[live-md-demo] Web server listening at http://127.0.0.1:${PORT}`)
  console.log(`[live-md-demo] Expect MDsrv at ${MDSRV_URL}`)
})