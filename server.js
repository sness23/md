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

// CORS headers for all requests
app.use((req, res, next) => {
  res.header('Access-Control-Allow-Origin', '*')
  res.header('Access-Control-Allow-Methods', 'GET, HEAD, OPTIONS')
  res.header('Access-Control-Allow-Headers', 'Content-Type')
  next()
})

// Serve data directory
app.use('/data', express.static(path.join(__dirname, DATA_DIR)))

// Static frontend
app.use(express.static(path.join(__dirname, 'public')))

// Simple config endpoint the frontend can read
app.get('/api/config', (req, res) => {
  res.json({
    mdsrvUrl: `http://127.0.0.1:${PORT}`,
    // We serve data files directly from Express instead of mdsrv
    topology: 'topology.pdb',
    trajectory: 'traj.dcd'
  })
})

app.listen(PORT, () => {
  console.log(`[live-md-demo] Web server listening at http://127.0.0.1:${PORT}`)
  console.log(`[live-md-demo] Expect MDsrv at ${MDSRV_URL}`)
})