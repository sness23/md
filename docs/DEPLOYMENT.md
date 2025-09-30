# Deployment Guide

## Overview

This guide covers deploying the live MD visualization system in various environments, from local development to cloud production setups.

## Deployment Architectures

### 1. Local Development (Current Setup)

**Components:**
- OpenMM simulation (Python)
- Express server (Node.js)
- Shared filesystem (`./data/`)

**Pros:**
- Simple setup
- No network configuration
- Direct file access

**Cons:**
- Single machine only
- No remote access
- Limited scalability

### 2. LAN Deployment

Share the visualization across local network:

```
┌─────────────┐
│ Workstation │ ← Running OpenMM + Express
│ 192.168.1.5 │
└──────┬──────┘
       │ LAN
       ↓
┌─────────────┐
│ Lab Computer│ ← Browser viewer
│ 192.168.1.X │
└─────────────┘
```

**Setup:**

1. **Configure Express to listen on all interfaces:**
   ```javascript
   // server.js
   app.listen(PORT, '0.0.0.0', () => {
       console.log(`Server listening on http://0.0.0.0:${PORT}`)
   })
   ```

2. **Open firewall:**
   ```bash
   # macOS
   # System Preferences → Security & Privacy → Firewall → Firewall Options
   # Allow incoming connections for Node

   # Linux (ufw)
   sudo ufw allow 5173/tcp

   # Linux (firewalld)
   sudo firewall-cmd --add-port=5173/tcp --permanent
   sudo firewall-cmd --reload
   ```

3. **Access from other machines:**
   ```
   http://192.168.1.5:5173
   ```

### 3. Compute Cluster Deployment

Run simulation on HPC cluster, serve visualization:

```
┌──────────────┐
│ Login Node   │ ← Express server
│ frontend.hpc │
└──────┬───────┘
       │ shared filesystem (/scratch/user/md/)
       ↓
┌──────────────┐
│ Compute Node │ ← OpenMM simulation
│ node0042     │
└──────────────┘
```

**Setup:**

1. **Run simulation on compute node:**
   ```bash
   # Submit SLURM job
   sbatch run_simulation.sh
   ```

   ```bash
   #!/bin/bash
   #SBATCH --job-name=md_sim
   #SBATCH --nodes=1
   #SBATCH --gres=gpu:1
   #SBATCH --time=48:00:00

   conda activate live-md
   cd /scratch/$USER/md
   python python/openmm_run.py
   ```

2. **Run Express on login node:**
   ```bash
   # On login node
   cd /scratch/$USER/md
   npm start
   ```

3. **Set up SSH tunnel from local machine:**
   ```bash
   ssh -L 5173:localhost:5173 user@frontend.hpc
   ```

4. **Access locally:**
   ```
   http://localhost:5173
   ```

### 4. Cloud Deployment (AWS/GCP/Azure)

Full cloud-hosted solution:

```
┌────────────────┐
│  Web Browser   │
└───────┬────────┘
        │ HTTPS
        ↓
┌────────────────┐
│ Load Balancer  │
│  (SSL/TLS)     │
└───────┬────────┘
        │
        ↓
┌────────────────┐
│  Web Server    │ ← Express + static files
│  (EC2/GCE)     │
└───────┬────────┘
        │ shared storage
        ↓
┌────────────────┐
│  EFS/GCS/Blob  │ ← data/topology.pdb, data/traj.dcd
└───────▲────────┘
        │
        │
┌───────┴────────┐
│  GPU Instance  │ ← OpenMM simulation
│  (p3/A100)     │
└────────────────┘
```

**AWS Example:**

1. **Launch GPU instance (p3.2xlarge):**
   ```bash
   aws ec2 run-instances \
       --image-id ami-0c55b159cbfafe1f0 \
       --instance-type p3.2xlarge \
       --key-name my-key \
       --security-group-ids sg-xxxxxx \
       --subnet-id subnet-xxxxxx
   ```

2. **Mount EFS:**
   ```bash
   sudo mount -t efs fs-xxxxxxxx:/ /mnt/md-data
   ```

3. **Run simulation:**
   ```bash
   conda activate live-md
   cd /mnt/md-data
   python openmm_run.py
   ```

4. **Launch web server (t3.medium):**
   ```bash
   # Mount same EFS
   sudo mount -t efs fs-xxxxxxxx:/ /mnt/md-data

   # Start server
   cd /mnt/md-data
   npm start
   ```

5. **Configure security groups:**
   - GPU instance: No inbound (private)
   - Web server: Port 5173 from load balancer only
   - Load balancer: Port 443 from 0.0.0.0/0

## Production Considerations

### Security

#### 1. HTTPS/TLS

Use reverse proxy (nginx) for HTTPS:

```nginx
server {
    listen 443 ssl http2;
    server_name md.example.com;

    ssl_certificate /etc/letsencrypt/live/md.example.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/md.example.com/privkey.pem;

    location / {
        proxy_pass http://localhost:5173;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }

    location /data/ {
        proxy_pass http://localhost:5173/data/;
        proxy_set_header Host $host;

        # CORS headers
        add_header 'Access-Control-Allow-Origin' '*';
        add_header 'Access-Control-Allow-Methods' 'GET, HEAD, OPTIONS';
    }
}
```

#### 2. Authentication

Add basic auth with nginx:

```nginx
location / {
    auth_basic "MD Viewer";
    auth_basic_user_file /etc/nginx/.htpasswd;

    proxy_pass http://localhost:5173;
}
```

Generate password file:
```bash
htpasswd -c /etc/nginx/.htpasswd username
```

Or implement JWT in Express:

```javascript
const jwt = require('jsonwebtoken')

app.use('/data', (req, res, next) => {
    const token = req.headers.authorization?.split(' ')[1]
    if (!token) return res.status(401).send('Unauthorized')

    try {
        jwt.verify(token, process.env.JWT_SECRET)
        next()
    } catch {
        res.status(403).send('Forbidden')
    }
})
```

#### 3. Rate Limiting

```javascript
const rateLimit = require('express-rate-limit')

const limiter = rateLimit({
    windowMs: 15 * 60 * 1000, // 15 minutes
    max: 100 // limit each IP to 100 requests per windowMs
})

app.use('/data/', limiter)
```

### Performance

#### 1. Caching

**Client-side:**
```javascript
// Add Cache-Control headers
app.use('/data', (req, res, next) => {
    if (req.path.endsWith('.pdb')) {
        res.setHeader('Cache-Control', 'public, max-age=3600')
    } else if (req.path.endsWith('.dcd')) {
        res.setHeader('Cache-Control', 'no-cache')
    }
    next()
})
```

**CDN:** Use CloudFlare/CloudFront for static assets:
```javascript
// Serve Babylon.js from CDN (already done)
<script src="https://cdn.babylonjs.com/babylon.js"></script>
```

#### 2. Compression

```javascript
const compression = require('compression')
app.use(compression())
```

#### 3. Load Balancing

Use PM2 for clustering:

```bash
npm install -g pm2

# Start with cluster mode
pm2 start server.js -i max

# Save configuration
pm2 save

# Auto-restart on reboot
pm2 startup
```

Or use nginx upstream:

```nginx
upstream md_backend {
    server 127.0.0.1:5173;
    server 127.0.0.1:5174;
    server 127.0.0.1:5175;
}

server {
    location / {
        proxy_pass http://md_backend;
    }
}
```

### Monitoring

#### 1. Application Monitoring

```javascript
// Add /health endpoint
app.get('/health', (req, res) => {
    res.json({
        status: 'ok',
        uptime: process.uptime(),
        timestamp: Date.now()
    })
})

// Add /metrics endpoint
app.get('/metrics', (req, res) => {
    res.json({
        memory: process.memoryUsage(),
        cpu: process.cpuUsage()
    })
})
```

#### 2. Log Management

```javascript
const winston = require('winston')

const logger = winston.createLogger({
    level: 'info',
    format: winston.format.json(),
    transports: [
        new winston.transports.File({ filename: 'error.log', level: 'error' }),
        new winston.transports.File({ filename: 'combined.log' })
    ]
})

app.use((req, res, next) => {
    logger.info(`${req.method} ${req.url}`)
    next()
})
```

#### 3. Uptime Monitoring

Use external services:
- UptimeRobot (free tier)
- Pingdom
- StatusCake

Or self-hosted:
```bash
# Simple cron job
*/5 * * * * curl -fsS -m 10 --retry 3 http://md.example.com/health > /dev/null || echo "Server down!" | mail -s "MD Viewer Alert" admin@example.com
```

### Backup & Recovery

#### 1. Trajectory Backup

```bash
#!/bin/bash
# backup_trajectory.sh

DATE=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="/backup/md-trajectories"

# Create backup
tar -czf "$BACKUP_DIR/traj_$DATE.tar.gz" data/

# Keep only last 7 days
find "$BACKUP_DIR" -name "traj_*.tar.gz" -mtime +7 -delete

# Upload to S3 (optional)
aws s3 cp "$BACKUP_DIR/traj_$DATE.tar.gz" s3://my-bucket/backups/
```

Run daily:
```bash
crontab -e
0 2 * * * /path/to/backup_trajectory.sh
```

#### 2. Database Backup (if using)

```bash
# MongoDB example
mongodump --out /backup/mongodb_$(date +%Y%m%d)

# PostgreSQL example
pg_dump mydb > /backup/mydb_$(date +%Y%m%d).sql
```

## Environment-Specific Configurations

### Development (.env.development)

```bash
NODE_ENV=development
PORT=5173
DATA_DIR=./data
LOG_LEVEL=debug
```

### Production (.env.production)

```bash
NODE_ENV=production
PORT=8080
DATA_DIR=/mnt/efs/md-data
LOG_LEVEL=info
JWT_SECRET=your-secret-key
ALLOWED_ORIGINS=https://md.example.com
```

Load with:
```javascript
require('dotenv').config({ path: `.env.${process.env.NODE_ENV}` })
```

## Containerization (Docker)

### Dockerfile

```dockerfile
FROM node:18-slim

WORKDIR /app

# Install dependencies
COPY package*.json ./
RUN npm ci --only=production

# Copy application
COPY server.js ./
COPY public ./public

# Expose port
EXPOSE 5173

# Run
CMD ["node", "server.js"]
```

### Docker Compose

```yaml
version: '3.8'

services:
  web:
    build: .
    ports:
      - "5173:5173"
    volumes:
      - ./data:/app/data:ro
    environment:
      - NODE_ENV=production
    restart: unless-stopped

  simulation:
    image: conda/miniconda3
    volumes:
      - ./data:/data
      - ./python:/app
    command: >
      bash -c "
      conda create -n md -c conda-forge openmm openmmtools -y &&
      conda run -n md python /app/openmm_run.py
      "
    restart: unless-stopped
```

Run:
```bash
docker-compose up -d
```

## Scaling Strategies

### Horizontal Scaling

**Multiple simulations:**
```
Simulation 1 → /data/sim1/
Simulation 2 → /data/sim2/
Simulation 3 → /data/sim3/
```

**Load balancer routing:**
```nginx
location /sim1/ {
    proxy_pass http://backend1:5173/;
}

location /sim2/ {
    proxy_pass http://backend2:5173/;
}
```

### Vertical Scaling

**GPU upgrades:**
- V100 → A100 (2.5x faster)
- Multi-GPU: Run multiple simulations in parallel

**Memory:**
- Larger systems require more RAM
- Rule of thumb: 1 GB per 100k atoms

## Troubleshooting Production Issues

### High Memory Usage

```bash
# Monitor memory
free -h
ps aux --sort=-%mem | head

# Limit Node.js heap
node --max-old-space-size=4096 server.js
```

### File Descriptor Limits

```bash
# Check limits
ulimit -n

# Increase
ulimit -n 65536

# Permanent (add to /etc/security/limits.conf)
* soft nofile 65536
* hard nofile 65536
```

### Slow File Serving

```bash
# Enable nginx sendfile
sendfile on;
tcp_nopush on;
tcp_nodelay on;

# Or use Express static options
app.use('/data', express.static('data', {
    maxAge: '1h',
    etag: true
}))
```

## Cost Estimation (AWS Example)

**Small deployment:**
- Web server: t3.medium ($30/month)
- GPU instance: g4dn.xlarge ($390/month)
- EFS: 100 GB ($30/month)
- Data transfer: 500 GB ($45/month)
- **Total: ~$495/month**

**Large deployment:**
- Web servers: 3× t3.large ($180/month)
- Load balancer: ALB ($23/month)
- GPU instances: 2× p3.2xlarge ($1,836/month)
- EFS: 1 TB ($300/month)
- Data transfer: 2 TB ($180/month)
- **Total: ~$2,519/month**

## Maintenance

### Regular Tasks

**Daily:**
- Check disk space
- Monitor logs for errors
- Verify simulations running

**Weekly:**
- Review performance metrics
- Update dependencies (dev)
- Backup trajectory data

**Monthly:**
- Security updates
- Cost analysis
- Archive old trajectories

### Update Procedure

```bash
# 1. Backup current version
cp -r /app /app.backup

# 2. Pull latest code
git pull origin main

# 3. Install dependencies
npm ci

# 4. Test in development
NODE_ENV=development npm start

# 5. Restart production (zero-downtime with PM2)
pm2 reload server

# 6. Monitor logs
pm2 logs

# 7. Rollback if needed
pm2 stop server
cp -r /app.backup/* /app/
pm2 start server
```

## References

- [PM2 Documentation](https://pm2.keymetrics.io/)
- [Nginx Configuration](https://nginx.org/en/docs/)
- [AWS EFS](https://docs.aws.amazon.com/efs/)
- [Docker Compose](https://docs.docker.com/compose/)
- [Let's Encrypt SSL](https://letsencrypt.org/)