#!/usr/bin/env python
"""
Verify trajectory file integrity and print stats
"""
import mdtraj as md
import sys

try:
    print("Loading trajectory...")
    traj = md.load('data/traj.dcd', top='data/topology.pdb')

    print(f"\n=== Trajectory Info ===")
    print(f"Frames: {traj.n_frames}")
    print(f"Atoms: {traj.n_atoms}")
    print(f"Time span: {traj.time[0]:.2f} - {traj.time[-1]:.2f} ps")
    print(f"Timestep: {traj.timestep:.2f} ps")

    print(f"\n=== First frame coordinates (first 3 atoms) ===")
    for i in range(min(3, traj.n_atoms)):
        x, y, z = traj.xyz[0, i, :]
        print(f"Atom {i}: ({x:.3f}, {y:.3f}, {z:.3f}) nm")

    print(f"\n=== Last frame coordinates (first 3 atoms) ===")
    for i in range(min(3, traj.n_atoms)):
        x, y, z = traj.xyz[-1, i, :]
        print(f"Atom {i}: ({x:.3f}, {y:.3f}, {z:.3f}) nm")

    print("\n✓ Trajectory is valid!")

except Exception as e:
    print(f"\n✗ Error loading trajectory: {e}", file=sys.stderr)
    sys.exit(1)