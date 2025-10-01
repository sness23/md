#!/usr/bin/env python3
"""
Check whether OpenMM can use your GPU (Metal/OpenCL/etc.) on this machine.
Works with both `openmm` (>=7.6) and legacy `simtk.openmm`.
"""

import sys
import platform as pyplat
import traceback
from time import perf_counter

# --- import OpenMM (new name first, then legacy) ---
try:
    import openmm as mm
    import openmm.unit as unit
    LIB_NAME = "openmm"
except Exception:
    try:
        from simtk import openmm as mm
        from simtk import unit
        LIB_NAME = "simtk.openmm"
    except Exception as e:
        print("❌ Could not import OpenMM. Install with: conda install -c conda-forge openmm")
        sys.exit(1)

from openmm import Platform  # works via either import path above

def list_platforms():
    names = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    return names

def try_platform(name, props=None):
    """Attempt to create a minimal Context on the given platform."""
    props = props or {}
    # Minimal system: 1 particle, no forces
    system = mm.System()
    system.addParticle(12.0)  # arbitrary mass in amu
    integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)

    try:
        plat = Platform.getPlatformByName(name)
        t0 = perf_counter()
        context = mm.Context(system, integrator, plat, props)
        t1 = perf_counter()
        which = context.getPlatform().getName()
        # Clean up
        del context
        del integrator
        return True, which, (t1 - t0), None
    except Exception as e:
        return False, None, None, e

def main():
    print("=== OpenMM GPU Capability Checker ===")
    print(f"- OpenMM Python package: {LIB_NAME}")
    print(f"- Python: {sys.version.split()[0]}")
    print(f"- OS: {pyplat.system()} {pyplat.release()} ({pyplat.machine()})")
    if pyplat.system() == "Darwin":
        print("  note: on macOS, CUDA is unavailable; GPU paths are typically Metal (plugin) or OpenCL.")

    platforms = list_platforms()
    if not platforms:
        print("❌ No OpenMM platforms found.")
        sys.exit(2)

    print(f"- Detected OpenMM platforms: {platforms}")

    # Platforms that may indicate GPU acceleration in various setups
    gpu_candidates = []
    for p in platforms:
        if p.lower() in {"metal", "opencl", "cuda", "hip"}:
            gpu_candidates.append(p)

    if not gpu_candidates:
        print("⚠️  No GPU-capable platforms detected (expected one of: Metal, OpenCL, CUDA, HIP).")
        print("    You can still run CPU via the 'Reference' or 'CPU' platform if present.")
        sys.exit(0)

    print(f"- Trying GPU platforms: {gpu_candidates}")

    any_ok = False
    for name in gpu_candidates:
        # Try with sensible defaults; for OpenCL also try device index 0 explicitly.
        attempts = [({}, "default props")]
        if name.lower() == "opencl":
            attempts.append( ({"OpenCLPlatformIndex":"0","OpenCLDeviceIndex":"0"}, "OpenCLDeviceIndex=0") )
        if name.lower() == "cuda":
            attempts.append( ({"CudaPrecision":"mixed"}, "CudaPrecision=mixed") )

        for props, label in attempts:
            ok, which, dt, err = try_platform(name, props)
            if ok:
                any_ok = True
                print(f"✅ {name} [{label}] : created Context on '{which}' in {dt*1000:.1f} ms")
                break  # success on this platform; move on
            else:
                # If first attempt fails, show a short reason; keep it readable
                msg = str(err).strip().split("\n")[0]
                print(f"❌ {name} [{label}] failed: {msg}")

    if any_ok:
        print("\n🎉 Result: OpenMM can run on your GPU on this machine.")
    else:
        print("\n🙁 Result: No GPU Context could be created.")
        print("   Tips:")
        print("   - On Apple Silicon: install the OpenMM Metal plugin or try the OpenCL platform.")
        print("   - On Intel macOS: ensure OpenCL drivers are present; try 'OpenCLDeviceIndex=0'.")
        print("   - No NVIDIA CUDA on macOS; use Linux/Windows with an NVIDIA GPU for CUDA.")

if __name__ == "__main__":
    main()

