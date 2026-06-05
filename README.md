<div align='center'>
	<h1>Parallel Particle Filter — 6-DOF LiDAR Localization</h1>
	<h3>Advanced HPC course, University of Trieste</h3>
	<h5>Author: Christian Faccio</h5>
</div>

A 6-DOF Monte-Carlo Localization (particle filter) for a UAV navigating with a
3D LiDAR, written in C. A simulated LiDAR scans a known map; the filter
estimates the full rigid-body pose `(x, y, z, roll, pitch, yaw)` of the
vehicle. This repository is the **serial correctness baseline**; the
SIMD / OpenMP / CUDA backends build on the hot kernel identified below.

## How it works

Each time step runs the classic predict → update → resample → estimate loop:

```
world (lidarsim scene)
  └─ build once: dense multi-view scan → reference point cloud
                 → 3D Euclidean distance field   (the "likelihood field" map)

per step t along a scripted ground-truth trajectory:
  live scan  = simulated LiDAR at the true pose (+ sensor noise)   [sensor frame]
  control    = noisy relative motion between consecutive GT poses

  pf_prediction     : compose the noisy motion onto every particle
  pf_update_weights : for each particle, for each scan point —
                        world = R(particle)·point + t(particle)
                        d     = distance_field_lookup(world)        ← O(1)
                        logw += -d² / 2σ²
  pf_resample       : multinomial (CDF scan)
  pf_estimate       : weighted-mean position + averaged orientation
  → accumulate ATE against ground truth
```

The `N_particles × M_points` loop in `pf_update_weights` is the **HPC hot
kernel** — arithmetic-intense and data-parallel. Runtime scales linearly in
both `N` and `M`, with accuracy independent of them.

## Layout

| Path | Role |
|------|------|
| `include/pose.h`, `src/pose.c` | 6-DOF pose; quaternion math (no gimbal lock) |
| `include/world.h`, `src/world.c` | wrapper over lidarsim: builds the scene + LiDAR, returns scans. **Only** TU that includes lidarsim |
| `include/lidar_map.h`, `src/lidar_map.c` | reference cloud + 3D distance field (separable Felzenszwalb–Huttenlocher EDT) |
| `include/particle_filter.h`, `src/particle_filter.c` | the filter (init / predict / update / resample / estimate) |
| `src/main.c` | world + map setup, trajectory, the time-step loop, ATE report |
| `third_party/lidarsim/` | vendored LiDAR simulator (see below) |

## Build & run

```sh
make            # build bin/pf_serial  (-O3 -march=native)
make run        # build then run with defaults
make asan       # AddressSanitizer + UBSan build (run leak checks on Linux)
make debug      # -O0 -g for gdb
make clean
```

Run with optional arguments `pf_serial [num_particles] [num_steps] [n_az] [n_el]`:

```sh
./bin/pf_serial                 # 3000 particles, 200 steps, 180×16 LiDAR
./bin/pf_serial 4000 200 360 32 # more particles and a denser scan
```

Everything is **simulated**, so no external dataset is needed and ground-truth
poses are known exactly. Output reports the mean Absolute Trajectory Error
(ATE) in translation and rotation; the baseline converges to ≈ 0.06 m / 0.15°.

## Visualization

Set `PF_DUMP=1` to write ASCII `.pcd` clouds for a snapshot step:

```sh
PF_DUMP=1 ./bin/pf_serial
# -> map.pcd  scan.pcd  particles.pcd  gt_traj.pcd  est_traj.pcd
```

Open them in any point-cloud viewer (CloudCompare, `pcl_viewer`, Open3D) — each
file gets its own colour. For a quick, dependency-free PNG preview:

```sh
python3 tools/render_scene.py scene.png   # clean top-down + side views (recommended)
python3 tools/render_pcd.py   preview.png # raw oblique 3D point dump
```

`render_scene.py` is the readable one: a bird's-eye floor plan and a side
elevation, with the room and feature spheres as outlines, one live scan in
cyan, and the estimated path (red) overlaying ground truth (green). The side
view shows the UAV's vertical bob. (`render_pcd.py` dumps every cloud in one
oblique 3D view — denser and harder to read.)

## Dependency: lidarsim

The LiDAR sensor and the ground-truth world are provided by **lidarsim**
(Ricardo Casimiro, MIT licence — <https://github.com/stm32f303ret6/lidarsim>),
vendored under `third_party/lidarsim/`. It is a header-with-definitions, so it
is compiled in exactly one translation unit (`src/world.c`); everything else
uses the opaque `World` / `ScanCloud` interface in `world.h`.

Two notes on using it here:
- **macOS:** lidarsim calls `qsort_r` with the glibc signature; `world.c`
  provides an `__APPLE__`-guarded shim so it also builds on macOS. On the
  project's Linux container the shim compiles out.
- **Geometry:** lidarsim's `Box`/`Plane` intersection orients a *transformed*
  primitive by the transpose of its rotation (correct only when axis-aligned).
  Because the scene is rotated into the sensor frame at every scan, the world
  is built only from `Mesh` (room shell — transforms vertices, then
  `ray_triangle`) and `Sphere` (rotation-invariant) primitives, both correct
  under arbitrary rotation.

## Status & next steps

- [x] Serial 6-DOF LiDAR MCL baseline — builds clean (`-Wall -Wextra`),
      converges, linear `N×M` scaling.
- [ ] SIMD vectorization of the `pf_update_weights` kernel.
- [ ] OpenMP over particles (needs a thread-safe, per-particle RNG to replace
      `rand()`).
- [ ] CUDA / OpenACC offload of the kernel and the distance-field build.
