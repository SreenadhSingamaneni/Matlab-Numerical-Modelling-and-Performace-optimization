# Matlab Numerical Modelling and Performance optimization 

Numerical simulation of the **1D heat equation** using **Forward Euler** and **classical 4th-order Runge-Kutta (RK4)**, validated against a closed-form analytical solution.

Includes stability analysis, temporal and spatial convergence benchmarking, and a vectorization performance benchmark.

---

## Problem

$$\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}, \quad x \in (0,\, L),\quad t \geq 0$$

**Boundary conditions (Dirichlet):**
$$u(0,t) = 0, \quad u(L,t) = 0$$

**Initial condition:**
$$u(x,0) = \sin\!\left(\tfrac{\pi x}{L}\right)$$

This admits the closed-form solution used for validation:
$$u(x,t) = \exp\!\left(-\alpha\!\left(\tfrac{\pi}{L}\right)^{\!2} t\right)\sin\!\left(\tfrac{\pi x}{L}\right)$$

---

## Method

**Spatial discretisation:** Centred finite differences on `Nx` interior points with spacing `dx = L/(Nx+1)`:

$$\frac{\partial^2 u}{\partial x^2}\bigg|_{x_i} \approx \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2}$$

This produces the ODE system `dU/dt = A·U`, where `A` is a sparse tridiagonal Laplacian scaled by `α/Δx²`.

**Time integration:** Forward Euler and classical RK4.

---

## Repository Structure

```
heat-conduction-simulation/
├── src/
│   ├── build_laplacian.m      # Sparse tridiagonal Laplacian (spdiags)
│   ├── exact_solution.m       # Closed-form analytical solution
│   ├── euler_solver.m         # Forward Euler integrator
│   ├── rk4_solver.m           # Classical RK4 integrator
│   ├── simulate_heat_1d.m     # Main simulation: Euler + RK4 vs exact
│   ├── stability_analysis.m   # dt sweep, blow-up detection, Euler vs RK4
│   ├── error_convergence.m    # Temporal + spatial convergence studies
│   ├── profiler_benchmark.m   # Loop vs vectorised RHS timing
│   └── utils_timestamp.m      # Filesystem-safe timestamp helper
├── results/                   # Auto-generated plots and logs (git-ignored)
├── run_all.m                  # Master runner — executes all modules
├── LICENSE
└── README.md
```

---

## Quick Start

Open MATLAB, `cd` to the repo root, and run:

```matlab
run('run_all.m')
```

Or run any module individually:

```matlab
addpath('src');
simulate_heat_1d;
```

All output is saved to `results/`.

---

## Numerical Schemes

### Forward Euler — 1st-order accurate

$$\mathbf{U}_{n+1} = \mathbf{U}_n + \Delta t\, A\mathbf{U}_n$$

**Stability condition** (von Neumann / CFL):

$$\Delta t \leq \frac{\Delta x^2}{2\alpha} \approx \frac{2}{\rho(A)}$$

Violating this causes exponential blow-up.

### Classical RK4 — 4th-order accurate

$$\mathbf{U}_{n+1} = \mathbf{U}_n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

> **Note on stability:** For the diffusion equation, RK4's stability region is only **~1.39× larger** than Euler's (`dt ≤ 2.785/ρ(A)` vs `dt ≤ 2.0/ρ(A)`). Both methods are **conditionally stable**. RK4's main advantage over Euler is **accuracy** (O(Δt⁴) vs O(Δt¹)), not stability.

---

## Analysis Modules

### `stability_analysis`

Sweeps `dt` from `0.05×` to `4×` the Euler CFL limit for **both** Euler and RK4:
- Detects blow-up via norm threshold
- Marks stable vs unstable runs on a log-log plot
- Highlights the CFL boundary

### `error_convergence`

Two independent studies:

| Study | Fixed | Swept | Expected rate |
|---|---|---|---|
| **Temporal** | Nx = 200 | dt (separate per method) | Euler: O(Δt¹), RK4: O(Δt⁴) |
| **Spatial**  | dt = 1e-5 | Nx ∈ {10,20,40,80,160} | Both: O(Δx²) |

Separate dt sweeps are used for Euler and RK4 in the temporal study so each method operates in its temporal-error-dominated regime. Convergence slopes are estimated by log-log linear regression.

### `profiler_benchmark`

Compares two implementations of the RHS `dU/dt`:
- **Loop-based**: explicit `for` loop over interior points
- **Vectorised**: sparse matrix-vector multiply `A * U`

Both produce numerically identical results (verified by assertion). Timing uses `timeit()` and the speedup is reported in `results/benchmark_*.txt`.

---

## Default Parameters

| Parameter | Value | Notes |
|---|---|---|
| `L` | 1 | Domain length |
| `alpha` | 1 | Thermal diffusivity |
| `Nx` | 50 | Interior spatial points |
| `T` | 0.1 | Final simulation time |
| `dt` | 1×10⁻⁴ | CFL limit for Nx=50 is ≈1.92×10⁻⁴; this is safely below |

Edit parameters directly in `simulate_heat_1d.m`.

---

## Output Files

| File | Content |
|---|---|
| `solution_comparison_*.png` | Euler + RK4 vs exact solution at time T |
| `stability_sweep_*.png` | Blow-up map: stable/unstable across dt values |
| `convergence_temporal_*.png` | Temporal convergence log-log plot |
| `convergence_spatial_*.png` | Spatial convergence log-log plot |
| `benchmark_*.txt` | Loop vs vectorised timing report |

---

## Requirements

- MATLAB R2018b or later
- No additional toolboxes required
- All functions used (`spdiags`, `timeit`, `datetime`, `xline`) are in MATLAB core

---
