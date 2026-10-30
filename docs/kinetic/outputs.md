# Kinetic outputs

What a kinetic run writes, and the input parameters that control it.

## Parthenon outputs (`<problem_id>.outN.*.phdf`, `.rhdf`)

These are set up by the `<parthenon/outputN>` blocks of the deck. Each file holds
the particle swarm variables named in `particles_variables`, and every parameter
of the `Deck` package.

### Field-plot grid (`Hermite_Field_Plot_*`)

The `Deck` package holds the interpolated fields sampled on a uniform (R, Z, t)
grid: B, J, dB/dR, dB/dZ, E, db/dt, plus psi, B² and |B|. The t samples span one
MHD step. It is a `Deck` parameter, so it is copied into every output file.
`PlotFieldsTime` refills it before each output in `hybrid`, `avalanche` and `mhd`.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Output/plot_fields` | 1 | 0: the grid is not created, filled or written. |
| `Output/NR_plot` | 10 | R samples |
| `Output/NZ_plot` | 20 | Z samples |
| `Output/Nt_plot` | 3 | t samples across one MHD step |

`Hermite_Field_Plot_data_F` has shape `NR_plot x NZ_plot x Nt_plot x 3 x 7`. At
the defaults it is 12,600 doubles, about 100 KB. Before 2026-10 the defaults were
400 x 800 x 10, which added 537 MB to every output file. In the reduced `hybrid`
test deck, a particle output is now 172 KB with the grid and 68 KB without it,
down from 537 MB. To get the old resolution back, set
`Output/NR_plot=400 Output/NZ_plot=800 Output/Nt_plot=10`.

The grid is diagnostics only. Two things used to depend on it and no longer do:

- The momentum boundary (`UpdateMomentumBoundary`, the max |E_par| search that
  sets `p_BC` and `p_RE`) samples its own points, `Numerical/NR_maxE` x
  `Numerical/NZ_maxE` (default 400 x 800).
- The psi and magnetic-axis refresh before each output (`RefreshFluxAndAxis`)
  runs whether or not the grid is written. The particle push reads that psi.

Results are identical with `plot_fields` 0 or 1, and for any plot size.

## Diagnostics file (`<problem_id>.diag.h5`)

One HDF5 file per run, written by MPI rank 0 (`src/kinetic/DiagnosticsFile.hpp`).
Every column is a 1-D float64 dataset with a `units` attribute.

| Group | Columns | Written | Notes |
|-------|---------|---------|-------|
| `/timeseries/conservation` | `t` [s], `p_phi_mean`, `mu_mean`, `weight_total` | after every large-angle sub-step, when `Simulation/EnableComputeConservedQuantities=1` | weight-averaged invariants of the alive particles |
| `/timeseries/currents` | `t` [s], `I_re_particles`, `I_re_grid`, `I_ohmic` | before every output, in `hybrid` | runaway current as a particle sum and as the deposited-grid integral; Ohmic current from the Hermite field |
| `/timeseries/push` | `t` [s], `push_s` [s], `launches`, `slots`, `ns_per_slot` [ns] | before every output, when `Simulation/push_timing=1` | push-kernel time since the previous output (slowest rank) |
| `/timeseries/feedback` | `t` [s], `a_max`, `R` [a], `Z` [a], `intervals` | before every output, in `hybrid` with `Hybrid/feedback_gain=1` | largest loop gain of the in-step runaway-current correction of E over the cells and deposition intervals since the previous output, and where it was (`src/kinetic/FeedbackGain.hpp`). The lagged correction is stable for a < 2 and rings for a > 0.34. |
| `/profiles/collision` | `p` [m_e c], `gamma`, `dtSA` [tau_c], `psi`, `CB`, `CF`, `CouLogee_ratio`, `probability` | once, at set-up | collision coefficients against momentum |
| `/snapshots/c<cycle>_cd<iCD>_la<iLA>` | `p`, `xi`, `R`, `phi`, `Z`, `weight`, `status` of every alive particle; attributes `cycle`, `iCD`, `iLA`, `t_in_step_tau_c`, `n_particles` | every `Output/substep_snapshots`-th call of `MakeOutputs` | in-step state, see below |
| `/meta` | attributes `problem_id`, `format` | at set-up | |

`t` is the simulation time in seconds. The legacy `conservation.dat` time column
was the time within the current MHD step in units of tau_c, so it restarted at
zero every step.

Rows are buffered and written to disk before every Parthenon output and at the
end of the run. A job that is killed loses at most one output interval of rows
and leaves a readable file.

On restart (`-r`) the file is reopened, and the rows of every series with
`t` after the restart time are dropped before the first step, so the series stay
continuous with no duplicates. With collisions off, a restarted run reproduces
the uninterrupted run's series exactly. With collisions on, the rows after the
restart differ, because the random-number state is not part of the restart file.

Read it with `scripts/hymera_plot` (`hymera_plot.io.find_diagnostics(run_dir)`) or
any HDF5 reader.

### In-step snapshots (`MakeOutputs`)

`Kinetic::MakeOutputs(pmesh, cycle, iCD, iLA, t_in_step)` is a task a driver can
put anywhere in its task list to record the particles in the middle of a time
step. It gathers the alive particles of all ranks onto rank 0 and writes them as
one `/snapshots/` group. It reads the particles and changes nothing; Parthenon's
own outputs and their numbering are untouched.

`Output/substep_snapshots=N` writes every Nth call (default 0: never). Only
`ConservationDriver` adds the task, after each large-angle sub-step, and only
when N > 0. To use it elsewhere, add the same `AddTask` line to that driver.
Each snapshot holds 7 doubles per particle, so pick N with the particle count
and the number of sub-steps (`nLA x nCD` per step) in mind.

`hymera_plot phase-space` also plots the first and last snapshot.

### Legacy text files

`Simulation/legacy_dat=1` also writes the files the code wrote before the
diagnostics file existed, with the same content as before: `conservation.dat`
(name: `Simulation/conservation_log`) and `collision_profiles.dat`. Default 0.

## Other files

| File | Written by | Control |
|------|------------|---------|
| `current.out` | truncated by `Initialize` | `Simulation/file_path` |
| `single_particle.dat` | `single_particle` | `SingleParticle/output_file` |
| `<basename>.<n>.rdat` | `WorkBeforeRestartOutput`: MHD state next to each Parthenon restart | `MHD_Config/file_extention` |
| `fields.h5` | `hybrid` with `Simulation/dump_raw_fields=1`, then exits | `Simulation/raw_fields_file` |
| `TotalCurrent.txt` | MHD monitor (`src/mhd`) | `MHD_Config/monitor` |

## Run log

Rank 0 writes tagged lines (`[info]`, `[debug]`, `[warn]`) to stdout.
`Simulation/verbosity` selects `quiet` (warnings only), `info` (default: set-up
milestones and once-per-run values) or `debug` (adds per-step and per-block
messages, and the current values written to `/timeseries/currents`).

## Plotting

```sh
tests/regression/setup_env.sh                     # numpy, matplotlib, h5py in ~/.venv-hymera
cd scripts && ~/.venv-hymera/bin/python -m hymera_plot all RUN_DIR
```

`all` writes every plot the run's files allow into `RUN_DIR/plots` (or `-o DIR`):
`conservation.png`, `currents.png`, `collision_profiles.png`, `phase_space.png`
and `fields.png` (needs `Output/plot_fields=1`). The single commands are
`timeseries`, `profiles`, `phase-space` and `fields`. Runs made before the
diagnostics file existed are read from their `.dat` files.

## Particle push: serial or parallel

`Simulation/parallel_push` (default 1) selects how the particle push loop runs.

| Value | Host (OpenMP) build | Reproducibility |
|-------|---------------------|-----------------|
| 1 | Kokkos parallel loop over `OMP_NUM_THREADS` threads | Collisions off: bit-identical to 0 at any thread count. Collisions on: the random draws each particle gets depend on the thread schedule, so runs agree statistically, not bit for bit. |
| 0 | one serial loop (Parthenon's `DEFAULT_LOOP_PATTERN` on host) | Bit-identical at any thread count, including with collisions. The regression decks use this. |

On CUDA both values give the same parallel kernel.

Measured on the 16,384-particle `conservation.input` deck run by `profile` (push
only), seconds per MHD step:

| Threads | `parallel_push=0` | `parallel_push=1` |
|---------|------|------|
| 1 | 31.6 | 29.5 |
| 4 | 31.6 | 8.1 |
| 6 | 31.6 | 5.8 |

The container used for these numbers has a 6-CPU quota (`cpu.max`), so more
threads than that get slower, not faster (16 threads: 103 s). Set
`OMP_NUM_THREADS` to the CPUs actually available.

With collisions on, the alive-particle mean momentum after 2 steps differed
between the two settings by 0.01 standard errors.

`scripts/bench_push.sh` reruns this benchmark (`--help` for the options). Bind
the MPI rank to as many cores as threads (`mpirun --bind-to none`, or
`srun --cpus-per-task=N --cpu-bind=cores`); a rank pinned to one core runs all
its threads on that core and shows no speed-up.

## Particle push: work array in registers or team scratch

`Simulation/push_scratch=1` (read only when present) keeps each particle's
Runge-Kutta work array (rk4: 5 states, dopri5: 10) in Parthenon team scratch
memory instead of a per-thread array in registers and local memory. The push
then runs as a team kernel of `Simulation/push_team_size` threads (default 128),
one particle slot per thread. The arithmetic is the same, so results are
bit-identical to the register push wherever the random draws come in the same
order: collisions off, or collisions on with one host thread. With collisions on
more than one host thread (or on a GPU) the draw order differs, as for
`parallel_push=1`.

