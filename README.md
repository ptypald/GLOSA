# GLOSA Research Code

This repository collects research code and supporting notebooks for Green Light Optimal Speed Advisory (GLOSA) methods developed around two publications:

1. **Vehicle Trajectory Specification in Presence of Traffic Lights with Known or Uncertain Switching Times**  
   Panagiotis Typaldos, Ioanna Kalogianni, Kyriakos Simon Mountakis, Ioannis Papamichail, Markos Papageorgiou. TRB 2020 manuscript.

2. **Modified dynamic programming algorithms for GLOSA systems with stochastic signal switching times**  
   Panagiotis Typaldos, Markos Papageorgiou. *Transportation Research Part C*, 157 (2023), 104364.

The code focuses on deterministic and stochastic trajectory generation for vehicles approaching signalized intersections, with implementations based on Stochastic Dynamic Programming (SDP), Discrete Differential Dynamic Programming (DDDP), and Differential Dynamic Programming (DDP).

![Illustrative traffic-signal timing example](./all%20signal%20phases/MFTS22/signal_example.png)

## Research focus

The repository studies how an approaching vehicle should shape its trajectory when the next traffic-light switching time is:

- known in advance, as in fixed-time signal plans;
- uncertain but described by a switching-time distribution;
- solved with faster approximate dynamic-programming variants suitable for more demanding online use cases.

The implementation is research-oriented rather than packaged as a single application: each major source file corresponds to a concrete algorithmic experiment or paper-era comparison.

## Repository map

```text
red2green/
  dp/
    main.c                   Baseline dynamic programming implementation
  dddp/
    *.c
    mathematica/             Early stochastic and iterative DP notebooks
  ddp/
    *.nb                     Mathematica notebooks for DDP-related analysis
  journal TRC/
    sdp.c
    dddp.c
    ddp.c
    coefficients.h
    mathematica/             Notebooks associated with the TRC paper

all signal phases/
  MFTS22/
    mfts.cpp
    final_times.c
    plots.py
    signal_example.png
    *.txt                    Precomputed tables and sample outputs for this experiment set
    pet/                     Companion variant of the same experiment family
```

## How to read the repository

- `red2green/dp/` is the clearest starting point for the baseline stochastic DP formulation.
- `red2green/journal TRC/` is the most publication-ready code area for the 2023 TRC paper and contains the SDP, DDDP, and DDP comparison implementations.
- `red2green/dddp/` and `red2green/ddp/` retain earlier development notebooks and supporting implementations that help trace the research progression.
- `all signal phases/MFTS22/` contains a separate experiment family with supporting precomputed tables, plotting scripts, and an illustrative signal figure.

## Building representative programs

The programs are standalone research executables configured primarily through constants and scenario arrays near the top of each source file.

### Build with `make`

```bash
make
```

This creates:

- `build/glosa_dp`
- `build/glosa_sdp`
- `build/glosa_dddp`
- `build/glosa_ddp`

You can remove generated binaries with:

```bash
make clean
```

### Baseline DP example

```bash
./build/glosa_dp
```

### TRC paper examples

```bash
./build/glosa_sdp
./build/glosa_dddp
./build/glosa_ddp
```

Some programs write output files into folders such as `outputs/` or `results_sdp/`. These directories are treated as generated artifacts and are intentionally excluded from version control. When a particular experiment expects one of those folders, create it before running the executable.

`glosa_dddp` first looks for a saved SDP initialization trajectory. If that file is not present, it reports the missing path and automatically falls back to the built-in Scenario 1 initial trajectory so the executable remains runnable.

## Data and generated artifacts

The large `j_opt_*.txt` and `te_opt_*.txt` files under `all signal phases/MFTS22/` are retained because the accompanying C++ experiments read them as precomputed inputs. Other files named like `output_*.txt` are treated as run outputs rather than source assets.

## Notes for reproducibility

- Scenario definitions, initial conditions, and mode switches are currently encoded in the source files rather than exposed through a unified command-line interface.
- Several notebooks document derivations, plots, and experiment analysis alongside the C/C++ implementations.
- The repository has been kept close to the paper-era folder layout so that code, notebooks, and publication artifacts remain easy to cross-reference.

## Suggested repository metadata

A concise GitHub description for this repository:

> Dynamic-programming-based trajectory optimization for GLOSA systems with deterministic and stochastic traffic-signal timing.

Suggested GitHub topics:

- `glosa`
- `optimal-control`
- `dynamic-programming`
- `trajectory-optimization`
- `autonomous-vehicles`
- `traffic-signals`

## Publications and citation

For citation support on GitHub, see [`CITATION.cff`](./CITATION.cff).
