# moordyn-syrope-tests

Small driver project to validate SYROPE-style working-curve behavior against MoorDyn outputs.

## Purpose

This repository runs a set of MoorDyn line cases (linear, quadratic, exponential working curves),
computes an analytical/reference tension response, and compares both using an L2 relative error metric.

It is intended for:

- quick regression checks of working-curve behavior,
- side-by-side comparison of simulated vs analytical tension,


## Repository layout

- `main.cpp`: C++ validation driver.
- `input/`: OWC table (`owc.dat`) and case input folders (`linear_wc`, `quadratic_wc`, `exp_wc`).
- `scripts/mean_tension_verification.py`: Python post-processing/verification script.
- `scripts/syropepy.py`: SYROPE constitutive helper class used by the Python verifier.
- `CMakeLists.txt`: build configuration linking against installed MoorDyn + Eigen.

## Usage

### 1) Prerequisites

- CMake 3.18+
- C++ compiler with C++17 support
- Eigen3
- MoorDyn installed with CMake package files available at:
	- `/opt/moordyn/lib/cmake/moordyn`

If your MoorDyn install is elsewhere, edit `MoorDyn_DIR` in `CMakeLists.txt`.

### 2) Build

From this directory:

```bash
cmake -S . -B build
cmake --build build -j
```

### 3) Run

```bash
./build/main [--owc path] [--case all|linear|quadratic|exponential] [--superimpose-fast 0|1]
```

Defaults:

- `--owc input/owc.dat`
- `--case all`
- `--superimpose-fast 1`

Examples:

```bash
./build/main
./build/main --case linear
./build/main --case exponential --superimpose-fast 0
./build/main --owc input/owc.dat --case all
```

Driver notes:

- `--case` values for C++ are `linear`, `quadratic`, `exponential`, or `all`.
- `--superimpose-fast` accepts `0|1|true|false`.
- The C++ driver runs MoorDyn for the selected case(s) using each case `line.txt` under `input/`.

## Mean tension verification script

`scripts/mean_tension_verification.py` compares MoorDyn output against the SYROPE-based reference mean-tension history.

### CLI

```bash
python3 scripts/mean_tension_verification.py <all|linear|exp|quadratic> <True|False>
```

- First argument selects case(s):
	- `linear`, `exp`, `quadratic`: run one case
	- `all`: run all three in sequence
- Second argument is `has_fast_loading` and must be exactly `True` or `False`.

Note: Python mode names use `exp`, while the C++ driver uses `exponential`.

Examples:

```bash
python3 scripts/mean_tension_verification.py linear False
python3 scripts/mean_tension_verification.py exp True
python3 scripts/mean_tension_verification.py all True
```

### Outputs

For each executed case, the script writes:

- `tension_strain.png` in the same folder as the input `line.txt`
	- e.g. `input/linear_wc/tension_strain.png`

Console output includes:

- parsed case summary
- per-case L2 error (`L2-error between tension and tension_ref`)
- final end-of-run summary with:
	- working-curve formula (`WCType`, `k1`, `k2`)
	- whether fast loading is applied
	- where fast loading is applied (plot right panel strain-tension traces)
	- L2 error

## Python dependencies

The verification script requires:

- Python 3
- NumPy
- Matplotlib

Install example:

```bash
python3 -m pip install numpy matplotlib
```

## License

MIT (see `LICENSE`).

## Citation
* Wei, Zhilong; Bingham, Harry B.; Shao, Yanlin (2026). ESOMOOR D5.1: Extended MoorDyn solver and validation report. Technical University of Denmark. Online resource. https://doi.org/10.11583/DTU.31408806

## References

* Falkenberg, E., Åhjem, V., & Yang, L. (2017). Best practice for analysis of polyester rope mooring systems. Offshore Technology Conference, D031S034R006. https://doi.org/10.4043/27761-ms
* Hall, M., Duong, B., & Lozon, E. (2023, December). Streamlined Loads Analysis of Floating Wind Turbines With Fiber Rope Mooring Lines. ASME 2023 5th International Offshore Wind Technical Conference. https://doi.org/10.1115/iowtc2023-119524