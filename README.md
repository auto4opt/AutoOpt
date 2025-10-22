## AutoOptLib – Repository Guide

This repository hosts both the original MATLAB/Octave AutoOptLib implementation and the new Python port. The MATLAB directories remain unchanged so existing scripts continue to work; this document focuses on the Python workflow while pointing to the MATLAB/Octave resources when relevant. Reviewers should be able to reproduce every experiment described in the manuscript using the instructions below.

---

### 1. Environment & Dependencies

| Component | Recommended Version | Notes |
|-----------|--------------------|-------|
| Python    | ≥ 3.10             | Use a virtual environment (venv/conda). |
| Mandatory packages | `numpy`, `scipy`, `pytest`, `pyyaml` | Install via `pip install -r requirements.txt` (to be generated). |
| Optional | `matplotlib` | Only needed for plotting convergence curves. |
| Data files | `shift_data.txt`, `M_D*.mat/txt`, etc. | Already provided under `AutoOptLib for Matlab/Problems/.../data`. |

> **MATLAB/Octave note**: The Python version does not rely on MATLAB/Octave. For the original MATLAB code (located in `AutoOptLib for Matlab/`), ensure that the Statistics and Machine Learning Toolbox is installed. For Octave (`AutoOptLib for Octave/`), install the `statistics`, `control`, `signal`, `communications`, and `io` packages and execute `addpath(genpath('.'))` before calling `AutoOpt`.

---

### 2. Installation

```bash
git clone https://github.com/auto4opt/AutoOpt.git
cd AutoOpt
python -m venv .venv
source .venv/bin/activate      # On Windows: .venv\Scripts\activate
pip install -e .
```

After installation the package can be imported as `autooptlib`.

---

### 3. Quick Verification

```bash
pytest tests/unit
```

The unit suite covers all components (selection, search, crossover, reinitialization, update, parameter updates) and key utilities (`Design`, `space`, `cec2013_f1`). Expected runtime is < 10 seconds.

---

### 4. Minimal Examples

#### 4.1 Design mode

```python
from autooptlib import autoopt

final_algs, alg_trace = autoopt(
    Mode="design",
    Problem="cec2013_f1",
    InstanceTrain=[30, 30, 30],
    InstanceTest=[30],
    AlgN=3,
    AlgFE=120,
    ProbN=20,
    ProbFE=4000,
    Evaluate="racing",
    Compare="statistic",
    archive=["archive_best"],
    Seed=42,                          # ensures reproducibility
)
```

Output files are written to the working directory:

| File | Description |
|------|-------------|
| `Algs.pkl` | Serialization of the designed algorithms and performance history. |
| `Algs_best_algs_iter.csv` | Pseudocode of the best algorithm at each iteration. |
| `Algs_perf_final_algs.csv` | Performance matrix (instances × runs) on the test set. |
| `ConvergenceCurve.csv` | Design-stage average fitness over iterations. |

#### 4.2 Solve mode

```python
best_solutions, all_solutions = autoopt(
    Mode="solve",
    Problem="cec2013_f1",
    InstanceSolve=[30],
    AlgName="ICA",        # or "Continuous Random Search", "Differential Evolution", etc.
    AlgRuns=2,
    ProbN=30,
    ProbFE=5000,
    Metric="quality",
    Seed=123,
)
```

Output files:

| File | Description |
|------|-------------|
| `Solutions.pkl` | Best solutions per run and the full history. |
| `Solutions.csv` | Tabulated decision variables for each recorded iteration. |
| `Fitness_all_runs.csv` | Mean/std and per-run objective values. |

---

### 5. Parameters & Evaluation Budget

| Parameter | Description | MATLAB default | Python recommendation |
|-----------|-------------|----------------|-----------------------|
| `ProbN`   | Population size per solve iteration | 20 | 20–50 |
| `ProbFE`  | Function evaluations per solve | 5000 | 4000–20000 |
| `AlgN`    | Number of candidate algorithms maintained | 10 | 3–10 |
| `AlgFE`   | Total algorithm evaluations during design | 5000 | 120–2000 (tunable) |
| `AlgRuns` | Repetitions per candidate | 5 | 1–3 (increase for final runs) |
| `Evaluate`| Evaluation policy (`exact`, `racing`, `approximate`, `intensification`) | exact | `racing` recommended to reduce cost |
| `Compare` | Selection rule (`average`, `statistic`) | average | pair with `racing` via `statistic` |
| `IncRate` | Minimum improvement rate for inner loop | 0.05 | 0.03–0.05 |

**Accounting for evaluation cost:**  
Reviewers requested explicit accounting of design-phase cost. AutoOptLib now exposes the iteration count and allows inspection of `AlgTrace`. Users should log:
- `ProbFE * AlgRuns` (per candidate)  
- `AlgFE` (total evaluations across all candidates)  
- CPU time (via the surrounding script)  
These numbers should be reported when comparing with tuned baselines.

---

### 6. Training/Test Split & Overfitting Checks

*AutoOptLib* expects users to explicitly provide `InstanceTrain` and `InstanceTest`. The design pipeline records training performance in `Algs.pkl`; the test results are stored in `Algs_perf_final_algs.csv`. To inspect overfitting:
1. Aggregate `performance` from training runs.
2. Compare with the test matrix.
3. Repeat the whole design process (changing `Seed`) to evaluate variability.

For CEC2013, the typical setup is each function per dimension treated as a separate target. Example: `[30, 30, 30]` for training and `[30]` for testing means all instances use 30 dimensions.

---

### 7. Random Seeds & Reproducibility

- All random numbers rely on `numpy.random` via `default_rng`.
- Pass `Seed=<int>` to `autoopt`, or call `np.random.seed` before invoking.
- For exact reproducibility across runs, also fix `Evaluate`, `Compare`, `AlgFE`, `ProbFE`, and `ProbN`.
- The unit/system tests pin specific seeds to guarantee deterministic results.

---

### 8. Intermediate Checkpoints

Design can be lengthy. To avoid losing progress (a reviewer concern):
- Each iteration’s best algorithm is already stored in `Algs_best_algs_iter.csv`.
- If additional checkpoints are required, wrap `autoopt` in a script that periodically copies `Algs.pkl` (or increments an iteration counter).
- Future releases may offer a built-in checkpoint interval; users can open an issue if this feature is desired.

---

### 9. Problem Library

| Problem | Status | Notes |
|---------|--------|-------|
| CEC2013 f1 | Implemented | Available in `autooptlib.problems.cec2013`. |
| CEC2013 f2–f28 | Planned | Translation requires additional helper functions (`computeTAsym`, `constructLambda`, etc.) and data files. |
| Beamforming, portfolio selection, etc. | Planned | Can be migrated from the MATLAB directory as needed. |

Adding a problem only requires implementing the three callable modes used by MATLAB (`construct`, `repair`, `evaluate`) or supplying a custom Python function.

---

### 10. Testing & Continuous Integration

| Command | Purpose |
|---------|---------|
| `pytest tests/unit` | Unit tests for all components and utilities. |
| `pytest tests/system` | (To be added) End-to-end regression tests. |

Setting up continuous integration (GitHub Actions/GitLab CI) is straightforward: install the package, run the test commands, and optionally publish the documentation. This addresses the reviewers’ concern about syntactically incorrect code snippets in PDFs.

---