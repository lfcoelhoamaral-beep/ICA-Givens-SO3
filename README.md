# ICA-Givens-SO3

**Hybrid Independent Component Analysis on SO(3) Using Givens Rotations**

MATLAB implementation and experimental results for a hybrid independent component analysis (ICA) method combining fourth-order statistics and temporal covariance information.

## 1. Overview

This repository provides the experimental implementation of an ICA method based on a three-angle Givens rotation parameterization of the special orthogonal group SO(3).

The hybrid objective combines squared marginal kurtosis with a penalty for off-diagonal entries of time-lagged covariance matrices.

The method is evaluated against three established ICA algorithms: FastICA, JADE, and SOBI.

The study focuses on numerical separation performance and computational characteristics. It does not claim that the combination of higher-order and temporal statistics, or the use of Givens rotations in ICA, is fundamentally new.

## 2. Repository Contents

- `ica_givens_experimento4_corrigido.m` — Main MATLAB script implementing the four experimental scenarios and comparing the algorithms.
- `ica_givens_experimento4_resumo.csv` — Experimental summary containing performance statistics.

## 3. Experimental Configuration

The experiments use the following configuration:

| Parameter | Value |
|---|---|
| Number of sources | 3 |
| Samples per trial | 6000 |
| Monte Carlo trials | 20 per scenario |
| Givens multistart initializations | 20 |
| Hybrid weight | rho = 10 |
| Temporal lags | [1, 2, 3, 5, 8] |
| Optimization method | Nelder-Mead (`fminsearch`) |
| MATLAB version | R2026a Update 2 |

Four synthetic scenarios are investigated:

- **Scenario A — Distinct temporal structures:** sources with distinguishable temporal characteristics.
- **Scenario B — Equal temporal structures:** sources with similar temporal information.
- **Scenario C — Nearly Gaussian sources:** evaluation under reduced non-Gaussianity.
- **Scenario D — Complementary information:** sources with partially shared temporal characteristics.

## 4. Algorithms

The experimental comparison includes:

1. Givens-Hybrid ICA
2. FastICA
3. JADE
4. SOBI

FastICA, JADE, and SOBI require external MATLAB implementations.

The functions `fastica.m`, `jadeR.m`, and `sobi.m` must be available on the MATLAB path.

These external implementations are not redistributed in this repository. Users should obtain them from their respective authorized sources and comply with applicable licenses.

Note that implementation details, normalization conventions, and temporal-lag configurations may differ among algorithms.

## 5. Running the Experiments

1. Download or clone this repository.
2. Open MATLAB.
3. Add the required external ICA implementations to the MATLAB path.
4. Navigate to the repository directory.
5. Execute:

```matlab
ica_givens_experimento4_corrigido
```

The script evaluates all four scenarios using 20 Monte Carlo trials per scenario.

The experimental summary is provided in `ica_givens_experimento4_resumo.csv`.

## 6. Reproduced Results

The following mean Amari indices were obtained in the validated MATLAB execution:

| Scenario | Givens-Hybrid | FastICA | JADE | SOBI |
|---|---:|---:|---:|---:|
| A | 0.01006 | 0.01504 | 0.01921 | 0.00951 |
| B | 0.03230 | 0.02424 | 0.03211 | 0.38827 |
| C | 0.01118 | 0.36908 | 0.36370 | 0.01039 |
| D | 0.01882 | 0.01746 | 0.02388 | 0.16264 |

All algorithms completed 20 valid trials in each scenario.

Lower Amari values indicate better separation performance.

The Givens-Hybrid method achieved competitive performance across the four scenarios, although it was not the best-performing algorithm in every case.

## 7. Computational Considerations

The Givens-Hybrid implementation uses multiple Nelder-Mead optimizations and consequently requires substantially more computational time than the baseline algorithms in the reported experiments.

Execution times depend on hardware, MATLAB version, numerical libraries, and implementation details.

The hybrid objective uses temporal lags [1, 2, 3, 5, 8], while the SOBI implementation may use a different lag configuration. This distinction should be considered when interpreting comparisons.

## 8. Reproducibility and Limitations

The reported experiments were reproduced using MATLAB R2026a Update 2.

The main separation metrics matched the experimental results used in the associated manuscript to the displayed precision.

The experiments use synthetic mixtures with three sources. The results do not establish general superiority, global convergence, or universal identifiability of the hybrid objective.

External implementations of FastICA, JADE, and SOBI may affect reproducibility if different versions are used.

## 9. Author

Luís Fernando Coelho Amaral  
Universidade Federal do Maranhão (UFMA), Brazil

## 10. License

No software license has been assigned yet. Reuse and redistribution permissions remain reserved until a license is explicitly provided.
