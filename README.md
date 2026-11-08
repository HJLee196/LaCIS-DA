# README

## Contents

- [Dataset A Preprocessing](Preprocess/README.md)
- [Datasets and Preprocessing](#datasets-and-preprocessing)
- [Section 5.2 Code](RealDataAnalysis/ReadMe_Section52_RealDataAnalysisCode.md)
- [Section 5.2 - Write All](RealDataAnalysis/ReadMe_Section52_RealDataAnalysisCode_WriteAll.md)
- [Section 5.4 Code](SimulationStudy/ReadMe_Section54.md)
- [Section 5.4 - Write All](SimulationStudy/ReadMe_Section54_WriteAll.md)
- [Section B2 Code](SimulationStudy/ReadMe_SectionB2.md)
- [Section B2 - Write All](SimulationStudy/ReadMe_SectionB2_WriteAll.md)
- [Section 5 Figure 1 and S4 Code](RealDataAnalysis/ReadMe_Section5_Figures1andS4.md)


## Datasets and Preprocessing

### Dataset A — Alsema et al. (2020)

For Dataset A, we provide preprocessing code to reconstruct a Seurat object from the publicly available GEO data (GSE146639) and published donor metadata.

The reconstructed count matrix differs slightly from the matrix used in our original simulation studies. Among entries with positive counts in the original simulation matrix, approximately 99.78% had identical counts in both matrices. Rerunning the simulations with the GEO-derived matrix yielded results nearly identical to those reported in the manuscript.

The preprocessing scripts, required inputs, execution instructions, and limitations are documented in [Dataset A Preprocessing](Preprocess/README.md).

### Dataset B — Yang et al. (2022)

For Dataset B, our analysis began with the nucleus-level gene counts and cell-type annotations provided by [Yang et al. (2022)](https://doi.org/10.1038/s41586-021-04369-3).

The original investigators performed quality control and assigned cell-type labels using their expertise in human brain tissue, as described in the Methods section of their paper. Their quality-control procedure included filtering low-quality nuclei, removing doublets and multiplets, and reviewing cell-type-specific marker expression.

Because our analysis began with these processed data, the supplementary code covers only the steps we performed. 

The preprocessed data are available upon request from [lil155@pitt.edu](mailto:lil155@pitt.edu).

For the analysis code and execution instructions, see [Section 5.2 Code](RealDataAnalysis/ReadMe_Section52_RealDataAnalysisCode.md) and [Section 5.4 Code](SimulationStudy/ReadMe_Section54.md).


## File Naming Convention

The R scripts follow a numbering convention to indicate the analysis setting and statistical method.

For most analysis scripts, the leading two-digit number has the following structure:
- **First digit**: The analysis setting or simulation scenario.
- **Second digit**: The knockoff construction or variable selection method.

For example, `43_44` in `43_44_mother_code_multiLR_realdata-AllPart12.R` indicates real-data analysis (first digit `4`) using multi-LR knockoffs (second digit `3`) and the corresponding e-LR procedure (second digit `4`).

### Analysis Settings

| First digit | Description | Corresponding section |
|---|---|---|
| `0` | Standard simulation | Appendix B.2 |
| `1` | Negative binomial simulation | Appendix B.5 |
| `3` | Simulation with confounding effects | Section 5.4 |
| `4` | Real-data analysis | Section 5.2 |
| `5` | Simulation with de-biased imputation | Appendix E.5 |
| `6` | Real-data analysis using covariance estimation without conditioning on covariates | Appendix A.3 |
| `7` | Simulation with knockoffs generated directly from observed data | Appendix B.3 |
| `8` | Simulation with confounding effects and de-biased imputation | Appendix E.5 |

### Statistical Methods

| Second digit | Method |
|---|---|
| `1` | ASDP knockoffs |
| `2` | LR knockoffs |
| `3` | Multi-LR knockoffs |
| `4` | e-LR procedure |
| `5` | Decomp knockoffs |
| `6` | Multi-decomp knockoffs |
| `7` | e-decomp procedure |
| `8` | BH procedure |
| `9` | Standard Gaussian knockoffs (naive) |

Some scripts implement both a multiple-knockoff method and its corresponding e-BH procedure. These scripts have two numbers separated by an underscore, such as `03_04` (multi-LR and e-LR) or `06_07` (multi-decomp and e-decomp). For the simulation scripts, e-BH results are computed in the corresponding summary scripts. For the real-data AllPart12 scripts, e-BH selections are computed and saved in Part 2, and the summary script reads these saved selections.

### WriteAll and SummaryAll Scripts

Scripts beginning with `99_` serve two purposes:
- `99_WriteAll_*`: Generate analysis scripts with different settings and their corresponding SLURM batch scripts.
- `99_SummaryAll_*`: Summarize analysis results and generate tables or figures.

The number following `WriteAll` or `SummaryAll` identifies the corresponding analysis or script group.

## Shared Code

The `Section_Appendix_A1_Rank_Estimation` folder contains
the code for rank estimation using Ke's BEMA method, as
described in Appendix A.1. This code is used in both
real-data analysis and simulation studies and is therefore
maintained separately to avoid duplication.

## License

This repository is licensed under the GNU General Public License v3.0 (GPL-3.0-only). See the [LICENSE](LICENSE) file for details.