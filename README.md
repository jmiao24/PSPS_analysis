# PSPS experiments

## Simulations
All simulations are done using `R` with version 4.2.1 (2022-06-23). The required packages are `POPInf`, `IPD`, `quantreg`, `doParallel`, `data.table`, `knockoff`, `AER`, `MASS`, `sandwich`, and `lmtest`.

### Tasks that have been implemented for ML-assisted inference
* `./old_task/Fun.R`: Functions for simulation and PSPS protocol
* `./old_task/mean.R`: mean estimation
* `./old_task/ols.R`: linear regression
* `./old_task/logistic.R`: logistic regression

### Tasks that have not been implemented for ML-assisted inference
* `./new_task/Fun.R`: Functions for simulation and PSPS protocol
* `./new_task/qr.R`: quantile regression
* `./new_task/ngbr.R`: negative binomial regression
* `./new_task/iv.R`: instrumental variables regression
* `./new_task/dlasso.R`: debiased Lasso
* `./new_task/ranksum_t1e.R`: Wilcoxon rank-sum test (type-I error)
* `./new_task/ranksum_power.R`: Wilcoxon rank-sum test (power)

### FDR control
* `./fdr/Fun.R`: Functions for simulation and PSPS protocol
* `./fdr/Fun_PSPS-knockoff.R`: Functions for PSPS-knockoff
* `./fdr/BH`: PSPS-BH simulation
* `./fdr/data_knockoff`: generate data for PSPS-knockoff simulations
* `./fdr/dlasso_knockoff`: fit PSPS-dlasso on the generated data
* `./fdr/summarize_knockoff`: fit PSPS-knockoff

## Real data analysis
The real data analysis is done in UK Biobank, which is avaiable by application [this link](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access). We used [plink2](https://www.cog-genomics.org/plink/2.0/) and [QUAIL](https://github.com/qlu-lab/QUAIL) for fitting quantile regression to identify vQTL.
* `./real/Fun.R`: Functions for PSPS protocol
* `./real/vQTL.sh`: run genome-wide vQTL analysis
* `./real/combine.R`: apply PSPS protocol
