# Root Causal Strength from Perturbations (RCSP)

This is an R package implementing RCSP, an algorithm for discovering root causal genes from a combination of bulk RNA-seq and Perturb-seq data, each derived from possibly independent studies. RCSP estimates the root causal strength of a variable $X_i$ on a target variable $Y$. The root causal strength is defined as $$|E(Y|\textnormal{Pa}(X_i),X_i) - E(Y|\textnormal{Pa}(X_i))| = |E(Y|\textnormal{Pa}(X_i),E_i) - E(Y|\textnormal{Pa}(X_i))|$$ and thus measures the absolute root causal effect of the root vertex error term $E_i$ on $Y$ given $\textnormal{Pa}(X_i)$.

The academic article describing RCSP in detail can be found [here](). Please cite the article if you use any of the code in this repository.

# Installation

> library(devtools)

> install_github("ericstrobl/RCSP")

> library(RCSP)

# Run RCSP on sythetic data

> DAG = generate_DAG_big_same4(p=100,en=2) # generate DAG over 100 variables with an expected neighborhood size of 2

> nsamps = 200 # number of samples for bulk RNA-seq and per perturbation for Perturb-seq

> save_samps_by_file_mult2(DAG, sc_batches=1:DAG$DAGs$n_batch, nsamps=200) # generate Perturb-seq data

> samps = sample_DAG_NB_linear(nsamps,DAG$DAGb) # generate bulk RNA-seq data

> out = RCSP(samps)  # run RCSP

# Run RCSP on real data

> load("MS_bulk_and_perturb.RData")

> desL = find_des2(samps)

> out = RCSP(samps,desL)  ## outputs RCS scores and gene names
