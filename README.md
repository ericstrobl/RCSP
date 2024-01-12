# Root Causal Strength from Perturbations (RCSP)

This is an R package implementing RCSP, an algorithm for discovering root causal genes from a combination of bulk RNA-seq and Perturb-seq. RCSP estimatrd root causal strength of a variable $X_i$ on a target variable $Y$. The root causal strength score is defined as $$|E(Y|\textnormal{Pa}(X_i),X_i) - E(Y|\textnormal{Pa}(X_i))| = |E(Y|\textnormal{Pa}(X_i),E_i) - E(Y|\textnormal{Pa}(X_i))|$$ and thus measures the root causal effect of the root vertex error term $E_i$ on $Y$ given $Pa(X_i)$.

The academic article describing DDR in detail can be found [here](). Please cite the article if you use any of the code in this repository.

# Installation

> library(devtools)

> install_github("ericstrobl/RCSP")

> library(RCSP)

# Run RCSP on sythetic data

> desL = find_des2(samps)

> out = RCSP(samps,desL)  ##

# Run RCSP on real data

> load("MS_bulk_and_perturb.RData")

> desL = find_des2(samps)

> out = RCSP(samps,desL)  ## outputs RCS scores and gene names
