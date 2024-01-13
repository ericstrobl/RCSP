# Root Causal Strength from Perturbations (RCSP)

This is an R package implementing RCSP, an algorithm for discovering root causal genes from a combination of bulk RNA-seq and Perturb-seq data, each derived from possibly independent studies. RCSP estimates the root causal strength of a variable $X_i$ on a target variable $Y$. The root causal strength is defined as $$\Phi_i = |E(Y|\textnormal{Pa}(X_i),X_i) - E(Y|\textnormal{Pa}(X_i))| = |E(Y|\textnormal{Pa}(X_i),E_i) - E(Y|\textnormal{Pa}(X_i))|$$ and thus measures the absolute root causal effect of the root vertex error term $E_i$ on $Y$ given $\textnormal{Pa}(X_i)$. We say that $X_i$ is a root causal gene if $\Phi_i > 0$.

The academic article describing RCSP in detail can be found [here](). Please cite the article if you use any of the code in this repository.

# Installation

> if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

> BiocManager::install(c("fgsea","qvalue","org.Hs.eg.db"))

> library(devtools)

> install_github("ericstrobl/RCSP")

> library(RCSP)

# Run RCSP on sythetic data
Generate DAG over 100 variables with an expected neighborhood size of 2:
> DAG = generate_DAG_big_same4(p=100,en=2)

Generate Perturb-seq data with 200 samples per perturbation:
> save_samps_by_file_mult3(DAG, nsamps=200)

Generate 200 samples of bulk RNA-seq data:
> samps = sample_DAG_NB_linear(200,DAG$DAGb)

Run RCSP:
> out = RCSP(samps, reg="MARS")

# Run RCSP on AMD data
Download the Data folder and place its contents into the working directory. Then load bulk RNA-seq data:
> load("samps_bulk_AMD2.RData") 

Load descendants of each variable precomputed from Perturb-seq:
> load("desL_AMD.RData") # processed Perturb-seq data will be uploaded on Zenodo (200 MB)

Run RCSP:
> out = RCSP(samps,desL) # takes about 8 hours on my machine (2.30 GHz CPU, 16 GB RAM)

# Run RCSP on MS data
Download the Data folder and place its contents into the working directory. Then load bulk RNA-seq data:
> load("samps_bulk_MS2.RData")

Load descendants of each variable precomputed from Perturb-seq:
> load("desL_nonsparse_MS.RData") # processed Perturb-seq data will be uploaded on Zenodo (4.5 GB)

Run RCSP:
> out = RCSP(samps,desL) # takes about 2 hours on my machine
