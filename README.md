# Root Causal Strength using Perturbations (RCSP)

This is an R package implementing RCSP, an algorithm for discovering root causal genes from a combination of bulk RNA-seq and Perturb-seq data, each derived from possibly independent studies. RCSP estimates the root causal strength of a gene $\widetilde{X}_i$ on a target variable $Y$. The root causal strength is defined as $$\Phi_i = |E(Y|\textnormal{Pa}(\widetilde{X}_i),\widetilde{X}_i) - E(Y|\textnormal{Pa}(\widetilde{X}_i))| = |E(Y|\textnormal{Pa}(\widetilde{X}_i),E_i) - E(Y|\textnormal{Pa}(\widetilde{X}_i))|$$ where we have subsituted $\widetilde{X}_i$ with $E_i$. The root causal strength thus measures the absolute root causal effect of the root vertex error term $E_i$ on $Y$ given $\textnormal{Pa}(\widetilde{X}_i)$. We say that $\widetilde{X}_i$ is a root causal gene if $\Phi_i > 0$.

The academic article describing RCSP in detail can be found [here](https://www.biorxiv.org/content/10.1101/2024.01.13.574491v3). Please cite the article if you use any of the code in this repository.

The Experiments folder contains any additional code needed to replicate the experimental results in the paper. All code was tested in R version 4.3.1.

# Installation

> if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

> BiocManager::install(c("fgsea","qvalue","org.Hs.eg.db","reactome.db"))

> library(devtools)

> install_github("ericstrobl/RCSP")

> library(RCSP)

# Run RCSP on synthetic data
Generate DAG over 100 variables with an expected neighborhood size of 2:
> DAG = generate_DAG_big_same4(p=100,en=2)

Generate Perturb-seq data with 200 samples per perturbation:
> save_samps_by_file_mult3(DAG, nsamps=200)

Generate 200 samples of bulk RNA-seq data:
> samps = sample_DAG_NB_linear(200,DAG$DAGb)

Run RCSP:
> out = RCSP(samps, reg="MARS")

Print the signed (before taking absolute value) RCS value of each gene:
> print(cbind(out$genes,out$RCS))

# Run RCSP on AMD data
Download Bulk_data.zip, unzip its contents and place them into the working directory. Then load bulk RNA-seq data:
> load("samps_bulk_AMD.RData") 

Download Alg_outputs_AMD.zip and load descendants of each variable precomputed from Perturb-seq:
> load("desL_AMD.RData") 

Run RCSP:
> out = RCSP(samps,desL) # takes about 8 hours on my machine (2.30 GHz CPU, 16 GB RAM)

Print the signed RCS value of each gene:
> print(cbind(out$genes,out$RCS))

# Run RCSP on MS data
Download Bulk_data.zip, unzip its contents and place them into the working directory. Then load bulk RNA-seq data:
> load("samps_bulk_MS.RData")

Download Alg_outputs_MS.zip and load descendants of each variable precomputed from Perturb-seq:
> load("desL_MS.RData")

Run RCSP:
> out = RCSP(samps,desL) # takes about 2 hours on my machine

Print the signed RCS value of each gene:
> print(cbind(out$genes,out$RCS))

# Compute the list of descendants ('desL' files)
The Perturb-seq datasets are large, so we provide the desL files pre-computed. However, if you would like to compute the files on your own, then download 'ReplogleWeissman2022_rpe1.h5ad' from https://zenodo.org/records/10044268 and divide the h5ad file into manageable chunks:
> load("samps_bulk_AMD.RData")

> save_data_chunks_RPE1(colnames(samps$data),"directory/to/ReplogleWeissman2022_rpe1.h5ad")

Then compute desL for AMD:
> desL_AMD = find_des_AMD_final()$desL

Similarly, download 'ReplogleWeissman2022_K562_gwps.h5ad' from https://zenodo.org/records/10044268 and divide the h5ad file into manageable chunks:
> load("samps_bulk_AMD.RData")

> save_data_chunks_K562(colnames(samps$data),"directory/to/ReplogleWeissman2022_K562_gwps.h5ad")

Compute desL for MS:
> desL_MS = find_des_MS_final()$desL

