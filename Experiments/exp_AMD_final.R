

# load bulk RNA-seq data
load("samps_bulk_AMD.RData")
genes = colnames(samps$data)

# load precomputed descendants
load("desL_AMD.RData")

# age histogram
RCS_age = get_ctrl_RCE_nonsparse(samps,2)
sqrt(mean(RCS_age^2))
hist(abs(RCS_age))

# Run RCSP
# alg_out = RCSP(samps,desL)  ## takes ~8 hours
load("RCSP_AMD.RData")

# Compute D-SD
# alg_out_DSD = DSDP(samps,desL)  ## takes ~8 hours
load("DSDP_AMD.RData")

# histogram of D-RCS vs D-SD
hist(sqrt(colMeans(alg_out$RCS^2)))
hist(sqrt(colMeans(alg_out_DSD$RCS^2)))

# D-RCS values of genes
ix = order(sqrt(colMeans(alg_out$RCS^2)),decreasing=TRUE)
cbind(alg_out$genes[ix],sqrt(colMeans(alg_out$RCS^2))[ix])

# UMAP embedding
require(uwot)
cov0 = cov(abs(cbind(alg_out$RCS,RCS_age)))
eig = eigen(cov0)
data0 = abs(cbind(alg_out$RCS,RCS_age)) %*% eig$vectors[,1:10]

data.umap = uwot::umap(data0)
plot(data.umap[,1],data.umap[,2])

# cluster SS plot
require(fastcluster)
cl <- hclust.vector(data.umap, method="ward")
plot(rev(cl$height)[1:20])
ix = cutree(cl, k = 4)

# UMAP embedding with clusters
plot(data.umap,col=ix)

# pathway enrichment analysis

# require(fastcluster)
# cl <- hclust.vector(data.umap, method="ward") # use this for cluster-specific enrichment analysis
# ix = cutree(cl, k = 4)

stats = c()
for (j in 1:ncol(alg_out$RCS)){
  stats = c(stats,sqrt(mean(alg_out$RCS[,j]^2)))
  # stats = c(stats,sqrt(mean(alg_out$RCS[ix==1,j]^2))) # use this for cluster-specific enrichment analysis, ix=1 means cluster 1
}
genesC = alg_out$genes

require(org.Hs.eg.db)
hs <- org.Hs.eg.db
entrez = AnnotationDbi::select(org.Hs.eg.db, keys = genesC, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ix = match(genesC,entrez$SYMBOL)
iR = which(is.na(entrez$ENTREZID[ix]))
stats = stats[-iR]
entrez = as.character(entrez$ENTREZID[ix[-iR]])
names(stats) = entrez

require(fgsea)
pathways <- reactomePathways(names(stats))

fgseaRes <- fgsea(pathways, stats, nPermSimple = 100000, scoreType="pos")

collapsedPathways <- collapsePathways(fgseaRes[order(pval)], pathways, stats)
print(fgseaRes[pathway %in% collapsedPathways$mainPathways][order(pval)][1:20,])

# drug enrichment analysis

# require(fastcluster)
# cl <- hclust.vector(data.umap, method="ward") # use this for cluster-specific enrichment analysis
# ix = cutree(cl, k = 4)

stats = c()
for (j in 1:ncol(alg_out$RCS)){
  stats = c(stats,sqrt(mean(alg_out$RCS[,j]^2)))
  # stats = c(stats,sqrt(mean(alg_out$RCS[ix==1,j]^2))) # use this for cluster-specific enrichment analysis, ix=1 means cluster 1
}
genesC = alg_out$genes

require(org.Hs.eg.db)
hs <- org.Hs.eg.db
entrez = AnnotationDbi::select(org.Hs.eg.db, keys = genesC, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ix = match(genesC,entrez$SYMBOL)
iR = which(is.na(entrez$ENTREZID[ix]))
stats = stats[-iR]
entrez = as.character(entrez$ENTREZID[ix[-iR]])
names(stats) = entrez

library(org.Hs.eg.db)
# download this file: https://dsigdb.tanlab.org/Downloads/DSigDB_All_detailed.txt
dsig <- readr::read_tsv("directory to above file")

drug2gene=dsig[, c("Drug", "Gene")]

drugs = unique(dsig$Drug)
drug2gene = vector("list", length(drugs))
entrez = AnnotationDbi::select(org.Hs.eg.db, keys = dsig$Gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
ix = match(dsig$Gene,entrez$SYMBOL)
entrez = entrez$ENTREZID[ix]

for (d in 1:length(drugs)){
  id = which(dsig$Drug == drugs[d])
  ee = unique(entrez[id])
  drug2gene[[d]] <- ee[!is.na(ee)]
  
}
names(drug2gene) = drugs

fgseaRes <- fgsea(drug2gene, stats, nPermSimple = 100000, scoreType="pos")

head(fgseaRes[order(pval), ])
-log(fgseaRes[order(pval)][1:6,padj])

# graded gene UMAP embedding
require(dplyr)
gene_name = "SLC7A5"
plot(data.umap[,1],data.umap[,2],col=make_colour_gradient(ntile(abs(alg_out$RCS[,alg_out$genes==gene_name]),4)))

# graded severity UMAP embedding
plot(data.umap[,1],data.umap[,2],col=make_colour_gradient(ntile(samps$data[,ncol(samps$data)], 4)))


### correlation with UMAP dimension
ix = order(sqrt(colMeans(alg_out$RCS^2)),decreasing=TRUE)

stats = c()
ps = c()
CIs = c()
dim = 2 # choose UMAP dimension
for (j in 1:30){
  stats = c(stats,cor.test(data.umap[,dim],abs(alg_out$RCS[,ix[j]]),method="spearman")$estimate)
  ps = c(ps,cor.test(data.umap[,dim],abs(alg_out$RCS[,ix[j]]),method="spearman")$p.value)
  CIs = rbind(CIs,spearman_CI(data.umap[,dim],abs(alg_out$RCS[,ix[j]])))
}
iy = order(abs(stats),decreasing=TRUE)
mat = cbind(alg_out$genes[ix[1:30]][iy],stats[iy],CIs[iy,])
require(qvalue)
qvalue(ps[iy],pi0=1)

