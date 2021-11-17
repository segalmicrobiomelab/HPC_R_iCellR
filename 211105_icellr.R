#install iCellR 
#run through tutorial 

sessionInfo()

setwd("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/")
getwd()

install.packages("iCellR",lib="~/R/lib", repos='http://cran.us.r-project.org')
install.packages("gridExtra",lib="~/R/lib", repos='http://cran.us.r-project.org')

# Tutorial (dont need)
# save the URL as an object
# sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/pbmc3k_filtered_gene_bc_matrices.tar.gz"

# Tutorial (dont need) 
# download the file
# download.file(url = sample.file.url, 
#     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", 
#     method = "auto")  

# Tutorial (dont need) 
# unzip the file. 
# untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")

library("iCellR",lib="~/R/lib")
library("gridExtra",lib="~/R/lib")

# Each sample needs to be loaded separately

# Don't use this 
my.data <- load10x("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/2021-04-30/cellranger/aggregated/outs/count/filtered_feature_bc_matrix/")

my.data.single <- load10x("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/2021-04-30/cellranger/count-1_MOC_LU/outs/filtered_feature_bc_matrix/")

my.data.single.extra <- load10x("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/2021-04-30/cellranger/count-2_MOC_LU/outs/filtered_feature_bc_matrix/")

my.data.single.chronic <- load10x("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/2021-04-30/cellranger/count-C_MOC_LU/outs/filtered_feature_bc_matrix/")

my.data.PBS <- load10x("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/2021-04-30/cellranger/count-PBS_LU/outs/filtered_feature_bc_matrix/")

# Tutorial 
# This directory includes; barcodes.tsv, genes.tsv/features.tsv and matrix.mtx files
# Data could be zipped or unzipped.
# if your data is in a csv or tsv format read it like this example
# my.data <- read.delim("CITE-Seq_sample_RNA.tsv.gz",header=TRUE)
# if your data is in a h5 format read it like this example
# data <- load.h5("filtered_feature_bc_matrix.h5")

?load10x

dim(my.data)
head(my.data)[1:2]
head(my.data)[1:100]

dim(my.data.single)
head(my.data.single)[1:2]
head(my.data.single)[1:100]

dim(my.data.single.extra)
head(my.data.single.extra)[1:2]
head(my.data.single.extra)[1:100]

dim(my.data.single.chronic)
head(my.data.single.chronic)[1:2]
head(my.data.single.chronic)[1:100]

dim(my.data.PBS)
head(my.data.PBS)[1:2]
head(my.data.PBS)[1:100]

# Merging data together

#merge all of your samples to make a single aggregated file   
my.data.all <- data.aggregation(samples = c("my.data.single","my.data.single.extra", "my.data.single.chronic", "my.data.PBS"),
    condition.names = c("single.MOC","single.MOC.extra","chronic.MOC","PBS"))
    
dim(my.data.all)
head(my.data.all)[1:2]
head(my.data.all)[1:100]
    
    
my.obj <- make.obj(my.data.all)
my.obj    
    
my.obj <- qc.stats(my.obj)

# plot UMIs, genes and percent mito all at once and in one plot. 
# you can make them individually as well, see the arguments ?stats.plot.

png('three.in.one.png', width=10, height=10, units = 'cm', res=1200)
stats.plot(my.obj,
	plot.type = "three.in.one",
	out.name = "UMI-plot",
	interactive = FALSE,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green")
dev.off()	

# Scatter plots - doesn't work well 
# png('point.mito.umi.png', width=8, height=8, units = 'cm', res=600)
# stats.plot(my.obj, plot.type = "point.mito.umi", out.name = "mito-umi-plot")
# dev.off()

# png('point.gene.umi.png', width=8, height=8, units = 'cm', res=600)
# stats.plot(my.obj, plot.type = "point.gene.umi", out.name = "gene-umi-plot")
# dev.off()

# Filtering to keep cells of interest
my.obj <- cell.filter(my.obj,
	min.mito = 0,
	max.mito = 0.05,
	min.genes = 200,
	max.genes = 2400,
	min.umis = 0,
	max.umis = Inf)
	
# Normalization of samples, Reza example is for top 500 seems reasonable
my.obj <- norm.data(my.obj, 
     norm.method = "ranked.glsf",
     top.rank = 500) # best for scRNA-Seq
     
# Another QC step [optional] but will see
my.obj <- qc.stats(my.obj,which.data = "main.data")

png('three.in.one.filter.png', width=20, height=10, units = 'cm', res=1200)
stats.plot(my.obj,
	plot.type = "all.in.one",
	out.name = "UMI-plot",
	interactive = F,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green",
	back.col = "white")  
dev.off()

# Brief Gene stat output 
my.obj <- gene.stats(my.obj, which.data = "main.data")

head(my.obj@gene.data[order(my.obj@gene.data$numberOfCells, decreasing = T),])

# See model plot 
make.gene.model(my.obj, my.out.put = "plot",
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	out.name = "gene.model")
	
# Write the gene model data into the object

my.obj <- make.gene.model(my.obj, my.out.put = "data",
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	out.name = "gene.model")

head(my.obj@gene.model)
# "ACTB"  "ACTG1" "ACTR3" "AES"   "AIF1"  "ALDOA"

# get html plot (optional)
make.gene.model(my.obj, my.out.put = "plot",
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = T,
	out.name = "plot4_gene.model")
	
# Principal component analysis 

# If you run PCA (run.pca) there would be no batch alignment but if you run CPCA (using iba function) this would perform batch alignment and PCA after batch alignment. Example for batch alignment using iba function: 
# my.obj <- iba(my.obj,dims = 1:30, k = 10,ba.method = "CPCA", method = "gene.model", gene.list = my.obj@gene.model)

# run PCA in case no batch alignment is necessary
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main")

png('opt.pcs.plot.png', width=30, height=30, units = 'cm', res=1200)
opt.pcs.plot(my.obj)
dev.off()

# tSNE
my.obj <- run.pc.tsne(my.obj, dims = 1:10)

# UMAP
my.obj <- run.umap(my.obj, dims = 1:10)

# KNetL (for lager than 5000 cell use a zoom of about 400) 
# Because knetl has a very high resolution it's best to use a dim of 20 (this usually works best for most data)
my.obj <- run.knetl(my.obj, dims = 1:20, zoom = 400, dim.redux = "umap") # (Important note!) don't forget to set the zoom in the right range   

library(gridExtra)
A= cluster.plot(my.obj,plot.type = "pca",interactive = F)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F)
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F)

library("gridExtra",lib="~/R/lib")

png('four.plots.grid.png', width=30, height=30, units = 'cm', res=1200)
grid.arrange(plot(A),plot(B),plot(C), plot(D))
dev.off()

png('A.png', width=10, height=10, units = 'cm', res=1200)
plot(A)
dev.off()

png('B.png', width=10, height=10, units = 'cm', res=1200)
plot(B)
dev.off()

png('C.png', width=10, height=10, units = 'cm', res=1200)
plot(C)
dev.off()

png('D.png', width=10, height=10, units = 'cm', res=1200)
plot(D)
dev.off()

### Clustering ### 

# clustering based on KNetL
my.obj <- iclust(my.obj, sensitivity = 150, data.type = "knetl") 

# clustering based on PCA
my.obj <- iclust(my.obj, sensitivity = 150, data.type = "pca", dims=1:10)

# play with k to get the clusters right. Usually 150 is good.

# clustering based on PCA
my.obj <- iclust(my.obj,
    dist.method = "euclidean",
    sensitivity = 100,
    dims = 1:10,
    data.type = "pca")

# plot clusters (in the figures below clustering is done based on KNetL) 
# example: # my.obj <- iclust(my.obj, k = 150, data.type = "knetl") 

A <- cluster.plot(my.obj,plot.type = "pca",interactive = F,cell.size = 0.5,cell.transparency = 1, anno.clust=T)
B <- cluster.plot(my.obj,plot.type = "umap",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
C <- cluster.plot(my.obj,plot.type = "tsne",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
D <- cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)

library(gridExtra)

grid.arrange(A,B,C,D)

png('four.plots.grid.icluster.plot.png', width=10, height=10, units = 'cm', res=1200)
grid.arrange(plot(A),plot(B),plot(C), plot(D))
dev.off()

png('four.plots.grid.icluster.noplot.png', width=10, height=10, units = 'cm', res=1200)
grid.arrange(A,B,C,D)
dev.off()

my.obj <- clust.ord(my.obj,top.rank = 500, how.to.order = "distance")
my.obj <- clust.ord(my.obj,top.rank = 500, how.to.order = "random")

A= cluster.plot(my.obj,plot.type = "pca",interactive = F,cell.size = 0.5,cell.transparency = 1, anno.clust=T)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)

png('four.plots.grid.icluster.ordination.plot.png', width=10, height=10, units = 'cm', res=1200)
grid.arrange(plot(A),plot(B),plot(C), plot(D))
dev.off()

png('four.plots.grid.icluster.ordination.noplot.png', width=10, height=10, units = 'cm', res=1200)
grid.arrange(A,B,C,D)
dev.off()

# conditions 
A <- cluster.plot(my.obj,plot.type = "pca",col.by = "conditions",interactive = F,cell.size = 0.5)
B <- cluster.plot(my.obj,plot.type = "umap",col.by = "conditions",interactive = F,cell.size = 0.5)
C <- cluster.plot(my.obj,plot.type = "tsne",col.by = "conditions",interactive = F,cell.size = 0.5)
D <- cluster.plot(my.obj,plot.type = "knetl",col.by = "conditions",interactive = F,cell.size = 0.5)

plot(A)
plot(B)
plot(C)
plot(D)

library(gridExtra)
png('four.plots.grid.icluster.ordination.conditions.plot.png', width=30, height=30, units = 'cm', res=1200)
grid.arrange(plot(A),plot(B),plot(C), plot(D))
dev.off()

png('four.plots.grid.icluster.ordination.conditions.noplot.png', width=30, height=30, units = 'cm', res=1200)
grid.arrange(A,B,C,D)
dev.off()

### or 

png('cluster.plot.knetl.png', width=30, height=30, units = 'cm', res=1200)
cluster.plot(my.obj,
             cell.size = 0.5,
             plot.type = "knetl",
             cell.color = "black",
             back.col = "white",
             col.by = "conditions",
             cell.transparency = 1,
             clust.dim = 2,
             interactive = F,cond.facet = T)
dev.off()
             
save.image(file = "practice.Rdata")
load(file = "practice.Rdata")

## with memberships 
png('pseudotime.knetl.png', width=30, height=30, units = 'cm', res=1200)
pseudotime.knetl(my.obj,interactive = F,cluster.membership = F,conds.to.plot = NULL)
dev.off()

## with memberships 
png('pseudotime.knetl.membership.png', width=30, height=30, units = 'cm', res=1200)
pseudotime.knetl(my.obj,interactive = F,cluster.membership = T,conds.to.plot = NULL)
dev.off()

# intractive plot
# pseudotime.knetl(my.obj,interactive = T)

# Average expression per cluster 

# for all conditions
my.obj <- clust.avg.exp(my.obj, conds.to.avg = NULL)

# for one condition
my.obj <- clust.avg.exp(my.obj, conds.to.avg = "single.MOC")

# for two condition
my.obj <- clust.avg.exp(my.obj, conds.to.avg = c("single.MOC","single.MOC.extra"))

# for three condition
my.obj <- clust.avg.exp(my.obj, conds.to.avg = c("single.MOC","single.MOC.extra", "chronic.MOC"))

# for four condition
my.obj <- clust.avg.exp(my.obj, conds.to.avg = c("single.MOC","single.MOC.extra", "chronic.MOC", "PBS"))

head(my.obj@clust.avg)


# Cell cycle prediction 

# old method 
# my.obj <- cc(my.obj, s.genes = s.phase, g2m.genes = g2m.phase)

# new method 

G0 <- readLines(system.file('extdata', 'G0.txt', package = 'iCellR'))
G1S <- readLines(system.file('extdata', 'G1S.txt', package = 'iCellR'))
G2M <- readLines(system.file('extdata', 'G2M.txt', package = 'iCellR'))
M <- readLines(system.file('extdata', 'M.txt', package = 'iCellR'))
MG1 <- readLines(system.file('extdata', 'MG1.txt', package = 'iCellR'))
S <- readLines(system.file('extdata', 'S.txt', package = 'iCellR'))

# Tirosh scoring method (recomanded)
my.obj <- cell.cycle(my.obj, scoring.List = c("G0","G1S","G2M","M","MG1","S"), scoring.method = "tirosh")

# Coverage scoring method (recomanded)
# my.obj <- cell.cycle(my.obj, scoring.List = c("G0","G1S","G2M","M","MG1","S"), scoring.method = "coverage")

# plot cell cycle

A= cluster.plot(my.obj,plot.type = "pca",interactive = F,cell.size = 0.5,cell.transparency = 1, anno.clust=T,col.by = "cc")
B= cluster.plot(my.obj,plot.type = "umap",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T, col.by = "cc")
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T, col.by = "cc")
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T, col.by = "cc")

library(gridExtra)
png('cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))
dev.off()

## or 
png('cluster.plot.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
cluster.plot(my.obj,
              cell.size = 0.5,
              plot.type = "knetl",
              col.by = "cc",
              cell.color = "black",
              back.col = "white",
              cell.transparency = 1,
              clust.dim = 2,
              interactive = F,cond.facet = T)
dev.off()

# Pie
png('pie.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "pie.cc", interactive = F, conds.to.plot = NULL)
dev.off()

# bar
png('bar.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "bar.cc", interactive = F, conds.to.plot = NULL)
dev.off()

# or per condition
png('single.MOC.pie.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "pie.cc", interactive = F, conds.to.plot = "single.MOC")
dev.off()

png('single.MOC.extra.pie.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "pie.cc", interactive = F, conds.to.plot = "single.MOC.extra")
dev.off()

png('chronic.MOC.pie.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "pie.cc", interactive = F, conds.to.plot = "chronic.MOC")
dev.off()

png('PBS.pie.cell.cycle.four.grid.plots.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "pie.cc", interactive = F, conds.to.plot = "PBS")
dev.off()

# Cell frequencies 

clust.cond.info(my.obj, plot.type = "pie", normalize.ncell = TRUE, my.out.put = "plot", normalize.by = "percentage")

clust.cond.info(my.obj, plot.type = "bar", normalize.ncell = TRUE,my.out.put = "plot", normalize.by = "percentage")

clust.cond.info(my.obj, plot.type = "pie.cond", normalize.ncell = T, my.out.put = "plot", normalize.by = "percentage")

clust.cond.info(my.obj, plot.type = "bar.cond", normalize.ncell = T,my.out.put = "plot", normalize.by = "percentage")

my.obj <- clust.cond.info(my.obj)
head(my.obj@my.freq)

png('cluster.box.mitos.png', width=30, height=30, units = 'cm', res=1200)
clust.stats.plot(my.obj, plot.type = "box.mito", interactive = F)
dev.off()



# Gene to Gene correlation 

# impute more cells by increasing nn for better resulst. 
my.obj <- run.impute(my.obj,dims = 1:10,data.type = "pca", nn = 50)

# main data
A <- gg.cor(my.obj, 
	interactive = F, 
	gene1 = "GNLY",
	gene2 = "NKG7", 
	conds = NULL,
	clusts = NULL,
	data.type = "main")

# imputed data 
B <- gg.cor(my.obj, 
	interactive = F, 
	gene1 = "GNLY",
	gene2 = "NKG7", 
	conds = NULL,
	clusts = NULL,
	data.type = "imputed")

C <- gg.cor(my.obj, 
	interactive = F, 
	gene1 = "GNLY",
	gene2 = "NKG7", 
	conds = NULL,
	clusts = c(3,2),
	data.type = "imputed")


# imputed data 
D <- gg.cor(my.obj, 
	interactive = F, 
	gene1 = "GNLY",
	gene2 = "NKG7", 
	conds = c("PBS"),
	clusts = NULL,
	data.type = "imputed")

png('gene.to.gene.PBS.png', width=30, height=30, units = 'cm', res=1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))
dev.off()


# Find marker genes 

marker.genes <- findMarkers(my.obj,
	fold.change = 2,
	padjval = 0.1)

dim(marker.genes)

head(marker.genes)

# Heat map 

# find top genes
MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.2,filt.ambig = F)
MyGenes <- unique(MyGenes)

# main data 
png('heatmap.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = NULL)
dev.off()

# imputed data 
png('heatmap.imputed.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", conds.to.plot = NULL)
dev.off()

# sort cells and plot only one condition
png('heatmap.imputed.PBS.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", cell.sort = TRUE, conds.to.plot = c("PBS"))
dev.off()

png('heatmap.imputed.chronic.MOC.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", cell.sort = TRUE, conds.to.plot = c("chronic.MOC"))
dev.off()

png('heatmap.imputed.single.MOC.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", cell.sort = TRUE, conds.to.plot = c("single.MOC"))
dev.off()

png('heatmap.imputed.single.MOC.extra.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", cell.sort = TRUE, conds.to.plot = c("single.MOC.extra"))
dev.off()

# Pseudotime stile
png('heatmap.imputed.pseudotime.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "none", data.type = "imputed", cell.sort = TRUE)
dev.off()

# intractive 
# heatmap.gg.plot(my.obj, gene = MyGenes, interactive = T, out.name = "heatmap_gg", cluster.by = # "clusters")














             