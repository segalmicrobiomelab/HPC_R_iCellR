#install iCellR 
#run through tutorial 

sessionInfo()

setwd("/gpfs/scratch/wub02/projects/msq.singlecell.chronicmurine.true/")
getwd()

install.packages("iCellR",lib="~/R/lib", repos='http://cran.us.r-project.org')
install.packages("gridExtra",lib="~/R/lib", repos='http://cran.us.r-project.org')
install.packages("cowplot",lib="~/R/lib", repos='http://cran.us.r-project.org')
install.packages("BiocManager",lib="~/R/lib", repos='http://cran.us.r-project.org')
install.packages("DDRTree",lib="~/R/lib", repos='http://cran.us.r-project.org')
install.packages("shiny",lib="~/R/lib", repos='http://cran.us.r-project.org')

library("BiocManager",lib="~/R/lib")
library("DDRTree",lib="~/R/lib")  
  
BiocManager::install("monocle",lib="~/R/lib")

library("iCellR",lib="~/R/lib")
library("gridExtra",lib="~/R/lib")
library("cowplot",lib="~/R/lib")
library("shiny",lib="~/R/lib")
library("monocle",lib="~/R/lib")

load(file = "practice.Rdata")

# Impute data/impute more cells by increasing nn for better resulst. 
my.obj <- run.impute(my.obj,dims = 1:10,data.type = "pca", nn = 50)

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
# my.obj <- run.impute(my.obj,dims = 1:10,data.type = "pca", nn = 50)

# main data
# A <- gg.cor(my.obj, 
# 	interactive = F, 
# 	gene1 = "GNLY",
# 	gene2 = "NKG7", 
# 	conds = NULL,
# 	clusts = NULL,
# 	data.type = "main")

# imputed data 
# B <- gg.cor(my.obj, 
# 	interactive = F, 
# 	gene1 = "GNLY",
# 	gene2 = "NKG7", 
# 	conds = NULL,
# 	clusts = NULL,
# 	data.type = "imputed")

# C <- gg.cor(my.obj, 
# 	interactive = F, 
# 	gene1 = "GNLY",
# 	gene2 = "NKG7", 
# 	conds = NULL,# 
# 	clusts = c(3,2),
# 	data.type = "imputed")


# imputed data 
# D <- gg.cor(my.obj, 
# 	interactive = F, 
# 	gene1 = "GNLY",
# 	gene2 = "NKG7", 
# 	conds = c("PBS"),
# 	clusts = NULL,
# 	data.type = "imputed")

# png('gene.to.gene.PBS.png', width=30, height=30, units = 'cm', res=1200)
# grid.arrange(plot(A),plot(B),plot(C),plot(D))
# dev.off()


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

write.table(MyGenes, file = "Mygenes.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))
            
# find top genes cluster 1
MyGenes.one <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 1)
MyGenes.one <- unique(MyGenes.one)

write.table(MyGenes.one, file = "MyGenes.one.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))
            
# find top genes cluster 1
MyGenes.one <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 1)
MyGenes.one <- unique(MyGenes.one)

write.table(MyGenes.one, file = "MyGenes.one.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))                           

# find top genes cluster 2
MyGenes.two <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 2)
MyGenes.two <- unique(MyGenes.two)

write.table(MyGenes.two, file = "MyGenes.two.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))
            
# find top genes cluster 3
MyGenes.three <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 3)
MyGenes.three <- unique(MyGenes.three)

write.table(MyGenes.three, file = "MyGenes.three.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))                  

# find top genes cluster 4
MyGenes.four <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 4)
MyGenes.four <- unique(MyGenes.four)

write.table(MyGenes.four, file = "MyGenes.four.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# find top genes cluster 5
MyGenes.five <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 4)
MyGenes.five <- unique(MyGenes.five)

write.table(MyGenes.five, file = "MyGenes.five.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double")) 

# find top genes cluster 5
MyGenes.five <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 5)
MyGenes.five <- unique(MyGenes.five)

write.table(MyGenes.five, file = "MyGenes.five.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# find top genes cluster 6
MyGenes.six <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 6)
MyGenes.six <- unique(MyGenes.six)

write.table(MyGenes.six, file = "MyGenes.six.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# find top genes cluster 7
MyGenes.seven <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 7)
MyGenes.seven <- unique(MyGenes.seven)

write.table(MyGenes.seven, file = "MyGenes.seven.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))     
     
# find top genes cluster 8
MyGenes.eight <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 8)
MyGenes.eight <- unique(MyGenes.eight)

write.table(MyGenes.eight, file = "MyGenes.eight.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))  

# find top genes cluster 9
MyGenes.nine <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 9)
MyGenes.nine <- unique(MyGenes.nine)

write.table(MyGenes.nine, file = "MyGenes.nine.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# find top genes cluster 10
MyGenes.ten <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 10)
MyGenes.ten <- unique(MyGenes.ten)

write.table(MyGenes.ten, file = "MyGenes.ten.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# find top genes cluster 11
MyGenes.eleven <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 11)
MyGenes.eleven <- unique(MyGenes.eleven)

write.table(MyGenes.eleven, file = "MyGenes.eleven.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# find top genes cluster 12
MyGenes.twelve <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 12)
MyGenes.twelve <- unique(MyGenes.twelve)

write.table(MyGenes.twelve, file = "MyGenes.twelve.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double")) 

# find top genes cluster 13
MyGenes.thirteen <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 13)
MyGenes.thirteen <- unique(MyGenes.thirteen)

write.table(MyGenes.thirteen, file = "MyGenes.thirteen.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double")) 

# find top genes cluster 14
MyGenes.fourteen <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2,filt.ambig = F, cluster = 14)
MyGenes.twelve <- unique(MyGenes.twelve)

write.table(MyGenes.fourteen, file = "MyGenes.fourteen.txt", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))

# main data 
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = NULL)

png('heatmap.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = NULL)
dev.off()

# imputed data 
# heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", conds.to.plot = NULL)

# png('heatmap.imputed.png', width=30, height=30, units = 'cm', res=1200)
# heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", conds.to.plot = NULL)
# dev.off()

# sort cells and plot only one condition, not imputed
png('heatmap.PBS.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = c("PBS"))
dev.off()

png('heatmap.imputed.chronic.MOC.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = c("chronic.MOC"))
dev.off()

png('heatmap.imputed.single.MOC.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = c("single.MOC"))
dev.off()

png('heatmap.imputed.single.MOC.extra.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = c("single.MOC.extra"))
dev.off()

# Pseudotime stile
png('heatmap.imputed.pseudotime.png', width=30, height=30, units = 'cm', res=1200)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "none")
dev.off()

# intractive 
# heatmap.gg.plot(my.obj, gene = MyGenes, interactive = T, out.name = "heatmap_gg", cluster.by = # "clusters")


# Plot Genes

A <- gene.plot(my.obj, gene = "Pdcd1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Pdcd1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Pdcd1", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Pdcd1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.Pdcd1.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Pf4", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Pf4", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Pf4", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Pf4", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster2.Pf4.platelets.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Hbb.bt", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Hbb.bt", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Hbb.bt", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Hbb.bt", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster12.Hbb.bt.Redcells.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Hey1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Hey1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Hey1", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Hey1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster3.Hey1.endothelial.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Dcn", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Dcn", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Dcn", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Dcn", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster13.Dcn.fibroblasts.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Ms4a6c", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Ms4a6c", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Ms4a6c", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Ms4a6c", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster10.Ms4a6c.macrophages.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Tcrg.C1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Tcrg.C1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Tcrg.C1", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Tcrg.C1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster6.Tcrg.C1.Tcells.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Cd79a", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Cd79a", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Cd79a", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Cd79a", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster7.Cd79a.Bcells.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Cyp4b1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Cyp4b1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Cyp4b1", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Cyp4b1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster1.Cyp4b1.clubcells.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Mfap4", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Mfap4", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Mfap4", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Mfap4", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster9.Mfap4.monocytes.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Gzma", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Gzma", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Gzma", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Gzma", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster14.Gzma.CD8Tcells.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Txk", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Txk", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Txk", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Txk", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster4.Txk.NKcells.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Mmp9", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Mmp9", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Mmp9", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Mmp9", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster5.Mmp9.classicmonocyte.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "C1qa", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "C1qa", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "C1qa", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "C1qa", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster11.C1qa.nonclassicmonocyte.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Plot Genes

A <- gene.plot(my.obj, gene = "Socs3", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Socs3", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Socs3", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Socs3", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.cluster5.Socs3.granulocyte.neutrophil.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

### same on imputed data 

A <- gene.plot(my.obj, gene = "Pdcd1", 
	plot.type = "scatterplot",
	interactive = F,
	data.type = "main",
	out.name = "scatter_plot")
	
# PCA 2D	
B <- gene.plot(my.obj, gene = "Pdcd1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	data.type = "main",
	plot.data.type = "umap")
	
# Box Plot
C <- gene.plot(my.obj, gene = "Pdcd1", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	data.type = "main",
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "Pdcd1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	data.type = "main",
	out.name = "bar_plot")
	
library(gridExtra)
png('gene.plots.Pdcd1_imputed.png', width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))	
dev.off()

# Multiple plots 
# 
# genelist = c("MS4A1","GNLY","FCGR3A","NKG7","CD14","CD3E","CD8A","CD4","GZMH","CCR7","CD68")
# 
# rm(list = ls(pattern="PL_"))
# for(i in genelist){
####
#     MyPlot <- gene.plot(my.obj, gene = i,
#         interactive = F,
#         cell.size = 0.1,
#         plot.data.type = "knetl",
#         data.type = "main",
#         scaleValue = T,
#         min.scale = 0,max.scale = 2.0,
#         cell.transparency = 1)
####
#     NameCol=paste("PL",i,sep="_")
#     eval(call("<-", as.name(NameCol), MyPlot))
# }

# library(cowplot)
# filenames <- ls(pattern="PL_")

# B <- cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.1,cell.transparency = 1,anno.clust=T)
# filenames <- c("B",filenames)

# png('genes_KNetL.png',width = 30, height = 24, units = 'cm', res = 1200)
# plot_grid(plotlist=mget(filenames))
# dev.off()

# or heatmap 
# png('genes_KNetL.heatmap.png',width = 30, height = 24, units = 'cm', res = 1200)
# heatmap.gg.plot(my.obj, gene = genelist, interactive = F, cluster.by = "clusters")
# dev.off()


# Customized plots 

gene.plot(my.obj, gene = "Pdcd1", write.data = T, scaleValue = F, data.type = "main")

# This would create a text file called "MS4A1.tsv".
 head(read.table("Pdcd1.tsv"))
 
 # Annotating clusters 
 
###### Labeling the clusters 
# CD3E: only in T Cells
# FCGR3A (CD16): in CD16+ monocytes and some expression NK cells
# GNLY: NK cells
# MS4A1: B cells
# GZMH: in GZMH+ T8 cells and some expression NK cells
# CD8A: in T8 cells
# CD4: in T4 and some myeloid cells
# CCR7: expressed more in memory cells 
# CD14: in CD14+ monocytes
# CD68: in monocytes/MF

my.obj <- change.clust(my.obj, change.clust = 1, to.clust = "001.MG")
my.obj <- change.clust(my.obj, change.clust = 2, to.clust = "002.NK")
my.obj <- change.clust(my.obj, change.clust = 3, to.clust = "003.CD16+.Mono")
my.obj <- change.clust(my.obj, change.clust = 4, to.clust = "004.MF")
my.obj <- change.clust(my.obj, change.clust = 5, to.clust = "005.CD14+.Mono")
my.obj <- change.clust(my.obj, change.clust = 6, to.clust = "006.Naive.T8")
my.obj <- change.clust(my.obj, change.clust = 7, to.clust = "007.GZMH+.T8")
my.obj <- change.clust(my.obj, change.clust = 8, to.clust = "008.B")
my.obj <- change.clust(my.obj, change.clust = 9, to.clust = "009.Memory.T4")
my.obj <- change.clust(my.obj, change.clust = 10, to.clust = "010.Naive.T4")

A= cluster.plot(my.obj,plot.type = "pca",interactive = F,cell.size = 0.5,cell.transparency = 1, anno.clust=T)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)

png('labeling.cluster.plots.png',width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))
dev.off() 


save.image(file = "practice.individual.genes.Rdata")
load(file = "practice.individual.genes.Rdata")


A <- gene.plot(my.obj, gene = "Pdcd1", 
   plot.type = "scatterplot",
   interactive = F,
   cell.transparency = 1,
   scaleValue = TRUE,
   min.scale = 0,
   max.scale = 2.5,
   back.col = "white",
   cond.shape = TRUE)
   
B <- gene.plot(my.obj, gene = "Pdcd1", 
   plot.type = "scatterplot",
   interactive = F,
   cell.transparency = 1,
   scaleValue = TRUE,
   min.scale = 0,
   max.scale = 2.5,
   back.col = "white",
   cond.shape = TRUE,
   conds.to.plot = c("PBS","chronic.MOC"))

C <- gene.plot(my.obj, gene = "Pdcd1", 
   plot.type = "boxplot",
   interactive = F,
   back.col = "white",
   cond.shape = TRUE,
   conds.to.plot = c("PBS"))

D <- gene.plot(my.obj, gene = "Pdcd1", 
   plot.type = "barplot",
   interactive = F,
   cell.transparency = 1,
   back.col = "white",
   cond.shape = TRUE,
   conds.to.plot = c("PBS","chronic.MOC"))

png('pdcd1.plot.four.box.png',width = 30, height = 30, units = 'cm', res = 1200)
grid.arrange(plot(A),plot(B),plot(C),plot(D))
dev.off() 


# Interactive plot

# example
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "umap",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	cond.shape = T,
	interactive = T,
	out.name = "2d_UMAP_clusters_conds")

# 2D
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "tsne",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = F)
	
# interactive 2D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 2,
	interactive = T,
	out.name = "tSNE_2D_clusters")

# interactive 3D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 3,
	interactive = T,
	out.name = "tSNE_3D_clusters")

# Density plot for clusters
png('diffusion.clusters.clusters.png',width = 30, height = 30, units = 'cm', res = 1200)		 
cluster.plot(my.obj,
	plot.type = "pca",
	col.by = "clusters",
	interactive = F,
	density=T)
dev.off()

# Density plot for conditions 
png('diffusion.clusters.conditions.png',width = 30, height = 30, units = 'cm', res = 1200)		 
cluster.plot(my.obj,
	plot.type = "pca",
	col.by = "conditions",
	interactive = F,
	density=T)
dev.off()

# png('diffusion.clusters.cluster2.png',width = 30, height = 30, units = 'cm', res = 1200)		
# cluster.plot(my.obj,
#	cell.size = 1,
#	plot.type = "diffusion",
#	cell.color = "black",
#	back.col = "white",
#	col.by = "clusters",
#	cell.transparency = 0.5,
#	clust.dim = 2,
#	interactive = F)
# dev.off()
	
# png('diffusion.clusters.cluster3.png',width = 30, height = 30, units = 'cm', res = 1200)	
# cluster.plot(my.obj,
#	cell.size = 1,
# 	plot.type = "diffusion",
#	cell.color = "black",
#	back.col = "white",
#	col.by = "clusters",
#	cell.transparency = 0.5,
#	clust.dim = 3,
#	interactive = F)
# dev.off()


# Differential Expression analysis 

# diff.res <- run.diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))
# diff.res1 <- as.data.frame(diff.res)
# diff.res1 <- subset(diff.res1, padj < 0.05)
# head(diff.res1)

# more examples 

# Comparing a condition/conditions with different condition/conditions (e.g. PBS vs chronic.MOC)
diff.res <- run.diff.exp(my.obj, de.by = "conditions", cond.1 = c("PBS"), cond.2 = c("chronic.MOC"))

# Comparing a condition/conditions with different condition/conditions (e.g. single.MOC vs chronic.MOC)
# diff.res <- run.diff.exp(my.obj, de.by = "conditions", cond.1 = c("single.MOC"), cond.2 = c("chronic.MOC"))

# Comparing a condition/conditions with different condition/conditions (e.g. single.MOC vs chronic.MOC)
# diff.res <- run.diff.exp(my.obj, de.by = "conditions", cond.1 = c("single.MOC.extra"), cond.2 = c("chronic.MOC"))

# Comparing a cluster/clusters with different cluster/clusters (e.g. cluster 1 and 2 vs. 4)
# diff.res <- run.diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))

# Comparing a condition/conditions with different condition/conditions only in one/more cluster/clusters (e.g. cluster 1 WT vs cluster 1 KO)
# diff.res <- run.diff.exp(my.obj, de.by = "clustBase.condComp", cond.1 = c("PBS"), cond.2 = c("chronic.MOC"), base.cond = 1)

# Comparing a cluster/clusters with different cluster/clusters only in one/more condition/conditions (e.g. cluster 1 vs cluster 2 but only the WT sample)
# diff.res <- run.diff.exp(my.obj, de.by = "condBase.clustComp", cond.1 = c(1), cond.2 = c(2), base.cond = "chronic.MOC")


# volcano / MA 

# Volcano Plot 
png('volcano.differential.png',width = 30, height = 30, units = 'cm', res = 1200)
volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "volcano",
	interactive = F)
	dev.off() 

# MA Plot
png('MA.differential.png',width = 30, height = 30, units = 'cm', res = 1200)
volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "ma",
	interactive = F)
	dev.off() 


# Merging, resetting, renaming, and removing clusters

# let's say you  want to merge cluster 3 and 2.
# my.obj <- change.clust(my.obj, change.clust = 3, to.clust = 2)

# to reset to the original clusters run this.
# my.obj <- change.clust(my.obj, clust.reset = T)

# you can also re-name the cluster numbers to cell types. Remember to reset after this so you can ran other analysis. 
# my.obj <- change.clust(my.obj, change.clust = 7, to.clust = "B Cell")

# Let's say for what ever reason you want to remove acluster, to do so run this.
# my.obj <- clust.rm(my.obj, clust.to.rm = 1)

# Remember that this would perminantly remove the data from all the slots in the object except frrom raw.data slot in the object. If you want to reset you need to start from the filtering cells step in the biginging of the analysis (using cell.filter function). 

# To re-position the cells run tSNE again 
# my.obj <- run.tsne(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")

# Use this for plotting as you make the changes

# png('change.cluster.test.plot.png',width = 30, height = 30, units = 'cm', res = 1200)
# cluster.plot(my.obj,
#   cell.size = 1,
#   plot.type = "tsne",
#   cell.color = "black",
#   back.col = "white",
#   col.by = "clusters",
#   cell.transparency = 0.5,
#   clust.dim = 2,
#   interactive = F)
# dev.off() 

save.image(file = "practice.individual.genes.before.gate.Rdata")
load(file = "practice.individual.genes.before.gate.Rdata")


# Cell Gating

# my.plot <- gene.plot(my.obj, gene = "Pdcd1", 
#   plot.type = "scatterplot",
#   clust.dim = 2,
#   interactive = F)

# png('pdcd1.cell.gate.png',width = 30, height = 30, units = 'cm', res = 1200)
# my.plot <- gene.plot(my.obj, gene = "Pdcd1", 
#   plot.type = "scatterplot",
#   clust.dim = 2,
#   interactive = F)
# dev.off() 

# cell.gating(my.obj, my.plot = my.plot, plot.type = "tsne")	


# Pseudotime analysis 

MyGenes <- top.markers(marker.genes, topde = 50, min.base.mean = 0.2)
MyGenes <- unique(MyGenes)

pseudotime.tree(my.obj,
   marker.genes = MyGenes,
   type = "unrooted",
   clust.method = "complete")

png('pseudotime.mygenes.png',width = 30, height = 30, units = 'cm', res = 1200)   
 pseudotime.tree(my.obj,
   marker.genes = MyGenes,
   type = "unrooted",
   clust.method = "complete") 
dev.off()    

# or 

# pseudotime.tree(my.obj,
#   marker.genes = MyGenes,
#   type = "classic",
#   clust.method = "complete")
   
# pseudotime.tree(my.obj,
#   marker.genes = MyGenes,
#   type = "jitter",
#   clust.method = "complete")	


# pseudotime analysis monocle 

library("monocle",lib="~/R/lib")

MyMTX <- my.obj@main.data
GeneAnno <- as.data.frame(row.names(MyMTX))
colnames(GeneAnno) <- "gene_short_name"
row.names(GeneAnno) <- GeneAnno$gene_short_name
cell.cluster <- (my.obj@best.clust)
Ha <- data.frame(do.call('rbind', strsplit(as.character(row.names(cell.cluster)),'_',fixed=TRUE)))[1]
clusts <- paste("cl.",as.character(cell.cluster$clusters),sep="")
cell.cluster <- cbind(cell.cluster,Ha,clusts)
colnames(cell.cluster) <- c("Clusts","iCellR.Conds","iCellR.Clusts")
Samp <- new("AnnotatedDataFrame", data = cell.cluster)
Anno <- new("AnnotatedDataFrame", data = GeneAnno)
my.monoc.obj <- newCellDataSet(as.matrix(MyMTX),phenoData = Samp, featureData = Anno)

## find disperesedgenes 
my.monoc.obj <- estimateSizeFactors(my.monoc.obj)
my.monoc.obj <- estimateDispersions(my.monoc.obj)
disp_table <- dispersionTable(my.monoc.obj)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my.monoc.obj <- setOrderingFilter(my.monoc.obj, unsup_clustering_genes$gene_id)

# tSNE
my.monoc.obj <- reduceDimension(my.monoc.obj, max_components = 2, num_dim = 10,reduction_method = 'tSNE', verbose = T)
# cluster 
my.monoc.obj <- clusterCells(my.monoc.obj, num_clusters = 10)

## plot conditions and clusters based on iCellR analysis 
A <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "iCellR.Conds")
B <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "iCellR.Clusts")

## plot clusters based monocle analysis 
C <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "Cluster")

# get marker genes from iCellR analysis
MyGenes <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2)
my.monoc.obj <- setOrderingFilter(my.monoc.obj, MyGenes)

my.monoc.obj <- reduceDimension(my.monoc.obj, max_components = 2,method = 'DDRTree')
# order cells 
my.monoc.obj <- orderCells(my.monoc.obj)

# plot based on iCellR analysis and marker genes from iCellR
D <- plot_cell_trajectory(my.monoc.obj, color_by = "iCellR.Clusts")

## heatmap genes from iCellR

plot_pseudotime_heatmap(my.monoc.obj[MyGenes,],
  cores = 1,
  cluster_rows = F,
  use_gene_short_name = T,
  show_rownames = T)


png('pseudotime.heatmap.png',width = 30, height = 30, units = 'cm', res = 1200)   
plot_pseudotime_heatmap(my.monoc.obj[MyGenes,],
  cores = 1,
  cluster_rows = F,
  use_gene_short_name = T,
  show_rownames = T)
dev.off()    

save.image(file = "practice.individual.genes.end.Rdata")
load(file = "practice.individual.genes.end.Rdata")










             