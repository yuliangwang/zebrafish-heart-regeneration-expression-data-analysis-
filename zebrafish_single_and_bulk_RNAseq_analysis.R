#################################################
########Analysis code for scRNA-seq##############
#scRNA-seq data include day 0 (pre-injury), day 3, and 7 post injury
setwd("./zebrafish_sc/")
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)
a1.data<-Read10X(data.dir ="./zf_0d/")
a1<-CreateSeuratObject(raw.data = a1.data, project = "day 0")
projects<-c("_3d","_7d")
project_names<-c("day 3","day 7")
for (i in seq_along(projects)){
  project_data<-Read10X(data.dir=paste("./zf",projects[i],sep=""))
  a1 <- AddSamples(object = a1, new.data = project_data, add.cell.id = project_names[i])
}

cm<-a1
rm(a1)


mito.genes <- grep(pattern = "^MT-", x = rownames(x = cm@data), value = TRUE)
percent.mito <- Matrix::colSums(cm@raw.data[mito.genes, ])/Matrix::colSums(cm@raw.data)
cm <- AddMetaData(object = cm, metadata = percent.mito, col.name = "percent.mito")

cm <- FilterCells(object = cm, subset.names = c("nGene", "percent.mito","nUMI"), 
                  low.thresholds = c(200, -Inf,2000), high.thresholds = c(Inf, 0.4,Inf))
                  
cm <- NormalizeData(object = cm, normalization.method = "LogNormalize", 
                    scale.factor = 10000)


cm <- FindVariableGenes(object = cm, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

cm <- ScaleData(object = cm, vars.to.regress = c("nUMI", "percent.mito"))

cm <- RunPCA(object = cm, pc.genes = cm@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5)

VizPCA(object = cm, pcs.use = 1:2)

PCAPlot(object = cm, dim.1 = 1, dim.2 = 2)


PCHeatmap(object = cm, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)



PCElbowPlot(object = cm)

cm <- FindClusters(object = cm, reduction.type = "pca", dims.use = 1:10, 
                   resolution = 0.3, print.output = 0, save.SNN = TRUE,force.recalc=T)

cm <- RunTSNE(object = cm, dims.use = 1:10, do.fast = TRUE)

##Clustering
pdf("tSNE_plot_zebrafish.pdf")
TSNEPlot(object = cm,do.label=T,label.size = 10,pt.size = 2)
dev.off()


pdf("tSNE_plot_zebrafish_coded_by_date.pdf")
TSNEPlot(object = cm,do.label=F,label.size = 10,pt.size = 2,group.by = "orig.ident")
dev.off()

#Plot marker genes
#cell type markers#
markers<-c('ENC1', 'DCN', 'COL5A1','THY1')
markers<-tolower(markers)
pdf("FB_marker_genes_zebrafish.pdf",height=3,width=15)
FeaturePlot(object = cm, features.plot = markers, cols.use = c("grey", "red"), reduction.use = "tsne",nCol = 5)
dev.off()

markers<-c('MYL7', 'GATA4',"MEF2AA")
markers<-tolower(markers)
pdf("cardiac_marker_genes_zebrafish.pdf",height=3,width=12)
FeaturePlot(object = cm, features.plot = markers, cols.use = c("grey", "red"), reduction.use = "tsne",nCol = 4)
dev.off()

#Find marker genes for each cluster
maturation.markers <- FindAllMarkers(object = cm, only.pos = TRUE, min.pct = 0.25, 
                                     thresh.use = 0.25)

maturation.markers<-maturation.markers[-grep("si:",maturation.markers$gene),]
top10 <- maturation.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf("Top_positive_marker_genes_by_cluster_zf.pdf",width=5,height=8)
DoHeatmap(object = cm,genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
dev.off()

#Cell cycle genes#
cc.genes <- readLines(con = "./regev_lab_cell_cycle_genes.txt")
s.genes <- tolower(cc.genes[1:43])
g2m.genes <- tolower(cc.genes[44:97])


pdf("S_phase_genes_zebrafish.pdf",width=5,height=6.5)
DoHeatmap(object = cm,genes.use = tolower(s.genes), slim.col.label = TRUE, remove.key = FALSE,use.scaled=TRUE)
dev.off()

pdf("G2M_phase_genes_zebrafish.pdf",width=5,height=6.5)
DoHeatmap(object = cm,genes.use = tolower(g2m.genes), slim.col.label = TRUE, remove.key = FALSE,use.scaled=TRUE)
dev.off()

########Trajectory analysis using Monocle#########

ind<-match(colnames(cm@data),colnames(cm@raw.data))
count_matrix<-as.matrix(cm@raw.data[,ind])
pd<-new("AnnotatedDataFrame", data = cm@meta.data)
gene_info<-data.frame(gene_short_name=rownames(count_matrix))
rownames(gene_info)<-rownames(count_matrix)
fd<-new("AnnotatedDataFrame", data = gene_info)

save(fd,pd,file="zebrafish_sc_data_for_monocle.RData")
load("zebrafish_sc_data_for_monocle.RData")
CM <- newCellDataSet(as(count_matrix, "sparseMatrix"),
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

CM <- estimateSizeFactors(CM)
CM <- estimateDispersions(CM)

CM <- detectGenes(CM, min_expr = 0.1)
print(head(fData(CM)))

expressed_genes <- row.names(subset(fData(CM),
                                    num_cells_expressed >= 10))

disp_table <- dispersionTable(CM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
CM <- setOrderingFilter(CM, unsup_clustering_genes$gene_id)
plot_ordering_genes(CM)

plot_pc_variance_explained(CM, return_all = F) # norm_method='log'

CM <- reduceDimension(CM, max_components = 2, num_dim = 10,
                      reduction_method = 'tSNE', verbose = T)
CM <- clusterCells(CM, num_clusters = 5)

pdf("monocle_clusters.pdf")
plot_cell_clusters(CM, 1, 2, color = "Cluster")
dev.off()

pdf("monocle_clusters_markers.pdf")
plot_cell_clusters(CM, 1, 2, color = "CellType",
                   markers = c("cmlc1", "myl6"))
dev.off()

#Wnt, endothelial marker etc.
pdf("monocole_cell_type_markers.pdf")
plot_cell_clusters(CM, 1, 2, color = "CellType",
                   markers = c("cmlc1", "myl6","wnt11r","thy1","col1a2","dcn","rgs5a"))
dev.off()

pdf("monocole_apoptosis_markers.pdf")
plot_cell_clusters(CM, 1, 2, color = "CellType",
                   markers = c("casp2","bcl2l10","fosl2","hsp90b1"))
dev.off()


pdf("monocle_clusters_by_time.pdf")
plot_cell_clusters(CM, 1, 2, color = "orig.ident")
dev.off()

pdf("cluster_day_relations.pdf",height=4,width=10)
plot_cell_clusters(CM, 1, 2, color = "Cluster") +
  facet_wrap(~orig.ident)
dev.off()


#trajectory
diff_test_res <- differentialGeneTest(CM[expressed_genes,],
                                      fullModelFormulaStr = "~orig.ident")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

CM <- setOrderingFilter(CM, ordering_genes)
plot_ordering_genes(CM)

CM <- reduceDimension(CM, max_components = 2,
                      method = 'DDRTree')

CM <- orderCells(CM)

pdf("trajectory_by_day.pdf")
plot_cell_trajectory(CM, color_by = "orig.ident")
dev.off()

pdf("trajectory_by_cluster.pdf")
plot_cell_trajectory(CM, color_by = "Cluster")
dev.off()

pdf("trajectory_by_state.pdf")
plot_cell_trajectory(CM, color_by = "State")
dev.off()

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$orig.ident)[,"day 0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

CM <- orderCells(CM, root_state = GM_state(CM))
save(CM,file="monocle_CM_cluster_and_trajectory.RData")

pdf("pseudo_time_plot.pdf")
plot_cell_trajectory(CM, color_by = "Pseudotime")
dev.off()

CM_expressed_genes <-  row.names(subset(fData(CM),
                                        num_cells_expressed >= 10))
CM_filtered <- CM[CM_expressed_genes,]
my_genes <- row.names(subset(fData(CM_filtered),
                             gene_short_name %in% c("cmlc1", "myl7", "myl6","myl9a")))
cds_subset <- CM_filtered[my_genes,]

pdf("marker_genes_pseudotime.pdf")
plot_genes_in_pseudotime(cds_subset, color_by = "orig.ident")
dev.off()

pdf("marker_genes_pseudotime_by_cluster.pdf")
plot_genes_in_pseudotime(cds_subset, color_by = "Cluster")
dev.off()
#####################################################
#####################################################


#########Analysis code for bulk RNA-seq##############
#Bulk RNA-seq samples include young day 3, adult control, and 3 and 7 days post MTZ-injury#
library(DESeq)

#Fig.S3B
load("./zf_count_mat.RData")
ind2rm<-grep("^ERR",colnames(count_mat))
count_mat<-count_mat[,-ind2rm]
cds<-newCountDataSet(count_mat[,c(1,2,5:10)],factor(rep(c("day3","adult-ours-ctrl","adult-Bednarek"),times=c(2,2,4)),levels=c("day3","adult-ours-ctrl","adult-Bednarek"),ordered=T))
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
rs<-rowSums(counts(cds,normalized=T))
cds<-cds[rs>10,]

vst_mat<-getVarianceStabilizedData(cds)
#vst_mat<-counts(cds,normalized=T)
#vst_mat<-log2(vst_mat+0.01)
zf<- t(scale(t(vst_mat),scale=F,center=T))
set.seed(10)
pca<-prcomp(t(zf))
sd<-pca$sdev
loadings<-pca$rotation
var<-sd^2
var.percent<-var/sum(var)*100
labels<-rep(c("day3","adult-ours-ctrl","adult-Bednarek"),times=c(2,2,4))
scores <- data.frame(sample.groups=labels, pca$x[,1:3],labs=rep(c("ours","Bednarek"),times=c(4,4)))
colnames(scores)[2:4]<-c("PC1","PC2","PC3")
#scores$PC1<- -(scores$PC1)
#scores$PC2<- -(scores$PC2)
min_axis<-min(scores[,2:4])
max_axis<-max(scores[,2:4])
library(grid)
cbPalette<-c("azure4","black","blue","brown","cadetblue","chartreuse","cyan",
             "darkorange","darkorchid","deeppink","gold","lightcoral","lightseagreen","magenta","red","lightsalmon","yellow","mediumorchid4","deepskyblue","mediumvioletred","olivedrab","cornsilk","lavender","navajowhite4")

#pie(rep(1,length(cbPalette)),col=cbPalette)
library(ggplot2)
p1=ggplot(data=scores, aes(x=PC1,y=PC2,color=sample.groups))+ 
  #geom_point(aes(shape=factor(labs)),size=0.7) +
  geom_point(size=0.7) +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.15,"cm"),legend.justification=c(0,1),legend.margin=unit(-0.1,"cm"),legend.position="top",legend.text=element_text(size=5),axis.text.x = element_text(size=6),axis.text.y=element_text(size=6),axis.title=element_text(size=6),legend.title=element_blank()) + guides(color=guide_legend(ncol=4)) +
  xlim(min_axis,max_axis) + ylim(-52,max_axis+2) + scale_color_manual(values=cbPalette[c(2,15,8,3,4)]) + xlab("PC 1 (75%)") + ylab("PC2 (13%)") 
png("./zebrafish_day3_adult_PCA_plot_all_genes_08182016.png",res=600,width=4,height=2,units="in")
p1
dev.off()


#Fig.S3C
load("./zf_count_mat.RData")
ind2rm<-grep("^ERR",colnames(count_mat))
cds<-newCountDataSet(count_mat[,-ind2rm],factor(rep(c("day3","adult-ours-mtz","adult-ours-ctrl","adult-Bednarek","adult-Bednarek-cryoinjury"),times=c(2,2,2,4,4)),levels=c("day3","adult-ours-mtz","adult-ours-ctrl","adult-Bednarek","adult-Bednarek-cryoinjury"),ordered=T))
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
rs<-rowSums(counts(cds,normalized=T))
cds<-cds[rs>10,]

vst_mat<-getVarianceStabilizedData(cds)


zf<- t(scale(t(vst_mat),scale=F,center=T))
set.seed(10)
pca<-prcomp(t(zf))
#pca<-princomp(zf)
sd<-pca$sdev
loadings<-pca$rotation
var<-sd^2
var.percent<-var/sum(var)*100
labels<-rep(c("day3","adult-ours-mtz","adult-ours-ctrl","adult-Bednarek","adult-Bednarek-cryoinjury"),times=c(2,2,2,4,4))
scores <- data.frame(sample.groups=labels, pca$x[,1:3],labs=rep(c("ours","Bednarek"),times=c(6,8)))


colnames(scores)[2:4]<-c("PC1","PC2","PC3")

min_axis<-min(scores[,2:4])
max_axis<-max(scores[,2:4])
library(grid)
cbPalette<-c("azure4","black","blue","brown","cadetblue","chartreuse","cyan",
             "darkorange","darkorchid","deeppink","gold","lightcoral","lightseagreen","magenta","red","lightsalmon","yellow","mediumorchid4","deepskyblue","mediumvioletred","olivedrab","cornsilk","lavender","navajowhite4")


library(ggplot2)
p1=ggplot(data=scores, aes(x=PC1,y=PC2,color=sample.groups))+ 
  #geom_point(aes(shape=factor(labs)),size=0.7) +
  geom_point(size=0.7) +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.15,"cm"),legend.justification=c(0,1),legend.margin=unit(-0.1,"cm"),legend.position="top",legend.text=element_text(size=5),axis.text.x = element_text(size=6),axis.text.y=element_text(size=6),axis.title=element_text(size=6),legend.title=element_blank()) + guides(color=guide_legend(ncol=4)) +
  xlim(min_axis,max_axis) + ylim(-100,max_axis+2) + scale_color_manual(values=cbPalette[c(2,15,8,3,4)]) + xlab("PC 1 (60%)") + ylab("PC2 (17%)") 
png("zebrafish_day3_adult_cryoinjury_mtz_PCA_plot_all_genes_08182016.png",res=600,width=4,height=2,units="in")
p1
dev.off()

