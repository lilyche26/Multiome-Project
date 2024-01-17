library(devtools)
library(Seurat) 
library(dplyr)
library(sctransform)
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
set.seed(1368)
library(ggplot2)

test<-readRDS("work.space.che.test.rds")

##different timepoint/Compartment/cluster:
TimePoint=c("E7.5","E8.5","E9.5","E10.5","E12.5","E14.5","E18.5","P1")
Compartment=c("Cardiomyocytes","Endothelium","Mesenchymal","Others")


# for (i in TimePoint){
for (i in Compartment){
Stage1 <- subset(x= test,subset = compartment==i)
# Stage1 <- subset(x= test,subset = time_point==i)
Idents(Stage1) <- paste("sub.clustet",i,sep="")
NewPath<-paste("/volume2/volume2/multiome/NewData18Samples/",i,sep="")
dir.create(NewPath)
setwd(NewPath)

##resolution:
RES=1.0
##umap reduction name
ATACname=paste("umap.atac.",i,sep="")
RNAname=paste("umap.rna.",i,sep="")
WNNname=paste("umap.multimodal.",i,sep="")


# RNA analysis
DefaultAssay(Stage1) <- "RNA"
Stage1_filter1 <- NormalizeData(Stage1)
Stage1_filter2 <- FindVariableFeatures(Stage1_filter1)
Stage1_filter2.list <- SplitObject(Stage1_filter2, split.by = "replicate")
Stage1_filter2.list <- lapply(X = Stage1_filter2.list, FUN = function(x) {
    x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
features <- SelectIntegrationFeatures(object.list = Stage1_filter2.list)
head(features, n=500)
Stage1_filter3.anchors <- FindIntegrationAnchors(object.list = Stage1_filter2.list, anchor.features = features)
Stage1_filter3.combined <- IntegrateData(anchorset = Stage1_filter3.anchors)
DefaultAssay(Stage1_filter3.combined) <- "integrated"
Stage1_filter4 <- ScaleData(Stage1_filter3.combined)
Stage1_filter5 <- RunPCA(Stage1_filter4)
Stage1_filter6 <- FindNeighbors(Stage1_filter5, dims = 1:40)
Stage1_filter7 <- FindClusters(Stage1_filter6, resolution = RES)
Stage1_filter8 = RunUMAP(Stage1_filter7, dims = 1:50, reduction.name = RNAname, reduction.key = paste('umaprna', gsub("\\.", "", i),"_",sep="")) 




# ATAC analysis
library(harmony)
DefaultAssay(Stage1_filter8) <- "ATAC"
Stage1_filter9 <- RunTFIDF(Stage1_filter8)
Stage1_filter10 <- FindTopFeatures(Stage1_filter9, min.cutoff = 'q0')
Stage1_filter11 <- RunSVD(Stage1_filter10, n = 50, reduction.name = 'lsi', reduction.key = 'LSI_')
Stage1_filter12.integrated <- RunHarmony(object = Stage1_filter11,group.by.vars = 'replicate',reduction = 'lsi',assay.use = 'ATAC',project.dim = FALSE)
Stage1_filter13 <- RunUMAP(Stage1_filter12.integrated, dims = 2:50, reduction = 'harmony',reduction.name = ATACname,reduction.key= paste("umapatac",gsub("\\.", "", i),"_",sep=""), assay = "ATAC")
saveRDS(Stage1_filter13, file = paste(i,"filter13.rds",sep="_"))
# Stage1_filter13 <- readRDS("E7.5_filter13.rds")



# Calculate a WNN graph
Stage1_filter14 <- FindMultiModalNeighbors(Stage1_filter13, reduction.list = list("pca", "harmony"), dims.list = list(1:50, 2:50),modality.weight.name = "RNA.weight", verbose = TRUE)
Stage1_filter15 <- RunUMAP(Stage1_filter14, nn.name = "weighted.nn", reduction.name = WNNname,verbose = TRUE)
Stage1_filter16 <- FindClusters(Stage1_filter15, graph.name = "wsnn", resolution = RES,algorithm = 3, verbose = FALSE)
test_final<-Stage1_filter16
table(test_final$seurat_clusters)
DefaultAssay(test_final) <- "RNA"
length(rownames(test_final@meta.data))##cell number
saveRDS(test_final, file = paste("work.space.test_final.",i,".res",RES,".final.rds",sep=""))
# dev.off()





# Visualization
pdf(file="./1.uMap.Timepoints.pdf",width = 30, height =10)
p1 <- DimPlot(test_final, reduction = RNAname, label = TRUE,label.size=5,group.by = "time_point")+ ggtitle("RNA")
p2 <- DimPlot(test_final, reduction = ATACname, label = TRUE, label.size=5,group.by = "time_point") + ggtitle("ATAC")
p3 <- DimPlot(test_final, reduction = WNNname, label = TRUE, label.size=5,group.by = "time_point")+ ggtitle("WNN")
print (p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

pdf(file=paste("./1.uMap.WNN.cluster.res.",RES,".pdf",sep=""),width = 30, height =10)
p1 <- DimPlot(test_final, reduction = RNAname, label = TRUE,label.size=5,group.by = "seurat_clusters")+ ggtitle("RNA")
p2 <- DimPlot(test_final, reduction = ATACname, label = TRUE, label.size=5,group.by = "seurat_clusters") + ggtitle("ATAC")
p3 <- DimPlot(test_final, reduction = WNNname, label = TRUE, label.size=5,group.by = "seurat_clusters")+ ggtitle("WNN")
print (p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

pdf(file="./1.uMap.All.Sample.pdf",width = 30, height =10)
p1 <- DimPlot(test_final, reduction = RNAname, label = TRUE,label.size=5,group.by = "sample")+ ggtitle("RNA")
p2 <- DimPlot(test_final, reduction = ATACname, label = TRUE, label.size=5,group.by = "sample") + ggtitle("ATAC")
p3 <- DimPlot(test_final, reduction = WNNname, label = TRUE, label.size=5,group.by = "sample")+ ggtitle("WNN")
print (p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

pdf(file="./1.uMap.All.Repli.pdf",width = 30, height =10)
p1 <- DimPlot(test_final, reduction = RNAname, label = TRUE,label.size=5,group.by = "replicate")+ ggtitle("RNA")
p2 <- DimPlot(test_final, reduction = ATACname, label = TRUE, label.size=5,group.by = "replicate") + ggtitle("ATAC")
p3 <- DimPlot(test_final, reduction = WNNname, label = TRUE, label.size=5,group.by = "replicate")+ ggtitle("WNN")
print (p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()





####check the data
library(gprofiler2)
DefaultAssay(test_final) <- "RNA"
log.UMI <- log(Matrix::colSums(test_final@assays$RNA[,]),10)
test_final <- AddMetaData(object = test_final, metadata = log.UMI, col.name = "log.UMI")
s.genes <- read.table(file='/volume2/volume2/sharedData/G2Mannotation/s_genes.txt', header=FALSE, stringsAsFactors=FALSE)$V1
g2m.genes <- read.table(file='/volume2/volume2/sharedData/G2Mannotation/g2m_genes.txt', header=FALSE, stringsAsFactors=FALSE)$V1
mmus_s = gorth(s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(g2m.genes , source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
test_final <- CellCycleScoring(test_final, mmus_g2m, mmus_s, set.ident = FALSE)
test_final@meta.data$mols.per.gene<-test_final@meta.data$nCount_RNA/(test_final@meta.data$nFeature_RNA)###add UMI.per.gene
png(file="./nGene.png"); print (FeaturePlot(test_final, features =c('nFeature_RNA'),reduction = WNNname, pt.size=0.5)); dev.off()
png(file="./per.Mito.png"); print (FeaturePlot(test_final, features =c('percent.mt'),reduction = WNNname, pt.size=0.5)); dev.off()
png(file="./nCount_RNA.png");print (FeaturePlot(test_final, features =c('nCount_RNA'),reduction = WNNname, pt.size=0.5)); dev.off()
png(file="./G2M.Score.png"); print (FeaturePlot(test_final, features =c('G2M.Score'), reduction = WNNname, pt.size=0.5)); dev.off()
png(file="./S.Score.png"); print (FeaturePlot(test_final, features =c('S.Score'), reduction =WNNname, pt.size=0.5)); dev.off()
png(file="./mols.per.gene.png"); print (FeaturePlot(test_final, features =c('mols.per.gene'), reduction = WNNname, pt.size=0.5)); dev.off()
png(file="./log.UMI.png"); print (FeaturePlot(test_final, features =c('log.UMI'), reduction =WNNname, pt.size=0.2)); dev.off()
png(file="./nucleosomePercentile.png"); print (FeaturePlot(test_final, features =c('nucleosome_percentile'), reduction =WNNname, pt.size=0.2)); dev.off()
png(file="./TSS.enrichment.png"); print (FeaturePlot(test_final, features =c('TSS.enrichment'), reduction =WNNname, pt.size=0.2)); dev.off()
png(file="./TSS.percentile.png"); print (FeaturePlot(test_final, features =c('TSS.percentile'), reduction =WNNname, pt.size=0.2)); dev.off()
png(file="./nCount_ATAC.png"); print (FeaturePlot(test_final, features =c('nCount_ATAC'),reduction = WNNname, pt.size=0.5)); dev.off()
png(file="./nFeature_ATAC.png"); print (FeaturePlot(test_final, features =c('nFeature_ATAC'),reduction = WNNname, pt.size=0.5)); dev.off()



# 6.0 Finding DE genes
DefaultAssay(test_final) <- "RNA"
markers_all.RNA <- FindAllMarkers(test_final, only.pos = TRUE,test.use = "bimod",thresh.use = 0.25, logfc.threshold=0.25,return.thresh=0.01,pseudocount.use=0,min.cells.gene=3,min.cells.group=3, min.pct = 0.25)
x = markers_all.RNA 
write.table(x,paste("DeGene.list.",i,".txt",sep=""))
save(markers_all.RNA,file=paste("work.space.18Samples.",i,".res",RES,".markers_all.RNA.Rdata",sep=""))





##top20 gene umap:
markers_all.RNA %>% group_by(cluster) %>% top_n(20, wt = avg_log2FC)->top20    ##show top20 genes each cluster
top20<-as.data.frame(top20)
NewPath_top20<-paste("/volume2/volume2/multiome/NewData18Samples/",i,"1","/Top20",sep="")
dir.create(NewPath_top20)
setwd(NewPath_top20)
candi <-top20$gene
##check gene
xx <-rownames(test_final@assays$RNA)
for (ii in candi ){if (ii%in%xx==TRUE) print ("YES") else print (paste(ii,"NO"))}
for (ii in candi ){if (ii%in%xx==TRUE) print ("YES") else candi<-candi[!candi%in%c(ii)]}
##plot
for (m in candi){
filename = paste(m,".png",sep="")
png(file= filename)
print (FeaturePlot(test_final, features=c(m), reduction = WNNname,pt.size=0.5,cols = c("grey", "blue")))
dev.off()
}

##Other gene umap:
Others <- c("Runx1t1","Isl1","Tbx5","Gata4","Hand1","Mab21l2","Tbx18","Tbx3","Wt1","Tcf21","Upk3b","Lum","Postn","Pdgfra","Col12a1","Shox2","Hcn4","Myl7","Tnnt2","Myh6","Myh7","Ttn","Nkx2-5","Runx1","Gata1","Gata5","Mef2c","Hand2","Fgf10","Fgf8","Nppa","Hoxb1","Twist1","Msx2","Bmp4","Bmp2","Wnt2","Wnt5a","Tbx2","Rspo3","Sfrp5","Lhx2","Afp","Foxa2","Foxc2","Foxc1","Myl2","Myh11","Tagln","Actn1","Ttr","Cdh1","Pecam1","Vwf","Egr1","Npr3","Tek","Flt1","Hba-a1","Hba-x","Cd68","Lyve1","Hey2","Cck","Sox2","Nr2f2","Cdh5","Tnni3","Vsnl1","Stard10","Nr2f1","Cav1","Fxyd5","Itm2a","Loxl2","Sema3a","Sema3c","Dlx5","Dlx2","Otx2","Crabp2","Fxyd5","Cxcl12","Tnc","Ptn","Pla2g7","Aldh1a2","Col1a2","Homer2","Pitx2","Nfatc1","Adamtsl1","Nr6a1","Prtg","Irx4","Mb","Myl3","Kcnh7","Fgf12","Rgs6","Smoc2","Tmem108","Nrxn3","Hmcn1","Dach1","Cd36","Prox1","Tbx20","Dcn","Robo2","Ebf2","Msln","Sox10","Foxd3","Tbx1","Sfrp1","Osr1","Adamts8","Bmp10","Angpt1")
xx <-rownames(test_final@assays$RNA)
for (ii in Others ){if (ii%in%xx==TRUE) print ("YES") else print (paste(ii,"NO"))}

for (ii in Others){
filename = paste(ii,".png",sep="")
png(file= filename)
print (FeaturePlot(test_final, features=c(ii), reduction = WNNname,pt.size=0.5,cols = c("grey", "blue")))
dev.off()
}
}
