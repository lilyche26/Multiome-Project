#! /usr/bin/Rscript
library(ArchR)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(motifmatchr)
library(GenomicRanges)
set.seed(4502)
#### Setting default number of Parallel threads to 16.
addArchRThreads(threads = 16) 
## Setting default genome to mm10 ## Setting a Genome and GeneAnnotation
addArchRGenome("mm10") 


# --------------------------------------------------------------------------------------
# 1.Creating Arrow Files
# --------------------------------------------------------------------------------------
dataset_loc <- "/volume2/volume2/multiome/cellrangeOutput"
ids<- c("E10.5_WT1","E10.5_WT2","E14.5_WT1","E14.5_WT2","E18.5_WT1","E18.5_WT2","E10.5_WT3","E8.5_WT1","E8.5_WT2","E9.5_WT1","E9.5_WT2","E12.5_WT1","E12.5_WT2","P1_WT1","P1_WT2","E7.5_WT1","E7.5_WT2")

inputFiles <- sapply (ids, function(i){fragpath <- file.path(dataset_loc,i,"outs/atac_fragments.tsv.gz")})
name <- c("E10.5_WT1","E10.5_WT2","E14.5_WT1","E14.5_WT2","E18.5_WT1","E18.5_WT2","E10.5_WT3","E8.5_WT1","E8.5_WT2","E9.5_WT1","E9.5_WT2","E12.5_WT1","E12.5_WT2","P1_WT1","P1_WT2","E7.5_WT1","E7.5_WT2")
setwd("/scATAC_17Samples_ArchR")
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,sampleNames = name,filterTSS = 4,filterFrags = 1000, addTileMat = TRUE,addGeneScoreMat = TRUE)


# --------------------------------------------------------------------------------------
# 3.Doublet Inference with ArchR
# --------------------------------------------------------------------------------------
#If these R2 values are much lower (i.e. less than 0.9), this often indicates that the cells within the Arrow file have very little heterogeneity. 
#This makes the accuracy of doublet calling worse
doubScores <- addDoubletScores(input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
    outDir = '/scATAC_17Samples_ArchR/QualityControl')



# --------------------------------------------------------------------------------------
# 4.Creating an ArchRProject and Filtering Doublets
# --------------------------------------------------------------------------------------
projMultiome_ATAC <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "Multiome_snATAC_ArchR2",copyArrows = FALSE)#TRUE )
paste0("Memory Size = ", round(object.size(projMultiome_ATAC) / 10^6, 3), " MB")
getAvailableMatrices(projMultiome_ATAC)
quantile(projMultiome_ATAC$TSSEnrichment)


#Plotting QC metrics
df <- getCellColData(projMultiome_ATAC, select = c("log10(nFrags)", "TSSEnrichment","Sample"))
ll <- c("E10.5_WT1","E10.5_WT2","E14.5_WT1","E14.5_WT2","E18.5_WT1","E18.5_WT2","E10.5_WT3","E8.5_WT1","E8.5_WT2","E9.5_WT1","E9.5_WT2","E12.5_WT1","E12.5_WT2","P1_WT1","P1_WT2","E7.5_WT1","E7.5_WT2")
for (i in ll){
    df1<-df[df$Sample==i,]
    p <- ggPoint(x = df1[,1], y = df1[,2], colorDensity = TRUE,continuousSet = "sambaNight",xlabel = "Log10 Unique Fragments",ylabel = "TSS Enrichment",xlim = c(log10(500), quantile(df[,1], probs = 0.99)),ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
    plotPDF(p, name = paste("TSS-vs-Frags.",i,".pdf",sep=""), ArchRProj = projMultiome_ATAC, addDOC = FALSE)}

#Make a ridge plot for each sample for the TSS enrichment scores.
p1 <- plotGroups(ArchRProj = projMultiome_ATAC,colorBy = "cellColData",name = "TSSEnrichment",plotAs = "ridges")
#plotPDF(p1, name = "RidgePlot.TSS.pdf", ArchRProj = projHeme1, addDOC = FALSE)
p3 <- plotGroups(ArchRProj = projMultiome_ATAC, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "ridges")
#plotPDF(p3, name = "RidgePlot.log10.nFrags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

#violin plot for each sample for the TSS enrichment scores.
p2 <- plotGroups(ArchRProj = projMultiome_ATAC, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin", alpha = 0.4,addBoxPlot = TRUE)
#plotPDF(p2, name = "ViolinPlot.TSS.pdf", ArchRProj = projHeme1, addDOC = FALSE)
p4 <- plotGroups(ArchRProj = projMultiome_ATAC, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
#plotPDF(p4, name = "ViolinPlot.log10.nFrags.pdf", ArchRProj = projHeme1, addDOC = FALSE)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projMultiome_ATAC, addDOC = FALSE, width = 4, height = 4)

#Fragment size distributions
p1 <- plotFragmentSizes(ArchRProj = projMultiome_ATAC)
#TSS enrichment profiles
p2 <- plotTSSEnrichment(ArchRProj = projMultiome_ATAC)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projMultiome_ATAC, addDOC = FALSE, width = 5, height = 5)

###Filtering Doublets
projMultiome_ATAC2 <- filterDoublets(projMultiome_ATAC, filterRatio = 1)
head(projMultiome_ATAC2@cellColData)

##Additional quality control removal of cell outliers.
projMultiome_ATAC2.1 <- projMultiome_ATAC2[which(projMultiome_ATAC2$TSSEnrichment > 6)]



# --------------------------------------------------------------------------------------
# 5.Dimensionality Reduction with ArchR and Batch Effect Correction
# --------------------------------------------------------------------------------------
projMultiome_ATAC3 <- addIterativeLSI(ArchRProj = projMultiome_ATAC2.1,useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10), varFeatures = 25000,dimsToUse = 1:30)
##TimePoints:
TimePoints <- gsub("_WT1|_WT2|_WT3|_MAB21L2","",projMultiome_ATAC3$Sample)
table(TimePoints)
projMultiome_ATAC3$TimePoints <- TimePoints


# --------------------------------------------------------------------------------------
# 6.Clustering in ArchR
# --------------------------------------------------------------------------------------
###1. Clustering using Seurat’s FindClusters
projMultiome_ATAC5 <- addClusters(input = projMultiome_ATAC4,reducedDims = "IterativeLSI",method = "Seurat",name = "Clusters",resolution = 0.8)
table(projMultiome_ATAC5$Clusters)


# --------------------------------------------------------------------------------------
# 6.Single-cell Embeddings:Dimensionality Reduction After Harmony
# --------------------------------------------------------------------------------------
projMultiome_ATAC6 <- addUMAP(ArchRProj = projMultiome_ATAC5, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")
p1 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projMultiome_ATAC6, addDOC = FALSE, width = 5, height = 5)

p3 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "TimePoints", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "Rep", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "Plot-UMAP2Harmony-TimePoints-Replicates.pdf", ArchRProj = projMultiome_ATAC6, addDOC = FALSE, width = 5, height = 5)



##TimePoint Highlight
##add columns to cellColData
HighlightE7.5 <- gsub("E8.5|E9.5|E10.5|E12.5|E14.5|E18.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightE8.5 <- gsub("E7.5|E9.5|E10.5|E12.5|E14.5|E18.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightE9.5 <- gsub("E7.5|E8.5|E10.5|E12.5|E14.5|E18.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightE10.5 <- gsub("E7.5|E8.5|E9.5|E12.5|E14.5|E18.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightE12.5 <- gsub("E7.5|E8.5|E9.5|E10.5|E14.5|E18.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightE14.5 <- gsub("E7.5|E8.5|E9.5|E10.5|E12.5|E18.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightE18.5 <- gsub("E7.5|E8.5|E9.5|E10.5|E12.5|E14.5|P1","Non.Highlighted",projMultiome_ATAC6$TimePoints)
HighlightP1 <- gsub("E7.5|E8.5|E9.5|E10.5|E12.5|E14.5|E18.5","Non.Highlighted",projMultiome_ATAC6$TimePoints)



p3 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightE8.5", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightE8.5=="E8.5",]),pal=c("E8.5" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightE9.5", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightE9.5=="E9.5",]),pal=c("E9.5" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "Plot-UMAP2Harmony-TimePoints-Highlight.E8.5.E9.5.pdf", ArchRProj = projMultiome_ATAC6, addDOC = FALSE, width = 5, height = 5)


p3 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightE10.5", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightE10.5=="E10.5",]),pal=c("E10.5" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
plotPDF(p3,name = "Plot-UMAP2Harmony-TimePoints-Highlight.E10.5.pdf", ArchRProj = projMultiome_ATAC6, addDOC = FALSE, width = 5, height = 5)



p3 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightE12.5", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightE12.5=="E12.5",]),pal=c("E12.5" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightE14.5", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightE14.5=="E14.5",]),pal=c("E14.5" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "Plot-UMAP2Harmony-TimePoints-Highlight.E12.5.E14.5.pdf", ArchRProj = projMultiome_ATAC6, addDOC = FALSE, width = 5, height = 5)


p3 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightE18.5", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightE18.5=="E18.5",]),pal=c("E18.5" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projMultiome_ATAC6, colorBy = "cellColData", name = "HighlightP1", highlightCells=rownames(projMultiome_ATAC6@cellColData[projMultiome_ATAC6$HighlightP1=="P1",]),pal=c("P1" = "red", "Non.Highlighted" = "lightgrey"),embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "Plot-UMAP2Harmony-TimePoints-Highlight.E18.5.P1.pdf", ArchRProj = projMultiome_ATAC6, addDOC = FALSE, width = 5, height = 5)


saveRDS(projMultiome_ATAC6,'/volume2/volume2/multiome/scATAC_17Samples_ArchR/projMultiome_ATAC6.rds')
saveArchRProject(ArchRProj = projMultiome_ATAC6, outputDirectory = "Save-projMultiome_ATAC", load = FALSE, overwrite=FALSE)
projMultiome_ATAC6 <- loadArchRProject(path = "Save-projMultiome_ATAC", force = FALSE, showLogo = TRUE)
ArchRBrowser(projMultiome_ATAC6)

##highlight timepoints
for (i in unique(projMultiome_ATAC6$TimePoints)){
filename = paste("Highlight.",i,".pdf",sep="")
p=plotEmbedding(projMultiome_ATAC6,embedding = "UMAPHarmony",colorBy = "cellColData",name = "TimePoints",size = 0.05,sampleCells = NULL,highlightCells = getCellNames(ArchRProj = projMultiome_ATAC6)[which(projMultiome_ATAC6@cellColData$TimePoints == i)],baseSize = 10,plotAs = "points")
plotPDF(p, name = filename, ArchRProj = projMultiome_ATAC6, addDOC = FALSE,width = 5, height = 5)
}





# --------------------------------------------------------------------------------------
# 6.1. Calculating Gene Scores and Marker Genes with ArchR //jump, we not annotate cluster by gene scores, we annotated with snRNA
# --------------------------------------------------------------------------------------
##use MAGIC to impute gene scores by smoothing signal across nearby cells,this greatly improves the visual interpretation of gene scores. 
#add impute weights to our ArchRProject
# projMultiome_ATAC7 <- addImputeWeights(projMultiome_ATAC6)

# markerGenes  <- c(
#     "CD34",  #Early Progenitor
#     "GATA1", #Erythroid
#     "PAX5", "MS4A1", "MME", #B-Cell Trajectory
#     "CD14", "MPO", #Monocytes
#     "CD3D", "CD8A"#TCells
#   )
# p <- plotEmbedding(
#     ArchRProj = projHeme2, 
#     colorBy = "GeneScoreMatrix", 
#     name = markerGenes, 
#     embedding = "UMAP",
#     imputeWeights = getImputeWeights(projHeme2))

# #Rearrange for grid plotting
# p2 <- lapply(p, function(x){
#     x + guides(color = FALSE, fill = FALSE) + 
#     theme_ArchR(baseSize = 6.5) +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#     theme(
#         axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(), 
#         axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank()
#     )
# })
# do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

# plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", ArchRProj = projHeme2,addDOC = FALSE, width = 5, height = 5)








