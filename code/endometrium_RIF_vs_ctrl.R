
###all packaegs we need

library(DoubletFinder)
library(Seurat) 
library(SeuratDisk)
library(SeuratWrappers)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


### (1)Load endometrium dataset and QC


endometrium_data<-readRDS(file="your pathway/Endometrium_demo.rds")

#QC
Endometriumlist<-list()
for (i in 1:3)
{
   temp_data<-endometrium_data[,which(endometrium_data@meta.data$orig.ident == unique(endometrium_data@meta.data$orig.ident)[i])]
   temp_data <- CreateSeuratObject(counts =temp_data,meta.data=temp_data@meta.data, min.cells = 3, min.features = 200)
   temp_data <- subset(temp_data , subset = nFeature_RNA > 200 & percent.mt < 20)
   temp_data <- DoubletFinder_function(temp_data,0.5)
   Endometriumlist[[i]]<- temp_data 
}


for (i in 2:length(Endometriumlist)){
  Endometrium <- merge(Endometrium,Endometriumlist[[i]])
}


remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
remove.packages("DoubletFinder")
remotes::install_github('https://github.com/ekernf01/DoubletFinder', force = T)


###(2)Remove doublets

DoubletFinder_function<-function(temp_data,resolution =resolution)
{

library(DoubletFinder)

# Pre-process Seurat object -------------------------------------------------------------------------------------------------
temp_data<- NormalizeData(temp_data)
temp_data<- FindVariableFeatures(temp_data,selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(temp_data)
temp_data<- ScaleData(temp_data, features = all.genes)
temp_data<- RunPCA(temp_data, features = VariableFeatures(object = temp_data))
temp_data<- FindNeighbors(temp_data, dims = 1:20)
temp_data <- FindClusters(temp_data, resolution = resolution,algorithm = 2)
temp_data<- RunUMAP(temp_data, dims = 1:20)

# pK Identification ----------------------------------------------------------
sweep.res.list <- paramSweep_v3(temp_data, PCs = 1:30, sct =FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) 
# Best parameters is mpK：
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))


#Homotypic Doublet Proportion Estimate 
annotations <- temp_data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(temp_data)*8*1e-6   #Calculate by increasing the ratio of doublets by 8 thousandths for every additional 1000 cells
nExp_poi <- round(DoubletRate*length(temp_data$seurat_clusters))
# Calculate the ratio of doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


#Identify doublets using established parameters
temp_data <- doubletFinder_v3(temp_data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
col_name<-colnames(temp_data)
temp_data@meta.data[,"DF_hi.lo"] <- temp_data@meta.data[,grep("DF.classifications",col_name)]
return(temp_data)

}



###(3)Intergate data form three samples and cluster

library(Seurat)
library(SeuratWrappers)


##Endometrium <- readRDS(file="your pathway/Endometrium_demo.rds")

#check doublets were delected
Endometrium <-Endometrium[,which(Endometrium@meta.data$DF_hi.lo == "Singlet")]

# Normalization
Endometriumlist <- SplitObject(Endometrium, split.by = "orig.ident")
for (i in 1:length(Endometriumlist )){
 temp <- Endometriumlist[[i]]
 temp  <- NormalizeData(temp , normalization.method = "LogNormalize")
 temp <- FindVariableFeatures(temp , selection.method = "vst", nfeatures = 3000)
 Endometriumlist[[i]]<-temp
}


#Intergation

Endometrium  <- NormalizeData(Endometrium, normalization.method = "LogNormalize")
Endometrium <- FindVariableFeatures(Endometrium, selection.method = "vst", nfeatures = 3000)
Endometrium <- RunFastMNN(object.list = c(Endometriumlist[[1]],Endometriumlist[[2]],Endometriumlist[[3]]))
Endometrium <- RunUMAP(Endometrium, reduction = "mnn", dims = 1:30)

#Find clusters

Endometrium <- FindNeighbors(Endometrium, reduction = "mnn", dims = 1:30)
Endometrium <- FindClusters(Endometrium, resolution =1.5,algorithm = 2)
Endometrium <- FindClusters(Endometrium, resolution =2.25,algorithm = 2)


#Find DEGs
Idents(Endometrium) <- 'celltype'
all_marker <- FindAllMarkers(Endometrium,only.pos = TRUE, min.pct = 0.01)



###(4)Visualization results

library(ggplot2)
library(RColorBrewer)


#Fig3.a
meta_data<-Endometrium@meta.data
meta_data<-cbind(meta_data,as.data.frame(Endometrium@reductions$umap@cell.embeddings))
stage_loc<-c("#b6766c","#4e6597","#b6a8cd","#fecea0")
pdf(file="your pathway/Endometrium_stage_mnn_normal.pdf",width=9,height=5)
for(i in 1:3)
{
  colors_grey<-rep("#eaeaea",3)
  colors_grey[i] <- stage_loc[i]
  stage<-c("Endometrium-RIF-1","Endometrium-RIF-2","Endometrium-RIF-Ct.")
  # Adjusting the order cell type appeared
  meta_data <- meta_data %>%
    arrange(desc(orig.ident != stage[i]))
  
  p<-ggplot(data =meta_data,aes(x=UMAP_1,y=UMAP_2,colour=as.factor(orig.ident)))+geom_point(size=0.2)+theme_bw()+
    scale_color_manual(values = colors_grey)+
    theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_blank(),legend.title = element_text(size = 16),
          legend.text = element_text(size = 12)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  print(p)
  
}
dev.off()


FeaturePlot(Endometrium, features = c( "CDH1","HMGCR","TPPP3","PIFO","CDH2","CD44","PLAT","MKI67","CCNB1"),col=c("lightgrey",brewer.pal(11, "Blues")[5:8]),ncol=3)

#Fig3.b

library(ComplexHeatmap)

marker_show<-c("FAU","RPL30","RPS19","HINT1","DKK1","SERPINE1","GREM1","MAGI2","HMGA2","COL1A1","EXT1","CCND1","IGFBP4","WNT5B","CDH13","COL8A1","TIMP1","TGFBI","CGA","SPON2","ISG15","PENK","HTRA1","PLAT",
 "MMP3","RGCC","LOXL2","COL3A1","CCNB1","MKI67","MAD2L1","SMC4","PLK1","EPCAM","CLDN3","EGLN3","CD24","FOS","PGK1","SPINT2","CXCL1","S100A14","ERBB4","ROBO2","PRKX","CALB1","LGR4","ESR1","CFAP47",
 "ROPN1L","DNAH11","EFHC1","CFAP43","PIERCE1","HYDIN")

all_marker %>%
  group_by(cluster) %>%
  slice_head(n = 50) %>%
  ungroup() -> top50

gene_cell_exp <- AverageExpression(Endometrium,
                                   features = unique(top50$gene),
                                   group.by = 'celltype',
                                   slot = 'data')
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[,c("Ciliated cell","Epithelium-1","Epithelium-2","Proliferating Epithelium","Proliferating Stromal cell","Stromal cell-1","Stromal cell-2","Stromal cell-3")]

marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))

pdf(file="E:/implantation_result/2014.02.22_version_2/0_RIF_WT/DEGs/new_version/top50_marker_hetamap.pdf",width=4.5,height=10)
gene<-rownames(marker_exp)
Heatmap(marker_exp,
            cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list( title=' '),
        col = colorRampPalette(rev(brewer.pal(11,"RdBu")[2:10]))(100),
        border = NA,row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize =7)
)+rowAnnotation(link = anno_mark(at = which(rownames(marker_exp) %in% marker_show), labels = gene[which(rownames(marker_exp) %in% marker_show)], labels_gp = gpar(fontsize = 6)))
dev.off()

#Fig3.c

library(scRNAtoolVis)

markers <- read.csv("your pathway/marker_rif_vs_wt_respectively.csv",head=T,sep=",")
markerVocalno(markers = markers,
              topn =20,
              labelCol = ggsci::pal_npg()(9))


#Extended data Fig6.e
markers_extended <- read.csv("your pathway/marker_rif_vs_wt_together.csv",head=T,sep=",")
for (i in 1:3)
{
temp <- markers_extended[which(markers_extended$cluster == unique(markers_extended$cluster)[i]),]
temp$log2_pct<- log2(temp$pct.1 /temp$pct.2)
gene_filter<-temp[which(abs(temp$avg_log2FC) >1 & abs(temp$log2_pct) >0.5 ),7]
pdf(file="your pathway/volcano.pdf",width=8,height=8)
ggplot(temp,aes(x=avg_log2FC, y=log2_pct))+
 geom_hline(yintercept = c(-0.5,0.5), linetype = "dashed", color = "black")+
 geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
 geom_point(aes(color =abs(avg_log2FC),size=2))+
 scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))+
 scale_size_continuous(range = c(0,2))+
theme_bw(base_size = 12)+
 theme(panel.grid = element_blank(),
      legend.position = 'right',
     legend.justification = c(0,1))+
  geom_text_repel(data = filter(temp, gene %in% gene_filter),
                  #max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  # filter the genes you want to label
                  aes(label = gene),
                  size = 4, 
                  color = 'black') 

dev.off()


###(5)Intergation and comparison with human endometrium data

NM_data<-readRDS(file="/share/home/yangyuan/implatation/2024.02.22version_2/NM_data/GSE111976_ct_endo_10x.rds")

NM_data_seurat <- NormalizeData(NM_data_seurat)
NM_data_seurat <- FindVariableFeatures(NM_data_seurat,selection.method = "vst", nfeatures = 3000)

anchors <- FindTransferAnchors(reference = Endometrium, query = NM_data_seurat, dims = 1:30,k.anchor = 10,reference.reduction = "pca")
predictions <- TransferData(anchorset =anchors, refdata = Endometrium$RNA_snn_res.1, dims = 1:30)
NM_data_seurat <- AddMetaData(NM_data_seurat, metadata = predictions)

NM_data_seurat <- MapQuery(anchorset = anchors, reference = Endometrium, query = NM_data_seurat,refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

#Visualization after stricter quality control

NM_data_seurat_1  <- subset(NM_data_seurat , subset = nFeature_RNA > 500 & percent.mt < 20 & nCount_RNA > 10000 )


umap_ref<-cbind(Endometrium@meta.data,Endometrium@reductions$umap@cell.embeddings)
umap_ref<-umap_ref[,c("UMAP_1","UMAP_2","celltype")]
umap_quary<-cbind(NM_data_seurat_temp_1@meta.data,NM_data_seurat_temp_1@reductions$ref.umap@cell.embeddings)
umap_quary<-umap_quary[,c("UMAP_1","UMAP_2","celltype")]

color<-c("Ciliated"="#FB8072","Unciliated epithelia 1"="#8DD3C7","Unciliated epithelia 2"="#80B1D3","Stromal fibroblasts"="#BC80BD","Singlet"="#D9D9D9")
umap<-rbind(umap_ref,umap_quary)
p<-ggplot(data =umap,aes(x=UMAP_1,y=UMAP_2,colour=as.factor(celltype)))+geom_point(size=0.5)+theme_bw()+
  scale_color_manual(values=color)+
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_blank(),legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

pdf(file="your pathway/ref_quary_all_umap.pdf",width=10,height=8)
p
dev.off()


###(5)Cell communication between endometrium and stromal cell

library(CellChat)

meta_data<- Endometrium@meta.data
meta_data$celltype_cellchat <- rep("Epithelium",dim(meta_data)[1])
meta_data[which(meta_data$celltype_adj %in% c("Stromal cell-1","Stromal cell-2","Stromal cell-3")),23]<-"Stromal cell"
new_clltype_all<-paste(meta_data$celltype_cellchat,meta_data$orig.ident,sep="_")
new_clltype_all<-gsub("Endometrium-","",new_clltype_all)
meta_data$celltype_cellchat <- new_clltype_all

Endometrium@meta.data$celltype_cellchat <- new_clltype_all
data.input <- Endometrium@assays$RNA@data
cellchat <- createCellChat(object = data.input, meta =Endometrium@meta.data, group.by = "celltype_cellchat")

cellchat <- addMeta(cellchat, meta = Endometrium@meta.data)
cellchat <- setIdent(cellchat, ident.use = "celltype_cellchat")  # set "labels" as default cell identity
levels(cellchat@idents)  # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#Set up ligand receptor interaction database

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use


#Preprocessing expression data for cell communication analysis

#Only take the genes in the database
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#Calculate communication probability and infer cellchat network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
LR.net<-subsetCommunication(cellchat)
LR.net$source_cell<-apply(matrix(as.character(LR.net$source),ncol=1),1,function(x){strsplit(x,"_")[[1]][1]})
LR.net$target_cell<-apply(matrix(as.character(LR.net$target),ncol=1),1,function(x){strsplit(x,"_")[[1]][1]})
LR.net_new<-LR.net[which(LR.net$source_cell != LR.net$target_cell),]
LR.net_new$S_T<- paste(LR.net_new$source,LR.net_new$target,sep=">")

LR.net_new$log<- -LR.net_new$pval 
pdf(file="your pathway/cellchat.pdf",width=4,height=10)
ggplot(LR.net_new, aes(x=S_T, y=interaction_name_2))+geom_point(aes(size=log,color=prob),shape=16)+
  scale_size_continuous(range=c(1.5,3.5)) +
  scale_color_gradientn(colours=rev(brewer.pal(11,"RdYlBu")[c(-6,-7)]),name="negLog_qvalue")+
theme_bw()+
theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90,size=5),axis.text.y = element_text(size=7),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
dev.off()
