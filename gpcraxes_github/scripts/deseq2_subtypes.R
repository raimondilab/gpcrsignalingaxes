#!/usr/bin/env Rscript

# Custom addition function to handle character concatenation
`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

# Retrieve the cancer type from command line arguments
args = commandArgs(trailingOnly = TRUE)
cancer = args[1]

# Set the output path based on the cancer type
out_path = '/home/carora/subtypes_TvG_all/' + cancer + '/'
dir.create(out_path)

# Load the merged dataset
setwd(out_path)
full <- readRDS(file = "/home/carora/TCGA_GTEX_meta_count_merged.rds")
full2 <- subset(full, ((X_sample_type != 'Solid Tissue Normal')))

# Subset data based on the provided cancer type
if (cancer != 'Thyroid') {
  sub_full1 = subset(full2, X_primary_site == cancer)
} else {
  sub_full1 = subset(full2, (X_primary_site == 'Thyroid') | (X_primary_site == 'Thyroid Gland'))
}

# Load TCGAbiolinks package and retrieve subtypes
library(TCGAbiolinks)
subtypes <- PanCancerAtlas_subtypes()

# ... Additional preprocessing on subtypes ...

# Loop through each subtype
a = unique(subtypes[c("Subtype_Selected")])
for (i in seq(from = 1, to = dim(a)[1], by = 1)) {
  tcga1 <- subset(sub_full1, X_study == 'TCGA')
  gtex <- subset(sub_full1, X_study == 'GTEX')
  
  sub_filt <- subset(subtypes, Subtype_Selected == array(a)[[1]][i])
  
  tcga <- tcga1[tcga1$sample %in% sub_filt$pan.samplesID, ]
  
  sub_full <- rbind(tcga, gtex)
  
  if ((dim(tcga)[1] > 0) & (dim(gtex)[1] > 0)) {
  
  mapfile<-read.csv(file='/data/probemap.csv',sep='\t')[1:2]
  metadata<-sub_full[1:7]
  count_data<-subset(sub_full, select = -c(2:7) )
  count_data <- as.data.frame(t(as.matrix(count_data)))
  names(count_data) <- as.matrix(count_data[1, ])
  count_data <- count_data[-1, ]
  count_data[] <- lapply(count_data, function(x) type.convert(as.character(x)))
  count_data <- cbind(newColName = rownames(count_data), count_data)
  rownames(count_data) <- 1:nrow(count_data)
  colnames(count_data)[1]='gene_id'
  
# Perform DE analysis using DESeq2
  suppressMessages(library( "DESeq2" ))
  suppressMessages(library(ggplot2))
  dds <- DESeqDataSetFromMatrix(countData=count_data, 
                                colData=metadata, 
                                design=~X_study, tidy=TRUE )
  
  dds <- DESeq(dds)
  res <- results(dds, c("X_study", "TCGA", "GTEX"))
  res <- results(dds)
  res <- res[order(res$padj),]
  resultsNames(dds)
  
# Annotate DE results with gene symbols
  write.csv(res,out_path+cancer+'_'+array(a)[[1]][i]+"_DE.csv")
  res2<-read.csv(out_path+cancer+'_'+array(a)[[1]][i]+"_DE.csv")
  colnames(res2)[1]="id"
  res2 = merge(x = res2, y = mapfile, by = "id")
  require(dplyr)
  res2 <- res2 %>% relocate(gene, .before = baseMean)
  write.csv(res2,out_path+cancer+'_'+array(a)[[1]][i]+"_DE_annotated.csv")
  
  
  #we can use plotCounts fxn to compare the normalized counts
  #between treated and control groups for our top 6 genes
  
  png(out_path+cancer+'_'+array(a)[[1]][i]+'_top6_count_comp.png',width = 10, height = 10, units = 'in', res = 600)
  par(mfrow=c(2,3))
  
  plotCounts(dds, gene=res@rownames[1], intgroup="X_study")
  plotCounts(dds, gene=res@rownames[2], intgroup="X_study")
  plotCounts(dds, gene=res@rownames[3], intgroup="X_study")
  plotCounts(dds, gene=res@rownames[4], intgroup="X_study")
  plotCounts(dds, gene=res@rownames[5], intgroup="X_study")
  plotCounts(dds, gene=res@rownames[6], intgroup="X_study")
  
  dev.off()
  
  #VOLCANO PLOT
  
  png(out_path+cancer+'_'+array(a)[[1]][i]+'_volcano_plot.png',width = 10, height = 10, units = 'in', res = 600)
  
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  dev.off()
  
  #PCA
  #First we need to transform the raw count data
  #vst function will perform variance stabilizing transformation
  png(out_path+cancer+'_'+array(a)[[1]][i]+'_PCA.png',width = 10, height = 10, units = 'in', res = 600)
  
  vsdata <- vst(dds, blind=FALSE)
  plotPCA(vsdata, intgroup="X_study") #using the DESEQ2 plotPCA fxn
  
  dev.off()
  
  
  gs<-as.character(subset(res2,padj<0.01)$gene)
  
  
  #map symbol to entrezid for ggea
  library(org.Hs.eg.db)
  hs <- org.Hs.eg.db
  gs2<-AnnotationDbi::select(hs, 
                             keys = gs,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
  
  
  ## GSEA using WEBGESTALT
  
  ## creating a df with two columns and condition logFC>=1 or logFC<=-1 with padj<0.01
  library(dplyr)
  #x<- dplyr::select(subset(res2,((log2FoldChange>=1)|(log2FoldChange<=-1))&((padj<0.01))), c('gene', 'log2FoldChange'))
  x<- dplyr::select(subset(res2,(padj<0.01)),c('gene', 'log2FoldChange'))
  
  #res2$id <- gsub("\\..*","",res2$id)
  #x<- dplyr::select(subset(res2,(padj<0.01)),c('id', 'log2FoldChange'))
  
  ## Run GSEA
  library(WebGestaltR)
  WebGestaltR(
    enrichMethod = "GSEA",
    organism = "hsapiens",
    enrichDatabase = "pathway_Reactome",
    enrichDatabaseFile = NULL,
    enrichDatabaseType = "genesymbol",
    enrichDatabaseDescriptionFile = NULL,
    interestGeneFile = NULL,
    interestGene =x ,
    fdrThr = 1,
    interestGeneType = "genesymbol",
    collapseMethod = "mean",
    referenceGeneFile = NULL,
    referenceGene = NULL,
    referenceGeneType = NULL,
    referenceSet = "genome",
    outputDirectory=out_path,
    projectName=cancer+'_'+array(a)[[1]][i]+'_WebGestalt-GSEA')
  ###############################################################################
  
  }
}


