################################################

CREATE_CTS = function(DATA_DIR){
  
  # DATA_DIR => dir containing STAR output
  
  
}
  
################################################

IMPORT_CTS = function(FILE){
  
  # file => table containing raw counts 
  
  cts_all = read.table(FILE)
  row.names(cts_all) = cts_all$ensembl
  cts_all = cts_all[,-1]
  return(cts_all)
  
}

################################################

IMPORT_COLDATA = function(FILE){
  
  coldata = read.table(FILE,sep = ";",header = T)
  
  row.names(coldata) = coldata$Run
  
  return(coldata)
  
}

################################################

FILTER_CTS_BY_BP = function(cts_all,coldata,BIOPROJECT,samples2exclude=NULL){
  
  # return a list containing coldata and cts filtered by BIOPROJECT
  
  # cts_all => full cts
  # coldata => full coldata
  # BIOPROJECT => vector containing 
  # samples2exclude => samples (Runs) to exclude
  
  suppressPackageStartupMessages(library())
  
  coldata_BP = coldata[coldata$BioProject %in% BIOPROJECT,]
  
  cts_BP = cts_all[,names(cts_all) %in% coldata_BP$Run]
  
  if(!is.null(samples2exclude)){ # run only if samples2exclude is not null
    `%notin%` <- Negate(`%in%`)
    coldata_BP = coldata_BP[rownames(coldata_BP) %notin% samples2exclude,]
    cts_BP = cts_BP[,colnames(cts_BP) %notin% samples2exclude]
  }

  # 
  if(!all(rownames(coldata_BP) == colnames(cts_BP))){
    cts_BP <- cts_BP[, rownames(coldata_BP)]
  }
  
  DATA_BP = list("coldata" = coldata_BP,"cts" = cts_BP)
  
  return(DATA_BP)
  
}

################################################

CREATE_DDS = function(DATA_BP){
  
  suppressPackageStartupMessages(library("DESeq2"))
  
  dds_BP <- DESeqDataSetFromMatrix(countData = DATA_BP[["cts"]],
                                   colData = DATA_BP[["coldata"]],
                                   design = ~ STATUS)

  # keep genes with more than 10 counts in all samples
  keep <- rowSums(counts(dds_BP)) >= 10
  dds_BP <- dds_BP[keep,]
  
  return(dds_BP)
  
}

################################################

PLOT_PCA = function(rld_BP,BIOPRJECT){
  
  suppressPackageStartupMessages(library(ggplot2))
  
  pcaData <- plotPCA(rld_BP, intgroup=c("STATUS"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  P = ggplot(pcaData, aes(PC1, PC2, color=STATUS, shape=STATUS)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()+ 
      ggtitle(BIOPRJECT)+
      geom_text(aes(label=name),hjust=0.5, vjust=1,color="black",size=2)
  
  return(P)
  
}

################################################

RUN_DESEQ2 = function(dds){
  
  library("DESeq2")
  library("BiocParallel")
  register(MulticoreParam(12))
  
  dds = DESeq(dds,
              parallel = T,
              quiet = T)
  
  return(dds)
  
}

################################################

CREATE_RESihw_OBJECT=function(dds,COND1,COND2,ALPHA=0.05,LFCth=0){
  
  # COND1='WT'
  # COND2='WT4h'  
  # ALPHA = 0.05 | 0.01 | 0.001 => corte do pvalue ajustado
  
  library('biomaRt')
  library('stats')
  library('DESeq2')
  library("IHW")
  library("BiocParallel")
  register(MulticoreParam(12))
  
  res <- results(dds, filterFun=ihw, contrast = c('STATUS',COND2,COND1),alpha = 0.05, 
                 parallel = TRUE)
  
  res$ensembl <- sapply( strsplit( rownames(res), split="nn+" ), "[", 1 )
  
  ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
  
  genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = res$ensembl,
                    mart = ensembl )
  
  idx <- match( res$ensembl, genemap$ensembl_gene_id )
  
  res$entrez <- genemap$entrezgene[ idx ]
  
  res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
  
  res$Sample_1=COND1
  res$Sample_2=COND2
  res=res[,c("ensembl","hgnc_symbol",
             "Sample_1","Sample_2",
             "baseMean","log2FoldChange",
             "lfcSE","stat",
             "pvalue","padj")]
  
  return(res)
  
}

################################################
