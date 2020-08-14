run_maSigPro <-
  function(Species, ## One scientific Species abbreviation 
           Meta,    ## Sample information including Run, ScientificAbbreviation, SampleType, RelativeTime_PostFertilization
           Outpath, ## The path of output files (including Gene_CPM.Rds and Gene_CPM.tsv)
           Count,   ## The matrix of gene expression counts
           TPM      ## The matrix of gene tpm is used to filter genes
  ){
    library(data.table)
    library(maSigPro)
    library(edgeR)
    
    ## prepare data 
    if(class(Meta)!= "character"){
      Sampfiles <- data.table(Meta)
    } else{
      Sampfiles <- fread(Meta)
    }
    
    Samp_sub <- Sampfiles[ScientificAbbreviation == Species,]
    tissue = unique(Samp_sub$SampleType)
    
    SampInfo <- 
      lapply(
        tissue,
        function(x) {
          SampInfo <- Samp_sub[SampleType==x,]
          SampInfo$RelativeTime_PostFertilization <- log2(1+SampInfo$RelativeTime_PostFertilization)
          return(SampInfo)
        })
    names(SampInfo) <- tissue
    
    edesign <- 
      lapply(
        SampInfo,
        function(meta) {
          times <- meta$RelativeTime_PostFertilization
          edesign <- 
            data.frame(
              Time = times,
              Replicate = rep(1:length(unique(times)),
                              time = as.vector(table(times))),
              Group = rep(1,nrow(meta)),
              row.names = meta$Run
            )
        }
      )
    
    ## loading data 
    
    if(class(Count) != "character"){
      Cnt_sp_g <- as.data.frame(Count)
    } else {
      Cnt_sp_g <- read.table(Count, header=T)
    }
    
    if(class(TPM) != "character"){
      TPM_sp_g <- as.data.frame(TPM)
    } else {
      TPM_sp_g <- read.table(TPM, header=T)
    }
    
    
    CPM_spG_fil <- 
      lapply(
        SampInfo,
        function(meta) {
          tmp <- Cnt_sp_g[,meta$Run]
          rownames(tmp) <- Cnt_sp_g$gene_id
          
          ## cpm
          DGE <- DGEList(counts = tmp, group = meta$RelativeTime_PostFertilization)
          tmp <- data.frame(cpm(DGE,normalized.lib.sizes = TRUE))
          rownames(tmp) <- Cnt_sp_g$gene_id
          
          TPM_fil <- TPM_sp_g[,meta$Run]
          tmp <- tmp[rowSums(TPM_fil >= 1) >= (.15 * ncol(TPM_fil)), ]
          
          return(tmp)
        })
    
    ## maSigPro
    NBp.cpm <- list()
    NBt.cpm <- list()
    for (ts in tissue) {
      d <- make.design.matrix(edesign[[ts]], degree = 3)
      NBp <- p.vector(CPM_spG_fil[[ts]] , d, counts = TRUE, Q = 0.05) #theta=10
      NBp$i
      NBt <- T.fit(NBp, step.method = "backward", alfa = 0.05)
      NBp.cpm[[ts]] <- NBp
      NBt.cpm[[ts]] <- NBt
    }
    
    sp_gene_cpm <- list(pvector=NBp.cpm, Tfit= NBt.cpm)
    
    ## save
    for (ts in tissue) {
      res <- NBt.cpm[[ts]]$sol[, c(1,2)]
      ts <- gsub(" ", "", ts)
      write.table(
        res,
        file = paste0(Outpath,'/',Species,"_",ts,'_Gene_CPM.tsv'),
        row.names = T, 
        col.names = T, 
        quote = F, 
        sep='\t')
    }
    
    saveRDS(sp_gene_cpm,file=paste0(Outpath, '/',Species,"_Gene_CPM.Rds"))
  }

