library(magrittr)
library(parallel)
library(tibble)
library(tidyr)
library(dplyr)
library(plyr)

ts_Tau <- 
  function(x, 
           na = "rm") {
    if (any(!is.na(x))) {
      if (na == "rm") {
        x <- x[!is.na(x)]
      } else if (na == "zero") {
        x[is.na(x)] <- 0
      } else {
        stop("'na' error!")
      }
      
      if (min(x, na.rm = TRUE) >= 0) {
        if(max(x) != 0) {
          x <- (1 - (x/max(x)))
          res <- sum(x, na.rm = TRUE)
          res <- res/(length(x) - 1)
        } else {
          res <- 0
        }
      } else {
        res <- NA
        #print("Expression values have to be positive!")
      }
    } else {
      res <- NA
      #print("No data for this gene avalable.")
    }
    return(res)
  }

sum_foo_psi <- 
  function(mat, ## matrix or df, with rownames (ExonicPartName) and colnames (sample ID)
           colDat, ## metainfo, with rownames (sample ID)
           col, ## colname of stages, ordered factor
           keep = 0.3 ## maximum fraction of NAs to keep when call correlation
  ) {
    require(magrittr)
    require(parallel)
    
    if (!all(rownames(colDat) %in% colnames(mat))) {
      stop("Number of samples in 'colDat' are less than that in 'mat'!")
    }
    
    if (!(col %in% colnames(colDat))) {
      stop("'col' is not in 'colDat'")
    }
    
    colDat[, col] <- droplevels(colDat[, col])
    designLvls <- levels(colDat[[col]])
    
    rownms <- rownames(mat)
    
    mat <- mat[, rownames(colDat), drop = FALSE]
    
    if (max(mat, na.rm = T) <= 1) {
      message("'mat' * 100 to be used!")
      mat <- mat * 100
    }
    
    #################################correlation####################################
    cor_table <- 
      suppressWarnings(
        mclapply(
          1:nrow(mat), 
          function(n) {
            dat <- 
              mat[n, ] %>% 
              gather(Run, PSI) %>% 
              merge(colDat)
            NAs <- is.na(dat$PSI)
            nonNAs <- which(!NAs)
            NAs_pct <- sum(NAs) / nrow(dat)
            if (NAs_pct > keep) {
              res_kw <- 
                tryCatch(
                  kruskal.test(
                    dat$PSI ~ dat[, col],
                    na.action = "na.omit"),
                  error = function(cond) {
                    res <- list(p.value = NA)
                    return(res)
                  }
                )
              res <- 
                c(
                  NAs_pct,
                  NA, 
                  NA,
                  res_kw$p.value)
            } else {
              res_spearman <- 
                cor.test(
                  x = as.integer(dat[, col]), 
                  y = dat$PSI, 
                  method = "spearman", 
                  exact = FALSE,
                  na.action = "na.omit")
              res_kw <- 
                kruskal.test(
                  dat$PSI ~ dat[, col],
                  na.action = "na.omit")
              
              res <- 
                c(
                  NAs_pct,
                  res_spearman$estimate, 
                  res_spearman$p.value,
                  res_kw$p.value)
            }
            return(res)
          },
          mc.cores = 20)
      ) %>%
      data.frame %>%
      do.call(rbind, .) %>%
      data.frame
    
    colnames(cor_table) <- 
      c(
        "NAs_fraction",
        "Spearman.cor",
        "Spearman.pvalue",
        "KW.pvalue")
    rownames(cor_table) <- rownms
    cor_table$KW.padj <- p.adjust(cor_table$KW.pvalue, method = "BH")
    cor_table$ExonicPartName <- rownms
    
    #################################basic stat####################################
    stat <- 
      mclapply(
        designLvls, 
        function(l) {
          num <- 
            grep(
              paste0("^", sub("\\)", ".", sub("\\(", ".", l)), "$"), 
              sub("\\)", ".", sub("\\(", ".", as.vector(colDat[[col]])))
              )
          basic_sum <- 
            apply(
              mat[, num, drop = FALSE], 
              1, 
              function(x) {
                res1 <- summary(x)
                res1[7] <- sd(x)
                res1 <- round(res1, 2)
                names(res1) <- 
                c("Min", "Q1st", "Median", "Mean", 
                  "Q3rd", "Max", "sd")
                return(res1)
              }
            ) %>%
            t %>%
            data.frame %>%
            add_column(
              ExonicPartName = rownms,
              StageAbbreviation = l
            )
          basic_sum$NumberOfSample <- length(num)
          
          return(basic_sum)
        },
        mc.cores = 20) %>%
      do.call(rbind, .)
    
    StageTau <- 
      dlply(
        stat, 
        "ExonicPartName",
        function(x) {
          res <- ts_Tau(x$Mean)
          return(res)
        }) %>%
      unlist %>% 
      round(3)
    
    stat$TauCrossStage <- StageTau[stat$ExonicPartName]
    
    res <- merge(stat, cor_table, by = "ExonicPartName", all =T)
    
    return(res)
  }

