##
# Calculate cell type-trait associations
#
#' @param sscore1(eclipser) A dgeMatrix of eclipser specificity scores where
#' each row is a cell_name, column names are gene identifiers, and the final columns are the "specificity scores" "pvalue" "zscore".
#
#' @param sscore2(scdrs) A data.frame ("*.full_score.gz") for scDRS
#' with columns: index is "cell_name", normol_score is the "specificity scores", z-scores is the "zscore".
#
#' @param sscore3(sc-linker) A data.frame for sc-linker
#' with columns: trait_name is "cell_name", e_score is the "specificity scores", z-score is the "zscore".
#
#' @return A data.frame of trait integrated scores for each cells.
#' @export
#' 
#' #load R packages
#' library(data.table)
#' library(dplyr)
#' library(tidyr)

integrated_scores <- function(sscore1, sscore2, sscore3) {
  pvalue <- FDR <- NULL # due to non-standard evaluation notes in R
  
  # if given file path to file, check that it is loadable
  scores <- c("eclipser","scdrs","sesmic")
  scores[1] <- fread(paste0(scores[1],".csv"))
  scores[2] <- fread(paste0(scores[2],".csv"))
  scores[3] <- fread(paste0(scores[3],".csv"))
  
 
  # clean data formatting
  for (i in c(1:3)) {
    sscore.exp <- as.data.table(as.matrix(scores[i]), keep.rownames = T)
    setnames(sscore.exp, "ids", "cell_name")
    sscore.exp <- melt(sscore.exp, id.vars = "cell_name", variable.name = "zscore", value.name = "specificity")
    sscore.exp$cell_name <- as.character(sscore.exp$cell_name)
    sscore.exp$zscore <- as.character(sscore.exp$zscore)
    sscore.exp$specificity <- as.numeric(sscore.exp$specificity)
    
    paste0("sscore",i) <- sscore.exp
  }
  
  # set column names of scoring data
  setnames(sscore1, c("zscore","specificity"),c("zscore1","specificity1"))
  setnames(sscore2, c("zscore","specificity"),c("zscore2","specificity2"))
  setnames(sscore3, c("zscore","specificity"),c("zscore3","specificity3"))
  
  # combine the data with cell_name, specificity score, and zscore
  dt <- merge(sscore1,sscore2,sscore3,by = "cell_name")
  
  # calculate the integrated scores for each cell
  res <- as.data.frame(dt)
  
  exp <- res[,c("zscore1","zscore2","zscore3")]
  exp$avg_zscore <- apply(exp,1,mean)
  res$avg_zscore <- exp$avg_zscore
  
  res$integrated_score <- NA
  res$integrated_pval <- NA
  for (j in 1:length(res$cell_name)){
    slm <- speedglm::speedlm(res$avg_zscore ~ res$specificity1+res$specificity2+res$specificity3) # fast lm
    slm_summ <- summary(slm)$coefficients
    
    # test p-value
    integrated_pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4]/2 else (1 - slm_summ[2, 4]/2)
  }
  
  res[, FDR := stats::p.adjust(integrated_pval, method = "fdr")]
  return(res)
}


