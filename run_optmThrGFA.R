message("Begin script: run_optmThrGFA.R")
message(Sys.time())

# I recommend running this script on the server, but
# for testing on your local computer
local.nodename <- "L00019154" # to figure out your computers nodename,
# type `Sys.info()["nodename"]` into the console. The returned string
# is the nodename of your computer.
local.ed <- "/Volumes/T1000/Analysis/kforthman/ABCD_GFA" # The path to the
# working directory on your computer.

# If you run the script locally, it will set the working directory.
if(Sys.info()["nodename"] == local.nodename){setwd(local.ed)}

## 0.1 Required packages
library(GFA)
library(knitr)
library(gplots)
library(ggplot2)
library(foreach)

## 0.2 Set GFA Options
opts <- getDefaultOpts()
opts$convergenceCheck <- T

## 0.3 Load optmThrGFA
library(devtools)
install_github("kforthman/optmThrGFA")
library(optmThrGFA)

# Set the number of replicates
R <- 10

# File name
mthd <- "1"
message(paste0("\nmethod ", mthd, "\n"))
folder <- paste0("res/mthd", mthd,"/")
data_dir <- "data/"
match.filename <- paste0(folder, "match_orig_mthd", mthd,".rda")


# Set to TRUE if you would like to overwrite the files
# if they have already been created.
overwrite <- F

# Checks if the match file has been created.
# If the file has not been created (or overwrite = TRUE),
# the match file is generated.
if(!file.exists(match.filename) || overwrite){

  rep.summ <- data.frame(Replicate=1:R, conv=rep(NA, R), K=rep(NA, R))
  gfaList_full <- list()

  # Load in the generated GFA results. Recorded in the data.frame `rep.summ` are the values
  # - conv: An estimate of the convergence of the model's reconstruction based on Geweke diagnostic.
  # Values significantly above 0.05 imply a non-converged model,
  # and hence the need for a longer sampling chain.
  # - K: The number of components inferred.
  for(r in 1:R){
    # nda18_2.01_COG_CBCL_SOC_SMA[1:10]_10.31.2019.RDATA
    message(paste0("Loading replicate ", r, " of ", R))
    load(paste0(data_dir, "nda18_2.01_COG_CBCL_SOC_SMA", r, "_10.31.2019.RData"))
    gfaList_full[[r]] <- myres
    rep.summ$conv[r] <- myres$conv
    rep.summ$K[r] <- myres$K
  }
  remove(myres)

  message("\nAll replicates loaded.\n")

  message("\nLoading original dataset.")
  load(data_dir, "nda18_2.01_COG_CBCL_SOC_SMAdata.RData")
  message("Original dataset loaded.\n")

  Y <- as.matrix(MYdf)

  # corGrids and matchGrids are the threshold parameters.
  # corThr: How close two components are required to be, in terms of correlation, in order to match them.
  # matchThr: The proportion of sampling chains that need to contain a component in order to include it in the robust components.
  corGrids   <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))
  matchGrids <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))


  # Method 1: Using a list of original GFA objects
    # Create an empty matrix.
    tmp <- matrix(NA,
                  nrow = length(corGrids),
                  ncol = length(matchGrids),
                  dimnames = list(corThr=corGrids, matchThr=matchGrids))

    # K.girds is a matrix with 2 dimensions: corThresh (rows, i) and matchThresh (columns, j). It will contain the number of robust components found with the given thresholds i,j.
    # indicies is a list with 2 dimensions: corThresh and matchThresh. Every element is an RxK matrix containing the corresponding component indices. Negative indices denote the closest component in the corresponding repetition, having no components above the threshold.
    match_orig <- list(K.girds = tmp, mse.m = tmp, indices = rep( list(list()), length(corGrids) ))

    message("\nTesting thresholds...")
    for(i in 1:length(corGrids)){
      for (j in 1:length(matchGrids)){
        this.filename <- paste0(folder, "rComp_corThr", corGrids[i], "_matchThr", matchGrids[j], ".rda")

        if(!file.exists(this.filename) || overwrite){
          message(paste0("cor thresh ", i, " of ", length(corGrids),": ", corGrids[i], ",",
                         " match thresh ", j, " of ", length(matchGrids),": ", matchGrids[j]))
          rcomp <- robustComponents(gfaList_full, corThr=corGrids[i], matchThr=matchGrids[j])
          save(rcomp, file = this.filename)

        }
        else{
          load(this.filename)
        }

        match_orig$K.girds[i,j] <- rcomp$Krobust
        match_orig$indices[[i]][[j]] <- rcomp$indices

        if(rcomp$Krobust>0){

          # rcomp$effect is the component effect in the data space; and array of size NxDxK
          # rows: observations, columns: variables
          # the following sums the effects over the first 2 dimensions: N and D. From this we
          # get expected observations. We can compare this to the true observations to get MSE.
          yhat <- apply(rcomp$effect, 1:2, sum, na.rm=T)
          match_orig$mse.m[i,j] <- mean((yhat - Y)^2, na.rm = T)

        }
      }
    }
    message("\nSaving match_orig file.\n")
    save(match_orig, file = match.filename)

}else{
  load(match.filename)
}

opt.par <- optimizeK(K.grids=match_orig$K.grid, mse.array=match_orig$mse$all)
message(paste0("The min. MSE = ", round(opt.par$mse.min, 3)))
message(paste0("The 1-SE MSE threshold = ", round(opt.par$mseThr, 3)))
message(paste0("min. MSE criterion gives ", opt.par$Krobust.min, " matched factors"))
message(paste0("1-SE MSE criterion gives ", opt.par$Krobust.1se, " matched factors"))
# opt.par_full$par.min

message("Done.")
message(Sys.time())
