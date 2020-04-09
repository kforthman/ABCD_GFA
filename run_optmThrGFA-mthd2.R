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

# Number of replicates
R <- 10

# file name
mthd <- "2"
message(paste0("\nmethod ", mthd, "\n"))
folder <- paste0("res/mthd", mthd,"/")
data_dir <- "data/"
match.filename <- paste0(folder, "match_orig_mthd", mthd,".rda")

overwrite_match <- T
overwrite_gfaList <- F
overwrite_xw <- T

if(!file.exists(match.filename) || overwrite_match){
  
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
  load(paste0(data_dir, "nda18_2.01_COG_CBCL_SOC_SMAdata.RData"))
  message("Original dataset loaded.\n")
  
  Y <- as.matrix(MYdf)
  
  # corGrids and matchGrids are the threshold parameters.
  # corThr: How close two components are required to be, in terms of correlation, in order to match them.
  # matchThr: The proportion of sampling chains that need to contain a component in order to include it in the robust components.
  corGrids   <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))
  matchGrids <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))
  
  
  
  # Method 2: Using a list of posterior means of reconstructed data, lower memory required
  gfaList.filename <- paste0(folder, "gfaList_p50.rda")
  if(!file.exists(gfaList.filename) | overwrite_gfaList){
    message("\nCreating GFA list.")
    gfaList_p50 <- list()
    for (r in 1:R){ 
      message(paste0("running replicate ", r, " of ", R))
      gfaList_p50[[r]] <- psSummary(gfa.res=gfaList_full[[r]], credible.lv=0.95)
    }
    save(gfaList_p50, file = gfaList.filename)
    message("GFA list created.")
  }else{
    message("\nLoading GFA list.")
    load(gfaList.filename)
    message("GFA list loaded.")
  }

  xw_by_rep_comp.filename <- paste0(folder, "xw_by_rep_comp.rda")
  if(!file.exists(xw_by_rep_comp.filename) | overwrite_xw){
    message("\nCreating xw_by_rep_comp.")
    xw <- list()
    for(r in 1:R){
      xw[[r]] <- gfaList_p50[[r]]$Yhat.p50
      names(xw[[r]]) <- paste0("K", 1:length(xw[[r]]))
    }
    names(xw) <- paste0("Rep", 1:R)
    # Sort replicates by convergence value
    rep.use <- sort(rep.summ$conv, index.return=T)
    
    # !!! Re-orders replicates but does not save new order in output file!
    xw <- xw[rep.use$ix[1:10]]
    save(xw, file = xw_by_rep_comp.filename)
    message("xw_by_rep_comp created.")
  }else{
    message("\nLoading xw_by_rep_comp.")
    load(xw_by_rep_comp.filename)
    message("xw_by_rep_comp loaded.")
  }
  
  match.mse <- MSE.Grids(
    Ymtx=Y,                       ## the observed (normalized) data matrix (N x D)
    maxK = max(rep.summ$K),       ## the maximal K among the GFA replicates
    comps=xw,                     ## a list of GFA replicates with posterior medians
    corGrids=corGrids,            ## the grids of corThr values to be assessed
    matchGrids=matchGrids)        ## the grids of matchThr values to be assessed
  save(match.mse, file = match.filename)
  
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
