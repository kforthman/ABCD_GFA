message("Begin script: run_optmThrGFA.R")
message(Sys.time())

# ---- Set working directory ----

# I recommend running this script on the server, but
# for testing on your local computer
local.nodename <- "L00019154" # to figure out your computers nodename,
# type `Sys.info()["nodename"]` into the console. The returned string
# is the nodename of your computer.
local.ed <- "/Volumes/T1000/Analysis/kforthman/ABCD_GFA" # The path to the
# working directory on your computer.

# If you run the script locally, it will set the working directory.
if(Sys.info()["nodename"] == local.nodename){setwd(local.ed)}


# ---- Load packages ----

## Required packages
library(GFA)
library(knitr)
library(gplots)
library(ggplot2)
library(foreach)
library(DMwR)

## Load optmThrGFA
library(devtools)
install_github("kforthman/optmThrGFA")
library(optmThrGFA)


# ---- Set GFA Options ----
opts <- getDefaultOpts()
opts$convergenceCheck <- T
# opts$verbose <- 0 # Set to 0 if you do not 
# want to print results of each gfa replicate


# ---- Set file paths ----
folder <- "res" # The folder where you want to put the results
if(!file.exists(folder)){dir.create(folder)} # Will create the results folder if it does not exist
data_dir <- "data" # The folder that contains your raw data
match.filename <- paste0(folder, "/match_orig.rda")


# ---- Load raw data ----
# Raw data should be a matrix or data frame with a row 
# for each observation/participant,
# and a column for each variable
message("\nLoading original dataset.")
load(paste0(data_dir, "/nda18_2.01_COG_CBCL_SOC_SMAdata.RData"))
message("Original dataset loaded.\n")

Y <- as.matrix(MYdf)
if(sum(is.na(Y)) > 0){
  Y <- knnImputation(Y)
}
Y_grouped <- list()

# ---- Setting up the blocks ----  
# Labels for different variable sets:
cbclabels <- c('cbcl_scr_syn_anxdep_t', 'cbcl_scr_syn_withdep_t', 'cbcl_scr_syn_somatic_t', 
               'cbcl_scr_syn_social_t', 'cbcl_scr_syn_thought_t', 'cbcl_scr_syn_attention_t', 
               'cbcl_scr_syn_rulebreak_t', 'cbcl_scr_syn_aggressive_t', 'cbcl_scr_syn_internal_t',
               'cbcl_scr_syn_external_t', 'cbcl_scr_syn_totprob_t')

coglabels <- c('nihtbx_picvocab_agecorrected', 'nihtbx_flanker_agecorrected', 'nihtbx_list_agecorrected', 
               'nihtbx_cardsort_agecorrected', 'nihtbx_pattern_agecorrected', 'nihtbx_picture_agecorrected', 
               'nihtbx_reading_agecorrected', 'nihtbx_fluidcomp_agecorrected', 'nihtbx_cryst_agecorrected', 
               'nihtbx_totalcomp_agecorrected', 'pea_ravlt_sd_trial_vi_tc', 'pea_ravlt_ld_trial_vii_tc', 
               'pea_wiscv_tss')

soclabels <- c('fes_ss_fc_p', 'fes_y_ss_fc', 'crpbi_ss_studycaregiver', 'crpbi_y_ss_caregiver', 
               'prosocial_ss_mean_p', 'prosocial_ss_mean_y', 'pmq_y_ss_mean')

smalabels  <- c('screentime_wkdy_1', 'screentime_wkdy_2', 'screentime_wkdy_3', 'screentime_wkdy_4', 
                'screentime_wkdy_5', 'screentime_wkdy_6', 'screentime_wknd_7', 'screentime_wknd_8', 
                'screentime_wknd_9', 'screentime_wknd_10', 'screentime_wknd_11', 'screentime_wknd_12')

Y_grouped$cbc <- Y[,cbclabels]
Y_grouped$cog <- Y[,coglabels]
Y_grouped$soc <- Y[,soclabels]
Y_grouped$sma <- Y[,smalabels]

# Normalize the data
mynorm <- normalizeData(Y_grouped, type="scaleFeatures")
Y_norm <- mynorm$train # pull out the normalized data
Y_norm_bound <- do.call(cbind, Y_norm) # bind the groups into matrix


block.names <- c("cbc", "cog", "soc", "sma")
varIdx.by.block <- list(which(is.element(colnames(Y_norm_bound), cbclabels)),
                        which(is.element(colnames(Y_norm_bound), coglabels)),
                        which(is.element(colnames(Y_norm_bound), soclabels)),
                        which(is.element(colnames(Y_norm_bound), smalabels)))

block.length <- unlist(lapply(varIdx.by.block, length))
n.blocks <- length(block.length)

block4vars <-unlist(
  lapply(1:n.blocks, function(x){
    rep(block.names[x], block.length[x])
  }))

allvarlabs <- c(cbclabels, coglabels, soclabels, smalabels)

save(block.names, file = paste0(folder, "/block.names.rda"))
save(varIdx.by.block, file = paste0(folder, "/varIdx.by.block.rda"))
save(block4vars, file = paste0(folder, "/block4vars.rda"))

# Set the default K (the initial guess to the number of factors) 
# to the number of variables in the dataset.
startK <- dim(Y)[2]

# ---- Overwrite  settings ---
# Indicate if you would like to overwrite files if they already exist.
# This is useful for testing, if you made a mistake and need to
# re-run a step in the process.
overwrite_rep <- F
overwrite_gfaList <- F
overwrite_xw <- F
overwrite_match <- F

# ---- Create and load replicates ---- 

# Number of GFA replicates
R <- 10

# Writes output from gfa function to a .rda file.
foreach (r = 1:R) %dopar% {
  this.filename <- paste0(data_dir, "/GFA_rep_", r, ".rda")
  if(!file.exists(this.filename) | overwrite_rep){
    message(paste0("Creating replicate ", r, " of ", R))
    set.seed(r)
    res <- gfa(Y_norm, opts=opts, K=startK)
    save(res, file = this.filename)
    remove(res)
  }
}

# Load in the generated GFA results. Recorded in the data.frame `rep.summ` are the values
# - conv: An estimate of the convergence of the model's reconstruction based on Geweke diagnostic.
# Values significantly above 0.05 imply a non-converged model,
# and hence the need for a longer sampling chain.
# - K: The number of components inferred.
rep.summ <- data.frame(Replicate=1:R, conv=rep(NA, R), K=rep(NA, R))
gfaList_full <- list()
for(r in 1:R){
  message(paste0("Loading replicate ", r, " of ", R))
  load(paste0(data_dir, "/GFA_rep_", r, ".rda"))
  gfaList_full[[r]] <- res
  rep.summ$conv[r] <- res$conv
  rep.summ$K[r] <- res$K
  message(paste0("Replicate ", r, " has ", rep.summ$K[r], " factors."))
}
remove(res)

message("\nAll replicates loaded.\n")


# ---- Analysis ---- 

# gfaList
gfaList.filename <- paste0(folder, "/gfaList_p50.rda")
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

# xw
xw_by_rep_comp.filename <- paste0(folder, "/xw_by_rep_comp.rda")
if(!file.exists(xw_by_rep_comp.filename) | overwrite_xw){
  message("\nCreating xw_by_rep_comp.")
  xw <- list()
  for(r in 1:R){
    message(paste(c("Loading|", rep("-", r), rep(" ", R-r), "|"), collapse = ""))
    load(paste0(data_dir, "/GFA_rep_", r, ".rda"))
    xw[[r]] <- pmXW_by_factor(res)
  }
  
  save(xw, file = xw_by_rep_comp.filename)
  message("xw_by_rep_comp created.")
}else{
  message("\nLoading xw_by_rep_comp.")
  load(xw_by_rep_comp.filename)
  message("xw_by_rep_comp loaded.")
}

# corGrids and matchGrids are the threshold parameters.
# corThr: How close two components are required to be, 
# in terms of correlation, in order to match them.
corGrids   <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))
# matchThr: The proportion of sampling chains that need to contain a component 
# in order to include it in the robust components.
matchGrids <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))

match.xw.filename <- paste0(folder, "/match.xw_pm.rda")
if(!file.exists(match.xw.filename) | overwrite_xw){
  message("\nCreating match.xw")
  
  match.xw <- match_dopar(comps=xw, corGrids, matchGrids)
  
  save(match.xw, file = match.xw.filename)
  message("match.xw created.")
}else{
  message("\nLoading match.xw.")
  load(match.xw.filename)
  message("match.xw loaded.")
}

if(!file.exists(match.filename) || overwrite_match){
  # 3. Run MSE.Grids
  match.mse <- MSE.Grids(
    Ymtx=Y_norm_bound,            ## the observed (normalized) data matrix (N x D)
    comps=xw,                     ## a list of GFA replicates with posterior medians
    match.res = match.xw,         ## use the matching results from match_dopar()
    corGrids=corGrids,            ## the grids of corThr values to be assessed
    matchGrids=matchGrids)        ## the grids of matchThr values to be assessed
  save(match.mse, file = match.filename)
  
}else{
  load(match.filename)
}

# ---- Print the results ---- 
tmp.filename <- paste0(folder, "/opt.par.rda")
message(tmp.filename)
opt.par <- optimizeK(K.grids=match.mse$K.grid, mse.array=match.mse$mse$all)
save(opt.par, file = tmp.filename)

opt.cor <- opt.par$par.1se[1, "opt.corThr"]
opt.match <- opt.par$par.1se[1, "opt.matchThr"]
opt.K <- opt.par$par.1se[1, "optK"]
opt.cor.idx <- which(corGrids == opt.cor)
opt.match.idx <- which(matchGrids == opt.match)

Kgrids.filename <- paste0(folder, "/Kgrids.rda")
message(Kgrids.filename)
Kgrids <- match.mse$K.grids
save(Kgrids, file = Kgrids.filename)

optParams.filename <- paste0(folder, "/optParams.rda")
message(optParams.filename)
optParams <- Kgrids
optParams[opt.par$mse.m > opt.par$mseThr | Kgrids != opt.par$Krobust.1se] <- NA
save(optParams, file = optParams.filename)

varexp.filename <- paste0(folder, "/varexp.rda")
message(varexp.filename)

rob.ind <- list(indices=match.xw$corThr_0.1$matchThr_0.5$indices) 
varexp <- rob.var.exp(models=gfaList_p50,
                      indices=rob.ind,
                      block.names=block.names,
                      varIdx.by.block=varIdx.by.block,
                      use.unmatched=T,
                      by.block=T)
save(varexp, file = varexp.filename)

robust.xw.filename <- paste0(folder, "/robust.xw.rda")
message(robust.xw.filename)
robust.xw <- rob_wx(models=gfaList_p50, indices=varexp$indices, block.labs=block4vars, var.labs=allvarlabs)
save(robust.xw, file=robust.xw.filename)
# ---- Robust Components ----

rcomp.filename <- paste0(folder, "/rcomp.rda")
message(rcomp.filename)
if(!file.exists(rcomp.filename)){
  rcomp <- robustComponents(gfaList_full, 
                            corThr=opt.par$par.1se[1, "opt.corThr"], 
                            matchThr=opt.par$par.1se[1, "opt.matchThr"]
  )
  save(rcomp, file = rcomp.filename)
}else{
  load(rcomp.filename)
}

message("Creating report.")
rmarkdown::render("createReport.Rmd", output_file = "Report.pdf", clean = T)



message("Done.")
message(Sys.time())
