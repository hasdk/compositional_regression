# Script to estimate and impute cell type proportions, and run compositional regression on example data 
# In this example, we estimate the expected change in placenta cell type composition associated with 1 week increase
# in gestational age, along with the corresponding bootstrapped confidence intervals and corrected p-values.

# Author: Hachem Saddiki

# import required libraries
library(planet)
library(EpiDISH)
library(compositions)
library(robCompositions)
library(data.table)
library(parallel)

set.seed(seed=9871)

# load example data from planet package
data(plBetas)
data(plPhenoData)

############################################################################
## (1) Estimate placenta cell proportions from example methylation data
############################################################################

# load third trimester methylation placenta reference from planet package
data('plCellCpGsThird') 

# estimate cell type proportions using Houseman's constrained projection (CP) 
cell_comp <- EpiDISH::epidish(beta.m=plBetas, ref.m = plCellCpGsThird,
                              method="CP")
cell_comp <- data.frame(cell_comp$estF)

# round proportions less than pre-specified cutoff to zero and 
# store detection limit for imputation
cutoff <- 10^-4
DL <- matrix(0,ncol=ncol(cell_comp),nrow=1) 
for(c in 1:ncol(cell_comp)){
  zeros_idx <- which(cell_comp[,c] < cutoff)
  if(length(zeros_idx) > 0){
    # round zeros, set DL to smallest non-zero value after rounding
    cell_comp[zeros_idx,c] <- 0
    DL[c] <- min(cell_comp[-zeros_idx,c])
  }
}

# impute rounded zeros using detection limit
imp <- robCompositions::imputeBDLs(cell_comp, dl=DL, 
                  maxit=50,eps=0.1,R=50,method="pls", variation=F)
cell_comp_imp <- imp$x
cell_comp_imp$sample_id <- rownames(cell_comp)

# construct data set for compositional regression
comp_data <- merge(plPhenoData, cell_comp_imp, by='sample_id')


######################################################################################################
## (2) Perform compositional regression using isometric log-ratio transformation of cell composition
# Outcome: ilr transformed placenta cell type proportions
# Covariates: infant sex and gestational weeks
######################################################################################################
cell_types <- colnames(cell_comp)
outcome <- compositions::acomp(comp_data[,cell_types])
comp_mod <- lm(ilr(outcome) ~ sex + gestation_wk, data=comp_data)

# extract coefficients in inverse-ILR scale
compositions::ilrInv(coef(comp_mod), orig=outcome)

# extract pvalue associated with gestational week
anova(comp_mod)['gestation_wk',6]

## NOTE: One of the limitations of the 'compositions' package is that the 
##       p-value calculation is only correct for the last covariate in the 
##       model formula.  
##       In this example, the p-value for "gestation_wk" is correct, but not 
##       for infant sex. In order to calculate the correct p-value for "sex", 
##       you need to refit the same model with "sex" at the end of the formula: 
##      --> "ilr(outcome) ~ gestation_wk + sex"
## Reference: for further details, see Section 5.3.5 of Van den Boogaart,
##            Tolsana and Delgado (2013).

# calculate expected cell composition for 1 week increase in gestational week
gest_week_coef <- compositions::ilrInv(coef(comp_mod), orig=outcome)['gestation_wk',]
exp_ccomp <- mean(compositions::acomp(gest_week_coef) + outcome)
exp_ccomp

# calculate expected change in cell composition for 1 week increase in gestational week
exp_change <- as.numeric(exp_ccomp) - as.numeric(mean(outcome))
names(exp_change) <- names(outcome)
exp_change

################################################################################
# (3) Bootstrap approach to estimate confidence intervals and corrected p-value
#     for the effect of increasing gestational age by 1 week on placenta cell
#     type composition.
################################################################################

nboot = 100 # number of bootstrap repetitions
no_cores = detectCores() - 1 # number of computer cores to use for parallel loop 
cl <- makeCluster(no_cores) # create parallel backend
clusterEvalQ(cl, {
  library(compositions)
  library(data.table)
  library(parallel)
})

clusterExport(cl, varlist = c("nboot","comp_data","cell_types"))

# parallel loop that runs compositional regressions for each bootstrap data set and stores results
tmp = parLapplyLB(cl,1:nboot,function(b){
  
  # sample from original data with replacement
  rep_idx <- sample(x=1:nrow(comp_data), size=nrow(comp_data), replace=T)
  boot_data <- comp_data[rep_idx,]
  
  # perform analysis on bootstrap data and save estimates
  outcome_b = acomp(boot_data[,cell_types])

  # fit compositional model on current bootstrap data
  comp_mod_b <- lm(ilr(outcome_b) ~ sex + gestation_wk, data=boot_data)
  
  # extract F statistic
  Fstat_b = anova(comp_mod_b)['gestation_wk', 3]
  
  # extract exposure coefficient in inverse ILR scale
  coef_b <- ilrInv(coef(comp_mod_b), orig=outcome_b)['gestation_wk',]
  
  # calculate the expected cell composition after 1 unit increase in exposure
  exp_outcome = mean(acomp(coef_b) + outcome_b)
  
  # expected change in cell composition after 1 unit increase in exposure (original scale)
  exp_outcome_change = as.numeric(exp_outcome) - as.numeric(mean(outcome_b))
  
  # store results
  dt <- data.table(t(exp_outcome_change), 
                   Fstat = Fstat_b)
  
})
stopCluster(cl) # stop parallel backend

# extract bootstrap results
bootstrap_results = rbindlist(tmp)
colnames(bootstrap_results)[1:length(cell_types)] <- cell_types

# calculate 95% bootstrap Confidence intervals
boot_CIs <- as.data.frame(apply(bootstrap_results[,-7], 2, 
                              quantile, probs=c(0.025,0.975)))

# calculate bootstrap corrected p-value using the median F-statistic 
boot_pvalue <- pf(median(bootstrap_results$Fstat), 5, 17, lower.tail=F)

