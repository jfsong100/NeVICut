

# load library ----
rm(list = ls())

library(vacalibration)


# path to data and output
current_path <- rstudioapi::getSourceEditorContext()$path
base_dir <- dirname(dirname(current_path))
data.path <- file.path(base_dir, "comsa_data")

output.path = file.path(base_dir, "downstream_results")
if(!dir.exists(output.path)){
  
  dir.create(output.path)
  
}



# neonate ----
cohort = 'neonate'

# unlabeled data
comsadeaths_wo <- read.table(file.path(data.path, paste0(cohort, '.csv')))
comsadeaths_wo


## eava ----
vaalgo = 'eava'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### multiple imputation ----
# misclassification sample size from module 1
# nSamp.Mmat = 10
(nSamp.Mmat = dim(Mmat_wo)[1])

# calibrated samples from module 2
# nSamp.pcalib = 10
# nSamp.pcalib = 100
nSamp.pcalib = 1000

pcalib.multimpute = array(dim = c(nSamp.Mmat, nSamp.pcalib, length(comsadeaths_wo)),
                          dimnames = list(NULL, NULL, names(comsadeaths_wo)))
ptm0 = proc.time()
pb = txtProgressBar(min = 1, max = nSamp.Mmat, style = 3)
for(r in 1:nSamp.Mmat){

  # r = 1
  
  multimpute.out_r = vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                                                   colnames(comsadeaths_wo))),
                                                                     list(vaalgo)),
                                                  missmat_type = "fixed",
                                                  missmat = setNames(list(Mmat_wo[r,,]),
                                                                  list(vaalgo)),
                                                  donotcalib_type = "fixed",
                                                  donotcalib = setNames(list(NULL),
                                                                        list(vaalgo)),
                                                  nMCMC = nSamp.pcalib, nBurn = 10000,
                                                  path_correction = F, pshrink_strength = 4,
                                                  verbose = F, plot_it = F)
  
  dim(multimpute.out_r$p_calib)
  pcalib.multimpute[r,,] = multimpute.out_r$p_calib[1,,]

  setTxtProgressBar(pb, r)

}
ptm1 = proc.time()

class(pcalib.multimpute)
dim(pcalib.multimpute)
dimnames(pcalib.multimpute)

# pcalib.multimpute_melt = reshape2::melt(pcalib.multimpute, varnames = 2:3)
# class(pcalib.multimpute_melt)
# dim(pcalib.multimpute_melt)
# dimnames(pcalib.multimpute_melt)
# vioplot::vioplot(pcalib.multimpute, ylim = c(0,1))
# points(deathcount/sum(deathcount), pch = 8, col = 4)

multimpute.out = list('pcalib' = pcalib.multimpute,
                      'puncalib' = comsadeaths_wo/sum(comsadeaths_wo),
                      'runtime' = ptm1 - ptm0)

saveRDS(multimpute.out,
        file.path(output.path,
                  paste0(cohort, '_', vaalgo, '_multimpute')))


## insilicova ----
vaalgo = 'insilicova'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### multiple imputation ----
# misclassification sample size from module 1
# nSamp.Mmat = 10
(nSamp.Mmat = dim(Mmat_wo)[1])

# calibrated samples from module 2
# nSamp.pcalib = 10
# nSamp.pcalib = 100
nSamp.pcalib = 1000

pcalib.multimpute = array(dim = c(nSamp.Mmat, nSamp.pcalib, length(comsadeaths_wo)),
                          dimnames = list(NULL, NULL, names(comsadeaths_wo)))
ptm0 = proc.time()
pb = txtProgressBar(min = 1, max = nSamp.Mmat, style = 3)
for(r in 1:nSamp.Mmat){
  
  # r = 1
  
  multimpute.out_r = vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                                                   colnames(comsadeaths_wo))),
                                                                     list(vaalgo)),
                                                  missmat_type = "fixed",
                                                  missmat = setNames(list(Mmat_wo[r,,]),
                                                                     list(vaalgo)),
                                                  donotcalib_type = "fixed",
                                                  donotcalib = setNames(list(NULL),
                                                                        list(vaalgo)),
                                                  nMCMC = nSamp.pcalib, nBurn = 10000,
                                                  path_correction = F, pshrink_strength = 4,
                                                  verbose = F, plot_it = F)
  
  dim(multimpute.out_r$p_calib)
  pcalib.multimpute[r,,] = multimpute.out_r$p_calib[1,,]
  
  setTxtProgressBar(pb, r)
  
}
ptm1 = proc.time()

class(pcalib.multimpute)
dim(pcalib.multimpute)
dimnames(pcalib.multimpute)

# pcalib.multimpute_melt = reshape2::melt(pcalib.multimpute, varnames = 2:3)
# class(pcalib.multimpute_melt)
# dim(pcalib.multimpute_melt)
# dimnames(pcalib.multimpute_melt)
# vioplot::vioplot(pcalib.multimpute, ylim = c(0,1))
# points(deathcount/sum(deathcount), pch = 8, col = 4)

multimpute.out = list('pcalib' = pcalib.multimpute,
                      'puncalib' = comsadeaths_wo/sum(comsadeaths_wo),
                      'runtime' = ptm1 - ptm0)

saveRDS(multimpute.out,
        file.path(output.path,
                  paste0(cohort, '_', vaalgo, '_multimpute')))


## interva ----
vaalgo = 'interva'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### multiple imputation ----
# misclassification sample size from module 1
# nSamp.Mmat = 10
(nSamp.Mmat = dim(Mmat_wo)[1])

# calibrated samples from module 2
# nSamp.pcalib = 10
# nSamp.pcalib = 100
nSamp.pcalib = 1000

pcalib.multimpute = array(dim = c(nSamp.Mmat, nSamp.pcalib, length(comsadeaths_wo)),
                          dimnames = list(NULL, NULL, names(comsadeaths_wo)))
ptm0 = proc.time()
pb = txtProgressBar(min = 1, max = nSamp.Mmat, style = 3)
for(r in 1:nSamp.Mmat){
  
  # r = 1
  
  multimpute.out_r = vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                                                   colnames(comsadeaths_wo))),
                                                                     list(vaalgo)),
                                                  missmat_type = "fixed",
                                                  missmat = setNames(list(Mmat_wo[r,,]),
                                                                     list(vaalgo)),
                                                  donotcalib_type = "fixed",
                                                  donotcalib = setNames(list(NULL),
                                                                        list(vaalgo)),
                                                  nMCMC = nSamp.pcalib, nBurn = 10000,
                                                  path_correction = F, pshrink_strength = 4,
                                                  verbose = F, plot_it = F)
  
  dim(multimpute.out_r$p_calib)
  pcalib.multimpute[r,,] = multimpute.out_r$p_calib[1,,]
  
  setTxtProgressBar(pb, r)
  
}
ptm1 = proc.time()

class(pcalib.multimpute)
dim(pcalib.multimpute)
dimnames(pcalib.multimpute)

# pcalib.multimpute_melt = reshape2::melt(pcalib.multimpute, varnames = 2:3)
# class(pcalib.multimpute_melt)
# dim(pcalib.multimpute_melt)
# dimnames(pcalib.multimpute_melt)
# vioplot::vioplot(pcalib.multimpute, ylim = c(0,1))
# points(deathcount/sum(deathcount), pch = 8, col = 4)

multimpute.out = list('pcalib' = pcalib.multimpute,
                      'puncalib' = comsadeaths_wo/sum(comsadeaths_wo),
                      'runtime' = ptm1 - ptm0)

saveRDS(multimpute.out,
        file.path(output.path,
                  paste0(cohort, '_', vaalgo, '_multimpute')))


# child ----
cohort = 'child'

# unlabeled data
comsadeaths_wo <- read.table(file.path(data.path, paste0(cohort, '.csv')))
comsadeaths_wo


## eava ----
vaalgo = 'eava'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### multiple imputation ----
# misclassification sample size from module 1
# nSamp.Mmat = 10
(nSamp.Mmat = dim(Mmat_wo)[1])

# calibrated samples from module 2
# nSamp.pcalib = 10
# nSamp.pcalib = 100
nSamp.pcalib = 1000

pcalib.multimpute = array(dim = c(nSamp.Mmat, nSamp.pcalib, length(comsadeaths_wo)),
                          dimnames = list(NULL, NULL, names(comsadeaths_wo)))
ptm0 = proc.time()
pb = txtProgressBar(min = 1, max = nSamp.Mmat, style = 3)
for(r in 1:nSamp.Mmat){
  
  # r = 1
  
  multimpute.out_r = vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                                                   colnames(comsadeaths_wo))),
                                                                     list(vaalgo)),
                                                  missmat_type = "fixed",
                                                  missmat = setNames(list(Mmat_wo[r,,]),
                                                                     list(vaalgo)),
                                                  donotcalib_type = "fixed",
                                                  donotcalib = setNames(list(NULL),
                                                                        list(vaalgo)),
                                                  nMCMC = nSamp.pcalib, nBurn = 10000,
                                                  path_correction = F, pshrink_strength = 4,
                                                  verbose = F, plot_it = F)
  
  dim(multimpute.out_r$p_calib)
  pcalib.multimpute[r,,] = multimpute.out_r$p_calib[1,,]
  
  setTxtProgressBar(pb, r)
  
}
ptm1 = proc.time()

class(pcalib.multimpute)
dim(pcalib.multimpute)
dimnames(pcalib.multimpute)

# pcalib.multimpute_melt = reshape2::melt(pcalib.multimpute, varnames = 2:3)
# class(pcalib.multimpute_melt)
# dim(pcalib.multimpute_melt)
# dimnames(pcalib.multimpute_melt)
# vioplot::vioplot(pcalib.multimpute, ylim = c(0,1))
# points(deathcount/sum(deathcount), pch = 8, col = 4)

multimpute.out = list('pcalib' = pcalib.multimpute,
                      'puncalib' = comsadeaths_wo/sum(comsadeaths_wo),
                      'runtime' = ptm1 - ptm0)

saveRDS(multimpute.out,
        file.path(output.path,
                  paste0(cohort, '_', vaalgo, '_multimpute')))


## insilicova ----
vaalgo = 'insilicova'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### multiple imputation ----
# misclassification sample size from module 1
# nSamp.Mmat = 10
(nSamp.Mmat = dim(Mmat_wo)[1])

# calibrated samples from module 2
# nSamp.pcalib = 10
# nSamp.pcalib = 100
nSamp.pcalib = 1000

pcalib.multimpute = array(dim = c(nSamp.Mmat, nSamp.pcalib, length(comsadeaths_wo)),
                          dimnames = list(NULL, NULL, names(comsadeaths_wo)))
ptm0 = proc.time()
pb = txtProgressBar(min = 1, max = nSamp.Mmat, style = 3)
for(r in 1:nSamp.Mmat){
  
  # r = 1
  
  multimpute.out_r = vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                                                   colnames(comsadeaths_wo))),
                                                                     list(vaalgo)),
                                                  missmat_type = "fixed",
                                                  missmat = setNames(list(Mmat_wo[r,,]),
                                                                     list(vaalgo)),
                                                  donotcalib_type = "fixed",
                                                  donotcalib = setNames(list(NULL),
                                                                        list(vaalgo)),
                                                  nMCMC = nSamp.pcalib, nBurn = 10000,
                                                  path_correction = F, pshrink_strength = 4,
                                                  verbose = F, plot_it = F)
  
  dim(multimpute.out_r$p_calib)
  pcalib.multimpute[r,,] = multimpute.out_r$p_calib[1,,]
  
  setTxtProgressBar(pb, r)
  
}
ptm1 = proc.time()

class(pcalib.multimpute)
dim(pcalib.multimpute)
dimnames(pcalib.multimpute)

# pcalib.multimpute_melt = reshape2::melt(pcalib.multimpute, varnames = 2:3)
# class(pcalib.multimpute_melt)
# dim(pcalib.multimpute_melt)
# dimnames(pcalib.multimpute_melt)
# vioplot::vioplot(pcalib.multimpute, ylim = c(0,1))
# points(deathcount/sum(deathcount), pch = 8, col = 4)

multimpute.out = list('pcalib' = pcalib.multimpute,
                      'puncalib' = comsadeaths_wo/sum(comsadeaths_wo),
                      'runtime' = ptm1 - ptm0)

saveRDS(multimpute.out,
        file.path(output.path,
                  paste0(cohort, '_', vaalgo, '_multimpute')))


## interva ----
vaalgo = 'interva'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### multiple imputation ----
# misclassification sample size from module 1
# nSamp.Mmat = 10
(nSamp.Mmat = dim(Mmat_wo)[1])

# calibrated samples from module 2
# nSamp.pcalib = 10
# nSamp.pcalib = 100
nSamp.pcalib = 1000

pcalib.multimpute = array(dim = c(nSamp.Mmat, nSamp.pcalib, length(comsadeaths_wo)),
                          dimnames = list(NULL, NULL, names(comsadeaths_wo)))
ptm0 = proc.time()
pb = txtProgressBar(min = 1, max = nSamp.Mmat, style = 3)
for(r in 1:nSamp.Mmat){
  
  # r = 1
  
  multimpute.out_r = vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                                                   colnames(comsadeaths_wo))),
                                                                     list(vaalgo)),
                                                  missmat_type = "fixed",
                                                  missmat = setNames(list(Mmat_wo[r,,]),
                                                                     list(vaalgo)),
                                                  donotcalib_type = "fixed",
                                                  donotcalib = setNames(list(NULL),
                                                                        list(vaalgo)),
                                                  nMCMC = nSamp.pcalib, nBurn = 10000,
                                                  path_correction = F, pshrink_strength = 4,
                                                  verbose = F, plot_it = F)
  
  dim(multimpute.out_r$p_calib)
  pcalib.multimpute[r,,] = multimpute.out_r$p_calib[1,,]
  
  setTxtProgressBar(pb, r)
  
}
ptm1 = proc.time()

class(pcalib.multimpute)
dim(pcalib.multimpute)
dimnames(pcalib.multimpute)

# pcalib.multimpute_melt = reshape2::melt(pcalib.multimpute, varnames = 2:3)
# class(pcalib.multimpute_melt)
# dim(pcalib.multimpute_melt)
# dimnames(pcalib.multimpute_melt)
# vioplot::vioplot(pcalib.multimpute, ylim = c(0,1))
# points(deathcount/sum(deathcount), pch = 8, col = 4)

multimpute.out = list('pcalib' = pcalib.multimpute,
                      'puncalib' = comsadeaths_wo/sum(comsadeaths_wo),
                      'runtime' = ptm1 - ptm0)

saveRDS(multimpute.out,
        file.path(output.path,
                  paste0(cohort, '_', vaalgo, '_multimpute')))



