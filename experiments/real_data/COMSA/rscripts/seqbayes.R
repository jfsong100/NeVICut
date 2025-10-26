

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


### modular calibration ----
vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                              colnames(comsadeaths_wo))),
                                                list(vaalgo)),
                             missmat_type = "samples",
                             missmat = setNames(list(Mmat_wo),
                                             list(vaalgo)),
                             donotcalib_type = "fixed",
                             donotcalib = setNames(list(NULL),
                                                   list(vaalgo)),
                             nMCMC = 1000000, nBurn = 10000,
                             path_correction = F, pshrink_strength = 4,
                             saveoutput = T,
                             output_dir = output.path,
                             output_filename = paste0(cohort, '_', vaalgo, '_seqbayes'))


## insilicova ----
vaalgo = 'insilicova'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### modular calibration ----
vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                              colnames(comsadeaths_wo))),
                                                list(vaalgo)),
                             missmat_type = "samples",
                             missmat = setNames(list(Mmat_wo),
                                                list(vaalgo)),
                             donotcalib_type = "fixed",
                             donotcalib = setNames(list(NULL),
                                                   list(vaalgo)),
                             nMCMC = 1000000, nBurn = 10000,
                             path_correction = F, pshrink_strength = 4,
                             saveoutput = T,
                             output_dir = output.path,
                             output_filename = paste0(cohort, '_', vaalgo, '_seqbayes'))


## interva ----
vaalgo = 'interva'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### modular calibration ----
vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                              colnames(comsadeaths_wo))),
                                                list(vaalgo)),
                             missmat_type = "samples",
                             missmat = setNames(list(Mmat_wo),
                                                list(vaalgo)),
                             donotcalib_type = "fixed",
                             donotcalib = setNames(list(NULL),
                                                   list(vaalgo)),
                             nMCMC = 1000000, nBurn = 10000,
                             path_correction = F, pshrink_strength = 4,
                             saveoutput = T,
                             output_dir = output.path,
                             output_filename = paste0(cohort, '_', vaalgo, '_seqbayes'))



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


### modular calibration ----
vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                              colnames(comsadeaths_wo))),
                                                list(vaalgo)),
                             missmat_type = "samples",
                             missmat = setNames(list(Mmat_wo),
                                                list(vaalgo)),
                             donotcalib_type = "fixed",
                             donotcalib = setNames(list(NULL),
                                                   list(vaalgo)),
                             nMCMC = 1000000, nBurn = 10000,
                             path_correction = F, pshrink_strength = 4,
                             saveoutput = T,
                             output_dir = output.path,
                             output_filename = paste0(cohort, '_', vaalgo, '_seqbayes'))


## insilicova ----
vaalgo = 'insilicova'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### modular calibration ----
vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                              colnames(comsadeaths_wo))),
                                                list(vaalgo)),
                             missmat_type = "samples",
                             missmat = setNames(list(Mmat_wo),
                                                list(vaalgo)),
                             donotcalib_type = "fixed",
                             donotcalib = setNames(list(NULL),
                                                   list(vaalgo)),
                             nMCMC = 1000000, nBurn = 10000,
                             path_correction = F, pshrink_strength = 4,
                             saveoutput = T,
                             output_dir = output.path,
                             output_filename = paste0(cohort, '_', vaalgo, '_seqbayes'))


## interva ----
vaalgo = 'interva'

### misclassification samples ----
Mmat_wo = readRDS(file.path(data.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))
dim(Mmat_wo)
dimnames(Mmat_wo)


### modular calibration ----
vacalibration::vacalibration(va_data = setNames(list(setNames(as.numeric(comsadeaths_wo[vaalgo,]),
                                                              colnames(comsadeaths_wo))),
                                                list(vaalgo)),
                             missmat_type = "samples",
                             missmat = setNames(list(Mmat_wo),
                                                list(vaalgo)),
                             donotcalib_type = "fixed",
                             donotcalib = setNames(list(NULL),
                                                   list(vaalgo)),
                             nMCMC = 1000000, nBurn = 10000,
                             path_correction = F, pshrink_strength = 4,
                             saveoutput = T,
                             output_dir = output.path,
                             output_filename = paste0(cohort, '_', vaalgo, '_seqbayes'))


