

rm(list = ls()) # clear environment


# comsa-mozambique data ----
# install vacalibration. the ccva outputs of the data are in this package
# remotes::install_github("https://github.com/sandy-pramanik/vacalibration")
library(vacalibration)  # load

## where to save the prepared data
current_path <- rstudioapi::getSourceEditorContext()$path
base_dir <- dirname(dirname(current_path))
output.path <- file.path(base_dir, "comsa_data")
if(!dir.exists(output.path)){
  
  dir.create(output.path)
  
}


## neonate ----
cohort = 'neonate'

comsadeaths_wo = NULL
for(vaalgo in c('eava', 'insilicova', 'interva')){
  
  # mapping to broad causes
  comsa_broad = vacalibration::cause_map(df = vacalibration::comsamoz_CCVAoutput[[cohort]][[vaalgo]],
                                         age_group = cohort)
  head(comsa_broad)
  
  # removing 'other' cause and calculating cause-specific death counts
  comsa_broad_wo = comsa_broad[,colnames(comsa_broad)!="other"]
  head(comsa_broad_wo)
  
  # removing 'other' cause and calculating cause-specific death counts
  comsadeaths_wo = rbind(comsadeaths_wo, colSums(comsa_broad_wo))
  
}

comsadeaths_wo = as.data.frame(comsadeaths_wo)
rownames(comsadeaths_wo) = c('eava', 'insilicova', 'interva')
comsadeaths_wo

# saving
write.table(comsadeaths_wo, file.path(output.path, paste0(cohort, '.csv')))


## child ----
cohort = 'child'

comsadeaths_wo = NULL
for(vaalgo in c('eava', 'insilicova', 'interva')){
  
  # mapping to broad causes
  comsa_broad = vacalibration::cause_map(df = vacalibration::comsamoz_CCVAoutput[[cohort]][[vaalgo]],
                                         age_group = cohort)
  head(comsa_broad)
  
  # removing 'other' cause and calculating cause-specific death counts
  comsa_broad_wo = comsa_broad[,colnames(comsa_broad)!="other"]
  head(comsa_broad_wo)
  
  # removing 'other' cause and calculating cause-specific death counts
  comsadeaths_wo = rbind(comsadeaths_wo, colSums(comsa_broad_wo))
  
}

comsadeaths_wo = as.data.frame(comsadeaths_wo)
rownames(comsadeaths_wo) = c('eava', 'insilicova', 'interva')
comsadeaths_wo

# saving
write.table(comsadeaths_wo, file.path(output.path, paste0(cohort, '.csv')))





# ccva misclassification matrix samples ----

# download from: https://github.com/sandy-pramanik/CCVA-Misclassification-Matrices/releases/download/20241004/CCVA_missmat
# put in the working directory

ccva_missmat = readRDS(file.path(output.path,"CCVA_missmat"))
country = "Mozambique"  # to extract misclassification matrix. since we're using data comsa-mozambique.
nSamp.mmat = 1000 # how many samples to use. repository provides 5000 samples. extracts last nSamp.mmat samples.


## neonate ----
cohort = 'neonate'


### eava ----
vaalgo = 'eava'

# extract
mmat = ccva_missmat[[cohort]][[vaalgo]]$postsamples[[country]]
class(mmat)
dim(mmat)
dimnames(mmat)

# causes to use
(causes_sub = dimnames(mmat)[[2]][dimnames(mmat)[[2]]!="other"])

# removing "other" cause
sample_id = tail(1:dim(mmat)[1], nSamp.mmat)  # subset nSamp.mmat misclassification samples
Mmat_wo = array(dim = c(nSamp.mmat, length(causes_sub), length(causes_sub)),
                dimnames = list(NULL, causes_sub, causes_sub))
for(r in 1:nSamp.mmat){
  
  Mmat_wo[r,,] = mmat[sample_id[r],causes_sub,causes_sub]
  Mmat_wo[r,,] = Mmat_wo[r,,]/rowSums(Mmat_wo[r,,])
  
}

# checking
range(apply(Mmat_wo, 1:2, sum))  # each row sum to 1 for each sample
dim(Mmat_wo)  # dimension
dimnames(Mmat_wo)  # dimension names

saveRDS(Mmat_wo, file.path(output.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))


### insilicova ----
vaalgo = 'insilicova'

# extract
mmat = ccva_missmat[[cohort]][[vaalgo]]$postsamples[[country]]
class(mmat)
dim(mmat)
dimnames(mmat)

# removing "other" cause
Mmat_wo = array(dim = c(nSamp.mmat, length(causes_sub), length(causes_sub)),
                dimnames = list(NULL, causes_sub, causes_sub))
for(r in 1:nSamp.mmat){
  
  Mmat_wo[r,,] = mmat[sample_id[r],causes_sub,causes_sub]
  Mmat_wo[r,,] = Mmat_wo[r,,]/rowSums(Mmat_wo[r,,])
  
}

# checking
range(apply(Mmat_wo, 1:2, sum))  # each row sum to 1 for each sample
dim(Mmat_wo)  # dimension
dimnames(Mmat_wo)  # dimension names

saveRDS(Mmat_wo, file.path(output.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))


### interva ----
vaalgo = 'interva'

# extract
mmat = ccva_missmat[[cohort]][[vaalgo]]$postsamples[[country]]
class(mmat)
dim(mmat)
dimnames(mmat)

# removing "other" cause
Mmat_wo = array(dim = c(nSamp.mmat, length(causes_sub), length(causes_sub)),
                dimnames = list(NULL, causes_sub, causes_sub))
for(r in 1:nSamp.mmat){
  
  Mmat_wo[r,,] = mmat[sample_id[r],causes_sub,causes_sub]
  Mmat_wo[r,,] = Mmat_wo[r,,]/rowSums(Mmat_wo[r,,])
  
}

# checking
range(apply(Mmat_wo, 1:2, sum))  # each row sum to 1 for each sample
dim(Mmat_wo)  # dimension
dimnames(Mmat_wo)  # dimension names

saveRDS(Mmat_wo, file.path(output.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))


## child ----
cohort = 'child'


### eava ----
vaalgo = 'eava'

# extract
mmat = ccva_missmat[[cohort]][[vaalgo]]$postsamples[[country]]
class(mmat)
dim(mmat)
dimnames(mmat)

# causes to use
(causes_sub = dimnames(mmat)[[2]][dimnames(mmat)[[2]]!="other"])

# removing "other" cause
Mmat_wo = array(dim = c(nSamp.mmat, length(causes_sub), length(causes_sub)),
                dimnames = list(NULL, causes_sub, causes_sub))
for(r in 1:nSamp.mmat){
  
  Mmat_wo[r,,] = mmat[sample_id[r],causes_sub,causes_sub]
  Mmat_wo[r,,] = Mmat_wo[r,,]/rowSums(Mmat_wo[r,,])
  
}

# checking
range(apply(Mmat_wo, 1:2, sum))  # each row sum to 1 for each sample
dim(Mmat_wo)  # dimension
dimnames(Mmat_wo)  # dimension names

saveRDS(Mmat_wo, file.path(output.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))


### insilicova ----
vaalgo = 'insilicova'

# extract
mmat = ccva_missmat[[cohort]][[vaalgo]]$postsamples[[country]]
class(mmat)
dim(mmat)
dimnames(mmat)

# removing "other" cause
Mmat_wo = array(dim = c(nSamp.mmat, length(causes_sub), length(causes_sub)),
                dimnames = list(NULL, causes_sub, causes_sub))
for(r in 1:nSamp.mmat){
  
  Mmat_wo[r,,] = mmat[sample_id[r],causes_sub,causes_sub]
  Mmat_wo[r,,] = Mmat_wo[r,,]/rowSums(Mmat_wo[r,,])
  
}

# checking
range(apply(Mmat_wo, 1:2, sum))  # each row sum to 1 for each sample
dim(Mmat_wo)  # dimension
dimnames(Mmat_wo)  # dimension names

saveRDS(Mmat_wo, file.path(output.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))


### interva ----
vaalgo = 'interva'

# extract
mmat = ccva_missmat[[cohort]][[vaalgo]]$postsamples[[country]]
class(mmat)
dim(mmat)
dimnames(mmat)

# removing "other" cause
Mmat_wo = array(dim = c(nSamp.mmat, length(causes_sub), length(causes_sub)),
                dimnames = list(NULL, causes_sub, causes_sub))
for(r in 1:nSamp.mmat){
  
  Mmat_wo[r,,] = mmat[sample_id[r],causes_sub,causes_sub]
  Mmat_wo[r,,] = Mmat_wo[r,,]/rowSums(Mmat_wo[r,,])
  
}

# checking
range(apply(Mmat_wo, 1:2, sum))  # each row sum to 1 for each sample
dim(Mmat_wo)  # dimension
dimnames(Mmat_wo)  # dimension names

saveRDS(Mmat_wo, file.path(output.path, paste0("Mmat_wo_", cohort, "_", vaalgo)))


