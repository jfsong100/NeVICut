rm(list = ls()) # clear environment
current_path <- rstudioapi::getSourceEditorContext()$path
base_dir <- dirname(current_path)
output.path <- file.path(base_dir, "upstream")


Y11 = c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
Y12 = c(111, 71, 162, 188, 145, 215, 166, 37, 173, 143, 229, 696, 93)
Y21 = c(16, 215, 362, 97, 76, 62, 710, 56, 133, 28, 62, 413, 194)
Y22 = c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218)
Y22 = Y22/1000

shape1 = 1 + Y11
shape2 = 1 + Y12 - Y11

n = length(Y11)
N_MC = 1000
sample_mat = vector()

set.seed(1)
for (i in 1:n) {
  sample_mat = cbind(sample_mat,
                     rbeta(N_MC, shape1 = shape1[i], shape2 = shape2[i]))
}

colnames(sample_mat) = paste0("country",1:n)
row.names(sample_mat) = paste0("sample",1:N_MC)

saveRDS(sample_mat, file = file.path(output.path,"theta1_samples.rds"))

library(R.matlab)
writeMat(file.path(output.path,"theta1_samples.mat"), theta1 = sample_mat)
