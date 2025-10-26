rm(list = ls())

library(dplyr)
library(rstudioapi)
current_path = rstudioapi::getSourceEditorContext()$path
base_dir = dirname(dirname(current_path))

# ------------------ Converters ------------------
convert_df_function = function(df){
  df[] = lapply(df, function(x) {
    if (is.factor(x)) as.numeric(as.character(x))
    else if (is.character(x)) as.numeric(x)
    else x
  })
  as.matrix(df)
}

array3_to_matrix = function(arr){
  stopifnot(length(dim(arr)) == 3)
  matrix(arr, nrow = prod(dim(arr)[1:2]), ncol = dim(arr)[3])
}

# -------- optional CLR transform for compositional data (probability vectors) --------
clr_transform <- function(M, eps = 1e-12) {
  M <- pmax(M, eps)                  
  gm <- exp(rowMeans(log(M)))        
  sweep(log(M), 1, log(gm), "-")     
}

# -------- Sliced W1 using quantile interpolation --------
sliced_w1_interp <- function(X, Y, n_projections = 256, grid_size = 1001,
                             compositional = TRUE, seed = 42, eps = 1e-12) {
  set.seed(seed)
  if (compositional) {
    X <- clr_transform(X, eps)
    Y <- clr_transform(Y, eps)
  }
  
  p <- ncol(X)
  ugrid <- seq(0, 1, length.out = grid_size)
  acc <- 0
  
  for (b in seq_len(n_projections)) {
    v <- rnorm(p)
    v <- v / sqrt(sum(v^2))          # unit direction
    
    xproj <- drop(X %*% v)
    yproj <- drop(Y %*% v)
    
    qx <- quantile(xproj, probs = ugrid, names = FALSE, type = 7)
    qy <- quantile(yproj, probs = ugrid, names = FALSE, type = 7)
    
    acc <- acc + mean(abs(qx - qy))  # W1 on the line
  }
  
  acc / n_projections                 # sliced W1 (average over directions)
}


# ------------------ Config ------------------
populations = c("neonate")
algorithms  = c("eava", "insilicova", "interva")
adjustments = "noadj"
methods     = c("NeVI-Cut","Multiple-Imputation",
                "Full-Bayes",
                "Dirichlet")

# ------------------ Loaders ------------------
load_MI = function(population, algorithm, adjustment){
  f = paste0(base_dir, "/downstream_results/", population, "_", algorithm,
             "_multimpute")
  obj = readRDS(f)
  arr = obj[["pcalib"]]
  if (is.null(arr)) stop("MI file lacks 'pcalib': ", f)
  array3_to_matrix(arr)
}

load_FB = function(population, algorithm){
  f = paste0(base_dir, "/downstream_results/", population, "_", algorithm,
             "_seqbayes")
  obj = readRDS(f)
  mat = obj$p_calib
  if (is.null(mat)) stop("FB file lacks 'MCMCout$p_calib': ", f)
  drop(mat) 
}

load_NeVI = function(population, algorithm, adjustment){
  f = paste0(base_dir, "/downstream_results/", "/NeVI_Cut_",
             population, "_", algorithm, ".rds")
  convert_df_function(readRDS(f))
}

load_Dir = function(population, algorithm, adjustment){
  f = paste0(base_dir, "/downstream_results/", "/Dir_",
             population, "_", algorithm, ".rds")
  convert_df_function(readRDS(f))
}

load_method = function(method, population, algorithm, adjustment){
  switch(method,
         "Multiple-Imputation" = load_MI(population, algorithm, adjustment),
         "NeVI-Cut"            = load_NeVI(population, algorithm, adjustment),
         "Dirichlet"           = load_Dir(population, algorithm, adjustment),
         "Full-Bayes"          = load_FB(population, algorithm),  
         stop("Unknown method: ", method)
  )
}

# ------------------ Main loop ------------------
run_all_wasserstein = function(n_projections = 256, grid_size = 1001, seed = 2025,
                               sliced_w_fun = sliced_w2_interp){
  out = vector("list", length = 0L)
  k = 0L
  
  for (pop in populations) {
    print(pop)
    for (alg in algorithms) {
      print(alg)
      for (adj in adjustments) {
        print(adj)
        # Reference: MI
        MI_mat = load_method("Multiple-Imputation", pop, alg, adj)
        p = ncol(MI_mat)
        
        for (m in methods) {
          print(m)
          if (m == "Multiple-Imputation") next
          
          tgt = load_method(m, pop, alg, adj)
          if (ncol(tgt) != p) {
            stop(sprintf("Column mismatch for %s (%s/%s/%s). MI has p=%d, %s has p=%d",
                         m, pop, alg, adj, p, m, ncol(tgt)))
          }
          
          d = sliced_w_fun(MI_mat, tgt,
                               n_projections = n_projections,
                               grid_size = grid_size,
                               seed = seed)
          
          k = k + 1L
          out[[k]] = data.frame(
            population = pop,
            algorithm  = alg,
            adjustment = adj,
            method     = m,
            p          = p,
            n_ref      = nrow(MI_mat),
            n_tgt      = nrow(tgt),
            sw2        = d,
            stringsAsFactors = FALSE
          )
          
          rm(tgt); gc()
        }
        
        rm(MI_mat); gc()
      }
    }
  }
  
  bind_rows(out)
}

# ------------------ Run ------------------
dist_tblw1 = run_all_wasserstein(n_projections = 256, grid_size = 1001, seed = 123, sliced_w_fun = sliced_w1_interp)

print(dist_tblw1)

write.csv(dist_tblw1, file = file.path(base_dir, "downstream_results", "wasserstein_results.csv"), row.names = FALSE)

# ------------------ Time ------------------
load_MI_time = function(population, algorithm, adjustment){
  f = paste0(base_dir, "/downstream_results/", population, "_", algorithm,
             "_multimpute")
  obj = readRDS(f)
  obj$runtime[3]
}

load_NeVI_time = function(population, algorithm, adjustment){
  f = paste0(base_dir, "/downstream_results/", "/NeVI_Cut_",
             population, "_", algorithm, "_time.txt")
  times = read.table(f, quote="\"", comment.char="")[3] %>% as.numeric()
  times
}

combinations <- expand.grid(population = populations, algorithm = algorithms)

# --- Apply both functions ---
results <- combinations %>%
  rowwise() %>%
  mutate(
    MI_time   = load_MI_time(population, algorithm),
    NeVI_time = load_NeVI_time(population, algorithm)
  ) %>%
  ungroup()

output_path <- file.path(base_dir, "downstream_results", "NeVI_times_table.csv")
write.csv(results, output_path, row.names = FALSE)


# ------------------ Plot ------------------
library(readr)
library(dplyr)
library(ggplot2)

df <- dist_tblw1
df <- df %>%
  mutate(
    scenario = paste0(algorithm, " (", adjustment, ")"),
    population = factor(population,
                        levels = c("neonate", "child"),
                        labels = c("Neonate", "Child"))   # Uppercase facet labels
  )

scenarios <- c("eava (noadj)", 
               "insilicova (noadj)", 
               "interva (noadj)")
df$scenario <- factor(df$scenario, levels = scenarios)

library(tidyverse)
dist_sub = df[,c("population","algorithm","method","sw2")]
time_sub = results[,c("population","algorithm","MI_time","NeVI_time")]
method_levels <- c("Full-Bayes","Dirichlet","NeVI-Cut","Nested MCMC")
method_labels <- c("Full Bayes","Parametric Cut","NeVI-Cut","Nested MCMC")

time_long <- time_sub %>%
  select(algorithm, MI_time, NeVI_time) %>%
  pivot_longer(cols = c(MI_time, NeVI_time),
               names_to = "method",
               values_to = "time") %>%
  mutate(method = recode(method,
                         MI_time = "Multiple-Imputation",
                         NeVI_time = "NeVI-Cut"),
         algorithm = factor(algorithm, levels = c("eava", "insilicova", "interva")))

time_long <- time_sub %>%
  pivot_longer(c(MI_time, NeVI_time), names_to="method", values_to="time") %>%
  mutate(method=recode(method, MI_time="Nested MCMC", NeVI_time="NeVI-Cut"),
         method=factor(method, levels=method_levels),
         algorithm=factor(algorithm, levels=c("eava","insilicova","interva")))

dist_sub  <- dist_sub  %>%
  mutate(method=factor(method, levels=method_levels, labels = method_labels),
         algorithm=factor(algorithm, levels=c("eava","insilicova","interva")))

neonate_dist_sub = dist_sub
neonate_time_sub = time_long
save(neonate_dist_sub, neonate_time_sub, file = file.path(base_dir,"downstream_results","neonate_distance_time_sub.RData"))
