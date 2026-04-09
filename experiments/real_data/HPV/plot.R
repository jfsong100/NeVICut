rm(list = ls()) 
current_path <- rstudioapi::getSourceEditorContext()$path
base_dir <- dirname(current_path)
output.path <- file.path(base_dir, "downstream")

theta2_samples_Yu <- readRDS(file.path(output.path,"theta2_samples_Gaussian_VA.rds"))
theta2_samples_NeVi_Cut <- readRDS(file.path(output.path,"NeVI_Cut_theta2_samples.rds"))
multiple = readRDS(file.path(output.path,"hpv_multimpute"))
theta2_samples_mi = multiple[["theta2"]]
full = readRDS(file.path(output.path,"hpv_fullbayes"))
theta2_samples_full = full[["MCMCout"]][["theta2"]]
theta2_samples_naive <- readRDS(file.path(output.path,"theta2_samples_naive_cut.rds"))
theta2_samples_plummer <- readRDS(file.path(output.path,"theta2_samples_plummer.rds"))


library(ggplot2)
library(dplyr)

# Combine all samples into a data frame
df <- bind_rows(
  tibble(Value = as.numeric(theta2_samples_NeVi_Cut$Dim1), Method = "Ne-VI Cut", Parameter = "eta21"),
  tibble(Value = as.numeric(theta2_samples_Yu$Dim1), Method = "Gaussian VA-Cut", Parameter = "eta21"),
  tibble(Value = as.numeric(theta2_samples_naive$theta21_samples$`adaptive metropolis 1D`[[1]]), Method = "OpenBUGS-Cut", Parameter = "eta21"),
  tibble(Value = as.numeric(theta2_samples_plummer$theta21_samples), Method = "Tempered Cut", Parameter = "eta21"),
  tibble(Value = as.numeric(theta2_samples_mi[,1]), Method = "Nested MCMC", Parameter = "eta21"),
  tibble(Value = as.numeric(theta2_samples_full[,1]), Method = "Full Bayes", Parameter = "eta21"),
  
  tibble(Value = as.numeric(theta2_samples_NeVi_Cut$Dim2), Method = "Ne-VI Cut", Parameter = "eta22"),
  tibble(Value = as.numeric(theta2_samples_Yu$Dim2), Method = "Gaussian VA-Cut", Parameter = "eta22"),
  tibble(Value = as.numeric(theta2_samples_naive$theta22_samples$`adaptive metropolis 1D`[[1]]), Method = "OpenBUGS-Cut", Parameter = "eta22"),
  tibble(Value = as.numeric(theta2_samples_plummer$theta22_samples), Method = "Tempered Cut", Parameter = "eta22"),
  tibble(Value = as.numeric(theta2_samples_mi[,2]), Method = "Nested MCMC", Parameter = "eta22"),
  tibble(Value = as.numeric(theta2_samples_full[,2]), Method = "Full Bayes", Parameter = "eta22")
)

# Set method colors
method_colors <- c(
  "Full Bayes" = "magenta4",              # Purple
  "Nested MCMC" = "chocolate3",     # Orange
  "Gaussian VA-Cut" = "#1F78B4",                 # Light Blue
  "Tempered Cut" = "#0000FF",                 # Blue
  "OpenBUGS-Cut" = "gold1",               # Pink
  "Ne-VI Cut" = "green4"                 # Green
)

df_trimmed <- df %>%
  mutate(Value = ifelse(Parameter == "eta21" & (Value < -3 | Value > -1), NA, Value))

df_trimmed <- df_trimmed %>%
  mutate(Value = ifelse(Parameter == "eta22" & (Value > 40), NA, Value))

df_trimmed$Parameter <- factor(df_trimmed$Parameter, levels = c("eta21", "eta22"),
                               labels = c("eta[21]", "eta[22]"))

density_marginal = ggplot(df_trimmed, aes(x = Value, color = Method)) +
  geom_density(linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~Parameter, scales = "free", ncol = 2, labeller = label_parsed) +
  scale_color_manual(values = method_colors) +
  labs(x = "Value", y = "Density", color = "Method") +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  )


density_marginal

# Combine theta21 and theta22 into one dataset for joint plot
df_joint <- bind_rows(
  tibble(theta21 = as.numeric(theta2_samples_NeVi_Cut$Dim1),
         theta22 = as.numeric(theta2_samples_NeVi_Cut$Dim2),
         Method = "Ne-VI Cut"),
  tibble(theta21 = as.numeric(theta2_samples_Yu$Dim1),
         theta22 = as.numeric(theta2_samples_Yu$Dim2),
         Method = "Gaussian VA-Cut"),
  tibble(theta21 = as.numeric(theta2_samples_naive$theta21_samples$`adaptive metropolis 1D`[[1]]),
         theta22 = as.numeric(theta2_samples_naive$theta22_samples$`adaptive metropolis 1D`[[1]]),
         Method = "OpenBUGS-Cut"),
  tibble(theta21 = as.numeric(theta2_samples_plummer$theta21_samples),
         theta22 = as.numeric(theta2_samples_plummer$theta22_samples),
         Method = "Tempered Cut"),
  tibble(theta21 = as.numeric(theta2_samples_mi[,1]),
         theta22 = as.numeric(theta2_samples_mi[,2]),
         Method = "Nested MCMC"),
  tibble(theta21 = as.numeric(theta2_samples_full[,1]),
         theta22 = as.numeric(theta2_samples_full[,2]),
         Method = "Full Bayes")
)

# Plot joint density
density_joint = ggplot(df_joint, aes(x = theta21, y = theta22, color = Method, fill = Method)) +
  stat_density_2d(geom = "polygon", alpha = 0.3) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(x = expression(eta[21]), y = expression(eta[22])) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

density_joint

library(purrr)
library(dplyr)
library(ggplot2)
library(patchwork)

method_colors <- c(
  "Full Bayes" = "magenta4",
  "Nested MCMC" = "chocolate3",
  "Gaussian VA-Cut" = "#1F78B4",
  "Tempered Cut" = "#0000FF",
  "OpenBUGS-Cut" = "gold1",
  "Ne-VI Cut" = "green4"
)

ref_name <- "Nested MCMC"

df2 <- df %>%
  filter(Method != "Full Bayes") %>%
  filter(!(Parameter == "eta21" & (Value < -3 | Value > -1)),
         !(Parameter == "eta22" & (Value > 40))) %>%
  mutate(
    Parameter = factor(Parameter, levels = c("eta21","eta22"),
                       labels = c("eta[21]", "eta[22]"))
  )

methods_to_show <- setdiff(unique(df2$Method), ref_name)

make_pair_panel <- function(m) {
  dat <- df2 %>%
    filter(Method %in% c(ref_name, m)) %>%
    droplevels() %>%
    # ensure MI always comes first in legend
    mutate(Method = factor(Method, levels = c(ref_name, m)))
  
  ggplot(dat, aes(x = Value, color = Method, linetype = Method)) +
    geom_density(linewidth = 0.9) +
    facet_wrap(~ Parameter, ncol = 2, scales = "free",
               labeller = labeller(Parameter = label_parsed)) +
    scale_color_manual(values = method_colors[c(ref_name, m)], drop = TRUE) +
    scale_linetype_manual(
      values = setNames(c("solid","dashed"), c(ref_name, m)),
      drop = TRUE
    ) +
    labs(x = "Value", y = "Density") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

pair_plots <- map(methods_to_show, make_pair_panel)

p_com = wrap_plots(pair_plots, ncol = 2)
p_com
