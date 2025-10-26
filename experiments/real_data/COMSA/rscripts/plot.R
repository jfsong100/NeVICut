

rm(list = ls())
current_path = rstudioapi::getSourceEditorContext()$path
base_dir = dirname(dirname(current_path))

modeltype.label = c(
  # 'yu' = "dodgerblue2",
  'Nested MCMC' = 'Nested MCMC',
  'NeVI-Cut' = 'NeVI-Cut',
  # 'seqbayes' = 'blue',
  'Full Bayes' = 'Full Bayes',
  # 'Parametric Cut' = expression("Dir("*M^T*alpha*")"),
  'Parametric Cut' = "Parametric Cut",
  # 'dirich_MTalpha' = "#A3A500",
  # 'dirich_vecM' = "#9590FF",
  # 'nevicut_mod' = 'blue',
  'Uncalibrated' = 'Uncalibrated'
)

modeltype.color = c(
  # 'yu' = "dodgerblue2",
  'Nested MCMC' = 'chocolate3',
  'NeVI-Cut' = 'green4',
  # 'seqbayes' = 'blue',
  'Full Bayes' = 'magenta4',
  'Parametric Cut' = "#F8766D",
  # 'dirich_MTalpha' = "#A3A500",
  # 'dirich_vecM' = "#9590FF",
  # 'nevicut_mod' = 'blue',
  'Uncalibrated' = 'black'
)

modeltype.linetype = c(
  # 'yu' = "dodgerblue2",
  'Nested MCMC' = 1,
  'NeVI-Cut' = 2,
  # 'seqbayes' = 3,
  'Full Bayes' = 3,
  'Parametric Cut' = 4,
  # 'dirich_MTalpha' = 6,
  # 'dirich_vecM' = 7,
  # 'nevicut_mod' = 'blue',
  'Uncalibrated' = 2
)

vaalgo.label = c(
  'eava' = 'EAVA',
  'insilicova' = 'InSilicoVA',
  'interva' = 'InterVA'
)

causes.label = c('congenital_malformation' = 'Con. Mal.',
                 'infection' = 'Infection', 'ipre' = 'IPRE',
                 'other'= 'Other', 'prematurity' = 'Prematurity',
                 "malaria" = 'Malaria', "pneumonia" = 'Pneumonia',
                 'sepsis' = 'Sepsis', 'sepsis_meningitis_inf' = 'Sep./Menin./Inf.',
                 "diarrhea" = 'Diarrhea', "severe_malnutrition" = 'Sev. Mal.',
                 "hiv" = 'HIV', "other_infections" = 'Oth. Inf.',
                 'sepsis_meningitis' = 'Sep./Menin.', 'injury' = "Injury",
                 "nn_causes" = "Neonatal causes")

nGrid = 512


# unadjusted -----
## neonate ----
cohort = 'neonate'

plotdf = NULL
runtime_multimpute_neonate = NULL
# for(vaalgo in c('eava', 'insilicova', 'interva')[1]){
for(vaalgo in c('eava', 'insilicova', 'interva')){
  
  
  # vaalgo = c('eava', 'insilicova', 'interva')[1]
  
  
  ### seq bayes ----
  seqbayes_out = readRDS(file.path(base_dir, "downstream_results",
                                   paste0(cohort, '_', vaalgo, '_seqbayes')))
  # head(seqbayes_out$MCMCout$p_calib)
  dim(seqbayes_out$p_calib[1,,])
  # class(seqbayes_out$MCMCout$p_calib)
  # vioplot::vioplot(seqbayes_out$MCMCout$p_calib, ylim = c(0,1))
  seqbayes_out_melt = reshape2::melt(seqbayes_out$p_calib[1,,])
  head(seqbayes_out_melt)
  dt_seqbayes = seqbayes_out_melt[,-1]
  colnames(dt_seqbayes) = c("causes", "value")
  head(dt_seqbayes)
  dt_seqbayes = data.table::as.data.table(dt_seqbayes)
  head(dt_seqbayes)
  dt_seqbayes$value = as.numeric(dt_seqbayes$value)
  head(dt_seqbayes)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_seqbayes = dt_seqbayes[, {
    d <- density(value, 
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_seqbayes)
  
  print(paste0(vaalgo, ", seqbayes"))
  
  
  ### multi impute: unadjusted ----
  multimpute_out = readRDS(file.path(base_dir, "downstream_results",
                                    paste0(cohort, '_', vaalgo, '_multimpute')))
  # head(multimpute_out$pcalib)
  dim(multimpute_out$pcalib)
  dimnames(multimpute_out$pcalib)
  # class(multimpute_out$pcalib)
  # vioplot::vioplot(multimpute_out$MCMCout$p_calib, ylim = c(0,1))
  multimpute_out_melt = reshape2::melt(multimpute_out$pcalib)
  head(multimpute_out_melt)
  dt_multimpute = multimpute_out_melt[,-(1:2)]
  colnames(dt_multimpute) = c("causes", "value")
  head(dt_multimpute)
  dt_multimpute = data.table::as.data.table(dt_multimpute)
  head(dt_multimpute)
  dt_multimpute$value = as.numeric(dt_multimpute$value)
  head(dt_multimpute)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_multimpute = dt_multimpute[, {
    d <- density(value, 
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_multimpute)
  
  runtime_multimpute_neonate = c(runtime_multimpute_neonate,
                                 multimpute_out$runtime["elapsed"])
  
  print(paste0(vaalgo, ", multimpute"))
  
  
  ### nevi cut ----
  nevicut_out = readRDS(file.path(base_dir, "downstream_results",
                                  paste0('NeVI_Cut_', cohort, '_', vaalgo, '.rds')))
  nevicut_out = nevicut_out[tail(1:nrow(nevicut_out), nrow(seqbayes_out$p_calib[1,,])),]
  dim(nevicut_out)
  head(nevicut_out)
  colnames(nevicut_out) = colnames(seqbayes_out$p_uncalib)
  head(nevicut_out)
  nevicut_out_melt = reshape2::melt(as.matrix(nevicut_out))
  head(nevicut_out_melt)
  dt_nevicut = nevicut_out_melt[,-1]
  colnames(dt_nevicut) = c("causes", "value")
  head(dt_nevicut)
  dt_nevicut = data.table::as.data.table(dt_nevicut)
  head(dt_nevicut)
  dt_nevicut$value = as.numeric(dt_nevicut$value)
  head(dt_nevicut)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_nevicut = dt_nevicut[, {
    d <- density(value,
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_nevicut)

  print(paste0(vaalgo, ", nevicut"))


  ### dirich ----
  dirich_out = readRDS(file.path(base_dir, "downstream_results",
                                     paste0('Dir_', cohort, '_', vaalgo, '.rds')))
  dim(dirich_out)
  head(dirich_out)
  colnames(dirich_out) = colnames(seqbayes_out$p_uncalib)
  head(dirich_out)
  dirich_out_melt = reshape2::melt(as.matrix(dirich_out))
  head(dirich_out_melt)
  dt_dirich = dirich_out_melt[,-1]
  colnames(dt_dirich) = c("causes", "value")
  head(dt_dirich)
  dt_dirich = data.table::as.data.table(dt_dirich)
  head(dt_dirich)
  dt_dirich$value = as.numeric(dt_dirich$value)
  head(dt_dirich)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_dirich = dt_dirich[, {
    d <- density(value,
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_dirich)

  print(paste0(vaalgo, ", dirich"))
  
  
  plotdf = rbind.data.frame(
    plotdf,
    data.frame(dens_seqbayes, 
               'modeltype' = 'Full Bayes', 
               'vaalgo' = vaalgo),
    data.frame(dens_multimpute, 
               'modeltype' = 'Nested MCMC',
               'vaalgo' = vaalgo),
    data.frame(dens_nevicut,
               'modeltype' = 'NeVI-Cut',
               'vaalgo' = vaalgo),
    data.frame(dens_dirich,
               'modeltype' = 'Parametric Cut',
               'vaalgo' = vaalgo),
    data.frame('causes' = colnames(seqbayes_out$p_uncalib),
               'x' = unname(seqbayes_out$p_uncalib[1,]), 
               'density' = NA,
               'modeltype' = 'Uncalibrated', 
               'vaalgo' = vaalgo)
  )
  
  print(vaalgo)
  
}

head(plotdf)
tail(plotdf)

plotdf$causes = factor(x = plotdf$causes, levels = colnames(seqbayes_out$p_uncalib))
# plotdf$modeltype = factor(x = plotdf$modeltype, levels = names(modeltype.label))
plotdf$modeltype = as.factor(x = plotdf$modeltype)
plotdf$vaalgo = factor(x = plotdf$vaalgo, levels = names(vaalgo.label))

head(plotdf)
tail(plotdf)


ggplot2::ggplot(data = plotdf) +
  ggh4x::facet_grid2(vaalgo~causes, 
                     # scales = 'free', independent = "y",
                     scales = 'free_y', independent = "y",
                     # nrow = 2, ncol = 3,
                     labeller = ggplot2::labeller(
                       causes = causes.label,
                       vaalgo = vaalgo.label
                     )
  ) +
  # ggplot2::facet_grid(vaalgo~Var2, scales = 'free_y',
  #                     # nrow = 2, ncol = 3,
  #                     labeller = ggplot2::labeller(Var2 = cause.label)
  # ) +
  # ggplot2::facet_wrap(vaalgo~causes, scales = 'free_y',
  #                     # nrow = 1, ncol = 6,
  #                     nrow = 3, ncol = 5,
  #                     labeller = ggplot2::labeller(causes = cause.label)
  # ) +
  ggplot2::coord_cartesian(
    xlim = c(0,1),
    ylim = c(0,NA),
    expand = TRUE, default = FALSE, clip = 'on'
  ) +
  ggplot2::geom_line(
    data = plotdf[plotdf$modeltype!='Uncalibrated',],
    ggplot2::aes(x = x, y = density, color = modeltype, linetype = modeltype),
    linewidth = 1.3
  ) +
  # ggplot2::geom_density(
  #   data = plotdf[plotdf$modeltype!='uncalib',],
  #   ggplot2::aes(x = x, y = density, color = modeltype, linetype = modeltype),
  #   # fill = "#619CFF", color = '#619CFF',
  #   position = "identity", trim = T, #adjust = 5,
  #   alpha = .7, linewidth = .8) +
  ggplot2::geom_vline(
    data = plotdf[plotdf$modeltype=='Uncalibrated',],
    ggplot2::aes(xintercept = x, color = modeltype, linetype = modeltype),
    # linetype = "dashed",
    linewidth = 1.3
  ) +
  # ggplot2::scale_x_discrete(labels = modeltype.label) +
  ggplot2::scale_color_manual(
    values = modeltype.color,
    labels = modeltype.label
  ) +
  ggplot2::scale_linetype_manual(
    values = modeltype.linetype,
    labels = modeltype.label
  ) +
  # ggplot2::scale_fill_manual(values = modeltype.colors) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 16),
    axis.title.x = ggplot2::element_text(size=20,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=20,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 16,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 16),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 18,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 18,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(2.5, "cm"),
    legend.key.height = ggplot2::unit(.95, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.key.spacing.x = ggplot2::unit(2, 'cm'),
    legend.text = ggplot2::element_text(size=20),
    legend.position = 'bottom'
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow=FALSE),
                  linetype = ggplot2::guide_legend(nrow = 2, byrow=FALSE)) +
  ggplot2::labs(
    # title = "Neonate", subtitle = "Unadjusted",
    x = 'Calibrated CSMFs', 
    y = 'Posterior Density', fill = NULL
  ) # 17 X 9.5



## child ----
cohort = 'child'

plotdf = NULL
runtime_multimpute_child = NULL
# for(vaalgo in c('eava', 'insilicova', 'interva')[1]){
for(vaalgo in c('eava', 'insilicova', 'interva')){
  
  
  # vaalgo = c('eava', 'insilicova', 'interva')[1]
  
  
  ### seq bayes ----
  seqbayes_out = readRDS(file.path(base_dir, "downstream_results",
                                   paste0(cohort, '_', vaalgo, '_seqbayes')))
  # head(seqbayes_out$MCMCout$p_calib)
  dim(seqbayes_out$p_calib[1,,])
  # class(seqbayes_out$MCMCout$p_calib)
  # vioplot::vioplot(seqbayes_out$MCMCout$p_calib, ylim = c(0,1))
  seqbayes_out_melt = reshape2::melt(seqbayes_out$p_calib[1,,])
  head(seqbayes_out_melt)
  dt_seqbayes = seqbayes_out_melt[,-1]
  colnames(dt_seqbayes) = c("causes", "value")
  head(dt_seqbayes)
  dt_seqbayes = data.table::as.data.table(dt_seqbayes)
  head(dt_seqbayes)
  dt_seqbayes$value = as.numeric(dt_seqbayes$value)
  head(dt_seqbayes)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_seqbayes = dt_seqbayes[, {
    d <- density(value, 
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_seqbayes)
  
  print(paste0(vaalgo, ", seqbayes"))
  
  
  ### multi impute: unadjusted ----
  multimpute_out = readRDS(file.path(base_dir, "downstream_results",
                                     paste0(cohort, '_', vaalgo, '_multimpute')))
  # head(multimpute_out$pcalib)
  dim(multimpute_out$pcalib)
  dimnames(multimpute_out$pcalib)
  # class(multimpute_out$pcalib)
  # vioplot::vioplot(multimpute_out$MCMCout$p_calib, ylim = c(0,1))
  multimpute_out_melt = reshape2::melt(multimpute_out$pcalib)
  head(multimpute_out_melt)
  dt_multimpute = multimpute_out_melt[,-(1:2)]
  colnames(dt_multimpute) = c("causes", "value")
  head(dt_multimpute)
  dt_multimpute = data.table::as.data.table(dt_multimpute)
  head(dt_multimpute)
  dt_multimpute$value = as.numeric(dt_multimpute$value)
  head(dt_multimpute)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_multimpute = dt_multimpute[, {
    d <- density(value, 
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_multimpute)
  
  runtime_multimpute_child = c(runtime_multimpute_child,
                                 multimpute_out$runtime["elapsed"])
  
  print(paste0(vaalgo, ", multimpute"))
  
  
  ### nevi cut ----
  nevicut_out = readRDS(file.path(base_dir, "downstream_results",
                                  paste0('NeVI_Cut_', cohort, '_', vaalgo, '.rds')))
  nevicut_out = nevicut_out[tail(1:nrow(nevicut_out), nrow(seqbayes_out$p_calib[1,,])),]
  dim(nevicut_out)
  head(nevicut_out)
  colnames(nevicut_out) = colnames(seqbayes_out$p_uncalib)
  head(nevicut_out)
  nevicut_out_melt = reshape2::melt(as.matrix(nevicut_out))
  head(nevicut_out_melt)
  dt_nevicut = nevicut_out_melt[,-1]
  colnames(dt_nevicut) = c("causes", "value")
  head(dt_nevicut)
  dt_nevicut = data.table::as.data.table(dt_nevicut)
  head(dt_nevicut)
  dt_nevicut$value = as.numeric(dt_nevicut$value)
  head(dt_nevicut)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_nevicut = dt_nevicut[, {
    d <- density(value,
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_nevicut)
  
  print(paste0(vaalgo, ", nevicut"))
  
  
  ### dirich ----
  dirich_out = readRDS(file.path(base_dir, "downstream_results",
                                 paste0('Dir_', cohort, '_', vaalgo, '.rds')))
  dim(dirich_out)
  head(dirich_out)
  colnames(dirich_out) = colnames(seqbayes_out$p_uncalib)
  head(dirich_out)
  dirich_out_melt = reshape2::melt(as.matrix(dirich_out))
  head(dirich_out_melt)
  dt_dirich = dirich_out_melt[,-1]
  colnames(dt_dirich) = c("causes", "value")
  head(dt_dirich)
  dt_dirich = data.table::as.data.table(dt_dirich)
  head(dt_dirich)
  dt_dirich$value = as.numeric(dt_dirich$value)
  head(dt_dirich)
  # grid = seq(min[1], rng[2], length.out = 512)
  dens_dirich = dt_dirich[, {
    d <- density(value,
                 from = 0, to = 1, #n = length(grid),
                 n = nGrid)
    .(x = d$x, density = d$y)
  }, by = causes]
  head(dens_dirich)
  
  print(paste0(vaalgo, ", dirich"))
  
  
  plotdf = rbind.data.frame(
    plotdf,
    data.frame(dens_seqbayes, 
               'modeltype' = 'Full Bayes', 
               'vaalgo' = vaalgo),
    data.frame(dens_multimpute, 
               'modeltype' = 'Nested MCMC',
               'vaalgo' = vaalgo),
    data.frame(dens_nevicut,
               'modeltype' = 'NeVI-Cut',
               'vaalgo' = vaalgo),
    data.frame(dens_dirich,
               'modeltype' = 'Parametric Cut',
               'vaalgo' = vaalgo),
    data.frame('causes' = colnames(seqbayes_out$p_uncalib),
               'x' = unname(seqbayes_out$p_uncalib[1,]), 
               'density' = NA,
               'modeltype' = 'Uncalibrated', 
               'vaalgo' = vaalgo)
  )
  
  print(vaalgo)
  
}

head(plotdf)
tail(plotdf)

plotdf$causes = factor(x = plotdf$causes, levels = colnames(seqbayes_out$p_uncalib))
# plotdf$modeltype = factor(x = plotdf$modeltype, levels = names(modeltype.label))
plotdf$modeltype = as.factor(x = plotdf$modeltype)
plotdf$vaalgo = factor(x = plotdf$vaalgo, levels = names(vaalgo.label))

head(plotdf)
tail(plotdf)


ggplot2::ggplot(data = plotdf) +
  ggh4x::facet_grid2(causes~vaalgo, 
                     # scales = 'free', independent = "y",
                     scales = 'free_y', independent = "y",
                     # nrow = 2, ncol = 3,
                     labeller = ggplot2::labeller(
                       causes = causes.label,
                       vaalgo = vaalgo.label
                     )
  ) +
  # ggplot2::facet_grid(vaalgo~Var2, scales = 'free_y',
  #                     # nrow = 2, ncol = 3,
  #                     labeller = ggplot2::labeller(Var2 = cause.label)
  # ) +
  # ggplot2::facet_wrap(vaalgo~causes, scales = 'free_y',
  #                     # nrow = 1, ncol = 6,
  #                     nrow = 3, ncol = 5,
  #                     labeller = ggplot2::labeller(causes = cause.label)
  # ) +
  ggplot2::coord_cartesian(
    xlim = c(0,1),
    ylim = c(0,NA),
    expand = TRUE, default = FALSE, clip = 'on'
  ) +
  ggplot2::geom_line(
    data = plotdf[plotdf$modeltype!='Uncalibrated',],
    ggplot2::aes(x = x, y = density, color = modeltype, linetype = modeltype),
    linewidth = 1.3
  ) +
  # ggplot2::geom_density(
  #   data = plotdf[plotdf$modeltype!='uncalib',],
  #   ggplot2::aes(x = x, y = density, color = modeltype, linetype = modeltype),
  #   # fill = "#619CFF", color = '#619CFF',
  #   position = "identity", trim = T, #adjust = 5,
  #   alpha = .7, linewidth = .8) +
  ggplot2::geom_vline(
    data = plotdf[plotdf$modeltype=='Uncalibrated',],
    ggplot2::aes(xintercept = x, color = modeltype, linetype = modeltype),
    # linetype = "dashed",
    linewidth = 1.3
  ) +
  # ggplot2::scale_x_discrete(labels = modeltype.label) +
  ggplot2::scale_color_manual(
    values = modeltype.color,
    labels = modeltype.label
  ) +
  ggplot2::scale_linetype_manual(
    values = modeltype.linetype,
    labels = modeltype.label
  ) +
  # ggplot2::scale_fill_manual(values = modeltype.colors) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 16),
    axis.title.x = ggplot2::element_text(size=20,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=20,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 18,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 18),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 18,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 18,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(2.5, "cm"),
    legend.key.height = ggplot2::unit(.95, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.key.spacing.x = ggplot2::unit(2, 'cm'),
    legend.text = ggplot2::element_text(size=20),
    legend.position = 'bottom'
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow=FALSE),
                  linetype = ggplot2::guide_legend(nrow = 2, byrow=FALSE)) +
  ggplot2::labs(
    # title = "child", subtitle = "Unadjusted",
    x = 'Calibrated CSMFs', 
    y = 'Posterior Density', fill = NULL
  )
# 20 X 12

runtime_multimpute = rbind.data.frame(runtime_multimpute_neonate,
                                      runtime_multimpute_child)
rownames(runtime_multimpute) = c('neonate', 'child')
colnames(runtime_multimpute) = c('eava', 'insilicova', 'interva')
runtime_multimpute
write.csv(x = runtime_multimpute, file = file.path(base_dir, "downstream_results",
                                                   "runtime_multimpute.csv"))



# wasserstein distance, run time ----
load(file.path(base_dir, "downstream_results", "neonate_distance_time_sub.RData"))

## distance ----
neonate_dist_sub

neonate_dist_sub = rbind.data.frame(neonate_dist_sub,
                                    data.frame("population" = "Neonate",
                                               "algorithm" = c('eava', 'insilicova', 'interva'),
                                               "method" = "Nested MCMC",
                                               "sw2" = NA))
neonate_dist_sub


## time ----
neonate_time_sub

neonate_time_sub = rbind.data.frame(neonate_time_sub,
                                    data.frame("population" = "neonate",
                                               "algorithm" = rep(c('eava', 'insilicova', 'interva'),
                                                                 each = 2),
                                               "method" = rep(c("Full Bayes", "Parametric Cut"), 3),
                                               "time" = NA))
neonate_time_sub

# neonate_time_sub$algorithm = factor(x = neonate_time_sub$algorithm,
#                                     levels = c('eava', 'insilicova', 'interva'))
# neonate_time_sub$method = factor(x = neonate_time_sub$method,
#                                  levels = c("Nested MCMC", "NeVI-Cut", "Full Bayes", "Parametric Cut"))
# 
# neonate_time_sub
# 
# ggplot_time = ggplot2::ggplot(data = neonate_time_sub) +
#   ggplot2::coord_cartesian(
#     # xlim = c(0,1),
#     ylim = c(0,NA),
#     expand = TRUE, default = FALSE, clip = 'on'
#   ) +
#   ggplot2::geom_bar(ggplot2::aes(x = algorithm, y = time, fill = method),
#                     stat = "identity", color = 'black', alpha = .7, linewidth = .7,
#                     width = .7,
#                     position = 'dodge' #ggplot2::position_dodge(width = .1)
#   ) +
#   ggplot2::scale_fill_manual(values = modeltype.color,
#                              labels = modeltype.label) +
#   ggplot2::scale_x_discrete(labels = vaalgo.label) +
#   ggplot2::theme(
#     plot.title = ggplot2::element_text(size = 18, face = "bold"),
#     plot.subtitle = ggplot2::element_text(size = 16),
#     axis.title.x = ggplot2::element_text(size=20,
#                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
#     # axis.title.x = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_text(size=20,
#                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
#     axis.text.x = ggplot2::element_text(color = "black", size = 16,
#                                         angle = 0, hjust = .5, vjust = 1),
#     axis.text.y = ggplot2::element_text(color = "black", size = 16),
#     axis.ticks.x = ggplot2::element_line(linewidth = .5),
#     axis.ticks.length.x = ggplot2::unit(.2, "cm"),
#     axis.ticks.y = ggplot2::element_line(linewidth = .5),
#     axis.ticks.length.y = ggplot2::unit(.2, "cm"),
#     panel.background = ggplot2::element_blank(),
#     panel.border = ggplot2::element_rect(color='black', linetype = "solid",
#                                          fill = NA, linewidth = 1),
#     panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
#                                              colour = "grey90"),
#     panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
#                                              colour = "grey90"),
#     strip.text.x = ggplot2::element_text(
#       size = 18,
#       face = "bold"
#     ),
#     strip.text.y = ggplot2::element_text(
#       size = 18,
#       face = "bold"
#     ),
#     strip.background = ggplot2::element_rect(color="black", linewidth=1),
#     legend.title = ggplot2::element_blank(),
#     legend.key.width = ggplot2::unit(1.5, "cm"),
#     legend.key.height = ggplot2::unit(1, "cm"),
#     # legend.key.size = ggplot2::unit(.5, "cm"),
#     legend.key.spacing.x = ggplot2::unit(2, 'cm'),
#     legend.text = ggplot2::element_text(size=20),
#     legend.position = 'bottom'
#   ) +
#   ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow=FALSE),
#                   linetype = ggplot2::guide_legend(nrow = 2, byrow=FALSE)) +
#   ggplot2::labs(
#     title = "Runtime (Seconds)",
#     x = 'CCVA Algorithm',
#     y = 'Value', fill = NULL
#     )
# ggplot_time
# 
# library(patchwork)
# 
# (ggplot_dist | ggplot_time) +
#   patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")
# 
# ggpubr::ggarrange(ggplot_dist, ggplot_time, nrow = 1, ncol = 2,
#                   legend = "bottom", common.legend = TRUE)


## combine ----
mydf_dist = neonate_dist_sub
mydf_time = neonate_time_sub
colnames(mydf_dist) = colnames(mydf_time) =
  c("population", "algorithm", "method", "value")
plotdf = rbind.data.frame(cbind.data.frame(mydf_dist, 
                                           "value_type" = "Wasserstein Distance"),
                          cbind.data.frame(mydf_time, 
                                           "value_type" = "Runtime (Seconds)"))
head(plotdf)
tail(plotdf)

plotdf$algorithm = factor(x = plotdf$algorithm,
                          levels = c('eava', 'insilicova', 'interva'))
plotdf$method = factor(x = plotdf$method,
                       levels = names(modeltype.label))
plotdf$value_type = factor(x = plotdf$value_type,
                       levels = c("Wasserstein Distance", "Runtime (Seconds)"))

head(plotdf)
tail(plotdf)

ggplot2::ggplot(data = plotdf) +
  ggh4x::facet_grid2(.~value_type, 
                     # scales = 'free', independent = "y",
                     scales = 'free_y', independent = "y"
  ) +
  ggplot2::coord_cartesian(
    # xlim = c(0,1),
    ylim = c(0,NA),
    expand = TRUE, default = FALSE, clip = 'on'
  ) +
  ggplot2::geom_bar(ggplot2::aes(x = algorithm, y = value, fill = method),
                    stat = "identity", color = 'black', alpha = .7, linewidth = .7,
                    width = .7, position = 'dodge') +
  ggplot2::scale_fill_manual(values = modeltype.color,
                             labels = modeltype.label) +
  ggplot2::scale_x_discrete(labels = vaalgo.label) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 16),
    axis.title.x = ggplot2::element_text(size=22,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=22,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 18,
                                        angle = 0, hjust = .5, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 18),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.x = ggplot2::unit(.2, "cm"),
    axis.ticks.y = ggplot2::element_line(linewidth = .5),
    axis.ticks.length.y = ggplot2::unit(.2, "cm"),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.2, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(
      size = 22,
      face = "bold"
    ),
    strip.text.y = ggplot2::element_text(
      size = 22,
      face = "bold"
    ),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(1.5, "cm"),
    legend.key.height = ggplot2::unit(1, "cm"),
    # legend.key.size = ggplot2::unit(.5, "cm"),
    legend.key.spacing.x = ggplot2::unit(2, 'cm'),
    legend.text = ggplot2::element_text(size=20),
    legend.position = 'bottom'
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow=FALSE),
                  linetype = ggplot2::guide_legend(nrow = 2, byrow=FALSE)) +
  ggplot2::labs(
    # title = "Neonate", subtitle = "Unadjusted",
    x = 'CCVA Algorithms', 
    y = 'Value', fill = NULL
  ) # 7 X 13


