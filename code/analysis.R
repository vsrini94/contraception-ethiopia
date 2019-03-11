# Description: Clean and tabulate input survey microdata for analysis
# Inputs: DHS and PMA 2020 microdata (prepped for Global Burden of Disease)
# Outputs: Dataset with admin 1-level tabulations for analysis
rm(list = ls())
if(!require(pacman)) install.packages('pacman')
pacman::p_load(data.table, magrittr, haven, survey, boot, sp, maptools, rgdal, ggplot2, RColorBrewer, gridExtra,
               INLA, spdep)

# Read in prepped dataset(s)

data_df <- fread('../data/prepped.csv')[order(ADM1_CODE)]
spobj <- readRDS('../data/eth_shp_prepped.rds')


dhs_sp <- merge(spobj, data_df[source == 'dhs'], by = 'ADM1_CODE')
pma_sp <- merge(spobj, data_df[source == 'pma2020'], by = 'ADM1_CODE')

# Plot raw estimates from surveys to show differences between DHS/PMA

color_scale <- colorRampPalette(brewer.pal(9, "Blues"))(11)
raw_scale <- seq(min(data_df$prev_met_mod_need)-.05, max(data_df$prev_met_mod_need) + .05, l = 10)

dhs_raw <- spplot(dhs_sp, 'prev_met_mod_need', main = 'DHS, 2016', col = 'transparent', col.regions = color_scale, at = raw_scale)
pma_raw <- spplot(pma_sp, 'prev_met_mod_need', main = 'PMA 2020, 2016', col = 'transparent', col.regions = color_scale, at = raw_scale)

dot_plot_data <- copy(data_df) %>% .[, plot_order := mean(prev_met_mod_need), by = admin_1]

dot_plot_data <- dot_plot_data[order(plot_order)]
dot_plot_data[, plot_order := factor(plot_order, levels = unique(dot_plot_data$plot_order), labels = unique(dot_plot_data$admin_1))]

dot_plot <- ggplot(data = dot_plot_data, aes(x = plot_order, y = prev_met_mod_need)) + 
  geom_pointrange(aes(ymin = prev_met_mod_need - 1.96*prev_met_mod_need_se,
                      ymax = prev_met_mod_need + 1.96*prev_met_mod_need_se,
                      color = source)) + scale_color_brewer(palette = "Set2") +
  theme_bw() + coord_flip() + labs(x = 'Location', y = 'Prevalence Met Need for Modern Contraception',
                                   color = 'Survey')
png(filename = '../raw_prev_map.png', width = 800, height = 400, units = 'px')
grid.arrange(dhs_raw, pma_raw, nrow = 1)
dev.off()

png(filename = '../raw_prev_dotplot.png', width = 600, height = 200, units = 'px')
print(dot_plot)
dev.off()

# Create adjacency matrix for INLA
adj_mat <- poly2nb(spobj, row.names = spobj$ADM1_NAME, queen = F)
temp <- tempdir()
graph <- nb2INLA(paste0(temp, '/eth-INLA.adj'), adj_mat)

# Create numbered variable to match graph
data_df[, ADM1_CODE_N := seq(.N), by = source]
dhs_df <- data_df[source == 'dhs']

# Fit ICAR Model on individual surveys
dhs_only_icar <- inla(formula = logit(prev_met_mod_need) ~ 1 + f(ADM1_CODE_N, model="bym2", graph=paste0(temp, '/eth-INLA.adj'), 
                                                            scale.model=T, constr=T),
                 family  = 'gaussian', data = dhs_df, control.predictor = list(compute = T))

pma_only_icar <- inla(formula = logit(prev_met_mod_need) ~ 1 + f(ADM1_CODE_N, model="bym2", graph=paste0(temp, '/eth-INLA.adj'), 
                                                                 scale.model=T, constr=T),
                      family  = 'gaussian', data = data_df[source == 'pma2020'], control.predictor = list(compute = T))

# Fit ICAR model on pooled surveys with and without source type random effect
pooled_icar <- inla(formula = logit(prev_met_mod_need) ~ 1 + f(ADM1_CODE_N, model="bym2", graph=paste0(temp, '/eth-INLA.adj'), 
                                                                 scale.model=T, constr=T),
                      family  = 'gaussian', data = data_df, control.predictor = list(compute = T))

pooled_icar_source_re <- inla(formula = logit(prev_met_mod_need) ~ 1 + f(source, model = 'iid') + f(ADM1_CODE_N, model="bym2", graph=paste0(temp, '/eth-INLA.adj'), 
                                                                 scale.model=T, constr=T),
                      family  = 'gaussian', data = data_df, control.predictor = list(compute = T))

# Source REs show almost no differences between DHS/PMA, so don't need RE, but weight reg by inverse variance
inla.setOption(enable.inla.argument.weights = TRUE)
pooled_icar_weighted <- inla(formula = logit(prev_met_mod_need) ~ 1 + f(ADM1_CODE_N, model="bym2", graph=paste0(temp, '/eth-INLA.adj'), 
                                                                                       scale.model=T, constr=T),
                                            family  = 'gaussian', weights = 1/(data_df$prev_met_mod_need_se)^2, data = data_df, control.predictor = list(compute = T))

# Extract pooled icar for weighted model
total_res <- pooled_icar_weighted$summary.random$ADM1_CODE_N$mean[1:11]
spatial_res <- pooled_icar_weighted$summary.random$ADM1_CODE_N$mean[12:22]

# Extract fixed effects
fes <- pooled_icar_weighted$summary.fixed$mean
smoothed_preds <- inv.logit(fes + total_res)

spobj$SMOOTHED <- smoothed_preds

# Plot smoothed predictions
pooled_weighted <- spplot(spobj, 'SMOOTHED', main = 'Smoothed Predictions', col = 'transparent', col.regions = color_scale, at = raw_scale)

png(filename = '../smoothed_prev_map.png', width = 400, height = 400, units = 'px')
print(pooled_weighted)
dev.off()

# Comment on spatial random effects

# Plot distribution, comparison to raw survey data
shrinkage_df <- rbind(data_df[, .(admin_1, source, prev_met_mod_need)],
                      data.table(admin_1 = unique(data_df$admin_1), source = 'Weighted ICAR Model',
                                 prev_met_mod_need = smoothed_preds), use.names = T, fill = T)

# See variance across sources
shrinkage_df[, var(prev_met_mod_need), by = .(source)]
shrinkage_df[, var(prev_met_mod_need), by = .(!source %like% 'Weighted ICAR')]

dist_plot <- ggplot(data = shrinkage_df, aes(x = prev_met_mod_need, y = source, color = source)) + 
  geom_point(shape = 1, size = 2, stroke = 2) + scale_color_brewer(palette = "Set2") + guides(color = F) + 
  theme_bw()  + labs(x = 'Prevalence Met Need for Modern Contraception', y = 'Source',
                                   color = '')

png(filename = '../dist_comparisons.png', width = 500, height = 150, units = 'px')
print(dist_plot)
dev.off()

# Total variance of REs
hyperpar_summary <- pooled_icar_weighted$summary.hyperpar

# Proportion of variance from spatial RE
spatial_var <- list(med = hyperpar_summary[3,4], lower = hyperpar_summary[3, 3], upper = hyperpar_summary[3,5])
