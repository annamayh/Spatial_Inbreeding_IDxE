load("Deer_spatial_variation_ID/plots/surv_plots.RData")
load("Deer_spatial_variation_ID/plots/bw_plots.RData")


all=(bw_reg+inla_bw_plot)/(surv_reg+inla_surv_plot)+plot_annotation(tag_levels = 'A')
all

ggsave(all,
       file = "Deer_spatial_variation_ID/plots/bw_surv_spatial_var.png",
       width = 11,
       height = 12,
       dpi=1000)
