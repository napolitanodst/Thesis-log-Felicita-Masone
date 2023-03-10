#--- Plotting Spatial transcriptome data from coronal mouse brain sections after
#         striatal injection of heme and heme-hemopexin---

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# Total UMI per Tissue Covered Spot
plots <- list()

for (i in 1:length(sample_names)) {
  
  plots[[i]] <- bcs_merge %>% 
    filter(sample ==sample_names[i]) %>% 
    ggplot(aes(x=imagecol,y=imagerow,fill=sum_umi)) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_gradientn(colours = myPalette(100))+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[i]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    labs(fill = "Total UMI")+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_grid(plotlist = plots)

# Total Genes per Tissue Covered Spot
plots <- list()

for (i in 1:length(sample_names)) {
  
  plots[[i]] <- bcs_merge %>% 
    filter(sample ==sample_names[i]) %>% 
    ggplot(aes(x=imagecol,y=imagerow,fill=sum_gene)) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_gradientn(colours = myPalette(100))+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[i]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    labs(fill = "Total Genes")+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_grid(plotlist = plots)

# Cluster Assignments per Tissue Covered Spot
plots <- list()

for (i in 1:length(sample_names)) {
  
  plots[[i]] <- bcs_merge %>% 
    filter(sample == sample_names[i]) %>%
    filter(tissue == "1") %>% 
    ggplot(aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold", "#a65628", "#999999", "black", "grey", "white", "purple"))+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[i]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    labs(fill = "Cluster")+
    guides(fill = guide_legend(override.aes = list(size=3)))+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_grid(plotlist = plots)

# Gene of intrest

plots <- list()

for (i in 1:length(sample_names)) {
  
  plots[[i]] <- bcs_merge %>% 
    filter(sample ==sample_names[i]) %>% 
    bind_cols(select(matrix[[i]], "Hpca")) %>% 
    ggplot(aes(x=imagecol,y=imagerow,fill=Hpca)) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_gradientn(colours = myPalette(100))+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[i]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_grid(plotlist = plots)