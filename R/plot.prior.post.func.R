# Script of functions to plot prior and posterior/MMAP.

# function to plot normal prior and posterior MMAP
plot_norm_prior_post_simu <- function(data, mean, sd, par, par_val = NA,
                                      limits = c(NA, NA), bins = 30) {
  
  samp_df <- tibble(
    samp = data
  )
  
  calc_prior <- function(data, mean, sd, limits = c(NA, NA)) {
    if (is.na(limits[1])) {
      limits[1] = min(data)
    }
    if (is.na(limits[2])) { 
      limits[2] = max(data)
    }
    prior <- tibble(
      x = seq(limits[1], limits[2], length = 200),
      y = dnorm(x, mean, sd) / dnorm(mean, mean, sd)
    )
    invisible(prior)
  }
  
  p <- ggplot(samp_df, aes(x = samp)) +
    geom_histogram(aes(y = after_stat(count / max(count)), fill = "samples"),
                   color = "white", 
                   bins = bins) +
    geom_line(aes(x = x, y = y, colour = "Prior", linetype = 'Prior'), data = calc_prior(data, mean, sd, limits), size = 1) +
    scale_fill_manual(name = "Histogram",  values = c("samples"="darkgrey"), labels = c("samples" = "MMAP Estimates")) +
    theme_bw() +
    ggtitle(paste("Prior and MMAP Estimates of", par)) +
    xlab(paste(par, "MMAP Estimates")) + ylab("Density") +
    theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.width= unit(1, 'cm'),
          axis.title.x = element_text(size = 15, face = 'bold'),
          axis.title.y = element_text(size = 15, face = 'bold'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15))
  
  if (!is.na(par_val)) {
    p <- p +  geom_vline(aes(xintercept = par_val, colour = "TRUE", linetype = 'TRUE'), size = 1, show.legend = F) +
      scale_color_manual(name = "Line", values = c('Prior' = 'darkblue', 'TRUE' = 'red')) +
      scale_linetype_manual(name = "Line", values=c('Prior' = 'solid', "TRUE" = "dashed")) + 
      guides(fill = guide_legend(order = 1))
  }
  p
}

# function to plot gamma prior and posterior MMAP
plot_gamma_prior_post_simu <- function(data, shape, scale, par, par_val = NA, 
                                       limits = c(NA, NA), bins = 30) {
  samp_df <- tibble(
    samp = data
  )
  
  calc_prior <- function(data, shape, scale, limits = c(NA, NA)) {
    if (is.na(limits[1])) {
      limits[1] = min(data)
    }
    if (is.na(limits[2])) { 
      limits[2] = max(data)
    }
    
    prior <- tibble(x = seq(max(0, limits[1], na.rm = TRUE), limits[2], length = 200))
    
    if (shape <= 1) {
      prior <- prior %>% filter(x > 0)
      mode = min(prior$x)
    } else {
      mode = (shape - 1) * scale
    }
    
    prior <- prior %>% 
      mutate(y = dgamma(x, shape = shape, scale = scale) / dgamma(mode, shape = shape, scale = scale))
    invisible(prior)
  }
  
  p <- ggplot(samp_df, aes(x = samp)) +
    geom_histogram(aes(y = after_stat(count / max(count)), fill = "samples"),
                   color = 'white',
                   bins = bins) +
    geom_line(aes(x = x, y = y, colour = "Prior", linetype = 'Prior'), data = calc_prior(data, shape, scale, limits), size = 1) +
    scale_fill_manual(name = "Histogram",  values = c("samples"="darkgrey"), labels = c("samples" = "MMAP Estimates")) +
    theme_bw() +
    ggtitle(paste("Prior and MMAP Estimates of", par)) +
    xlab(paste(par, "MMAP Estimates")) + ylab("Density") +
    theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.width= unit(1, 'cm'),
          axis.title.x = element_text(size = 15, face = 'bold'),
          axis.title.y = element_text(size = 15, face = 'bold'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15))
  
  if (! is.na(par_val)) {
    p <- p +  geom_vline(aes(xintercept = par_val, colour = "TRUE", linetype = 'TRUE'), size = 1, show.legend = F) +
      scale_linetype_manual(name = "Line", values = c('Prior' = 'solid', "TRUE" = "dashed"))+
      scale_color_manual(name = "Line", values = c('Prior' = 'darkblue', 'TRUE' = 'red')) + 
      guides(fill = guide_legend(order = 1))
  }
  p
}

# function to plot normal prior and posterior MMAP
plot_norm_prior_post <- function(data, mean, sd, par, par_val = NA,
                                 limits = c(NA, NA), bins = 30) {
  
  samp_df <- tibble(
    samp = data
  )
  
  calc_prior <- function(data, mean, sd, limits = c(NA, NA)) {
    if (is.na(limits[1])) {
      limits[1] = min(data)
    }
    if (is.na(limits[2])) {
      limits[2] = max(data)
    }
    prior <- tibble(
      x = seq(limits[1], limits[2], length = 200),
      y = dnorm(x, mean, sd) / dnorm(mean, mean, sd)
    )
    invisible(prior)
  }
  
  p <- ggplot(samp_df, aes(x = samp)) +
    geom_histogram(aes(y = after_stat(count / max(count)), fill = "samples"),
                   color = "white",
                   bins = bins) +
    geom_line(aes(x = x, y = y, colour = "Prior", linetype = 'Prior'), data = calc_prior(data, mean, sd, limits), size = 1) +
    scale_fill_manual(name = "Histogram",  values = c("samples" = "darkgrey"), labels = c('samples' = 'Posterior Samples')) +
    theme_bw() +
    ggtitle(paste("Prior and Posterior Samples of", par)) +
    xlab(paste(par, "Posterior Samples")) + ylab("Density") +
    theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.width= unit(1, 'cm'),
          axis.title.x = element_text(size = 15, face = 'bold'),
          axis.title.y = element_text(size = 15, face = 'bold'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15))
  
  if (!is.na(par_val)) {
    p <- p +  geom_vline(aes(xintercept = par_val, colour = "MMAP", linetype = 'MMAP'), size = 1, show.legend = F) +
      scale_color_manual(name = "Line", values = c('Prior' = 'darkblue', 'MMAP' = 'red')) +
      scale_linetype_manual(name = "Line", values=c('Prior' = 'solid', "MMAP" = "dashed")) +
      guides(fill = guide_legend(order = 1))
  }
  p
}

# function to plot gamma prior and posterior MMAP
plot_gamma_prior_post <- function(data, shape, scale, par, par_val = NA,
                                  limits = c(NA, NA), bins = 30) {
  samp_df <- tibble(
    samp = data
  )
  
  calc_prior <- function(data, shape, scale, limits = c(NA, NA)) {
    if (is.na(limits[1])) {
      limits[1] = min(data)
    }
    if (is.na(limits[2])) {
      limits[2] = max(data)
    }
    
    prior <- tibble(x = seq(max(0, limits[1], na.rm = TRUE), limits[2], length = 200))
    
    if (shape <= 1) {
      prior <- prior %>% filter(x > 0)
      mode = min(prior$x)
    } else {
      mode = (shape - 1) * scale
    }
    
    prior <- prior %>%
      mutate(y = dgamma(x, shape = shape, scale = scale) / dgamma(mode, shape = shape, scale = scale))
    invisible(prior)
  }
  
  p <- ggplot(samp_df, aes(x = samp)) +
    geom_histogram(aes(y = after_stat(count / max(count)), fill = "samples"),
                   color = 'white',
                   bins = bins) +
    geom_line(aes(x = x, y = y, color = "Prior", linetype = 'Prior'), data = calc_prior(data, shape, scale, limits), size = 1) +
    scale_fill_manual(name = "Histogram",  values = c("samples" = "darkgrey"), labels = c('samples' = 'Posterior Samples')) +
    theme_bw() +
    ggtitle(paste("Prior and Posterior Samples of", par)) +
    xlab(paste(par, "Posterior Samples")) + ylab("Density") +
    theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.width= unit(1, 'cm'),
          axis.title.x = element_text(size = 15, face = 'bold'),
          axis.title.y = element_text(size = 15, face = 'bold'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15))
  
  if (! is.na(par_val)) {
    p <- p +  geom_vline(aes(xintercept = par_val, colour = "MMAP", linetype = 'MMAP'), size = 1, show.legend = F) +
      scale_color_manual(name = "Line", values = c('Prior' = 'darkblue', 'MMAP' = 'red')) +
      scale_linetype_manual(name = "Line", values=c('Prior' = 'solid', "MMAP" = "dashed")) +
      guides(fill = guide_legend(order = 1))
  }
  p
}