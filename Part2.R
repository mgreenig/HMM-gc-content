library(stringr)
library(dplyr)
library(ggridges)
library(ggpubr)
library(reshape2)

source('HMMfunctions.R')

if(!interactive()){
  parameter_file <- commandArgs(trailingOnly = T)
} else {
  parameter_file <- readline(prompt = 'Enter the filepath for the parameter file to be used: ')
}

# load in params
params <- readParamFile(parameter_file)

# load in sequence data
chr3 <- readLines('data/SCerevisiaeChr3.fasta')[-1]
chr3 <- paste(chr3, collapse = '')
# get 100bp non-overlapping windows from chr3
chr3_chunks <- substring(chr3, seq(1, nchar(chr3)-1, 100), seq(100, nchar(chr3), 100))
# only keep chunks that are 100bp
chr3_chunks <- chr3_chunks[nchar(chr3_chunks) == 100]
gc_contents <- str_count(chr3_chunks, 'G|C') / 100
gc_quantiles <- quantile(gc_contents, probs = seq(0, 1, 1 / (ncol(params$eMat))))

# get emission values 
emissions <- as.integer(cut(gc_contents, gc_quantiles, labels = params$E, include.lowest = T))

# put gc content and emissions into data frame
gc_content_density <- density(gc_contents)
gc_content_df <- data.frame(x = gc_content_density$x, y = gc_content_density$y)
gc_content_df$quantile <- factor(findInterval(gc_content_df$x, gc_quantiles, all.inside = T))

# make GC content histogram
gc_content_density <- ggplot(gc_content_df, aes(x, y)) +
  geom_line() + geom_ribbon(aes(ymin = 0, ymax = y, fill = quantile)) + theme_minimal() +
  labs(x = '\nGC content', y = 'Density\n', 
       title = 'Estimated density of GC content values',
       subtitle = 'Calculated from 100bp chunks of S. cerevisiae Chr3\n',
       fill = 'GC content bin\n(emission state)') +
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 14, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave('figures/fig2.png', gc_content_density, width = 7, height = 4)

# calculate likelihood with forward algorithm
logL <- scaledForwardAlgorithm(params[['E']],
                               params[['init']],
                               params[['tMat']],
                               params[['eMat']],
                               emissions)

# estimate parameters via Baum-Welch
theta_hat <- baumWelch(params[['E']],
                       params[['init']],
                       params[['tMat']],
                       params[['eMat']],
                       emissions,
                       L_threshold = 0.01)

baumWelch_plot <- plotBaumWelch(theta_hat$logLs, theta_hat$deltas)

ggsave('figures/fig3.png', baumWelch_plot, width = 16, height = 8)

# plot new emissions distribution compared to old one
old_emission_dist_plot <- plotEmissionDistribution(params[['S']], 
                                                   params[['E']], 
                                                   params[['eMat']]) +
  labs(title = 'Before Baum-Welch training', x = '') + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))
new_emission_dist_plot <- plotEmissionDistribution(params[['S']], 
                                                   params[['E']], 
                                                   theta_hat[['eMat']]) +
  labs(title = 'After Baum-Welch training') + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))
combined_emission_dists_plot <- ggarrange(old_emission_dist_plot, 
                                          new_emission_dist_plot,
                                          common.legend = T,
                                          nrow = 2, ncol = 1)
emission_plot_title <- text_grob('Emission distributions before and after Baum-Welch training',
                                size = 14)
combined_emission_dists_plot <- annotate_figure(combined_emission_dists_plot, 
                                                top = emission_plot_title)
ggsave('figures/fig4.png', combined_emission_dists_plot, width = 6, height = 6)

# calculate most probable state path using the Viterbi algorithm
states <- viterbiAlgorithm(params[['S']],
                           params[['E']],
                           theta_hat[['init']],
                           theta_hat[['tMat']],
                           theta_hat[['eMat']],
                           emissions)

# plot viterbi sequence alongside actual gc contents and model emissions
GC_content_plot <- plotViterbi(states, gc_contents) + 
  labs(x = '', y = 'GC content\n', title = 'Actual GC content values') + 
  theme(plot.title = element_text(hjust = 0.5))
emissions_plot <- plotViterbi(states, emissions) + 
  labs(title = 'Model emissions') + 
  theme(plot.title = element_text(hjust = 0.5))

# combine plots
combined_viterbi_plot <- ggarrange(GC_content_plot, emissions_plot, ncol = 1, nrow = 2, 
                                   common.legend = T, legend = 'top')
viterbi_plot_title <- text_grob('Most likely hidden state sequence for the first 120 segments of chromosome III, S. cerevisiae\n', size = 26)
combined_viterbi_plot <- annotate_figure(combined_viterbi_plot, 
                                         top = viterbi_plot_title)
ggsave('figures/fig5.png', combined_viterbi_plot, width = 16, height = 14)

# initialise emission parameter matrix for normal distribution
normal_eMat <- matrix(c(0.49, 0.05, 0.51, 0.05), 2, 2, byrow = T)

# estimate parameters via Baum-Welch
theta_hat_normal <- baumWelch(params[['E']],
                              params[['init']],
                              params[['tMat']],
                              normal_eMat,
                              gc_contents,
                              eFunc = dnorm, 
                              paramFunc = normalParams,
                              L_threshold = 0.01)

baumWelch_plot_normal <- plotBaumWelch(theta_hat_normal$logLs, theta_hat_normal$deltas,
                                       title = 'Change in log-likelihood during Baum-Welch training process, gaussian emissions', ylab1 = "Log likelihood (logL)\n", ylab2 = "\u0394 logL\n")

ggsave('figures/fig6.png', baumWelch_plot_normal, width = 16, height = 8)

# calculate most probable state path using the Viterbi algorithm
states_normal <- viterbiAlgorithm(params[['S']],
                                  params[['E']],
                                  theta_hat_normal[['init']],
                                  theta_hat_normal[['tMat']],
                                  theta_hat_normal[['eMat']],
                                  gc_contents, 
                                  eFunc = dnorm)

# plot gaussian emission distribution after Baum-welch estimation
old_emission_density_plot <- plotEmissionDensity(normal_eMat, dnorm) +
  labs(title = 'Before Baum-Welch training', x = '') + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))
new_emission_density_plot <- plotEmissionDensity(theta_hat_normal[['eMat']], dnorm) +
  labs(title = 'After Baum-Welch training') + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))

combined_emission_density_plot <- ggarrange(old_emission_density_plot, new_emission_density_plot,
                                             common.legend = T, nrow = 2, ncol = 1)
combined_emission_density_plot_title <- text_grob('Gaussian emission distributions before and after Baum-Welch training', size = 14)
combined_emission_density_plot <- annotate_figure(combined_emission_density_plot, 
                                                  top = combined_emission_density_plot_title)

ggsave('figures/fig7.png', combined_emission_density_plot, width = 7, height = 6)

# put viterbi results into a data frame
viterbi_plot_normal <- plotViterbi(states_normal, gc_contents) + 
  labs(title = 'Continuous emission distribution (gaussian)', y = 'GC content\n') +
  theme(plot.title = element_text(hjust = 0.5))

# modify title of previous gc content plot
GC_content_plot_to_compare <- GC_content_plot + 
  labs(title = 'Discrete emission distribution') + 
  theme(plot.title = element_text(hjust = 0.5))

# combine plots
combined_viterbi_plot_normal <- ggarrange(GC_content_plot_to_compare, viterbi_plot_normal, 
                                   ncol = 1, nrow = 2, 
                                   common.legend = T, legend = 'top')
viterbi_plot_title_normal <- text_grob('Most likely hidden state sequence in discrete vs. continuous emission model\n', size = 26)
combined_viterbi_plot_normal <- annotate_figure(combined_viterbi_plot_normal, 
                                         top = viterbi_plot_title_normal)

ggsave('figures/fig8.png', combined_viterbi_plot_normal, width = 16, height = 14)
