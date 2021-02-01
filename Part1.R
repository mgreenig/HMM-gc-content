source('HMMfunctions.R')

if(!interactive()){
  args <- commandArgs(trailingOnly = T)
  args <- unlist(strsplit(args, '\\s+'))
  parameter_file <- args[1]; observation_file <- args[2]
} else {
  parameter_file <- readline(prompt = 'Enter the filepath for the parameter file to be used: ')
  observation_file <- readline(prompt = 'Enter the filepath for the file of observations to be used: ')
}

params <- readParamFile(parameter_file)

simulated_seq <- simulateSequence(params[['S']],
                               params[['E']],
                               params[['init']],
                               params[['tMat']],
                               params[['eMat']],
                               n = 115)

simulation <- as.data.frame(simulated_seq)
simulation$states <- as.factor(simulation$states)
simulation$sequence <- 1:nrow(simulation)

simulation_plot <- ggplot(simulation, aes(x = sequence, y = emissions)) + 
  geom_line(size = 0.5, alpha = 0.25) + geom_point(aes(colour = states), size = 3) + 
  theme_minimal() + 
  theme(text = element_text(size = 18), plot.title = element_text(size = 20, hjust = 0.5)) + 
  labs(x = '\nSequence position', y = 'Emission\n', colour = 'Hidden state', 
       title = 'Emission sequence simulated from HMM\n') +
  xlim(0, 115)
ggsave('figures/fig1.png', simulation_plot, width = 14, height = 7)

print(paste('Simulated sequence:', paste(simulated_seq$emissions, collapse = '')))

observations <- readSequenceFile(observation_file)

logL <- scaledForwardAlgorithm(params[['E']],
                               params[['init']],
                               params[['tMat']],
                               params[['eMat']],
                               observations)

print(paste('Likelihood of the inputted sequence, given the model:', exp(logL)))
