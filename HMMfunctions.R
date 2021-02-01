# function for reading in a parameter file
readParamFile <- function(filepath){
  
  options(warn = -1)
  
  params <- readLines(filepath)
  params <- strsplit(params, '\\s+')
  params <- lapply(params, as.numeric)
  # check the dimensions of the transition/emission matrices are correct
  if(length(params[[4]]) != length(params[[1]]) ^ 2){
    stop(paste('Transition matrix is not the correct size, should be', 
               length(params[[1]]) ^ 2, 'elements'))
  }
  if(length(params[[5]]) != length(params[[1]]) * length(params[[2]])){
    stop(paste('Emission matrix is not the correct size, should be', 
               length(params[[1]]) * length(params[[2]]), 'elements'))
  }
  # set transition and emission matrices
  params[[4]] <- matrix(params[[4]], 
                        nrow = length(params[[1]]), 
                        ncol = length(params[[1]]), 
                        byrow = TRUE)
  params[[5]] <- matrix(params[[5]], 
                        nrow = length(params[[1]]), 
                        ncol = length(params[[2]]), 
                        byrow = TRUE)
  if(!all(apply(params[[4]], 1, sum) == 1)){
    stop('Rows of the transition matrix do not sum to 1. Ensure line 4 is formatted correctly.')
  }
  if(!all(apply(params[[5]], 1, sum) == 1)){
    stop('Rows of the emission matrix do not sum to 1. Ensure line 5 is formatted correctly.')
  }
  names(params) <- c('S', 'E', 'init', 'tMat', 'eMat')
  return(params)
}

# function for reading in a sequence file
readSequenceFile <- function(filepath){
  
  options(warn = -1)
  
  sequence <- readLines(filepath)
  sequence <- gsub('\\s+', '', sequence)
  sequence <- paste(sequence, collapse = '')
  sequence <- unlist(strsplit(sequence, ''))
  sequence <- as.numeric(sequence)
  if(any(!(sequence %in% params[['E']]))){
    stop('Some observations in the sequence file are not in the emission space. Please try again.')
  }
  return(sequence)
}

simulateSequence <- function(S, E, init, tMat, eMat, n, random_state = 123){
  
  set.seed(random_state)
  emissions <- vector(mode = typeof(E), length = n)
  states <- vector(mode = typeof(S), length = n)
  # get first state from initial distribution
  state <- sample(S, 1, prob = init)
  for(i in 1:n){
    # save state
    states[i] <- state
    # get new emission probs from emission matrix
    emission_probs <- eMat[which(S == state),]
    # sample from emission distribution to get emissions
    emissions[i] <- sample(E, 1, prob = emission_probs)
    # get transition probabilities from transition matrix
    transition_probs <- tMat[which(S == state),]
    # sample from states to get next state
    state <- sample(S, 1, prob = transition_probs)
  }
  results <- list('emissions' = emissions, 'states' = states)
  return(results)
  
}

# calculating a_hat in the scaled forward algorithm
calculate_a <- function(tMat, emission_probs, a, c){
  state_probs <- a %*% tMat
  new_a <- 1/c * state_probs * emission_probs
  return(new_a)
}

# calculating b_hat in the scaled backward algorithm
calculate_b <- function(tMat, emission_probs, b, c){
  weighted_bs <- emission_probs * b
  state_probs <- tMat %*% weighted_bs
  new_b <- as.numeric(1/c * state_probs)
  return(new_b)
}

# calculating c in the forward algorithm
calculate_c <- function(tMat, emission_probs, a){
  c <- as.numeric((a %*% tMat) %*% emission_probs)
  return(c)
}

scaledForwardAlgorithm <- function(E, init, tMat, eMat, obs, eFunc = NULL){
  
  # if an emission function provided, assume emissions matrix rows contain parameters for each state 
  if(!is.null(eFunc)){
    getEmissionProbs <- function(obs, eFunc, eMat){
      probs <- apply(eMat, 1, function(params) do.call(eFunc, as.list(c(obs, params))))
      return(probs)
    }
    # otherwise use the columns of the emission matrix
  } else {
    getEmissionProbs <- function(obs, eFunc, eMat) eMat[,which(E == obs)]
  }
  
  # initialise vector of c values
  c_values <- numeric(length = length(obs))
  # get c0 and a0 values
  emission_probs <- getEmissionProbs(obs[1], eFunc, eMat)
  c <- sum(emission_probs * init)
  c_values[1] <- c
  a <- 1/c * emission_probs * init
  # iterate through the sequence
  for(i in 2:length(obs)){
    emission_probs <- getEmissionProbs(obs[i], eFunc, eMat)
    c <- calculate_c(tMat, emission_probs, a)
    c_values[i] <- c
    a <- calculate_a(tMat, emission_probs, a, c)
  }
  # take logarithm to make computation easier
  c_values <- log(c_values)
  # log likelihood is the sum of log(c) values
  logL <- sum(c_values)
  return(logL)
}

# viterbi algorithm implementation
viterbiAlgorithm <- function(S, E, init, tMat, eMat, obs, eFunc = NULL){
  
  # if an emission function provided, assume emissions matrix rows contain parameters for each state 
  if(!is.null(eFunc)){
    getEmissionProbs <- function(obs, eFunc, eMat){
      probs <- apply(eMat, 1, function(params) do.call(eFunc, as.list(c(obs, params))))
      return(probs)
    }
    # otherwise use the columns of the emission matrix
  } else {
    getEmissionProbs <- function(obs, eFunc, eMat) eMat[,which(E == obs)]
  }
  
  # initialise objects for tracking viterbi variables
  phi_values <- matrix(nrow = length(obs), ncol = nrow(tMat))
  psi_values <- matrix(nrow = length(obs) - 1, ncol = nrow(tMat))
  
  # get phi_0
  phi <- log(init) + log(getEmissionProbs(obs[1], eFunc, eMat))
  phi_values[1,] <- phi

  # possible states
  states <- 1:nrow(tMat)
  
  # fill in the matrix of phi and psi values
  for(i in 2:length(obs)){

    # get the most likely previous state for each state s
    prev_states <- apply(tMat, 2, function(s) which.max(log(s) + phi_values[i-1,]) )
    psi_values[i-1,] <- prev_states
    
    # add phi values to the data frame
    phi_values[i,] <- mapply(function(s_prev, s){
      emission_prob <- getEmissionProbs(obs[i], eFunc, eMat)[s]
      phi <- log(emission_prob) + log(tMat[s_prev,s]) + phi_values[i-1,s_prev]
      return(phi)
    }, prev_states, states)
  }
  
  # get end state before traceback
  end_state <- which.max(phi_values[nrow(phi_values),])
  
  # initialise state path variable and add end state
  state_path <- integer(length = length(obs))
  state_path[length(obs)] <- end_state
  
  # traceback to find optimal state path
  for(i in (length(obs)-1):1){
    state_path[i] <- psi_values[i, state_path[i+1]]
  }
  
  results <- S[state_path]
  
  return(results)
  
}

# function for getting the parameters from a normal distribution
normalParams <- function(mean, var){
  return(c(mean, sqrt(var)))
}

# Baum-Welch parameter estimation
baumWelch <- function(E, init, tMat, eMat, obs, L_threshold = 0.01, eFunc = NULL, paramFunc = NULL){
  
  # if an emission function provided, assume emissions matrix rows contain parameters for each state 
  if(!is.null(eFunc)){
    getEmissionProbs <- function(obs, eFunc, eMat){
      probs <- apply(eMat, 1, function(params) do.call(eFunc, as.list(c(obs, params))))
      return(probs)
    }
  # otherwise use the columns of the emission matrix as emission probabilities
  } else {
    getEmissionProbs <- function(obs, eFunc, eMat) eMat[,which(E == obs)]
  }
  
  # set initial delta log-likelihood 
  deltaL <- Inf
  # vectors for tracking log-likelihoods and deltas
  logLs <- c()
  deltas <- c()
  # calculate log likelihood of initial model and append to vector
  logL <- scaledForwardAlgorithm(E, init, tMat, eMat, obs, eFunc = eFunc)
  logLs <- c(logLs, logL)
  count <- 0
  # iterate while the change in likelihood is above the specified threshold
  while(deltaL > L_threshold){
    
    # initialise vectors of a, b, and c values
    a_values <- matrix(nrow = length(obs), ncol = nrow(tMat))
    b_values <- matrix(nrow = length(obs), ncol = nrow(tMat))
    c_values <- numeric(length = length(obs))
    
    # get c0 and a0 values
    c <- sum(getEmissionProbs(obs[1], eFunc, eMat) * init)
    c_values[1] <- c
    a <- 1/c * getEmissionProbs(obs[1], eFunc, eMat) * init
    a_values[1,] <- a
    
    # forward algorithm
    for(i in 2:length(obs)){
      emission_probs <- getEmissionProbs(obs[i], eFunc, eMat)
      c <- calculate_c(tMat, emission_probs, a)
      c_values[i] <- c
      a <- calculate_a(tMat, emission_probs, a, c)
      a_values[i,] <- a
    }
    
    # get bn value
    b <- 1 / c_values[length(c_values)]
    b_values[length(obs),] <- b
    
    # backward algorithm
    for(i in (length(obs)-1):1){
      emission_probs <- getEmissionProbs(obs[i+1], eFunc, eMat)
      b <- calculate_b(tMat, emission_probs, b, c_values[i])
      b_values[i,] <- b
    }
    
    # update initial state probabilities
    init <- as.numeric(a_values[1,] * (tMat %*% (getEmissionProbs(obs[2], eFunc, eMat) * b_values[2,])))
    
    # update transition matrix
    new_tMat <- matrix(nrow = nrow(tMat), ncol = ncol(tMat))
    for(i in 1:nrow(tMat)){
      for(j in 1:ncol(tMat)){
        emission_probs <- sapply(obs[-1], function(ob) getEmissionProbs(ob, eFunc, eMat)[j])
        new_tMat[i,j] <- tMat[i,j] * sum(a_values[-nrow(a_values),i] * b_values[-1,j] * emission_probs)
      }
    }
    
    # update emission matrix
    new_eMat <- matrix(nrow = nrow(eMat), ncol = ncol(eMat))
    for(i in 1:nrow(eMat)){
      emission_probs <- t(sapply(obs[-1], function(ob) getEmissionProbs(ob, eFunc, eMat)))
      state_i_probs <- ((emission_probs * b_values[-1,]) %*% tMat[i,]) * a_values[-nrow(a_values),i]
      if(is.null(eFunc)){
        for(k in 1:ncol(eMat)){
          sum_over_vk <- sum(state_i_probs[which(obs[-length(obs)] == E[k])])
          normalised <- sum_over_vk / sum(new_tMat[i,])
          new_eMat[i,k] <- normalised
        }
      } else {
        sum_over_vk <- sum(state_i_probs * obs[-length(obs)])
        mean_estimate <- sum_over_vk / sum(new_tMat[i,])
        old_mean <- do.call(paramFunc, as.list(eMat[i,]))[1]
        sum_over_vk_sq <- sum(state_i_probs * (obs[-length(obs)] - old_mean)^2)
        variance_estimate <- sum_over_vk_sq / sum(new_tMat[i,])
        new_params <- paramFunc(mean_estimate, variance_estimate)
        new_eMat[i,] <- new_params
      }
    }
    
    # normalise transition matrix rows so they sum to 1
    tMat <- sweep(new_tMat, 1, rowSums(new_tMat), FUN = '/')
    # set new emission matrix
    eMat <- new_eMat

    # get likelihood of updated parameters and change in likelihood
    logL_prime <- scaledForwardAlgorithm(E, init, tMat, eMat, obs, eFunc = eFunc)
    deltaL <- logL_prime - logL
    logL <- logL_prime
    
    # save delta and logL
    deltas <- c(deltas, deltaL)
    logLs <- c(logLs, logL)
    
    # track progression
    count <- count + 1
    print(paste('Change in log likelihood of', round(deltaL, 2), 'in iteration', count))
  }
  
  print('Baum-Welch training finished.')
  
  results <- list('init' = init,
                  'tMat' = tMat,
                  'eMat' = eMat,
                  'deltas' = deltas,
                  'logLs' = logLs)
  
  return(results)
}

# function for plotting a discrete emission distribution
plotEmissionDistribution <- function(S, E, eMat){
  
  emission_dist <- data.frame(eMat)
  colnames(emission_dist) <- E
  emission_dist$state <- S
  emission_dist$state <- as.factor(emission_dist$state)
  emission_dist <- melt(emission_dist, id.vars = 'state', 
                        variable.name = 'emission')
  
  ggplot(emission_dist, aes(x = emission, y = value, fill = state)) + 
    geom_col(position = 'dodge', colour = 'black', width = 0.5) + 
    labs(x = '\nEmission', y = 'Probability\n', fill = 'Hidden state') +
    theme_minimal() + theme(text = element_text(size = 10))
  
}

# function for plotting likelihood evolution during Baum Welch training
plotBaumWelch <- function(logLs, deltas, 
                          title = 'Change in log-likelihood during Baum-Welch training process\n',
                          ylab1 = 'Log likelihood (logL)\n',
                          ylab2 = '\u0394 logL\n'){
  
  # plot change in likelihood
  baumWelch_results <- data.frame(logL = logLs, delta = c(0, deltas))
  baumWelch_results$iteration <- 1:nrow(baumWelch_results)
  
  logL_plot <- ggplot(baumWelch_results, aes(x = iteration, y = logL)) + 
    geom_point(size = 3) + geom_line() + 
    labs(x = '\nIteration', y = ylab1) + 
    theme_minimal() + 
    theme(text = element_text(size = 20))
  
  delta_logL_plot <- ggplot(baumWelch_results[-1,], aes(x = iteration, y = delta)) + 
    geom_point(size = 3) + geom_line() + 
    labs(x = '\nIteration', y = ylab2) + 
    theme_minimal() + 
    theme(text = element_text(size = 20))
  
  baumWelch_plot_title <- text_grob(title, size = 24)
  combined_baumWelch_plot <- ggarrange(logL_plot, delta_logL_plot, nrow = 1)
  combined_baumWelch_plot <- annotate_figure(combined_baumWelch_plot, top = baumWelch_plot_title)
  
  return(combined_baumWelch_plot)
}

# function for plotting observations with Viterbi sequence
plotViterbi <- function(states, observations, n = 120){
  
  
  results <- data.frame(state = as.factor(states), emissions = observations,
                        sequence = 1:length(states))
  # take first n emissions to plot
  results_to_plot <- results[1:n,]
  
  viterbi_plot <- ggplot(results_to_plot, aes(x = sequence, y = emissions)) +
    geom_point(aes(colour = state), size = 3) + geom_line(size = 0.5, alpha = 0.25) + 
    theme_minimal() + 
    theme(text = element_text(size = 20),
          plot.subtitle = element_text(hjust = 0.5),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 22),
          plot.margin = unit(c(1,1,1,1), "lines")) + 
    labs(x = '\nChromosome III position', y = 'Emission\n', colour = 'Hidden state') +
    scale_x_continuous(breaks = seq(20, 120, 20), labels = seq(20, 120, 20))
  
  return(viterbi_plot)
  
}

# plot densities for all the hidden states, using parameters in eMat in the eFunc
plotEmissionDensity <- function(eMat, eFunc, n_points = 1000){
  
  density_functions = vector(mode = 'list', length = nrow(eMat))
  densities <- data.frame(x = seq(0, 1, length.out = n_points))
  
  for(i in 1:nrow(eMat)){
    dFunc <- function(x) do.call(eFunc, as.list(c(x, eMat[i,])))
    colname <- as.character(i-1)
    densities[[colname]] <- sapply(densities$x, dFunc)
  }
  
  densities <- reshape2::melt(densities, id.vars = 'x')
  densities$variable <- as.factor(densities$variable)
  baumWelch_distribution_plot <- ggplot(densities, aes(x = x, y = value, colour = variable)) + 
    geom_line() + theme_minimal() +
    labs(x = '\nGC content', colour = 'Hidden state', y = 'Density\n') +
    theme(text = element_text(size = 10))
  
  return(baumWelch_distribution_plot)
  
}
