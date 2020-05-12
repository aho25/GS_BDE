#optimization after filter method
BDE_optim <- function(PHENO, MARKERS.pool, feature_pool.names, p_of_feature, D, CROSSVAL, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS) {
  ### Create initial population - Generation == 1
  Population <- list()
  Population$G1$X <- matrix(nrow = NP, ncol = D)
  colnames(Population$G1$X) <- feature_pool.names
  Population$G1$X[1,] <- rep.int(0, D)
  for (i in 2:NP) {
    for (j in 1:D) {
      if (p_of_feature[j] < runif(1)) {
        Population$G1$X[i,j] <- 1
      } else {
        Population$G1$X[i,j] <- 0
      }
    }
  }
  
  ### Fitness evaluation ###
  Population$G1$Fitness <- vector(length = NP)
  Population$G1$Fitness <- sapply(1:NP, function(i) {
    #	cat(i, 'iteration', '\n')
    feature_fit.idx <- which(Population$G1$X[i,] == 0)# column index of features used in algorithm 
    MARKERS.individual <- MARKERS.pool[,feature_fit.idx]
    Population$G1$Fitness[i] <- OBJFUNC(PHENO, MARKERS.individual, OBJFUNC.ARGS, CROSSVAL, SEEDRNG)
    return(Population$G1$Fitness[i])
  })
  
  ### Find the best individual from the population
  Population$G1$x_best <- which.max(Population$G1$Fitness)
  cat('1 generation best fitness =', Population$G1$Fitness[Population$G1$x_best], '\n')
  cat('1 generation best fitness =', Population$G1$Fitness[Population$G1$x_best], '\n', 
      file = paste0("results/", AnalyseName, "_Fithess_by_Generation.txt"), append = T)
  
  ### Count hamming based population diversity for initial population
  gtype <- vector()
  l <- 2
  for (i in 1:(NP-1)) {
    for (j in l:NP) {
      ham <- sum(abs(Population$G1$X[i,] - Population$G1$X[j,]))
      gtype <- sum(gtype, ham)
    }
    l <- l + 1
  }
  Population$G1$gtype <- (gtype/D) / ((NP*(NP-1))/2)
  
  ### Go on with BDE ###
  
  for (G in 2:GENERATION) {
    cat(G, 'generation ')
    cat(G, 'generation ', file = paste0("results/", AnalyseName, "_Fithess_by_Generation.txt"), append = T)
    Population[[paste0('G',G)]]$X <- matrix(nrow = NP, ncol = D)
    colnames(Population[[paste0('G',G)]]$X) <- feature_pool.names
    Population[[paste0('G',G)]]$Fitness <- vector(length = NP)
    
    ### Count minimum change value - C_min
    C_min <- ceiling(13*(1 - (G-1)/GENERATION)) + 4
    
    X_and_Fitness <- lapply(1:NP, function(i) {
      target <- Population[[paste0('G',G-1)]]$X[i,]
      ### Create binary mutation operator
      r1 <- sample(c(1:NP)[-i], 1)
      r2 <- sample(c(1:NP)[-c(i,r1)], 1)
      dif <- vector(length = D)
      for (j in 1:D) {
        if (Population[[paste0('G',G-1)]]$X[r1,j] == Population[[paste0('G',G-1)]]$X[r2,j]) {
          dif[j] <- 0
        } else {
          dif[j] <- Population$G1$X[r1,j]
        }
      }
      ### Create mutant vector
      mutant <- vector(length = D)
      for (j in 1:D) {
        if (dif[j] == 1 && runif(1) < MUTFACTOR) {
          mutant[j] <- 1
        } else {
          mutant[j] <- Population[[paste0('G',G-1)]]$X[Population[[paste0('G',G-1)]]$x_best,j]
        }
      }
      ### Create trial vector
      trial <- vector(length = D)
      for (j in 1:D) {
        if (runif(1) < CR) {
          trial[j] <- mutant[j]
        } else {
          trial[j] <- target[j]
        }
      }
      ### Modification of the trial vector with C_min
      modif <- sample(1:D, C_min)
      trial[modif[j]] <- abs(trial[modif[j]] - 1)
      
      ### Fitness evaluation ###
      feature_fit.idx <- which(Population[[paste0('G',G-1)]]$X[i,] == 0) # Features acting in fitness evaluation
      MARKERS.individual <- MARKERS.pool[,feature_fit.idx]
      trial.fit <- OBJFUNC(PHENO, MARKERS.individual, OBJFUNC.ARGS, CROSSVAL, SEEDRNG)
      
      ### Write best (target or trial) in the generation
      if (trial.fit > Population[[paste0('G',G-1)]]$Fitness[i]) {
        return(list(trial, trial.fit))
      } else {
        return(list(target, Population[[paste0('G',G-1)]]$Fitness[i]))
      }
    })
    for (i in 1:NP) {
      Population[[paste0('G',G)]]$X[i,] <- X_and_Fitness[[i]][[1]]
      Population[[paste0('G',G)]]$Fitness[i] <- X_and_Fitness[[i]][[2]]
    }
    
    ### Find the best individual from the population
    Population[[paste0('G',G)]]$x_best <- which.max(Population[[paste0('G',G)]]$Fitness)
    
    ### Count hamming based population diversity & CR
    gtype <- vector()
    l <- 2
    for (i in 1:(NP-1)) {
      for (j in l:NP) {
        ham <- sum(abs(Population[[paste0('G',G)]]$X[i,] - Population[[paste0('G',G)]]$X[j,]))
        gtype <- sum(gtype, ham)
      }
      l <- l + 1
    }
    Population[[paste0('G',G)]]$gtype <- (gtype/D) / ((NP*(NP-1))/2)
    CR <- Population[[paste0('G',G)]]$gtype / Population$G1$gtype
    CR <- 1 - min(0.9, max(0.1, CR))
    cat('best fitness =', Population[[paste0('G',G)]]$Fitness[Population[[paste0('G',G)]]$x_best], '\n')
    cat('best fitness =', Population[[paste0('G',G)]]$Fitness[Population[[paste0('G',G)]]$x_best], '\n', 
        file = paste0("results/", AnalyseName, "_Fithess_by_Generation.txt"), append = T)
  }
  return(Population)
}
