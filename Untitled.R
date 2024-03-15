library(rstan)
library(ggplot2)

poisson_model <- stan_model("/Users/yunliqi/EpiOS/poisson_model.stan")

theta <- read.csv("/Users/yunliqi/EpiOS/y_interpolated.csv")$y_interpolated
N <- length(theta)
omega <- read.csv("/Users/yunliqi/EpiOS/omega.csv")$omega
S <- length(omega)

data <- list(N=N, S=S, Theta=theta, revOmega=rev(omega))

# Function to return best optimized fit from Stan

n_opt <- 1

multiple_opt <- function(stan_data, n_opt, init_list){
  
  n_opt <- n_opt
  
  # empty lists to store multiple fits and log posteriors
  fit_temp <- vector(mode = "list", length = length(n_opt))
  lp_temp <- as.data.frame(matrix(ncol = 1, nrow = n_opt))
  
  # runnning optimize for n_opt times for a given LTLA
  for(i in 1:n_opt){
    fit_temp[[i]] <- optimizing(poisson_model, data=stan_data, as_vector = FALSE,
                                init='random', algorithm = 'Newton')
    
    lp_temp[i,] <- fit_temp[[i]]$value
    
    # extract optimized fit which has the highest log posterior probability
    index_opt <- which.max(lp_temp$V1)
    fit_optimized_output <- fit_temp[[index_opt]]
  }
  
  return(fit_optimized_output)
  
}

init_list <- 100

result <- multiple_opt(data, n_opt, init_list)

est_theta <- data.frame(x = 1:86, y = as.vector(result$theta_tilde))

ggplot(data = est_theta) + geom_line(aes(x = x, y = y))


