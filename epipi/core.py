import stan
import arviz as az
from epipi import Prior
import nest_asyncio
# import rpy2.robjects as robjects
nest_asyncio.apply()


class InferenceController():

    def __init__(self, theta, omega, kernel=[], prior=None, random_seed=None, **kwargs):
        self.theta = theta
        self.omega = omega
        self.model = self._setup_model(prior)
        self.random_seed = random_seed
        self.params_dict = kwargs
        if len(kernel) > len(theta):
            raise ValueError("Kernel size cannot be greater than the length of the local prevalence list")
        self.kernel = kernel

    
    def _setup_model(self, prior_list=None):
        model_pre_string = """
functions {
    real effective_no_infectives (
        int N, int S, int t, int k, array [] real aK, array [] real aI, array [] real aOmega) {
            real mean;
            if (t < N - k + 2) {
                if(t > S) {
                    mean = (
                        dot_product(aI[(t-S):(t-1)], aOmega));
                }
                else {
                    mean = (
                        dot_product(aI[:(t-1)], aOmega[(S-t+2):]));
                }
            }
            else {
                if(t > S) {
                    mean = (
                        dot_product(aI[(t-S):(N-k)], aOmega[1:(S-t+N-k + 1)]) + dot_product(aK[1:t-N+k-1], aOmega[(S-t+N-k+2):]));
                }
                else {
                    mean = (
                    dot_product(aI[:(N-k)], aOmega[(S-t+2):(S-t+N-k + 1)]) + dot_product(aK[1:t-N+k-1], aOmega[(S-t+N-k+2):]));
                }
            }
            return mean;
    }
}
data {
    int N; // number of days
    int S; // length omega
    int k; // kenel size
    array [N] real Theta; // local prevalence for N days
    array [S] real revOmega; // reversed PCR positive probability
    array [k] real K; // data in the kernel
}
parameters {
    array [N - k] real<lower=0> I; // local incidence for N days
    real<lower=0> sigma; // variance of the log-normal distribution
}
model {
    for(t in 1:N) {
        Theta[t] ~ normal (
            effective_no_infectives(
                N, S, t, k, K, I, revOmega), sigma); // likelihood
    }
    """
        model_post_string = """}
"""
        if prior_list is None:
            model_prior_string = """sigma ~ normal(0, 1); // prior of sigma
    for(t in 1:N - k) {
        I[t] ~ normal(5000, 1700); // prior of I
    }"""
        else:
            for prior in prior_list:
                if not isinstance(prior, Prior):
                    raise ValueError("prior_list must contain only Prior objects")
            I_prior_list = prior_list[:-1]
            sigma_prior = prior_list[-1]
            model_prior_string_list = []
            for ind, prior in enumerate(I_prior_list):
                prefix = "        I[{}] ~ ".format(ind+1)
                model_prior_string_list.append(prefix + prior.string() + ";")
            model_prior_string_list.append("        sigma ~ " + sigma_prior.string() + ";")
            model_prior_string = "\n".join(model_prior_string_list)
        model_string = model_pre_string + model_prior_string + "\n" + model_post_string
        return model_string
    
    def _compile_data(self, theta, omega, kernel=None):
        data = {
            "N": len(theta),
            "S": len(omega),
            "k": len(kernel),
            "Theta": theta,
            "revOmega": omega[::-1],
            "K": kernel
        }
        return data
    
    def _compile_model(self, model, data):
        if self.random_seed is not None:
            model = stan.build(model, data, random_seed=self.random_seed)
        else:
            model = stan.build(model, data)
        return model
    
    def run(self):
        data = self._compile_data(self.theta, self.omega, self.kernel)
        model = self._compile_model(self.model, data)
        fit = model.sample(**self.params_dict)
        samples = az.from_pystan(
            posterior=fit,
            observed_data='Theta',  # Use 'Theta' as the observed data
            coords={'observation': list(range(len(self.theta) - len(self.kernel)))},  # 'covariate' removed
            dims={'Theta': ['observation'], 'I': ['observation']}  # Adjusted dims
        )
        return fit, samples


def predict_incidences(theta, omega, prior=None, kernel=[], random_seed=None, **kwargs):
    inference_controller = InferenceController(
        theta=theta, omega=omega, kernel=kernel, prior=prior, random_seed=random_seed, **kwargs)
    _, samples = inference_controller.run()
    df = az.summary(samples)
    predicted_mean = df['mean'].values
    predicted_std = df['sd'].values
    return predicted_mean, predicted_std



# class InferenceController():

#     def __init__(self, theta, omega, prior=None, random_seed=None, **kwargs):
#         self.theta = theta
#         self.omega = omega
#         self.model = self._setup_model(prior)
#         self.random_seed = random_seed
#         self.params_dict = kwargs
    
#     def _setup_model(self, prior_list=None):
#         model_pre_string = """
# functions {
#     real effective_no_infectives (
#         int N, int S, int t, array [] real aI, array [] real aOmega) {
#             real mean;
#             if(t > S) {
#                 mean = (
#                     dot_product(aI[(t-S):(t-1)], aOmega));
#             }
#             else {
#                 mean = (
#                     dot_product(aI[:(t-1)], aOmega[(S-t+2):]));
#             }
#             return mean;
#     }
# }
# data {
#     int N; // number of days
#     int S; // length omega
#     array [N] real Theta; // local prevalence for N days
#     array [S] real revOmega; // reversed PCR positive probability
# }
# parameters {
#     array [N] real<lower=0> I; // local incidence for N days
#     real<lower=0> sigma; // variance of the log-normal distribution
# }
# model {
#     for(t in 1:N) {
#         Theta[t] ~ normal (
#             effective_no_infectives(
#                 N, S, t, I, revOmega), sigma); // likelihood
#     }
#     """
#         model_post_string = """}
# """
#         if prior_list is None:
#             model_prior_string = """sigma ~ normal(0, 1); // prior of sigma
#     for(t in 1:N) {
#         I[t] ~ normal(5000, 1700); // prior of I
#     }"""
#         else:
#             for prior in prior_list:
#                 if not isinstance(prior, Prior):
#                     raise ValueError("prior_list must contain only Prior objects")
#             I_prior_list = prior_list[:-1]
#             sigma_prior = prior_list[-1]
#             model_prior_string_list = []
#             for ind, prior in enumerate(I_prior_list):
#                 prefix = "        I[{}] ~ ".format(ind+1)
#                 model_prior_string_list.append(prefix + prior.string() + ";")
#             model_prior_string_list.append("        sigma ~ " + sigma_prior.string() + ";")
#             model_prior_string = "\n".join(model_prior_string_list)
#         model_string = model_pre_string + model_prior_string + "\n" + model_post_string
#         return model_string
    
#     def _compile_data(self, theta, omega):
#         data = {
#             "N": len(theta),
#             "S": len(omega),
#             "Theta": theta,
#             "revOmega": omega[::-1]
#         }
#         return data
    
#     def _compile_model(self, model, data):
#         if self.random_seed is not None:
#             model = stan.build(model, data, random_seed=self.random_seed)
#         else:
#             model = stan.build(model, data)
#         return model
    
#     def run(self):
#         data = self._compile_data(self.theta, self.omega)
#         model = self._compile_model(self.model, data)
#         fit = model.sample(**self.params_dict)
#         samples = az.from_pystan(
#             posterior=fit,
#             observed_data='Theta',  # Use 'Theta' as the observed data
#             coords={'observation': list(range(len(self.theta)))},  # 'covariate' removed
#             dims={'Theta': ['observation'], 'I': ['observation']}  # Adjusted dims
#         )
#         return fit, samples





# class OptimizerController(InferenceController):

#     def _r_code(self):
#         r_code = """
# library(rstan)
# library(ggplot2)

# model_code <- '%s'

# poisson_model <- stan_model(model_code = model_code)

# theta <- %s
# N <- length(theta)
# omega <- %s
# S <- length(omega)

# data <- list(N=N, S=S, Theta=theta, revOmega=rev(omega))

# # Function to return best optimized fit from Stan

# n_opt <- 1

# multiple_opt <- function(stan_data, n_opt){
  
#   n_opt <- n_opt
  
#   # empty lists to store multiple fits and log posteriors
#   fit_temp <- vector(mode = "list", length = length(n_opt))
#   lp_temp <- as.data.frame(matrix(ncol = 1, nrow = n_opt))
  
#   # runnning optimize for n_opt times for a given LTLA
#   for(i in 1:n_opt){
#     fit_temp[[i]] <- optimizing(poisson_model, data=stan_data, as_vector = FALSE,
#                                 init='random', algorithm = 'Newton')
    
#     lp_temp[i,] <- fit_temp[[i]]$value
    
#     # extract optimized fit which has the highest log posterior probability
#     index_opt <- which.max(lp_temp$V1)
#     fit_optimized_output <- fit_temp[[index_opt]]
#   }
  
#   return(fit_optimized_output)
  
# }

# result <- multiple_opt(data, n_opt)

# result
# """ % (self.model, robjects.FloatVector(self.theta).r_repr(), robjects.FloatVector(self.omega).r_repr())
#         return r_code

#     def run(self):
#         r_code = self._r_code()
#         # Execute the R code from Python
#         result = robjects.r(r_code)
#         return result
