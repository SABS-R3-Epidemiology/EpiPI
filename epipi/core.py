import stan
import arviz as az
from epipi import Prior
import nest_asyncio
# import rpy2.robjects as robjects
nest_asyncio.apply()


class InferenceController():

    def __init__(self, theta, omega, init_infection=None, kernel=[], prior=None, random_seed=None, **kwargs):
        if init_infection is not None:
            self.model, self.data = self._setup_init_infection_model(theta, omega, init_infection, kernel, prior)
            self._n_params = len(theta) - len(kernel) - 1
        else:
            self.model, self.data = self._setup_model(theta, omega, kernel, prior)
            self._n_params = len(theta) - len(kernel)
        self.random_seed = random_seed
        self.params_dict = kwargs
        if len(kernel) > len(theta):
            raise ValueError("Kernel size cannot be greater than the length of the local prevalence list")

    def _setup_init_infection_model(self, theta, omega, I0, kernel, prior_list=None):
        model_pre_string = """
functions {
    real effective_no_infectives (
        int N, int S, int t, int k, real I0, array [] real aK, array [] real aI, array [] real aOmega) {
            real mean;
            if (t > 1) {
                if (t < N - k + 2) {
                    if(t > S + 1) {
                        mean = (
                            dot_product(aI[(t-S-1):(t-2)], aOmega));
                    }
                    else {
                        mean = (
                            dot_product(aI[:(t-2)], aOmega[(S-t+3):]) + I0 * 0);
                    }
                }
                else {
                    if(t > S + 1) {
                        mean = (
                            dot_product(aI[(t-S-1):(N-k-1)], aOmega[1:(S-t+N-k + 1)]) + dot_product(aK[1:t-N+k-1], aOmega[(S-t+N-k+2):]));
                    }
                    else {
                        mean = (
                        dot_product(aI[:(N-k-1)], aOmega[(S-t+3):(S-t+N-k + 1)]) + I0 * 0 + dot_product(aK[1:t-N+k-1], aOmega[(S-t+N-k+2):]));
                    }
                }
            }
            else {
                mean = I0;
            }
            return mean;
    }
}
data {
    int N; // number of days
    int S; // length omega
    int k; // kenel size
    real I0; // initial infection
    array [N] real Theta; // local prevalence for N days
    array [S] real revOmega; // reversed PCR positive probability
    array [k] real K; // data in the kernel
}
parameters {
    array [N - k - 1] real<lower=0> I; // local incidence for N days
    real<lower=0> sigma; // variance of the log-normal distribution
}
model {
    for(t in 1:N) {
        Theta[t] ~ normal (
            effective_no_infectives(
                N, S, t, k, I0, K, I, revOmega), sigma); // likelihood
    }
    """
        model_post_string = """}
"""
        if prior_list is None:
            model_prior_string = """sigma ~ normal(0, 1); // prior of sigma
    for(t in 1:N - k - 1) {
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

        data = {
            "N": len(theta),
            "S": len(omega),
            "k": len(kernel),
            "Theta": theta,
            "I0": I0,
            "revOmega": omega[::-1],
            "K": kernel
        }
        return model_string, data
    
    def _setup_model(self, theta, omega, kernel, prior_list=None):
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

        data = {
            "N": len(theta),
            "S": len(omega),
            "k": len(kernel),
            "Theta": theta,
            "revOmega": omega[::-1],
            "K": kernel
        }
        return model_string, data
    
    def _compile_model(self, model, data):
        if self.random_seed is not None:
            model = stan.build(model, data, random_seed=self.random_seed)
        else:
            model = stan.build(model, data)
        return model
    
    def run(self):
        model = self._compile_model(self.model, self.data)
        fit = model.sample(**self.params_dict)
        samples = az.from_pystan(
            posterior=fit,
            observed_data='Theta',  # Use 'Theta' as the observed data
            coords={'observation': list(range(self._n_params))},  # 'covariate' removed
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
