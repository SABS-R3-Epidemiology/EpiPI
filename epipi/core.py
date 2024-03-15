import stan
import arviz as az
from epipi import Prior
import nest_asyncio
nest_asyncio.apply()


class InferenceController():
    '''
    A class to contain the inference process

    Parameters:
    -----------

    theta : list
        A list containing the prevalence data
    omega : list
        A list containing the PCR positive probability
    kernel : list
        A list containing the kernel data

        Kernel is defined since the model will not have enough information to infer the data at the end of time series.
        The values in the kernel is used to overwrite the poorly predicted values at the end of the time series.
    prior : list of Prior objects
        A list containing customised prior for parameters in inference
    random_seed : int
        A random seed for reproducibility
    kwargs : dict
        Keyword arguments for the stan model

    '''

    def __init__(self, theta, omega, kernel=[], prior=None, random_seed=None, **kwargs):
        # When all data provided, directly define the model needed
        self.model, self.data = self._setup_model(theta, omega, kernel, prior)
        self._n_params = len(theta) - len(kernel)
        self.random_seed = random_seed
        self.param_dict = kwargs
        if len(kernel) > len(theta):
            raise ValueError("Kernel size cannot be greater than the length of the local prevalence list")

    def _setup_model(self, theta, omega, kernel, prior_list=None):
        '''
        A private method to setup the stan model string, put the cumsomised prior into the model string.

        '''
        # The following is the main model string
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
            # If no prior provided, use the default ones
            model_prior_string = """sigma ~ normal(0, 1); // prior of sigma
    for(t in 1:N - k) {
        I[t] ~ normal(5000, 1700); // prior of I
    }"""
        else:
            for prior in prior_list:
                # Verify the prior is the right type
                if not isinstance(prior, Prior):
                    raise ValueError("prior_list must contain only Prior objects")
            # Separate the I priors and the sigma prior
            I_prior_list = prior_list[:-1]
            sigma_prior = prior_list[-1]
            model_prior_string_list = []
            # Formulate the prior string
            for ind, prior in enumerate(I_prior_list):
                prefix = "        I[{}] ~ ".format(ind + 1)
                model_prior_string_list.append(prefix + prior.string() + ";")
            # Add the sigma prior string
            model_prior_string_list.append("        sigma ~ " + sigma_prior.string() + ";")
            model_prior_string = "\n".join(model_prior_string_list)
        model_string = model_pre_string + model_prior_string + "\n" + model_post_string

        # Define the data dictionary that will be used in compiling the model
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
        '''
        A private method to compile the stan model, which put the model and data together
        '''
        # Compile the model, depending on whether a random seed is provided
        if self.random_seed is not None:
            model = stan.build(model, data, random_seed=self.random_seed)
        else:
            model = stan.build(model, data)
        return model

    def run(self):
        '''
        Method to run the MCMC sampling and return the fit object.

        Output:
        -------

        fit : stan fit object
            The fit object returned by the MCMC sampling, how to use this object can be found
            in the documentation of stan
        samples : arviz InferenceData object
            The samples obejct returned by arviz package, detailed documentation can be found
            in the documentation of arviz
        '''
        # Use the above methods to compile the model, then do sampling
        model = self._compile_model(self.model, self.data)
        fit = model.sample(**self.param_dict)
        samples = az.from_pystan(
            posterior=fit,
            observed_data='Theta',  # Use 'Theta' as the observed data
            coords={'observation': list(range(self._n_params))},  # 'covariate' removed
            dims={'Theta': ['observation'], 'I': ['observation']}  # Adjusted dims
        )
        return fit, samples


def predict_incidences(theta, omega, prior=None, kernel=[], random_seed=None, **kwargs):
    '''
    A function that can be run directly to generate inference results using the InferenceController class

    Parameters:
    -----------

    theta : list
        A list containing the prevalence data
    omega : list
        A list containing the PCR positive probability
    kernel : list
        A list containing the kernel data

        Kernel is defined since the model will not have enough information to
        infer the data at the end of time series.

        The values in the kernel is used to overwrite the poorly predicted
        values at the end of the time series.
    prior : list of Prior objects
        A list containing customised prior for parameters in inference
    random_seed : int
        A random seed for reproducibility
    kwargs : dict
        Keyword arguments for the stan model
    '''
    inference_controller = InferenceController(
        theta=theta, omega=omega, kernel=kernel, prior=prior, random_seed=random_seed, **kwargs)
    _, samples = inference_controller.run()
    df = az.summary(samples)
    predicted_mean = df['mean'].values
    predicted_std = df['sd'].values
    return predicted_mean, predicted_std
