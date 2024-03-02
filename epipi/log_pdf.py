import pints
import warnings
import numpy as np


class PredictIncidenceLogPDF(pints.LogPDF):
    def __init__(self, theta, omega):
        self.theta = np.asarray(theta)
        self.omega = np.asarray(omega)
        self.rev_omega = np.asarray(omega)[::-1]
        self.N = len(theta)
        self.S = len(omega)
        self._n_parameters = self.N + 1

    def _likelihood_precomp(self, I, t):
        '''
        Precompute quantities that are used in the likelihood function
        and its gradient.

        Parameters
        ----------
        I : np.ndarray
            The incidence vector.
        t : int
            The time point at which to compute the likelihood.

        '''
        if t >= self.S:
            return self.theta[t] - (self.rev_omega * I[(t - self.S+1):(t+1)]).sum()
        else:
            return self.theta[t] - (self.rev_omega[-(t+1):] * I[:(t + 1)]).sum()

    def __call__(self, parameters):
        '''
        Sigma is the variance
        '''
        try:
            I = parameters[:self.N]
            sigma = parameters[-1]
            log_likelihood = 0
            for t in range(self.N):
                log_likelihood += -0.5 * np.log(sigma)
                log_likelihood += -0.5 * (self._likelihood_precomp(I, t))**2 / sigma
            return log_likelihood
    
        except ValueError:  # pragma: no cover
            warnings.warn('RuntimeWarning: for x={}, the likelhood retrurned -infinity'.format(parameters))
            return -np.inf
    
    def evaluateS1(self, parameters):
        try:
            sigma = parameters[-1]
            I = parameters[:self.N]
            log_likelihood = self(parameters)
            dlog_likelihood = np.zeros(self._n_parameters)
            pre_comps = [self._likelihood_precomp(I, t) for t in range(self.N)]
            for t in range(self.N):
                pre_comp = pre_comps[t]
                for k in range(self.N):
                    if t - k >= 0 and t - k < self.S:
                        dlog_likelihood[k] += self.omega[t - k] * pre_comp / sigma
            
            dlog_likelihood[-1] = -0.5 * self.N / sigma + 0.5 * np.sum((pre_comps[t])**2 for t in range(self.N)) / sigma**2
            return log_likelihood, dlog_likelihood
    
        except ValueError:  # pragma: no cover
            warnings.warn('RuntimeWarning: for x={}, the likelhood retrurned -infinity'.format(parameters))
            return -np.inf, [-np.inf] * self.n_parameters()

    def n_parameters(self):
        return self._n_parameters
