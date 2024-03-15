class Prior():
    '''
    The abstract class for prior

    Parameters:
    -----------

    parameters : any
        The parameters of the prior distribution
    '''
    def __init__(self, parameters):
        self.parameters = parameters

    def string(self):
        return "dist({})".format(self.parameters)


class NormalPrior(Prior):
    '''
    A class for normal prior, can generate the corresponding normal prior string for the stan model

    Parameters:
    -----------

    mu : float
        The mean of the normal distribution
    sigma : float
        The standard deviation of the normal distribution
    '''
    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma

    def string(self):
        '''
        This is to output the normal prior string for the stan model
        '''
        return "normal({}, {})".format(self.mu, self.sigma)


class UniformPrior(Prior):
    '''
    A class for uniform prior, can generate the corresponding normal prior string for the stan model

    Parameters:
    -----------

    lower : float
        The lower bound of the uniform distribution
    upper : float
        The upper bound of the uniform distribution
    '''
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def string(self):
        '''
        This is to output the normal prior string for the stan model
        '''
        return "uniform({}, {})".format(self.lower, self.upper)
