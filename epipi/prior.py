class Prior():
    def __init__(self, parameters):
        self.parameters = parameters

    def string(self):
        return "dist({})".format(self.parameters)

class NormalPrior(Prior):
    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma

    def string(self):
        return "normal({}, {})".format(self.mu, self.sigma)

class UniformPrior(Prior):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def string(self):
        return "uniform({}, {})".format(self.lower, self.upper)