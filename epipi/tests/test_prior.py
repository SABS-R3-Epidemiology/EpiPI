import unittest
from prior import Prior, NormalPrior, UniformPrior

class TestPriors(unittest.TestCase):

    def test_normal_prior_initialization(self):
        mu, sigma = 0, 1
        normal_prior = NormalPrior(mu, sigma)
        self.assertEqual(normal_prior.mu, mu)
        self.assertEqual(normal_prior.sigma, sigma)

    def test_normal_prior_string(self):
        normal_prior = NormalPrior(0, 1)
        expected_string = "normal(0, 1)"
        self.assertEqual(normal_prior.string(), expected_string)

    def test_uniform_prior_initialization(self):
        lower, upper = 0, 10
        uniform_prior = UniformPrior(lower, upper)
        self.assertEqual(uniform_prior.lower, lower)
        self.assertEqual(uniform_prior.upper, upper)

    def test_uniform_prior_string(self):
        uniform_prior = UniformPrior(0, 10)
        expected_string = "uniform(0, 10)"
        self.assertEqual(uniform_prior.string(), expected_string)

    def test_prior_abstract_class(self):
        parameters = [1, 2, 3]
        prior = Prior(parameters)
        self.assertEqual(prior.parameters, parameters)

    def test_prior_string(self):
        parameters = [1, 2, 3]
        prior = Prior(parameters)
        expected_string = "dist([1, 2, 3])"
        self.assertEqual(prior.string(), expected_string)

if __name__ == '__main__':
    unittest.main()
