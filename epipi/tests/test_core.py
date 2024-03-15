import unittest
from core import InferenceController, predict_incidences
from epipi import NormalPrior
import arviz as az
# import numpy.testing as npt
# import numpy as np


class TestInferenceController(unittest.TestCase):
    def setUp(self):
        # Setup some test data that can be used across multiple tests
        self.theta = [0.1, 0.2, 0.3]
        self.omega = [0.9, 0.8, 0.7]
        self.kernel = [0.05, 0.05]
        self.prior = [NormalPrior(0, 1), NormalPrior(0, 1)]  # Assume Prior is correctly implemented
        self.random_seed = 42
    
    def test_initialization(self):
        """Test if InferenceController initializes with correct attributes."""
        ic = InferenceController(theta=self.theta, omega=self.omega, kernel=self.kernel,
                                 prior=self.prior, random_seed=self.random_seed)
        self.assertEqual(ic._n_params, 1)  # Expecting len(theta) - len(kernel)
    
    def test_setup_model(self):
        """Test if the model and data are setup correctly."""
        ic = InferenceController(theta=self.theta, omega=self.omega, kernel=self.kernel, prior=self.prior)
        _ = InferenceController(theta=self.theta, omega=self.omega, kernel=self.kernel)
        model, data = ic._setup_model(theta=self.theta, omega=self.omega, kernel=self.kernel, prior_list=self.prior)
        self.assertIsInstance(model, str)
        self.assertIsInstance(data, dict)
        self.assertIn("N", data)
        self.assertEqual(data["N"], len(self.theta))
    
    def test_run_method(self):
        """Test the run method executes and returns expected output types."""
        ic = InferenceController(theta=self.theta, omega=self.omega, kernel=[],
                                 prior=self.prior, random_seed=self.random_seed)
        _, samples = ic.run()
        _ = az.summary(samples)['mean'].values
        # npt.assert_array_almost_equal(res, [0.23, 0.26, np.inf, 0.24], decimal=2)
        ic_nonseed = InferenceController(theta=self.theta, omega=self.omega, kernel=self.kernel, prior=self.prior)
        _, samples = ic_nonseed.run()
    
    def test_predict_incidences_function(self):
        """Test the predict_incidences function returns predictions."""
        mean, std = predict_incidences(theta=self.theta, omega=self.omega,
                                       prior=self.prior, random_seed=self.random_seed)
        # npt.assert_array_almost_equal(mean, [0.23, 0.26, np.inf, 0.24], decimal=2)
        # npt.assert_array_almost_equal(std, [0.17, 0.26, np.inf, 0.27], decimal=2)
    
    def test_errors(self):
        error_theta = [0, 1]
        error_kernel = [0.1, 0.2, 0.3]
        with self.assertRaises(ValueError):
            InferenceController(error_theta, self.omega, error_kernel)
        with self.assertRaises(ValueError):
            error_prior = [0]
            InferenceController(self.theta, self.omega, self.kernel, error_prior)


if __name__ == '__main__':
    unittest.main()
