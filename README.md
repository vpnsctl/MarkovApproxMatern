# MarkovApproxMatern

Markov approximation of one-dimensional stationary Gaussian processes with Matérn covariance.

## Shiny Application

To run the Shiny app:

```r
# install.packages(c("shiny", "shinythemes", "shinyWidgets", "plotly"))
shiny::runApp("shiny_app/")
```

## Codes for the Methods

The distance computations can be found in the following files:

```bash
examples/compute_covariance_errors.R
python_codes/get_fourier_cov_errors.py
```

The prediction error computations were obtained using these files:

```bash
# Generate samples and obtain predictions based on the exact covariance:
python_codes/get_true_pred.py
python_codes/gen_samples_and_true_pred.py

# Compute the prediction errors for PCA, Fourier, and State-Space methods:
python_codes/get_pca_pred.py
python_codes/get_fourier_pred.py
python_codes/get_statespace_pred.py

# Compute prediction errors for Rational and nnGP methods:
examples/compute_rational_nngp_errors.R
```

The posterior probability errors were obtained using:

```bash
examples/probability_calculations_updated.R
```

The obtained quantities are stored in:

```bash
distance_tables/full_dists.RDS
pred_tables/pred_error.RDS
prob_tables/prob_errors.RDS
```

## Codes for Calibrations

The calibrations were performed using these files:

```bash
aux_functions/calibration_functions.R
aux_functions/create_calibration_fem.R
examples/example_calibration.R
examples/example_calibration_nngp.R
```

We calibrated the methods to ensure equal total runtime for assembling matrices (construction cost) and computing the posterior mean (prediction cost). This ensures fair comparisons based on actual performance rather than theoretical complexity, which can be influenced by constants or the study size. For a given set of parameters $(\kappa, \sigma, \nu)$ and a fixed value of $m$ for the proposed method, we calibrated the values of $m$ for the other methods to ensure equal computation time. Prediction times were averaged over 100 samples to obtain the calibration.

Some methods have costs dependent on the smoothness parameter $\nu$. Therefore, calibration was performed separately for these ranges:

- $0 < \nu < 0.5$  
- $0.5 < \nu < 1.5$  
- $1.5 < \nu < 2.5$

For $\nu < 0.5$, calibration was not feasible for nnGP and taper methods, as the rational approximation remained faster even with $m = 1$. In such cases, we set $m = 1$. In the Shiny app, users can either fix $m = 1$ for $\nu < 0.5$ or allow $m$ to vary between 1 and 6, though this may result in unfavorable comparisons with the rational approximation.

For the taper method, the taper range was selected such that each observation, on average, had $m$ neighbors within the range. The value of $m$ was adjusted to match the total computational cost of the rational approximation.

For PCA and Fourier methods, the prediction costs are identical since both depend on the number of basis functions. Thus, the value of $m$ for the Fourier method was set to match the PCA method. PCA was calibrated without considering the construction cost, assuming the eigenvectors of the covariance matrix were known, providing a theoretical lower bound for any low-rank method.

The state-space method offers an alternative Markov representation, enabling the use of the same computational strategies as the proposed method. Its value of $m$ was chosen as $m - \lfloor \alpha \rfloor$.

To minimize boundary effects in the FEM method, the original domain $[0, 50]$ was extended to $[-4\rho, 50 + 4\rho]$, where $\rho$ is the practical correlation range. This extension ensures accurate approximation of the Matérn covariance at target locations. Since the FEM method uses the same rational approximation as the proposed method, the value of $m$ was fixed to match the proposed method. Calibration for the FEM method was performed on the finite element mesh.

Specifically, a mesh with  $N = kn + 500 - (k + 3)$ locations was used, where 500 locations were in the extensions and $kn$ in the interior $[0, 50]$. The term $-(k + 3)$ ensures that the regular mesh contains the observation locations. These $kn$ locations were evenly spaced to include the observations, and $k \in \mathbb{N}$ was calibrated to match the computational cost of the proposed method for $\nu < 1.5$. For $1.5 < \nu < 2.5$, $k$ was chosen as the largest value yielding stable predictions, as matching the computational cost resulted in instability.
