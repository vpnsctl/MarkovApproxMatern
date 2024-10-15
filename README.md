# MarkovApproxMatern

Markov approximation of one-dimensional stationary Gaussian processes with Matérn covariance

## Shiny Application

To run the shiny app:

```{r}
# install.packages(c("shiny", "shinythemes", "shinyWidgets", "plotly"))
shiny::runApp("shiny_app/")
```

## Codes for the methods

The distance computations can be found on the files:
```{bash, eval=FALSE}
examples/compute_covariance_errors.R
python_codes/get_fourier_cov_errors.py
```

The prediction error computations were obtained using the following files:

```{bash, eval=FALSE}
# The following codes are used to generate samples and 
# obtain the predictions based on the exact covariance:
python_codes/get_true_pred.py
python_codes/gen_samples_and_true_pred.py

# The following codes compute the prediction errors for PCA, 
# Fourier and State-Space methods, respectively:
python_codes/get_pca_pred.py
python_codes/get_fourier_pred.py
python_codes/get_statespace_pred.py

# The following code computes the prediction errors for 
# Rational and nnGP (using the predictions and samples 
# from the previous files):
examples/compute_rational_nngp_errors.R
```

The posterior probability errors were obtained using the following file:

```{bash, eval=FALSE}
examples/probability_calculations_updated.R
```

The obtained quantities are stored in:

```{bash, eval=FALSE}
distance_tables/full_dists.RDS
pred_tables/pred_error.RDS
prob_tables/prob_errors.RDS
```

## Codes for calibrations

The calibrations were performed using the codes in the following files:

```{bash, eval=FALSE}
aux_functions/calibration_functions.R
aux_functions/create_calibration_fem.R
examples/example_calibration.R
examples/example_calibration_nngp.R
```

We calibrated the methods to have the same total runtime for assembling all matrices (the construction cost) and computing the posterior mean (the prediction cost). This ensures fair comparisons based on actual performance rather than theoretical complexity, which can be affected by the constants in the cost expression and the size of the study. Specifically, for a  given set of parameters \( (\kappa, \sigma, \nu) \), and a fixed value of $m$ for the proposed method, we calibrate the values of \( m \) for the other methods to ensure that the total computation time is the same. The total prediction times were averaged over 100 samples to obtain the calibration. 

Some of the methods have costs which depend on the smoothness parameter $\nu$. The calibration is therefore done separately for the ranges \( 0 < \nu < 0.5 \), \( 0.5 < \nu < 1.5 \), and \( 1.5 < \nu < 2.5 \). 

For \( \nu < 0.5 \), the calibration was not possible for nnGP and taper, because the rational approximation remained faster even with \( m = 1 \). For these cases, we set \( m = 1 \). In the Shiny application, users can either fix \( m = 1 \) for \( \nu < 0.5 \) or allow \( m \) to vary between 1 and 6, though the latter may result in an unfavorable comparison with the rational approximation.

For the taper method, the taper range was chosen so that each observation, on average, had \( m \) neighbors within the range, and the value of $m$ was then chosen to ensure that the total computational cost matches that of the rational approximation.

Given the number of basis functions, the PCA and Fourier methods have the same computational cost for prediction. Thus, the value of \( m \) for the Fourier method (the number of basis functions) was set to match the value obtained for the PCA method. The PCA method was calibrated disregarding the construction cost, which is equivalent to assuming that we know the eigenvectors of the covariance matrix. This is not realistic in practise, but makes the method act as a theoretical lower bound for any low-rank method.

As the state-space method provides an alternative Markov representation for which we could use the same computational methods as for the proposed method. Its value of \( m \) was therefore chosen as \( m - \lfloor \alpha \rfloor \).

To minimize boundary effects of the FEM method, we extended the original domain \( [0, 50] \) to \( [-4\rho, 50 + 4\rho] \), where \( \rho \) is the practical correlation range. This extension ensures accurate approximations of the Matérn covariance at the target locations. Because the FEM method uses the same type of rational approximation as the propose method, we fixed the value of $m$ for the FEM method to be equal to $m$ for the proposed method. The calibration was then instead performed on the finite element mesh. Specifically, a mesh with $N = kn + 500 - (k+3)$ locations in the extended domain used, where $500$ locations were in the extensions and $kn$ locations in the interior $[0,50]$, and the $-(k+3)$ term appears to ensure that the regular mesh contains the observation locations. These $kn$ locations were chosen equally spaced to include the observation locations and $k\in\mathbb{N}$ was calibrated to make the total computational cost match that of the proposed method for $\nu < 1.5$, and for  $1.5 < \nu < 2.5$ it was chosen as the largest values that yielded a stable prediction, as the value which matched the computational cost yielded unstable predictions.