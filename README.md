# MarkovApproxMatern

Markov approximation of one-dimensional stationary Gaussian processes with Matérn covariance

To run the shiny app:

```{r}
# install.packages(c("shiny", "shinythemes", "shinyWidgets", "plotly"))
shiny::runApp("shiny_app/")
```

The distance computations can be found on the files:
```{bash, eval=FALSE}
examples/compute_covariance_errors.R
python_codes/get_fourier_cov_errors.py
```

The prediction error computations were obtained using the following files:

```{bash, eval=FALSE}
# These codes are used to generate samples and obtain the predictions based on the exact covariance
python_codes/get_true_pred.py
python_codes/gen_samples_and_true_pred.py
# These codes compute the prediction errors for PCA, Fourier and State-Space methods, respectively
python_codes/get_pca_pred.py
python_codes/get_fourier_pred.py
python_codes/get_statespace_pred.py
# The following code computes the prediction errors for Rational and nnGP (using the predictions and samples from the previous files):
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