
<!-- README.md is generated from README.Rmd. Please edit that file -->

# affinitymatrix

<!-- badges: start -->

<!-- badges: end -->

The goal of `affinitymatrix` is to provide a set of tools to estimate
matching models without frictions and with Transferable Utility starting
from matched data. The package contains a set of functions to implement
the tools developed by Dupuy and Galichon (2014), Dupuy, Galichon and
Sun (2019) and Ciscato, Galichon and Gousse (2020).

## Installation

<!-- You can install the released version of affinitymatrix from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("affinitymatrix") 
```-->

You can install the released version of affinitymatrix directly from
[Github](https://github.com/edoardociscato/affinitymatrix) with:

``` r
devtools::install_github("edoardociscato/affinitymatrix")
```

## Example

The example below shows how to use the main function of the package,
`estimate.affinity.matrix`, and how to summarize its findings. We first
generate a random sample of *matches*: the values of the normalized
variance-covariance matrix used in the data generating process are taken
from Chiappori, Ciscato and Guerriero (2020). For the sake of clarity,
consider the marriage market example: every row of our data frame
corresponds to a couple. The first `Kx` columns correspond to the
husbands’ characteristics, or *matching variables*, while the last `Ky`
columns correspond to the wives’. The key inputs to feed to
`estimate.affinity.matrix` are two matrices `X` and `Y` corresponding to
the husbands’ and wives’ traits and sorted so that the i-th man in `X`
is married to the i-th woman in `Y`.

``` r
# Load
library(affinitymatrix)

# Parameters
Kx = 4; Ky = 4; # number of matching variables on both sides of the market
N = 500 # sample size
mu = rep(0, Kx+Ky) # means of the data generating process
Sigma = matrix(c(1, 0.326, 0.1446, -0.0668, 0.5712, 0.4277, 0.1847, -0.2883, 0.326, 1, -0.0372, 0.0215, 0.2795, 0.8471, 0.1211, -0.0902, 0.1446, -0.0372, 1, -0.0244, 0.2186, 0.0636, 0.1489, -0.1301, -0.0668, 0.0215, -0.0244, 1, 0.0192, 0.0452, -0.0553, 0.2717, 0.5712, 0.2795, 0.2186, 0.0192, 1, 0.3309, 0.1324, -0.1896, 0.4277, 0.8471, 0.0636, 0.0452, 0.3309, 1, 0.0915, -0.1299, 0.1847, 0.1211, 0.1489, -0.0553, 0.1324, 0.0915, 1, -0.1959, -0.2883, -0.0902, -0.1301, 0.2717, -0.1896, -0.1299, -0.1959, 1),
               nrow=Kx+Ky) # (normalized) variance-covariance matrix of the data generating process
labels_x = c("Educ.", "Age", "Height", "BMI") # labels for men's matching variables
labels_y = c("Educ.", "Age", "Height", "BMI") # labels for women's matching variables

# Sample
data = MASS::mvrnorm(N, mu, Sigma) # generating sample
X = data[,1:Kx]; Y = data[,Kx+1:Ky] # men's and women's sample data
w = sort(runif(N-1)); w = c(w,1) - c(0,w) # sample weights

# Main estimation
res = estimate.affinity.matrix(X, Y, w = w, nB = 500)
#> Setup...
#> Main estimation...
#> Saliency analysis...
#> Rank tests...
#> Saliency analysis (bootstrap)...
```

The output of `estimate.affinity.matrix` is a list of objects that
constitute the estimation results. The estimated affinity matrix is
stored under `Aopt`, while the vector of loadings describing men’s and
women’s matching factors are stored under `U` and `V` respectively. The
following functions can be used to produce summaries of the different
findings. For further details, Chiappori, Ciscato and Guerriero (2020)
contain tables and plots that are generated with these functions.

``` r
# Print affinity matrix with standard errors
show.affinity.matrix(res, labels_x = labels_x, labels_y = labels_y)

# Print diagonal elements of the affinity matrix with standard errors
show.diagonal(res, labels = labels_x)

# Print rank test summary
show.test(res)

# Print results from saliency analysis
show.saliency(res, labels_x = labels_x, labels_y = labels_y, ncol_x = 2, ncol_y = 2)

# Show correlation between observed variables and matching factors
plots = show.correlations(res, labels_x = labels_x, labels_y = labels_y,
                          label_x_axis = "Husband", label_y_axis = "Wife")
```

## Literature

Ciscato, Edoardo, Alfred Galichon, and Marion Gousse. “Like attract
like? a structural comparison of homogamy across same-sex and
different-sex households.” *Journal of Political Economy* 128, no. 2
(2020): 740-781.

Chiappori, Pierre-André, Edoardo Ciscato, and Carla Guerriero.
“Analyzing matching patterns in marriage: theory and application to
Italian data.” *HCEO Working Paper* no. 2020-080 (2020).

Dupuy, Arnaud, and Alfred Galichon. “Personality traits and the marriage
market.” *Journal of Political Economy* 122, no. 6 (2014): 1271-1319.

Dupuy, Arnaud, Alfred Galichon, and Yifei Sun. “Estimating matching
affinity matrices under low-rank constraints.” *Information and
Inference: A Journal of the IMA* 8, no. 4 (2019): 677-689.
