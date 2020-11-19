
<!-- README.md is generated from README.Rmd. Please edit that file -->

# affinitymatrix

<!-- badges: start -->

<!-- badges: end -->

The goal of affinitymatrix is to provide a set of tools to estimate
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
generate a random sample of *matches*. For the sake of clarity, consider
the marriage market example: every row of our data frame corresponds to
a couple. The first `Kx` columns correspond to the husbands’
characteristics, or *matching variables*, while the last `Ky` columns
correspond to the wives’. The key inputs to feed to
`estimate.affinity.matrix` are two matrices `X` and `Y` corresponding to
the husbands’ and wives’ traits and sorted so that the i-th man in `X`
is married to the i-th woman in `Y`.

## Literature

Ciscato, Edoardo, Alfred Galichon, and Marion Gousse. “Like attract
like? a structural comparison of homogamy across same-sex and
different-sex households.” *Journal of Political Economy* 128, no. 2
(2020): 740-781.

Chiappori, Pierre-André, Edoardo Ciscato, and Carla Guerriero}.
Analyzing matching patterns in marriage: theory and application to
Italian data. *HCEO Working Paper* no. 2020-080 (2020).

Dupuy, Arnaud, and Alfred Galichon. “Personality traits and the marriage
market.” *Journal of Political Economy* 122, no. 6 (2014): 1271-1319.

Dupuy, Arnaud, Alfred Galichon, and Yifei Sun. “Estimating matching
affinity matrices under low-rank constraints.” *Information and
Inference: A Journal of the IMA* 8, no. 4 (2019): 677-689.
