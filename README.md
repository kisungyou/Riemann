
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Learning on Riemannian Manifolds <a href='https://kyoustat.com/Riemann/'><img src='man/figures/logo.png' align="right" height="150" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/Riemann)](https://CRAN.R-project.org/package=Riemann)
[![Travis build
status](https://travis-ci.com/kyoustat/Riemann.svg?branch=master)](https://travis-ci.com/kyoustat/Riemann)
<!-- badges: end -->

**Riemann** is an R package for learning with data on **Riemannian
manifolds**. In statistics and machine learning, the term *manifold*
appears in two realms; one is *dimensionality reduction* where we assume
that low-dimensional data manifold is embedded in high-dimensional
Euclidean space. The other is *statistics on manifolds* - data lie on
some Riemannian manifolds that we are already well aware of. **Riemann**
aims at the latter. If you are interested in dimension reduction, please
check another R package [Rdimtools](https://kyoustat.com/Rdimtools/).

### Installation

  - Option 1 : **released** version from
    [CRAN](https://CRAN.R-project.org).

<!-- end list -->

``` r
install.packages("Riemann")
```

  - Option 2 : **development** version from
    [GitHub](https://github.com/).

<!-- end list -->

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("kyoustat/Riemann")
```
