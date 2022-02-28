
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Riemann <img src='man/figures/logo.png' alt="" align="right" height="150" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/Riemann)](https://CRAN.R-project.org/package=Riemann)
[![R-CMD-check](https://github.com/kisungyou/Riemann/workflows/R-CMD-check/badge.svg)](https://github.com/kisungyou/Riemann/actions)
<!-- badges: end -->

**Riemann** is an R package for learning with data on **Riemannian
manifolds**. In statistics and machine learning, the term *manifold*
appears in two realms; one is *dimensionality reduction* where we assume
that low-dimensional data manifold is embedded in high-dimensional
Euclidean space. The other is *statistics on manifolds* - data lie on
some Riemannian manifolds that we are already well aware of. **Riemann**
aims to achieve the latter. If you are interested in dimension
reduction, please check another R package
[Rdimtools](https://kisungyou.com/Rdimtools/).

### Installation

-   Option 1 : **released** version from
    [CRAN](https://CRAN.R-project.org).

``` r
install.packages("Riemann")
```

-   Option 2 : **development** version from
    [GitHub](https://github.com/).

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("kisungyou/Riemann")
```
