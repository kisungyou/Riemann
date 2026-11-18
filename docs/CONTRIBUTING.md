# Contributing to Riemann

Kisung You maintains Riemann, reviews changes, and prepares releases.
Questions, examples that reveal unexpected results, corrections to
documentation, and code contributions are welcome. Support is provided
as time permits.

## Ask a question or report a problem

Search the [issue tracker](https://github.com/kisungyou/Riemann/issues)
first. If an existing discussion does not address your question, open an
issue. Include the smallest example that reproduces the behavior, the
result you expected, the result you obtained, and the output of
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). For a
GitHub installation, include the commit if known. You can use simulated
data to reproduce a problem with a private dataset.

For a numerical issue, specify the manifold, geometry, and algorithm.
Include the random seed, the termination message or diagnostics when
available, and a reference value or calculation if you have one. For
installation problems, include the installation command and the first
error in its output.

## Propose and test a change

Fork the repository, clone your fork, and create a branch from `master`.
Use R 4.4 or later and a working compiler for R packages. The
[installation
instructions](https://github.com/kisungyou/Riemann#installation) link to
tools for each operating system. Install the development dependencies
from an R session at the repository root:

``` r

install.packages("devtools")
devtools::install_deps(dependencies = TRUE)
devtools::test()
devtools::check(cran = TRUE, run_dont_test = TRUE, error_on = "warning")
```

The check builds the package and its vignettes and runs the extended
examples. Rendering the vignettes requires
[Pandoc](https://pandoc.org/installing.html), also supplied by RStudio.
For work on a single test file, for example
`tests/testthat/test-kmeans.R`, use `devtools::test(filter = "kmeans")`
during development. Include a regression case for a numerical fix. A
result calculated independently, a closed-form example, or an
appropriate geometric invariance can make the test informative across
platforms.

Edit function help in the roxygen comments in `R/` and regenerate it
with
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html).
For changes to the public method inventory, follow the [maintenance
instructions](https://github.com/kisungyou/Riemann/blob/master/tools/maintenance/README.md).
To rebuild the website, install `pkgdown`, install the current checkout
with
[`devtools::install()`](https://devtools.r-lib.org/reference/install.html),
and run `Rscript pkgdown/build.R` from its root.

Open a pull request that explains the user-visible change, links any
relevant issue, and states how you checked the result. If you propose
changing a geometric convention or a public interface, start an issue to
discuss the effect on existing analyses. The maintainer reviews
proposals and uses the package checks on Linux, macOS, and Windows to
inform merge decisions.
