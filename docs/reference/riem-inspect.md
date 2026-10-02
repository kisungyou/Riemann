# Inspect Wrapped Observations and Geometric Summaries

Standard methods report the observation space, selected geometry, and
actual numerical termination. A finite estimate is not by itself proof
of convergence. Summary plots display matrix entries or ambient
coordinates; they do not imply an isometric view of the manifold.

## Usage

``` r
# S3 method for class 'riemdata'
print(x, ...)

# S3 method for class 'riemdata'
summary(object, ...)

# S3 method for class 'riem_summary'
print(x, ...)

# S3 method for class 'riem_summary'
summary(object, ...)

# S3 method for class 'riem_summary'
plot(x, ...)
```

## Arguments

- x, object:

  A wrapped dataset or fitted geometric summary.

- ...:

  Additional arguments passed to the plotting function where relevant.

## Value

Print methods return their object invisibly. Summary methods return
structured metadata and numerical diagnostics. Plot methods return
invisibly.
