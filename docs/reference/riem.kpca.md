# Kernel Principal Component Analysis

Although the method of Kernel Principal Component Analysis (KPCA) was
originally developed to visualize non-linearly distributed data in
Euclidean space, we graft this to the case for manifolds where extrinsic
geometry is explicitly available. The algorithm uses Gaussian kernel
with \$\$K(X_i, X_j) = \exp\left( - \frac{d^2 (X_i, X_j)}{2 \sigma^2}
\right )\$\$ where \\\sigma\\ is a bandwidth parameter and \\d(\cdot,
\cdot)\\ is an extrinsic distance defined on a specific manifold.

## Usage

``` r
riem.kpca(riemobj, ndim = 2, sigma = 1, negative = c("truncate", "error"))
```

## Arguments

- riemobj:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- ndim:

  an integer-valued target dimension (default: 2).

- sigma:

  A finite, strictly positive Gaussian bandwidth.

- negative:

  Treatment of an empirically indefinite kernel: `"truncate"` warns and
  returns a positive-spectrum embedding; `"error"` stops.

## Value

a named list containing

- embed:

  an \\(N\times ndim)\\ matrix whose rows are embedded observations.

- vars:

  a length-\\N\\ vector of eigenvalues from kernelized covariance
  matrix.

## Details

Scores are eigenvectors multiplied by square roots of the centered Gram
eigenvalues. They are ordered by decreasing eigenvalue. A Gaussian of an
arbitrary manifold dissimilarity is not automatically a
positive-definite kernel. For an indefinite Gram matrix, truncated
coordinates are an approximation, not exact kernel principal component
scores. The full spectrum, negative eigenvalues, and empirical kernel
status are returned. This function remains a training-data embedding,
without a predict method.

## References

Schölkopf B, Smola A, Müller K (1997). “Kernel principal component
analysis.” In Goos G, Hartmanis J, van Leeuwen J, Gerstner W, Germond A,
Hasler M, Nicoud J (eds.), *Artificial Neural Networks — ICANN'97*,
volume 1327, 583–588. Springer Berlin Heidelberg, Berlin, Heidelberg.
ISBN 978-3-540-63631-1 978-3-540-69620-9.

## Examples

``` r
#-------------------------------------------------------------------
#          Example for Gorilla Skull Data : 'gorilla'
#-------------------------------------------------------------------
## PREPARE THE DATA
#  Aggregate two classes into one set
data(gorilla)

mygorilla = array(0,c(8,2,59))
for (i in 1:29){
  mygorilla[,,i] = gorilla$male[,,i]
}
for (i in 30:59){
  mygorilla[,,i] = gorilla$female[,,i-29]
}

gor.riem = wrap.landmark(mygorilla)
gor.labs = c(rep("red",29), rep("blue",30))

## APPLY KPCA WITH DIFFERENT KERNEL BANDWIDTHS
kpca1 = riem.kpca(gor.riem, sigma=0.01)
kpca2 = riem.kpca(gor.riem, sigma=1)
#> Warning: The Gaussian dissimilarity kernel is indefinite; returning an explicitly truncated positive-spectrum embedding.
kpca3 = riem.kpca(gor.riem, sigma=100)
#> Warning: The Gaussian dissimilarity kernel is indefinite; returning an explicitly truncated positive-spectrum embedding.
## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3), pty="s")
plot(kpca1$embed, pch=19, col=gor.labs, main="sigma=1/100")
plot(kpca2$embed, pch=19, col=gor.labs, main="sigma=1")
plot(kpca3$embed, pch=19, col=gor.labs, main="sigma=100")

par(opar)
```
