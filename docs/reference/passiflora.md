# Data : Passiflora Leaves

Passiflora is a genus of about 550 species of flowering plants. This
dataset contains 15 landmarks in 2 dimension of 3319 leaves of 40
species. Papers listed in the reference section analyzed the data and
found 7 clusters.

## Usage

``` r
data(passiflora)
```

## Format

a named list containing

- data:

  a 3d array of size \\(15\times 2\times 3319)\\.

- species:

  a length-\\3319\\ vector of 40 species factors.

- class:

  a length-\\3319\\ vector of 7 cluster factors.

## References

Chitwood DH, Otoni WC (2017). “Divergent leaf shapes among Passiflora
species arise from a shared juvenile morphology.” *Plant Direct*,
**1**(5), e00028. ISSN 24754455.

Chitwood DH, Otoni WC (2017). “Morphometric analysis of Passiflora
leaves: the relationship between landmarks of the vasculature and
elliptical Fourier descriptors of the blade.” *GigaScience*, **6**(1).
ISSN 2047-217X.

## See also

[`wrap.landmark`](https://www.kisungyou.com/Riemann/reference/wrap.landmark.md)

## Examples

``` r
data(passiflora)                         # load the data
riemobj = wrap.landmark(passiflora$data) # wrap as RIEMOBJ
pga2d   = riem.pga(riemobj)$embed        # embedding via PGA

opar <- par(no.readonly=TRUE)            # visualize
plot(pga2d, col=passiflora$class, pch=19, cex=0.7,
     main="PGA Embedding of Passiflora Leaves",
     xlab="dimension 1", ylab="dimension 2")

par(opar)
```
