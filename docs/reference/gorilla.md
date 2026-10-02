# Data : Gorilla Skull

This is 29 male and 30 female gorillas' skull landmark data where each
individual is represented as 8-ad/landmarks in 2 dimensions. This is a
re-arranged version of the data from shapes package.

## Usage

``` r
data(gorilla)
```

## Format

a named list containing

- male:

  a 3d array of size \\(8\times 2\times 29)\\

- female:

  a 3d array of size \\(8\times 2\times 30)\\

## Source

The arrays match
[`shapes::gorm.dat`](https://rdrr.io/pkg/shapes/man/gorm.dat.html) and
[`shapes::gorf.dat`](https://rdrr.io/pkg/shapes/man/gorf.dat.html) after
reordering landmarks by `c(1, 5, 4, 3, 2, 8, 7, 6)`. This exact
correspondence was verified against shapes 1.2.8; coordinates and
individual ordering are otherwise unchanged. The upstream manual
attributes the data to Paul O'Higgins and cites O'Higgins and Dryden
(1993). <https://cran.r-project.org/package=shapes>

## Details

These are raw landmark configurations. `wrap.landmark` removes
translation and scale and uses orthogonal alignment, including
reflections. The resulting shape analysis is different from analysis of
the original sizes or an orientation-preserving quotient. Preserve the
documented landmark order in comparisons with shapes.

## References

Dryden IL, Mardia KV (2016). *Statistical shape analysis with
applications in R*, Wiley series in probability and statistics, Second
edition edition. John Wiley and Sons, Chichester, UK ; Hoboken, NJ. ISBN
978-1-119-07251-5 978-1-119-07250-8.

O'Higgins P, Dryden IL (1993). "Sexual dimorphism in hominoids: further
studies of craniofacial shape differences in Pan, Gorilla, Pongo."
*Journal of Human Evolution*, 24:183–205.

## See also

[`wrap.landmark`](https://www.kisungyou.com/Riemann/reference/wrap.landmark.md)

## Examples

``` r
data(gorilla)                               # load the data
riem.female = wrap.landmark(gorilla$female) # wrap as RIEMOBJ
opar <- par(no.readonly=TRUE)
for (i in 1:30){
  if (i < 2){
    plot(riem.female$data[[i]], cex=0.5, 
         xlim=c(-1,1)/2, ylim=c(-2,2)/5,
         main="30 female gorilla skull preshapes",
         xlab="dimension 1", ylab="dimension 2")
    lines(riem.female$data[[i]])
  } else {
    points(riem.female$data[[i]], cex=0.5)
    lines(riem.female$data[[i]])
  }
}

par(opar)
```
