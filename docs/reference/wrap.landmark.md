# Wrap Landmark Data on Shape Space

One of the frameworks used in shape space is to represent the data as
landmarks. Each shape is a point set of \\k\\ points in \\\mathbf{R}^p\\
where each point is a labeled object. We consider general landmarks in
\\p=2,3,\ldots\\. Note that when \\p \> 2\\, it is stratified space but
we assume singularities do not exist or are omitted. The wrapper takes
translation and scaling out from the data to make it *preshape*
(centered, unit-norm). Also, for convenience, orthogonal Procrustes
analysis is applied with the first observation being the reference so
that all the other data are aligned to that reference. The full
orthogonal group O(p), including reflections, is used. This is not an
SO(p)-only orientation-preserving analysis. Singular configurations are
rejected. Fitted tangent models retain and apply their own training
reference.

## Usage

``` r
wrap.landmark(input, reference = NULL)
```

## Arguments

- input:

  data matrices to be wrapped as `riemdata` class. Following inputs are
  considered,

  array

  :   a \\(k\times p\times n)\\ array where each slice along 3rd
      dimension is a \\k\\-ad in \\\mathbf{R}^p\\.

  list

  :   a length-\\n\\ list whose elements are \\k\\-ads.

- reference:

  optional reference configuration for alignment. If NULL, the first
  observation is used.

## Value

a named `riemdata` S3 object containing

- data:

  a list of preshapes in \\\mathbf{R}^p\\.

- size:

  size of each preshape.

- name:

  name of the manifold of interests, *"landmark"*

## References

Dryden IL, Mardia KV (2016). *Statistical shape analysis with
applications in R*, Wiley series in probability and statistics, Second
edition edition. John Wiley and Sons, Chichester, UK ; Hoboken, NJ. ISBN
978-1-119-07251-5 978-1-119-07250-8.

## Examples

``` r
## USE 'GORILLA' DATA
data(gorilla)
riemobj = wrap.landmark(gorilla$male)
```
