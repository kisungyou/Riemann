# Riemann software paper

The manuscript describes Riemann 0.2.0. The [Markdown source](paper.md) and
[bibliography](paper.bib) are the journal inputs. A [rendered draft](paper.pdf)
is available for reading. Open Journals supplies the draft watermark and
placeholder publication metadata.

## Reviewing the software

| Review item | Where to check it |
| --- | --- |
| Research need, comparison with related tools, design choices and documented research use | [Manuscript](paper.md) |
| Authorship, license and attribution | [Package metadata](../DESCRIPTION), [software license](../COPYING), [third-party notices](../inst/COPYRIGHTS) and [data provenance](../inst/provenance/README.md) |
| Installation, dependencies and first example | [Installation instructions and runnable example](../README.md) |
| Scientific workflows and limits of supported methods | [Geometry tutorial](../vignettes/geometry-workflows.Rmd), [data tutorial](../vignettes/data-workflows.Rmd) and [method contracts](../inst/method-contracts.csv) |
| Correctness and portability | [Mathematical regression tests](../tests/testthat/), [package check workflow](../.github/workflows/R-CMD-check.yaml) and [check results](https://github.com/kisungyou/Riemann/actions/workflows/R-CMD-check.yaml) |
| Development history, contributions and support | [Release notes](../NEWS.md), [commit history](https://github.com/kisungyou/Riemann/commits/master/), [public issues](https://github.com/kisungyou/Riemann/issues) and [contribution guide](../CONTRIBUTING.md) |

For a reproducible check, select a commit and follow the contribution guide
from that checkout. The [paper workflow](../.github/workflows/draft-pdf.yml)
uses the Open Journals renderer with `journal: joss` and
`paper-path: paper/paper.md`. Its **riemann-joss-draft** artifact contains
the PDF for the commit that triggered the build.

## Submission location

Use the [public software repository](https://github.com/kisungyou/Riemann),
branch `master`, manuscript path `paper/paper.md`. This directory keeps the
paper and its references under version control alongside the package, as
specified in the [JOSS submission instructions](https://joss.readthedocs.io/en/latest/submitting.html).
After successful review, tag the reviewed version and deposit it in a
preservation archive for its software DOI.
