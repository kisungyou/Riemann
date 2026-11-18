# Package maintenance

These developer tools are self-contained within the Riemann checkout and are
excluded from the R source package. Run the following commands from its root.

## Method contracts

```sh
Rscript tools/maintenance/generate-inventory.R
```

This statically resolves public functions and their direct test references,
then uses the explicit dispositions in `generate-method-contracts.R` to update
`inst/method-contracts.csv`. It evaluates no package functions and produces no
separate inventory files. A test reference is not evidence that a test passed.
Review the resulting diff after API, source, test, or contract changes. Repeated
runs on unchanged source produce identical output. An optional argument selects
another Riemann checkout; without it, paths resolve relative to this script,
independently of the working directory.

## Installed-package tests

First install the source being tested into a dedicated library, then run:

```sh
Rscript tools/maintenance/verify-package.R . /path/to/library /path/to/results
```

The verifier requires the explicitly selected installation to match the source
version, runs the complete testthat suite against that installation, and writes
session information and per-test results. It fails on any failure, error,
warning, or skip. Matching version numbers do not establish matching source:
install the exact checkout before each release validation and record its commit.

## Rendered website

After running the website build documented in `pkgdown/build.R`:

```sh
python3 tools/maintenance/audit-website.py
```

The read-only verifier checks local links and anchors, search and sitemap
targets, current version labels, reference topics, mathematical rendering, and
the documented limitations of stationary summary estimates. It does not test
external services. An optional path selects another website output directory;
use `--version` when its parent does not contain the package `DESCRIPTION`.

Data-provenance replay tools live in `provenance/`; their retained sources and
verification scope are documented separately.
