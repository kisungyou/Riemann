---
title: 'Riemann: Statistical learning with manifold-valued data in R'
tags:
  - R
  - statistical learning
  - Riemannian geometry
  - directional statistics
  - covariance matrices
  - shape analysis
authors:
  - name: Kisung You
    orcid: 0000-0002-8584-459X
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Mathematics, Baruch College, City University of New York, United States
    index: 1
date: 18 November 2026
bibliography: paper.bib
software_repository_url: https://github.com/kisungyou/Riemann
---

# Summary

Two animals on opposite sides of the international date line can be close
neighbors even when subtracting their recorded longitudes gives nearly 360
degrees. The question of how to measure distance also arises when an
observation is a covariance matrix or a set of anatomical landmarks. A geometry
specifies which differences between observations matter. `Riemann` brings
established learning algorithms for these data into R. Researchers can find
cluster centers, extract principal components, and predict a numerical
response while retaining the geometry used to fit a model. The package also
provides descriptive and inferential methods and probability models for
directional data. Applications include summarizing animal tracks and studying
variation in shape.

# Statement of need

Consider a collection of sensor recordings, each summarized by a
positive-definite covariance matrix. To cluster the recordings, an analyst must
decide how to measure the distance between two matrices and how to average the
matrices assigned to a cluster. For example, the log-Euclidean distance measures
differences after taking the matrix logarithm [@arsigny2006]. Under this
geometry, a cluster center is obtained by averaging the logarithms and taking
the matrix exponential. Averaging the original matrix entries solves a
different problem. Adapting k-means to these data therefore requires choosing
both a distance and a compatible update for the center.

The intended use of a new recording also matters. A tree built
by average linkage describes the recordings already in the sample. Its branches
can reveal smaller groups within a larger cluster. A fitted k-means model
answers a different question about a new recording by assigning it to the
nearest trained center. For a scalar response, kernel regression uses distances
to the training observations to average their responses, with larger weights
for nearer observations. An R user comparing these approaches needs to be able
to change the statistical task while keeping the definition of distance fixed.

# State of the field

In R [@rbase], `shapes` offers
Procrustes analysis and tools for studying variation in landmark
configurations, including permutation tests and graphical comparisons
[@shapes]. The `manifold` package supplies
geometric operations, including exponential and logarithm maps, for users who
want to construct their own estimators [@manifold]. In Python, `pyRiemann`
puts covariance-based estimation and classification into scikit-learn
pipelines [@pyriemann]. `Geomstats` offers a wider collection of manifolds
and metrics, together with learning algorithms expressed in terms of those
geometries [@geomstats].

The case for a separate R package rests on what these learning algorithms
share. A researcher may need the same scoring and prediction conventions
for directions, covariance matrices, and landmark shapes, alongside the usual
R tools for resampling and plotting. Adding a single estimator to a package
for one data type would address only one of those cases. In `Riemann`, the
shared work includes checking new observations against the training
representation, computing distances under the fitted geometry, and carrying
the reference into prediction. Existing R routines can still do the parts of
an analysis that depend only on dissimilarities. For instance, hierarchical
clustering passes the selected manifold distances to the standard R linkage
routine and returns a tree that can be cut or plotted in the usual way.

# Software design

Building on earlier CRAN releases [@riemann], [version 0.2.0](https://github.com/kisungyou/Riemann)
records the choice of geometry when a core model is fitted. A collection of symmetric positive-definite matrices can
be analyzed under either the affine-invariant or the log-Euclidean metric.
The affine-invariant distance is unchanged by a common nonsingular linear
transformation of the underlying measurements. The log-Euclidean choice makes
ordinary Euclidean calculations on matrix logarithms available. Both accept
the same data, so the matrices and their dimensions alone cannot identify how
a saved model was fitted. Recording the metric lets the prediction method
apply it to subsequent observations.
R handles model objects and their `predict` methods, with numerical operations
implemented in R, native C++, and selected R-package dependencies.

Tangent principal component analysis illustrates how much information a fit
must retain [@fletcher2004]. To score a new observation, the procedure needs
the training reference, the coordinates in which variation was measured, and
the fitted component directions. Even the choice of coordinates for a
symmetric matrix has a statistical consequence. An off-diagonal entry occurs
twice in the full matrix but once in its upper triangle. The package scales
that entry by the square root of two when storing it once, preserving its
contribution to the squared Frobenius norm. The fitted transformation can then
score one new matrix or a batch on the same axes. Mapping a score back to the
manifold shows the part of an observation represented by the retained
components.

For landmarks, the reference determines how a new specimen is aligned. In a
[worked example](https://www.kisungyou.com/Riemann/articles/data-workflows.html),
each gorilla skull is represented by corresponding landmarks
from a morphometric study [@ohiggins1993; @shapes]. The analysis removes
location and scale and permits rotations and reflections, so a specimen and its
mirror image receive the same shape scores. That identification is appropriate
only when handedness can be discarded. The tutorial uses the trained frame
when scoring a single specimen and reconstructs landmark positions along a
component so that the reader can see what a change in its score means.

Different algorithms require different geometric operations. Average linkage
needs pairwise distances, whereas a tangent analysis also needs coordinates
at a reference point. On the sphere, for example, distance to the antipode is
well defined, but there are many equally short paths to it. The round-sphere
tangent routine reports that ambiguity. The package documents supported
operations for each geometry and identifies experimental routines separately.
Tests include agreement with ordinary principal component analysis in
Euclidean space and recovery of the same scores after a fitted model is saved
and read back into R. The [geometry tutorial](https://www.kisungyou.com/Riemann/articles/geometry-workflows.html)
shows how to save a fit and apply it to a new covariance matrix.

# Research impact statement

@gilmour2025 used `Riemann` 0.1.4 to study habitat use by animals tracked around
a marine protected area in the central Pacific. The tags recorded positions at
different intervals, making individual tracks difficult to compare directly.
After processing the tracks, the researchers used the package to estimate a
geodesic median location for each animal on each day. Those daily locations
were inputs to kernel density estimates of habitat use and to some of the
species distribution models.

The package also distributes estimators developed in research on directional
distributions. @youSuh2022 describe likelihood estimation and mixture
clustering for the spherical normal distribution and identify `Riemann` as
their implementation. @youShung2023 developed related methods for the
spherical Laplace distribution. The two models give different roles to
distant observations. For a fixed concentration or scale, the location
criterion involves squared geodesic distance in the normal model and
geodesic distance itself in the Laplace model. Providing both in R allows a
researcher to fit either model to the same sample of directions.

In a comparison of software for brain connectivity analysis, @escobar
reported numerical problems with the affine-invariant mean in `Riemann`
0.1.6. This calculation matters beyond descriptive statistics because the
mean also supplies the reference for tangent principal components. Version
0.2.0 uses a common mean solver for these tasks and reports the objective
value and convergence status. Its tests compare the estimate for three
noncommuting covariance matrices with an independent optimization in base R,
checking both the matrix itself and the value of the objective.

# AI usage disclosure

OpenAI Codex assisted with numerical code revisions, regression tests, and
package documentation during development of version 0.2.0. Exact model
versions were not recorded for all earlier development sessions. Codex, using
GPT-6.2 Astra for this manuscript, also assisted with drafting, checking
citations against their sources, and preparing the PDF. Checks on the
software included comparisons with known mathematical results and execution
of package tests and worked examples.

<!-- Before submission, the author must verify the complete tool and model
history, review, edit and validate the assisted work, and confirm responsibility
for the core design decisions in the visible disclosure. -->

# Acknowledgements

This work received no specific funding. The author declares no competing interests.

# References
