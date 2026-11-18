# Packaged data provenance

This record documents the origins and interpretation of Riemann's packaged
data. The ERP numerical reconstruction and source-license inspection were
performed on **September 21, 2026**; relocation of the evidence into this
package does not constitute a new reconstruction. Dataset values are unchanged.
Paths beginning with `data/` refer to the source package. Replay helpers are
available in the repository at `tools/maintenance/provenance/`.

## Gorilla configurations

The male and female arrays derive from `gorm.dat` and `gorf.dat` in
[`shapes` 1.2.8](https://cran.r-project.org/package=shapes). Reorder the
landmark dimension by **`c(1, 5, 4, 3, 2, 8, 7, 6)`** to obtain exact equality
with `gorilla$male` (8 x 2 x 29) and `gorilla$female` (8 x 2 x 30), respectively.
No coordinate transformation, rescaling, rounding, or individual reordering is
needed. The package file `data/gorilla.rda` has MD5
`62a4f252f74c05a15bff33e860b795bd`.

The upstream documentation attributes the data to Paul O'Higgins and cites
O'Higgins and Dryden (1993), *Sexual dimorphism in hominoids: further studies
of craniofacial shape differences in Pan, Gorilla, Pongo*, Journal of Human
Evolution 24, 183–205, and Dryden and Mardia (2016), Chapter 1. The inspected
`shapes_1.2.8.tar.gz` archive had 369,109 bytes and SHA256
`2b1df1abbe8f4836e20bb23844fabe8acb3f7ca09b9fc9fbe9cd55d972463333`.
It declares GPL-2, contains these coordinates as data and R-source literals,
and had no observed data-specific exception. Preserve the upstream attribution,
GPL-2 terms and landmark permutation; the package's author-owned GPL-3 software
notice does not relicense these separate materials.

Run `Rscript tools/maintenance/provenance/verify-gorilla-provenance.R` from a
source checkout with `shapes` installed. The script reports its actual upstream
version, checks both arrays with `identical()`, and prints the package-data MD5.

## ERP covariance matrices

The supplied `data/ERP.rda` is byte-identical to its addition in Riemann commit
`1fe5f4695ef5e60b714be0a34848295f7525e6dc` on November 2, 2021. Its SHA256 is
`787bd1bf3e593dfb823440f23abbfe427c41f1b57f6afe79f8f8f938e9e9ae49`.
All 216 matrices and their complete label sequence numerically correspond to
the [historical pyRiemann v0.2.7 example](https://github.com/pyRiemann/pyRiemann/blob/176e766f540bd4c7846f38573165fc3d27fc69ca/examples/ERP/plot_embedding_EEG.py).
That example's SHA256 is
`73353b2764c5f28bf6a2580a257a4914e869508baeefba7870483d1d7c582e32`.

The historical replay's maximum unadjusted relative Frobenius discrepancy was
`5.556594353951526e-12` (median `9.508958511327822e-13`), below its `1e-10`
acceptance threshold. The maximum absolute entry discrepancy was
`9.049928968889244e-37`; this absolute quantity is not the acceptance criterion.
All labels matched in their stored order. No scaling, observation reordering,
ridge, or sign correction was used for acceptance. This establishes numerical
correspondence, not byte identity or recovery of the unknown 2021 environment.
The original exploratory run had no prespecified tolerance; the portable
follow-up explicitly tested `1e-10`.

### Construction and meaning of an observation

The example reads the MNE sample MEG recording after historical 40-Hz low-pass
filtering and decimation by four, then applies a 2-Hz IIR high-pass filter.
Despite its historical EEG filename, it selects MEG and excludes EEG. After
excluding `MEG 2443`, it constructs 288 epochs with 305 channels and 151 samples
over 0–1 seconds, without baseline correction or projection.

`train_test_split(train_size=0.25, random_state=42)` assigns 72 epochs to training
and 216 to testing. Four supervised Xdawn filters per class yield 16 filtered
trial coordinates; 16 fitted class-prototype coordinates are prepended before
forming centered empirical covariances with denominator 151. Thus each 32 x 32
matrix describes augmented trial features from **one sample participant**, not
32 EEG sensors or a separate participant. Stored label counts are LA = 57,
RA = 53, LV = 57 and RV = 49. The upper-left 16 x 16 prototype block is identical
across all matrices. Shared prototypes, repeated measurements and potentially
overlapping one-second windows preclude assuming independent observations from
the array count. Since Xdawn uses training labels, subsequent separation is not
an entirely unsupervised discovery from raw recordings.

### Source terms and bounded recording linkage

The [OpenNeuro ds000248 version 1.2.4 release](https://openneuro.org/datasets/ds000248/versions/1.2.4),
DOI `10.18112/openneuro.ds000248.v1.2.4`, identifies the MNE sample data and
declares CC0. It names Alexandre Gramfort and Matti S. Hämäläinen. The retained
DataCite metadata records publication year 2020. Pinned metadata is preserved
under [erp/dataset-metadata/](erp/dataset-metadata/); these dataset terms are
separate from the package's GPL-3 and pyRiemann's BSD-3-Clause code licenses.
The November 2010 MNE-C manual (printed page 307) restricted redistribution of
its then-supplied sample data. The present provenance relies on the affirmative
current CC0 release and the documented derivative correspondence; it does not
certify the terms of every historical distribution.

The raw OpenNeuro recording and processed input have different checksums.
The [recording-link result](erp/recording-link.json) establishes agreement of
376 ordered channels, kinds and units, factor-of-four sampling metadata, and
all 41,700 stimulus-channel blocks. All 319 supplied processed events are
recreated and each maps uniquely to a same-code raw event. One of 320 raw rising
events is omitted by the increasing-event extraction after aggregation; a
320-to-319 bijection is not claimed. The check does not compare filtered MEG
waveforms or replace historical MNE-C filtering with a modern approximation.

Please retain the MNE acknowledgements: Gramfort et al. (2014), *MNE software
for processing MEG and EEG data*, NeuroImage 86, 446–460,
[doi:10.1016/j.neuroimage.2013.10.027](https://doi.org/10.1016/j.neuroimage.2013.10.027),
and Gramfort et al. (2013), *MEG and EEG data analysis with MNE-Python*,
Frontiers in Neuroscience 7, 267,
[doi:10.3389/fnins.2013.00267](https://doi.org/10.3389/fnins.2013.00267).
The retained upstream metadata has a different NeuroImage DOI; its original
contents have not been silently corrected.

### Evidence and replay

- [comparison.json](erp/comparison.json): historical numerical results, input
  and source hashes, full class counts and numerical-library versions.
- [recording-link.json](erp/recording-link.json): bounded recording checks.
- [data-downloads.json](erp/data-downloads.json): recording URLs and checksums;
  the large recordings are not distributed here.
- [mne-generation-command.json](erp/mne-generation-command.json): pinned
  historical filter/decimation commands, recorded rather than re-executed.
- [R-session-info.txt](erp/R-session-info.txt): historical export environment.

The numerical replay used Python 3.12.14, MNE 1.8.0, NumPy 2.2.6, SciPy 1.15.3,
scikit-learn 1.6.1, modern pyRiemann 0.9 for an additional comparison, and R 4.5.2
for exporting the supplied array. The repository's provenance helpers retain
the dependency pins, unchanged upstream numerical sources and BSD license.
They check source hashes before extracting the historical classes/functions,
avoiding unrelated legacy imports. Historical JSON retains its original
temporary input paths and script hashes as evidence, not current dependencies.
The relocated helpers differ only in package-license wording and path setup;
new executions record their own script hashes and output paths. Large generated
arrays can be recreated and are not retained as package assets.
