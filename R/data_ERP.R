#' Data: xDAWN Covariances of MNE Sample MEG Epochs
#'
#' A collection of 216 symmetric positive-definite 32-by-32 matrices with
#' four event-label groups. The matrices represent test epochs from one
#' participant in the MNE sample audiovisual experiment, after supervised
#' xDAWN spatial filtering and prototype augmentation. The 32 coordinates
#' comprise 16 class-prototype features and 16 filtered trial features,
#' rather than 32 original EEG or MEG sensors.
#'
#' @details A September 21, 2026 reconstruction used the historical
#'   \pkg{pyRiemann} 0.2.7 example and numerical routines at commit
#'   \code{176e766f540bd4c7846f38573165fc3d27fc69ca}. All 216 stored matrices
#'   matched with maximum relative Frobenius error \eqn{5.56\times10^{-12}},
#'   and the complete label sequence matched exactly. No rescaling,
#'   reordering, or eigenvector-sign correction was needed. The reconstruction
#'   used recorded contemporary dependencies; it does not claim byte identity
#'   or recovery of the original 2021 runtime.
#'
#'   The recipe reads MNE's prefiltered 0--40 Hz sample recording, sampled
#'   at approximately 150.15 Hz after fourfold decimation, applies a 2 Hz
#'   IIR high-pass filter, excludes \code{MEG 2443}, and selects 305 MEG
#'   channels, with EEG excluded. Epochs cover 0--1 second without baseline
#'   correction or projection, giving 288 epochs with 151 samples each.
#'   A random 25 percent training split (seed 42, without stratification)
#'   supplies 72 epochs for supervised xDAWN fitting with four filters per
#'   class. The remaining 216 epochs supply the packaged matrices.
#'   Covariance estimation uses centered empirical covariance with denominator
#'   151, on the combined prototype and filtered-trial features.
#'
#'   These are repeated epochs from one participant, not independent
#'   participants. The representation was learned with labels from the
#'   separate 72-epoch training set, so it must not be described as an
#'   unsupervised transformation of the original sensors. A later split of
#'   the packaged matrices does not recreate raw-data preprocessing or
#'   establish performance for new participants.
#'
#'   The matrices retain their original numerical scale; small absolute
#'   eigenvalues alone do not imply poor conditioning. No unit conversion,
#'   regularization, data values, or labels were changed in this revision.
#'   The MNE sample data's public OpenNeuro record, \code{ds000248} version
#'   1.2.4, declares CC0. See the source record for attribution and terms.
#'
#' @source \href{https://github.com/pyRiemann/pyRiemann/blob/176e766f540bd4c7846f38573165fc3d27fc69ca/examples/ERP/plot_embedding_EEG.py}{Historical pyRiemann ERP recipe},
#'   by Pedro Rodrigues (the historical filename says EEG, but the recipe
#'   explicitly selects MEG and excludes EEG).
#'
#'   \href{https://doi.org/10.18112/openneuro.ds000248.v1.2.4}{MNE sample data on OpenNeuro, version 1.2.4};
#'   \href{https://raw.githubusercontent.com/OpenNeuroDatasets/ds000248/0a357f1276dbf7fcc3c3fc35f74cb24ffe4a68a8/dataset_description.json}{pinned source metadata and CC0 declaration}.
#'   The \href{https://mne.tools/stable/documentation/datasets.html#sample}{MNE sample description}
#'   identifies the audiovisual experiment. The separate writing workspace
#'   records the reconstruction in
#'   \code{maintenance/reports/erp-provenance-followup-2026-09-21.md}.
#'
#' @usage data(ERP)
#' 
#' @examples
#' \donttest{
#' ## LOAD THE DATA AND WRAP AS RIEMOBJ
#' data(ERP)
#' myriem = wrap.spd(ERP$covariance)
#' }
#' 
#' @format a named list containing\describe{
#' \item{covariance}{an \eqn{(32\times 32\times 216)} array of covariance matrices.}
#' \item{label}{a length-\eqn{216} factor with levels LA, LV, RA, RV and counts 57, 57, 53, 49.}
#' }
#' 
#' @references
#' Gramfort A et al. (2013). "MEG and EEG data analysis with MNE-Python."
#' \emph{Frontiers in Neuroscience}, 7, 267.
#' \doi{10.3389/fnins.2013.00267}.
#'
#' Gramfort A et al. (2014). "MNE software for processing MEG and EEG data."
#' \emph{NeuroImage}, 86, 446--460.
#' \doi{10.1016/j.neuroimage.2013.10.027}.
#'
#' @seealso \code{\link{wrap.spd}}
#' @concept data
"ERP"

# data(ERP)
# wrapERP = wrap.spd(ERP$covariance)
# embed2d = riem.mds(wrapERP)
# https://pyriemann.readthedocs.io/en/latest/auto_examples/ERP/plot_embedding_EEG.html#sphx-glr-auto-examples-erp-plot-embedding-eeg-py
# https://mne.tools/stable/generated/mne.datasets.sample.data_path.html#mne.datasets.sample.data_path
