# Data provenance replay helpers

The package's data-source record and historical results are in
`inst/provenance/README.md`. These scripts operate independently of any writing
workspace. They do not change packaged data. Author-owned helpers follow the
package's GPL-3 terms; unchanged files in `upstream/` retain the enclosed
pyRiemann BSD-3-Clause license.

From the repository root, check the gorilla permutation against an installed
`shapes` package:

```sh
Rscript --vanilla tools/maintenance/provenance/verify-gorilla-provenance.R
```

The September 21, 2026 ERP replay used Python 3.12.14 and the pinned packages
in `requirements.txt`. These pins document that runtime; they do not claim to
recover the original 2021 environment. Create a separate environment and output
directory, leaving the historical results unchanged:

```sh
erp_helpers="$PWD/tools/maintenance/provenance"
erp_work="$(mktemp -d)"
python3.12 -m venv "$erp_work/python"
"$erp_work/python/bin/python" -m pip install -r "$erp_helpers/requirements.txt"
curl -L --fail --retry 3 'https://osf.io/szh6d/download' \
  -o "$erp_work/sample_audvis_filt-0-40_raw.fif"
curl -L --fail --retry 3 'https://osf.io/xbvt2/download' \
  -o "$erp_work/sample_audvis_filt-0-40_raw-eve.fif"
curl -L --fail --retry 3 \
  'https://s3.amazonaws.com/openneuro.org/ds000248/sub-01/meg/sub-01_task-audiovisual_run-01_meg.fif?versionId=plteJGmBaF0Nlquy40KlEI2ZoHPcBowD' \
  -o "$erp_work/sub-01_task-audiovisual_run-01_meg.fif"
"$erp_work/python/bin/python" "$erp_helpers/check-recording-link.py" \
  --cc0-raw "$erp_work/sub-01_task-audiovisual_run-01_meg.fif" \
  --filtered-raw "$erp_work/sample_audvis_filt-0-40_raw.fif" \
  --events "$erp_work/sample_audvis_filt-0-40_raw-eve.fif" \
  --out "$erp_work/recording-link.json"
Rscript --vanilla "$erp_helpers/export-reference.R" \
  "$PWD/data/ERP.rda" "$erp_work/reference"
"$erp_work/python/bin/python" "$erp_helpers/check-historical-erp.py" \
  --raw "$erp_work/sample_audvis_filt-0-40_raw.fif" \
  --events "$erp_work/sample_audvis_filt-0-40_raw-eve.fif" \
  --historic "$erp_helpers/upstream" \
  --reference-array "$erp_work/reference/package-erp.float64" \
  --reference-labels "$erp_work/reference/package-erp-labels.txt" \
  --relative-tolerance 1e-10 --out "$erp_work/reconstruction"
```

`check-recording-link.py` verifies all three recording checksums before reading
them. `check-historical-erp.py` verifies the pinned upstream numerical-source
hashes, then checks all 216 covariance matrices and the complete label sequence.
The acceptance criterion is maximum **unadjusted** relative Frobenius error
below `1e-10`. Optional sign-adjustment diagnostics do not determine acceptance.
Sources and hashes are listed in `upstream-sources.json`; recording download
URLs and hashes are also retained in `inst/provenance/erp/data-downloads.json`.

The `--historic` path contains pinned sources whose selected classes and
functions are extracted unchanged through Python's abstract syntax tree. This
avoids unrelated legacy imports, without claiming that the entire old package
is installed. The historical source example and license are included alongside
the three numerical source files. The raw recordings, generated NPZ arrays and
reference exports are deliberately generated outside the repository.
