"""Read-only provenance test of pyRiemann v0.2.7's historical ERP recipe.

Numerical classes/functions are AST-extracted without changes from upstream
commit 176e766f540bd4c7846f38573165fc3d27fc69ca. Extraction avoids importing
unrelated legacy modules requiring distutils. NumPy/SciPy/sklearn/MNE remain
the explicitly recorded current investigation versions, not a claimed 2021
environment. Files are written only to the explicitly selected output folder.
This audit helper is supplied under the Riemann package's GPL-3 terms;
the unchanged upstream excerpts retain their separate BSD-3-Clause license.
"""
import argparse, ast, hashlib, json, os, pathlib, platform, sys, time

parser = argparse.ArgumentParser()
parser.add_argument('--raw', required=True)
parser.add_argument('--events', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--historic', required=True, help='Root containing pinned pyriemann/ sources')
parser.add_argument('--reference-array', required=True, help='R-exported little-endian float64 array')
parser.add_argument('--reference-labels', required=True, help='One stored class label per line')
parser.add_argument('--relative-tolerance', type=float, default=1e-10)
args = parser.parse_args()
if not 0 < args.relative_tolerance < 1:
    parser.error('--relative-tolerance must lie strictly between zero and one')
started = time.perf_counter()
out = pathlib.Path(args.out); out.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(out/'mpl'))
import numpy as np
import scipy, sklearn, mne, pyriemann
from scipy.linalg import eigh, inv
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.covariance import empirical_covariance, oas, ledoit_wolf, fast_mcd
from sklearn.model_selection import train_test_split

def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

upstream = pathlib.Path(args.historic)/'pyriemann'
ns = dict(np=np, eigh=eigh, inv=inv, BaseEstimator=BaseEstimator,
          TransformerMixin=TransformerMixin, empirical_covariance=empirical_covariance,
          oas=oas, ledoit_wolf=ledoit_wolf, fast_mcd=fast_mcd)
selected = {
    'utils/covariance.py': ['_lwf', '_oas', '_scm', '_mcd', '_check_est', 'covariances_EP'],
    'spatialfilters.py': ['Xdawn'],
    'estimation.py': ['XdawnCovariances'],
}
expected_sources = {
    'utils/covariance.py':'58568c6c081dfdcaf1e495dc61ae1a234916b1846609271f5a5909a933d3f27d',
    'spatialfilters.py':'14f07c2cebaa4799e6cb2280da9a1bfa78a12e9d2e9855643091dbf51bb05332',
    'estimation.py':'4b746651a51b960c6dbefe70dcd7562c9b2eb3c28519f18516d40c62035f5fd3',
}
source_hashes = {}
for file, names in selected.items():
    path = upstream/file; source_hashes[file] = sha(path)
    if source_hashes[file] != expected_sources[file]:
        raise ValueError(f'Pinned historical source checksum mismatch: {path}')
    tree = ast.parse(path.read_text())
    nodes = [x for x in tree.body if isinstance(x, (ast.ClassDef,ast.FunctionDef)) and x.name in names]
    assert {x.name for x in nodes} == set(names)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)

raw = mne.io.read_raw_fif(args.raw, preload=True, verbose=False)
raw_before = {'sfreq':float(raw.info['sfreq']), 'n_times':int(raw.n_times),
              'first_samp':int(raw.first_samp), 'last_samp':int(raw.last_samp),
              'lowpass':float(raw.info['lowpass']), 'highpass':float(raw.info['highpass'])}
raw.filter(2,None,method='iir')
events = mne.read_events(args.events)
raw.info['bads'] = ['MEG 2443']
picks = mne.pick_types(raw.info,meg=True,eeg=False,stim=False,eog=False,exclude='bads')
event_id = dict(aud_l=1,aud_r=2,vis_l=3,vis_r=4)
epochs = mne.Epochs(raw,events,event_id,0.,1.,proj=False,picks=picks,
                    baseline=None,preload=True,verbose=False)
X = epochs.get_data(copy=True); y = epochs.events[:,-1]
indices = np.arange(len(y))
train,test = train_test_split(indices,train_size=.25,random_state=42)
old = ns['XdawnCovariances'](estimator='scm',nfilter=4)
covs = old.fit(X[train],y[train]).transform(X[test])
modern = pyriemann.estimation.XdawnCovariances(estimator='scm',nfilter=4)
modern_covs = modern.fit(X[train],y[train]).transform(X[test])
reference = np.fromfile(args.reference_array,dtype='<f8').reshape((32,32,216),order='F').transpose(2,0,1)
labels = np.array(pathlib.Path(args.reference_labels).read_text().splitlines())
if labels.shape != (216,) or not np.isfinite(reference).all():
    raise ValueError('Reference must contain 216 labels and finite covariance matrices')
label_names = np.array(['','LA','RA','LV','RV'])

def metrics(a,b):
    if a.shape != b.shape:return {'shape_matches':False,'shapes':[list(a.shape),list(b.shape)]}
    diff=a-b; norms=np.linalg.norm(b,axis=(1,2)); rel=np.linalg.norm(diff,axis=(1,2))/norms
    return {'shape_matches':True,'byte_equal':bool(np.array_equal(a,b)),
            'maximum_absolute_element_error':float(np.max(np.abs(diff))),
            'maximum_relative_frobenius_error':float(rel.max()),
            'median_relative_frobenius_error':float(np.median(rel))}

result = {'versions':{k:v.__version__ for k,v in [('numpy',np),('scipy',scipy),('sklearn',sklearn),('mne',mne),('pyriemann_modern',pyriemann)]},
 'runtime':{'python':sys.version,'platform':platform.platform()},
 'script_sha256':sha(__file__),
 'historical_commit':'176e766f540bd4c7846f38573165fc3d27fc69ca',
 'historical_source_sha256':source_hashes,
 'inputs':{str(p):sha(p) for p in (args.raw,args.events,args.reference_array,args.reference_labels)},
 'raw':raw_before, 'epoch_shape':list(X.shape),'picked_channel_names':[raw.ch_names[i] for i in picks],
 'epoch_counts':{str(k):int(np.sum(y==k)) for k in event_id.values()},
 'train_counts':{str(k):int(np.sum(y[train]==k)) for k in event_id.values()},
 'test_counts':{str(k):int(np.sum(y[test]==k)) for k in event_id.values()},
 'test_label_sequence_matches':bool(np.array_equal(label_names[y[test]],labels)),
 'historical_vs_package':metrics(covs,reference),
 'modern_vs_historical':metrics(modern_covs,covs),
 'package_prototype_block_identical':bool(np.all(reference[:,:16,:16]==reference[0,:16,:16]))}

# Eigenvector signs may differ by numerical library. This diagnostic permits
# only one global repeated 16-feature sign vector, inferred from the fixed
# prototype block and then evaluated on all full matrices. It is not a refit,
# arbitrary congruence, scale correction, or evidence of byte identity.
if covs.shape==reference.shape:
    weight=covs[0,:16,:16]*reference[0,:16,:16]
    signs=np.ones(16); visited={0}
    while len(visited)<16:
        i,j=max(((i,j) for i in visited for j in range(16) if j not in visited),key=lambda ij:abs(weight[ij]))
        signs[j]=signs[i]*np.sign(weight[i,j]);visited.add(j)
    ds=np.tile(signs,2)
    adjusted=covs*ds[None,:,None]*ds[None,None,:]
    result['prototype_inferred_sign_vector']=signs.tolist()
    result['sign_adjusted_vs_package']=metrics(adjusted,reference)
    result['minimum_eigenvalues']={'historical':float(np.linalg.eigvalsh(covs).min()),'package':float(np.linalg.eigvalsh(reference).min())}
np.savez(out/'historical-results.npz',covariances=covs,modern_covariances=modern_covs,
         labels=y[test],train_indices=train,test_indices=test,events=epochs.events,
         filters=old.Xd_.filters_,prototypes=old.P_)
comparison = result['historical_vs_package']
result['acceptance'] = {
    'metric':'maximum unadjusted relative Frobenius norm over all 216 matrices',
    'relative_tolerance':args.relative_tolerance,
    'passed':bool(comparison['shape_matches'] and result['test_label_sequence_matches']
                  and comparison['maximum_relative_frobenius_error'] < args.relative_tolerance),
    'sign_adjusted_result_used_for_acceptance':False,
}
result['elapsed_seconds'] = time.perf_counter()-started
(out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='picked_channel_names'},indent=2))
if not result['acceptance']['passed']:
    raise SystemExit('ERP reconstruction failed the relative-error or label-order check')
