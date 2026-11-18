"""Bounded metadata, trigger-signal, and event correspondence check.

This does not compare filtered MEG waveforms or claim file-byte identity.
Supplied under the Riemann package's GPL-3 terms.
"""
import argparse, hashlib, json, os, pathlib, platform, sys

p = argparse.ArgumentParser()
p.add_argument('--cc0-raw', required=True)
p.add_argument('--filtered-raw', required=True)
p.add_argument('--events', required=True)
p.add_argument('--out', required=True)
a = p.parse_args()
output = pathlib.Path(a.out)
output.parent.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(output.parent/'mpl'))
import mne
import numpy as np

def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

expected = {
    'cc0_raw':'d92a44c6e6e62503ec75e9e07e70e9b650cf444bc0f0e8ff37e2bed162548966',
    'filtered_raw':'327e163c9d4ee6855d41d887d86cec3fa7d6dea3c870e8263edf802889b4f844',
    'events':'63b6cf49907f136c4ff0ad76a19db8db2ae3a404848922a04f9be4e700091b1f',
}
hashes = {key:sha(getattr(a,key)) for key in expected}
if hashes != expected:
    raise ValueError('Input checksum disagrees with the documented recording')
raw = mne.io.read_raw_fif(a.cc0_raw, preload=False, verbose=False)
filtered = mne.io.read_raw_fif(a.filtered_raw, preload=False, verbose=False)
events = mne.read_events(a.events)
raw_events = mne.find_events(raw, stim_channel='STI 014', shortest_event=1, verbose=False)
filtered_events = mne.find_events(filtered, stim_channel='STI 014', shortest_event=1, verbose=False)
raw_stim = raw.get_data(picks=['STI 014'])[0]
filtered_stim = filtered.get_data(picks=['STI 014'])[0]
block_max = raw_stim.reshape(-1,4).max(axis=1)
matched = []
for event in events:
    candidates = np.flatnonzero((raw_events[:,2] == event[2]) &
                                (raw_events[:,0]//4 == event[0]))
    if len(candidates) != 1:
        raise ValueError('Supplied filtered event has no unique raw-event match')
    matched.append(int(candidates[0]))
unmatched = np.setdiff1d(np.arange(len(raw_events)),matched)
checks = {
    'both_have_376_channels':len(raw.ch_names)==len(filtered.ch_names)==376,
    'ordered_channel_names_equal':raw.ch_names==filtered.ch_names,
    'channel_kinds_equal':all(x['kind']==y['kind'] for x,y in zip(raw.info['chs'],filtered.info['chs'])),
    'channel_units_equal':all(x['unit']==y['unit'] for x,y in zip(raw.info['chs'],filtered.info['chs'])),
    'sampling_frequency_exact_factor_four':raw.info['sfreq']==4*filtered.info['sfreq'],
    'sample_count_exact_factor_four':raw.n_times==4*filtered.n_times,
    'first_sample_exact_factor_four':raw.first_samp==4*filtered.first_samp,
    'half_open_last_sample_exact_factor_four':raw.last_samp+1==4*(filtered.last_samp+1),
    'all_41700_stimulus_blocks_match_maximum':np.array_equal(block_max,filtered_stim),
    'supplied_319_events_recreated_exactly':np.array_equal(events,filtered_events) and len(events)==319,
    'each_supplied_event_has_unique_code_and_floor_sample_match':len(matched)==len(events),
    'one_expected_raw_rising_event_omitted':np.array_equal(raw_events[unmatched],np.array([[142826,0,4]])),
}
def info(r):
    return {'n_channels':len(r.ch_names),'sampling_frequency':float(r.info['sfreq']),
            'n_samples':int(r.n_times),'first_sample':int(r.first_samp),
            'last_sample':int(r.last_samp),'bad_channels':r.info['bads']}
transitions = mne.find_events(filtered,stim_channel='STI 014',shortest_event=1,
                             consecutive=True,verbose=False)
result = {'input_sha256':hashes,'checks':{k:bool(v) for k,v in checks.items()},
    'passed':all(checks.values()),'cc0_raw':info(raw),'filtered_raw':info(filtered),
    'raw_rising_event_count':len(raw_events),'supplied_filtered_event_count':len(events),
    'unmatched_raw_rising_events':raw_events[unmatched].tolist(),
    'filtered_transitions_near_omitted_event':transitions[(transitions[:,0]>35700)&(transitions[:,0]<35710)].tolist(),
    'event_sample_mapping':'filtered_sample = floor(raw_sample / 4) for all 319 supplied events',
    'stimulus_sample_mapping':'filtered_STI014[j] = max(raw_STI014[4*j:4*j+4]) for all 41700 samples',
    'scope':'Metadata, trigger signal, and event linkage only; no filtered MEG waveform or file-byte identity claim.',
    'runtime':{'python':sys.version,'platform':platform.platform(),'mne':mne.__version__,'numpy':np.__version__},
    'script_sha256':sha(__file__)}
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if not result['passed']:
    raise SystemExit('Recording-link check failed')
