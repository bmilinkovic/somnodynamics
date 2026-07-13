# bandpower_calculation.py

import os
import glob
import mne
import numpy as np
import pandas as pd

# Paths (adjust these paths as needed)
preprocessed_dir = '/Users/borjan/code/healthy-sleep-project/data/preprocessed_sensor_space'
bandpower_dir = '/Users/borjan/code/healthy-sleep-project/results/bandpower'

# Create the bandpower directory if it doesn't exist
os.makedirs(bandpower_dir, exist_ok=True)

# Parameters for PSD calculation
freqs = np.linspace(0.5, 40, 159)
method = "welch"
fmin = 0.5
fmax = 40
n_fft = 1024
n_per_seg = n_fft
n_overlap = int(n_per_seg / 2)
window = "hamming"

freq_bands = {
    'delta': (0.5, 4),
    'theta': (4, 8),
    'alpha': (8, 12),
    'sigma': (12, 16),
    'beta': (16, 30),
    'gamma': (30, 40)
}

threshold = dict(eeg=400e-6)

# Columns of interest
coi = ['sub_id', 'sleep_state', 'channel',
       'abs_delta', 'abs_theta', 'abs_alpha', 'abs_sigma', 'abs_beta', 'abs_gamma',
       'rel_delta', 'rel_theta', 'rel_alpha', 'rel_sigma', 'rel_beta', 'rel_gamma']

cols_power = [
    'abs_delta', 'abs_theta', 'abs_alpha', 'abs_sigma', 'abs_beta', 'abs_gamma',
    'rel_delta', 'rel_theta', 'rel_alpha', 'rel_sigma', 'rel_beta', 'rel_gamma'
]

# Sleep states
sleep_states = ["W", "N1", "N2", "N3", "R"]

# Get list of epochs files
files = glob.glob(os.path.join(preprocessed_dir, "*_epochs.fif"))

bigdic = {f: [] for f in coi}

for i_file, file in enumerate(files):
    sub_id = os.path.basename(file).split('_')[0]
    print(f"...Processing {sub_id}")

    # Load epochs
    epochs = mne.read_epochs(file, preload=True)
    epochs.drop_bad(threshold)

    # Ensure metadata is available
    metadata = epochs.metadata
    if metadata is None or 'sleep_state' not in metadata.columns:
        print(f"No sleep_state metadata found in {file}, skipping.")
        continue

    for ss in sleep_states:
        print(f'Processing sleep state: {ss}')
        if ss not in metadata['sleep_state'].unique():
            continue

        # Select epochs for this sleep state
        epochs_ss = epochs[epochs.metadata['sleep_state'] == ss]

        if len(epochs_ss) == 0:
            continue

        # Compute PSD
        psds, freqs = epochs_ss.compute_psd(
            method=method,
            fmin=fmin,
            fmax=fmax,
            n_fft=n_fft,
            n_overlap=n_overlap,
            n_per_seg=n_per_seg,
            window=window,
            picks='eeg'
        ).get_data(return_freqs=True)

        # Loop over epochs
        for i_epoch in range(len(epochs_ss)):
            this_power = np.squeeze(psds[i_epoch])

            abs_bandpower_ch = {
                f"abs_{band}": np.nanmean(
                    this_power[:, np.logical_and(freqs >= borders[0], freqs <= borders[1])],
                    axis=1
                )
                for band, borders in freq_bands.items()
            }

            total_power = np.sum(
                [abs_bandpower_ch[f"abs_{band}"] for band in freq_bands.keys()],
                axis=0  # Sum over frequency bands
            )

            rel_bandpower_ch = {
                f"rel_{band}": abs_bandpower_ch[f"abs_{band}"] / total_power
                for band in freq_bands.keys()
            }

            for i_ch, channel in enumerate(epochs_ss.ch_names):
                bigdic['sub_id'].append(sub_id)
                bigdic['sleep_state'].append(ss)
                bigdic['channel'].append(channel)
                for col in cols_power:
                    if col.startswith('abs'):
                        bigdic[col].append(abs_bandpower_ch[col][i_ch])
                    elif col.startswith('rel'):
                        bigdic[col].append(rel_bandpower_ch[col][i_ch])

# Convert dictionary to DataFrame and save as CSV
df = pd.DataFrame.from_dict(bigdic)
df.to_csv(os.path.join(bandpower_dir, "bandpower.csv"), index=False)
print(f"Bandpower data saved to {os.path.join(bandpower_dir, 'bandpower.csv')}")
