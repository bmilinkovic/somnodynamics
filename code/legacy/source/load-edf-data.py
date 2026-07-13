#%%
import mne
import os
import shutil
import tempfile
import time
import pandas as pd
import numpy as np
import scipy.io as sio  # For saving as .mat files
import matplotlib.pyplot as plt
from mne import Epochs
from mne.beamformer import make_lcmv, apply_lcmv_epochs

# Paths
# Input data root (contains "raw/" and "Co-registered_average_positions.pos").
# Override with APHNY_DATA_DIR if needed.
data_dir = os.environ.get(
    'APHNY_DATA_DIR',
    '/Users/borjan/Library/CloudStorage/OneDrive-ICM/APHNY data',
)
raw_dir = os.path.join(data_dir, 'raw')

# Output directory for source-reconstructed time series.
# Override with HS_PREPROC_SOURCE_DIR if needed.
preprocessed_dir = os.environ.get(
    'HS_PREPROC_SOURCE_DIR',
    '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data',
)
results_dir = os.path.join(preprocessed_dir, '_legacy_aux_results')

# Create necessary directories if they don't exist
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(preprocessed_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# Helper function to convert time to sample indices
def time_to_sample(raw, time_in_seconds):
    return int(time_in_seconds * raw.info['sfreq'])


def get_local_montage_file(data_dir):
    env_montage = os.environ.get('APHNY_MONTAGE_FILE', '')
    candidates = []
    if env_montage:
        candidates.append(env_montage)
    candidates.extend([
        os.path.join(data_dir, 'Co-registered_average_positions.pos'),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Co-registered_average_positions.pos'),
    ])
    for path in candidates:
        if not path:
            continue
        if not os.path.exists(path):
            continue
        temp_path = os.path.join(tempfile.gettempdir(), 'aphny_montage_cache.pos')
        last_error = None
        for attempt in range(5):
            try:
                shutil.copyfile(path, temp_path)
                if file_is_all_zero_bytes(temp_path):
                    raise RuntimeError(f'{path} appears to be an all-zero cloud placeholder')
                return temp_path
            except Exception as err:
                last_error = err
                if attempt < 4:
                    time.sleep(2.0)
        raise RuntimeError(
            f'Could not copy montage file from {path}. '
            "If using OneDrive, mark it 'Always keep on this device'."
        ) from last_error
    raise FileNotFoundError('Montage file not found. Set APHNY_MONTAGE_FILE to a local copy.')


def file_is_all_zero_bytes(path, sample_size=512):
    with open(path, 'rb') as handle:
        sample = handle.read(sample_size)
    return bool(sample) and all(byte == 0 for byte in sample)


def read_montage_table(montage_path):
    try:
        return pd.read_csv(montage_path, sep='\t', header=None)
    except Exception:
        return pd.read_csv(montage_path, sep=r'\s+', header=None, engine='python')

#%% Loop over all subjects in the data directory
for subject_folder in os.listdir(raw_dir):
    subject_path = os.path.join(raw_dir, subject_folder)
    if not os.path.isdir(subject_path):
        continue  # Skip non-directories

    print(f"Processing subject: {subject_folder}")

    # Define the subject EDF and TXT files
    edf_file = os.path.join(subject_path, f"{subject_folder}.edf")
    sleep_stage_file = os.path.join(subject_path, f"{subject_folder}.txt")

    if not os.path.exists(edf_file) or not os.path.exists(sleep_stage_file):
        print(f"Missing data for subject {subject_folder}, skipping.")
        continue

    # Load the raw EEG data
    raw = mne.io.read_raw_edf(edf_file, preload=True)
    raw_filtered = raw.notch_filter(freqs=[50, 100], method='spectrum_fit')  # Notch filter 50 and 100 Hz
    raw = raw.filter(1, 45, fir_design='firwin')  # Band-pass filter between 1-45 Hz
    raw = raw.resample(256, npad='auto')  # Resample to 256 Hz

    # Load the sleep stage annotations
    sleep_stages_df = pd.read_csv(sleep_stage_file, sep='\t', header=None)
    sleep_stages_df.columns = ['label', 'start_time', 'duration']
    sleep_stages_mapping = {'W': 1, 'N1': 2, 'N2': 3, 'N3': 4, 'R': 5}

    # Define the montage
    coregistration_file = get_local_montage_file(data_dir)
    montage_data = read_montage_table(coregistration_file)
    electrode_names = montage_data[1].str.upper().tolist()
    coordinates = montage_data[[2, 3, 4]].values
    coordinates[:, [0, 1]] = coordinates[:, [1, 0]]  # Swap X and Y
    montage = mne.channels.make_dig_montage(ch_pos=dict(zip(electrode_names, coordinates)), coord_frame='head')

    # Set montage and clean data
    raw.rename_channels({ch_name: ch_name.replace('-Ref', '').upper() for ch_name in raw.ch_names if '-Ref' in ch_name})
    raw.set_montage(montage, on_missing='ignore')

    exclude_channels = ['ChEMG1', 'ChEMG2', 'RLEG-', 'RLEG+', 'LLEG-', 'LLEG+', 'EOG1', 'EOG2', 'ECG1', 'ECG2', 'SO1', 'SO2', 'ZY1', 'ZY2']

    # Identify bad channels that are missing location info
    bad_channels = [ch['ch_name'] for ch in raw.info['chs'] if not np.isfinite(ch['loc'][:3]).all()]

    # Exclude channels that are missing locations and non-EEG channels
    raw_eeg_cleaned = raw.copy().drop_channels(bad_channels + exclude_channels)

    # Apply average reference to the cleaned EEG data and apply projection
    raw_eeg_cleaned.set_eeg_reference('average', projection=True)
    raw_eeg_cleaned.apply_proj()

    # Create 2-second epochs directly for each sleep stage
    epochs_by_stage = {}
    for stage_label, events_df in sleep_stages_df.groupby('label'):
        if stage_label not in sleep_stages_mapping:
            continue  # Skip unwanted labels like 'L'
        
        stage_data = events_df.copy()
        events = []
        for _, row in stage_data.iterrows():
            onset_sample = time_to_sample(raw_eeg_cleaned, float(row['start_time']))
            event_id = sleep_stages_mapping[stage_label]
            events.append([onset_sample, 0, event_id])

        if events:
            # Create 2-second epochs directly
            epochs = mne.Epochs(raw_eeg_cleaned, np.array(events), event_id=sleep_stages_mapping[stage_label], 
                                tmin=0, tmax=2.0, baseline=None, preload=True)
            epochs_by_stage[stage_label] = epochs
            print(f"Created {len(epochs)} epochs for {stage_label}")

    # Compute noise covariance
    noise_cov = mne.compute_raw_covariance(raw_eeg_cleaned, method='empirical')

    # Create source space and forward solution
    subjects_dir = '/Users/borjan/mne_data/MNE-fsaverage-data'
    src = mne.setup_source_space('fsaverage', spacing='ico4', subjects_dir=subjects_dir, add_dist=False)
    bem_model = mne.make_bem_model(subject='fsaverage', ico=4, subjects_dir=subjects_dir)
    bem_solution = mne.make_bem_solution(bem_model)
    fwd = mne.make_forward_solution(raw_eeg_cleaned.info, trans='fsaverage', src=src, bem=bem_solution, eeg=True, mindist=5.0)

    # Compute LCMV inverse operator using cleaned EEG data
    inverse_operator = make_lcmv(raw_eeg_cleaned.info, fwd, noise_cov, reg=0.05, pick_ori='max-power', weight_norm='nai')

    # Save source time-series
    labels = mne.read_labels_from_annot('fsaverage', 'aparc', subjects_dir=subjects_dir)

    subject_preprocessed_dir = os.path.join(preprocessed_dir, subject_folder)
    os.makedirs(os.path.join(subject_preprocessed_dir, 'csv'), exist_ok=True)
    os.makedirs(os.path.join(subject_preprocessed_dir, 'mat'), exist_ok=True)

    for stage_label, epochs in epochs_by_stage.items():
        stcs = apply_lcmv_epochs(epochs, inverse_operator)

        # Parcellation and saving
        epoch_parcellation_list = []
        for stc in stcs:
            label_ts = mne.extract_label_time_course(stc, labels, fwd['src'], mode='mean_flip', return_generator=False, allow_empty=True)
            epoch_parcellation_list.append(label_ts)

        epoch_parcellation_array = np.transpose(np.array(epoch_parcellation_list), (1, 2, 0))

        # **Exclude rows (regions) where all time points are zero**
        non_zero_regions = np.any(epoch_parcellation_array != 0, axis=(1, 2))  # Check for regions that have non-zero time points
        epoch_parcellation_array = epoch_parcellation_array[non_zero_regions]  # Only keep non-zero regions
        valid_labels = [label.name for label, is_non_zero in zip(labels, non_zero_regions) if is_non_zero]  # Filter labels

        # Save CSV
        reshaped_array = epoch_parcellation_array.transpose(2, 1, 0)  # epochs x time points x regions
        flattened_data = reshaped_array.reshape(-1, reshaped_array.shape[-1])  # Flatten to (epochs * time points) x regions
        stage_df = pd.DataFrame(flattened_data, columns=valid_labels)
        csv_file_path = os.path.join(subject_preprocessed_dir, 'csv', f"{subject_folder}_{stage_label}_source_time_series.csv")
        stage_df.to_csv(csv_file_path, index=False)

        # Save as MAT
        mat_file_path = os.path.join(subject_preprocessed_dir, 'mat', f"{subject_folder}_{stage_label}_source_time_series.mat")
        sio.savemat(mat_file_path, {"source_ts": epoch_parcellation_array, "labels": valid_labels})

        print(f"Saved source time-series for {stage_label} in subject {subject_folder}")


# %%
