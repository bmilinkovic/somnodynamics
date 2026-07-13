#%%
import mne
import os
import shutil
import tempfile
import time
import pandas as pd
import numpy as np
import scipy.io as sio  # For saving as .mat files (optional)
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

# Output directory for sensor-space preprocessing.
# Override with HS_PREPROC_SENSOR_DIR if needed.
preprocessed_dir = os.environ.get(
    'HS_PREPROC_SENSOR_DIR',
    '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_sensor_space',
)

# Create necessary directories if they don't exist
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(preprocessed_dir, exist_ok=True)

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

    # Apply notch filter at 50 and 100 Hz
    raw.notch_filter(freqs=[50, 100], method='spectrum_fit')

    # Band-pass filter between 1-45 Hz
    raw.filter(1, 45, fir_design='firwin')

    # Resample to 256 Hz
    raw.resample(256, npad='auto')

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

    # Save the cleaned data to preprocessed_dir as .csv file
    cleaned_data = raw_eeg_cleaned.get_data()
    channel_names = raw_eeg_cleaned.ch_names
    times = raw_eeg_cleaned.times

    # Create a DataFrame with channels as columns and times as index
    df_cleaned = pd.DataFrame(data=cleaned_data.T, index=times, columns=channel_names)

    # Save the DataFrame to CSV
    csv_filename = os.path.join(preprocessed_dir, f"{subject_folder}_cleaned.csv")
    df_cleaned.to_csv(csv_filename)
    print(f"Cleaned data saved to {csv_filename}")

    # Create 10-second epochs for each sleep stage
    epochs_by_stage = {}
    for stage_label, events_df in sleep_stages_df.groupby('label'):
        if stage_label not in sleep_stages_mapping:
            continue  # Skip unwanted labels like 'L'

        stage_data = events_df.copy()
        events = []
        for _, row in stage_data.iterrows():
            onset_sample = time_to_sample(raw_eeg_cleaned, float(row['start_time']))
            duration_samples = int(float(row['duration']) * raw_eeg_cleaned.info['sfreq'])

            # Calculate the number of 10-second epochs in this duration
            epoch_length_samples = int(10 * raw_eeg_cleaned.info['sfreq'])
            num_epochs = duration_samples // epoch_length_samples

            for i in range(num_epochs):
                event_sample = onset_sample + i * epoch_length_samples
                event_id = sleep_stages_mapping[stage_label]
                events.append([event_sample, 0, event_id])

        if events:
            # Create 10-second epochs
            epochs = mne.Epochs(raw_eeg_cleaned, np.array(events), event_id=sleep_stages_mapping[stage_label],
                                tmin=0, tmax=10.0, baseline=None, preload=True)
            epochs_by_stage[stage_label] = epochs
            print(f"Created {len(epochs)} epochs of 10 seconds for {stage_label}")

    # Optionally, you can save the epochs to files if needed
    # For example, save epochs for each sleep stage as .fif files
    for stage_label, epochs in epochs_by_stage.items():
        epochs_filename = os.path.join(preprocessed_dir, f"{subject_folder}_{stage_label}_epochs.fif")
        epochs.save(epochs_filename, overwrite=True)
        print(f"Epochs for stage {stage_label} saved to {epochs_filename}")

# %%
