#%%
import mne
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Paths
raw_dir = '/Users/borjan/code/healthy-sleep-project/data/raw'
results_dir = os.path.join('/Users/borjan/code/healthy-sleep-project/data', 'results', 'PSD-raw')

# Create necessary directories if they don't exist
os.makedirs(results_dir, exist_ok=True)

#%% Loop over all subjects in the raw data directory
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

    # Load the raw EEG data (no filtering, no preprocessing)
    raw = mne.io.read_raw_edf(edf_file, preload=True)

    # Load the sleep stage annotations
    sleep_stages_df = np.loadtxt(sleep_stage_file, delimiter='\t', dtype=str)
    sleep_stages_df = pd.DataFrame(sleep_stages_df, columns=['label', 'start_time', 'duration'])

    # Define the sleep stage mapping, excluding unwanted labels like 'L'
    sleep_stages_mapping = {'W': 1, 'N1': 2, 'N2': 3, 'N3': 4, 'R': 5}

    # Helper function to convert time to sample indices
    def time_to_sample(raw, time_in_seconds):
        return int(time_in_seconds * raw.info['sfreq'])

    # Create epochs for each sleep stage using raw data (no filtering)
    epochs_by_stage = {}
    for stage_label, events_df in sleep_stages_df.groupby('label'):
        if stage_label not in sleep_stages_mapping:
            continue  # Skip unwanted labels like 'L'

        stage_data = events_df.copy()
        events = []
        for _, row in stage_data.iterrows():
            onset_sample = time_to_sample(raw, float(row['start_time']))
            event_id = sleep_stages_mapping[stage_label]
            events.append([onset_sample, 0, event_id])

        if events:
            # Create epochs
            epochs = mne.Epochs(raw, np.array(events), event_id=sleep_stages_mapping[stage_label], 
                                tmin=0, tmax=30.0, baseline=None, preload=True)

            epochs_by_stage[stage_label] = epochs
            print(f"Created {len(epochs)} epochs for {stage_label}.")

    # Power Spectral Density (PSD) Analysis and Plotting
    for stage_label, epochs in epochs_by_stage.items():
        # Concatenate all epochs into a continuous segment
        concatenated_data = np.concatenate(epochs.get_data(), axis=1)  # (n_channels, total_time_points)
        concatenated_raw = mne.io.RawArray(concatenated_data, epochs.info)  # Convert back to Raw format

        # Ensure the results directory for the subject exists
        subject_results_dir = os.path.join(results_dir, subject_folder)
        os.makedirs(subject_results_dir, exist_ok=True)

        # Plot the PSD for the entire sleep stage (concatenated data)
        concatenated_raw.plot_psd(fmax=45, average=True, dB=True)
        plt.title(f"PSD for {stage_label}")
        plt.savefig(os.path.join(subject_results_dir, f"{subject_folder}_{stage_label}_PSD_avg.png"))
        plt.close()

        print(f"Saved PSD plot for {stage_label} in subject {subject_folder}")

# %%
