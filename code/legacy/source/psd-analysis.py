#%%
import mne
import os
import numpy as np
import scipy.io as sio  # For saving as .mat files
import matplotlib.pyplot as plt
import pandas as pd
import h5py  # New import to handle MATLAB v7.3 files

# Paths
raw_dir = '/Users/borjan/code/healthy-sleep-project/data/raw'
artifact_dir = '/Users/borjan/code/healthy-sleep-project/data/Artifact matrix'
results_dir = os.path.join('/Users/borjan/code/healthy-sleep-project/data', 'results', 'PSD')

# Create necessary directories if they don't exist
os.makedirs(results_dir, exist_ok=True)

# Helper function to convert time to sample indices
def time_to_sample(raw, time_in_seconds):
    return int(time_in_seconds * raw.info['sfreq'])

#%% Loop over all subjects in the raw data directory
for subject_folder in os.listdir(raw_dir):
    subject_path = os.path.join(raw_dir, subject_folder)
    if not os.path.isdir(subject_path):
        continue  # Skip non-directories

    print(f"Processing subject: {subject_folder}")

    # Define the subject EDF and TXT files
    edf_file = os.path.join(subject_path, f"{subject_folder}.edf")
    sleep_stage_file = os.path.join(subject_path, f"{subject_folder}.txt")
    artifact_file = os.path.join(artifact_dir, f"{subject_folder}_artndxn.mat")

    if not os.path.exists(edf_file) or not os.path.exists(sleep_stage_file) or not os.path.exists(artifact_file):
        print(f"Missing data for subject {subject_folder}, skipping.")
        continue

    # Load the raw EEG data
    raw = mne.io.read_raw_edf(edf_file, preload=True)
    raw_filtered = raw.notch_filter(freqs=[50, 100], method='spectrum_fit')  # Notch filter 50 and 100 Hz
    raw_filtered = raw_filtered.filter(1, 45, fir_design='firwin')  # Band-pass filter between 1-45 Hz
    raw_filtered = raw_filtered.resample(256, npad='auto')  # Resample to 256 Hz

    # Load the sleep stage annotations
    sleep_stages_df = np.loadtxt(sleep_stage_file, delimiter='\t', dtype=str)
    sleep_stages_df = pd.DataFrame(sleep_stages_df, columns=['label', 'start_time', 'duration'])

    sleep_stages_mapping = {'W': 1, 'N1': 2, 'N2': 3, 'N3': 4, 'R': 5}

    # Load the artifact matrix using h5py (MATLAB v7.3 files)
    with h5py.File(artifact_file, 'r') as f:
        artifact_data = np.array(f['artndxn'])  # Load the artifact matrix (channels x epochs)

    n_epochs = artifact_data.shape[1]

    # Identify bad 30-second epochs (where any channel in that epoch contains a 0)
    bad_epochs = np.where(np.any(artifact_data == 0, axis=0))[0]  # Find columns with any 0s (bad epochs)

    # Create 30-second epochs for each sleep stage
    epochs_by_stage = {}
    for stage_label, events_df in sleep_stages_df.groupby('label'):
        if stage_label not in sleep_stages_mapping:
            continue  # Skip unwanted labels like 'L'
        
        stage_data = events_df.copy()
        events = []
        for _, row in stage_data.iterrows():
            onset_sample = time_to_sample(raw_filtered, float(row['start_time']))
            event_id = sleep_stages_mapping[stage_label]
            events.append([onset_sample, 0, event_id])

        if events:
            # Create 30s epochs
            epochs = mne.Epochs(raw_filtered, np.array(events), event_id=sleep_stages_mapping[stage_label], 
                                tmin=0, tmax=30.0, baseline=None, preload=True)

            # Drop the bad 30s epochs based on artifact detection
            epochs.drop(bad_epochs, reason='artifact')

            # Split the remaining clean epochs into smaller 2-second epochs
            epochs_resampled = epochs.crop(tmin=0, tmax=2.0)  # Adjust this if needed for smaller epoch lengths
            epochs_by_stage[stage_label] = epochs_resampled
            print(f"Created {len(epochs_resampled)} epochs for {stage_label} after artifact rejection.")

    # Power Spectral Density (PSD) Analysis and Plotting
    for stage_label, epochs in epochs_by_stage.items():
        # Concatenate all epochs into a continuous segment
        concatenated_data = np.concatenate(epochs.get_data(), axis=1)  # (n_channels, total_time_points)
        concatenated_raw = mne.io.RawArray(concatenated_data, epochs.info)  # Convert back to Raw format

        # Ensure the results directory for the subject exists
        subject_results_dir = os.path.join(results_dir, subject_folder)
        os.makedirs(subject_results_dir, exist_ok=True)

        # Plot the PSD for the entire sleep stage (concatenated data)
        concatenated_raw.plot_psd(method='welch', fmax=45, average=True, estimate='power', dB=True, n_jobs=-1)
        plt.title(f"PSD for {stage_label}")
        plt.savefig(os.path.join(subject_results_dir, f"{subject_folder}_{stage_label}_PSD_avg.png"))
        plt.close()

        # Compute PSD for the entire concatenated data
        power_freq = concatenated_raw.compute_psd(method='welch', fmax=45, n_jobs=-1)

        # Get the frequency values and PSD data
        psd, freqs = power_freq.get_data(return_freqs=True)
        mean_psd_per_stage = np.mean(psd, axis=0)

        # Save raw PSD data
        sio.savemat(os.path.join(subject_results_dir, f"{subject_folder}_{stage_label}_PSD.mat"),
                    {"psd": mean_psd_per_stage, "freqs": freqs})

        print(f"Saved PSD plot and data for {stage_label} in subject {subject_folder}")

        #### New Code for Additional Plots ####

        # 1. PSD per Epoch (each epoch's PSD plotted individually with the average)
        fig, ax = plt.subplots()
        psds_per_epoch = epochs.compute_psd(method='welch', fmax=45, n_jobs=-1)
        for i, psd_epoch in enumerate(psds_per_epoch.get_data()):
            ax.plot(freqs, 10 * np.log10(psd_epoch.mean(axis=0)), label=f'Epoch {i+1}', alpha=0.3)
        
        # Plot average across epochs
        mean_psd_all_epochs = 10 * np.log10(np.mean(psds_per_epoch.get_data(), axis=0).mean(axis=0))
        ax.plot(freqs, mean_psd_all_epochs, color='black', label='Average PSD', linewidth=2)
        ax.set_title(f"PSD per Epoch for {stage_label}")
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power (dB)')
        ax.legend(loc='upper right')

        # Save the PSD per epoch plot
        plt.savefig(os.path.join(subject_results_dir, f"{subject_folder}_{stage_label}_PSD_per_epoch.png"))
        plt.close()

        # 2. PSD per Channel (each channel's PSD plotted individually)
        fig, ax = plt.subplots()
        for ch_idx in range(len(epochs.info['ch_names'])):
            ax.plot(freqs, 10 * np.log10(psds_per_epoch.get_data().mean(axis=0)[ch_idx]), label=epochs.info['ch_names'][ch_idx], alpha=0.3)
        
        # Plot the average across channels
        mean_psd_all_channels = 10 * np.log10(np.mean(psds_per_epoch.get_data().mean(axis=0), axis=0))
        ax.plot(freqs, mean_psd_all_channels, color='black', label='Average PSD', linewidth=2)
        ax.set_title(f"PSD per Channel for {stage_label}")
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power (dB)')
        ax.legend(loc='upper right')

        # Save the PSD per channel plot
        plt.savefig(os.path.join(subject_results_dir, f"{subject_folder}_{stage_label}_PSD_per_channel.png"))
        plt.close()

        print(f"Saved additional PSD plots for {stage_label} in subject {subject_folder}")


# %%
