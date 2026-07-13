# %% Import Libraries
import os
import numpy as np
import pandas as pd
import mne
from glob import glob
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from mne.viz import plot_topomap

# %% Paths
raw_dir = '/Users/borjan/code/healthy-sleep-project/data/raw'
data_dir = '/Users/borjan/code/healthy-sleep-project/data'
preprocessed_dir = os.path.join(data_dir, 'preprocessed_data')
results_dir = os.path.join(data_dir, 'results', 'PSD')

# Create necessary directories if they don't exist
os.makedirs(preprocessed_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# %% Configuration Parameters
# Update with your actual EEG channel names
eeg_channels = [
    'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4',
    'O1', 'O2', 'F7', 'F8', 'T3', 'T4', 'T5', 'T6',
    'Fz', 'Cz', 'Pz', 'Oz', 'F9', 'F10',
    'T9', 'T10', 'P9', 'P10', 'AF7',
    'AF3', 'F11', 'F5', 'F1', 'FT11', 'FT9', 'FT7',
    'FC5', 'FC3', 'FC1', 'FCz', 'C5', 'C1', 'TP11',
    'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'P11', 'P5',
    'P1', 'PO7', 'PO3', 'POz', 'Fpz', 'AFz', 'AF4',
    'AF8', 'F2', 'F6', 'F12', 'FC2', 'FC4', 'FC6',
    'FT8', 'FT10', 'FT12', 'C6', 'C2', 'CPz', 'CP2',
    'CP4', 'CP6', 'TP8', 'TP10', 'TP12', 'P2', 'P6',
    'P12', 'PO4', 'PO8'
]

exclude_channels = [
    'ChEMG1', 'ChEMG2', 'RLEG-', 'RLEG+', 'LLEG-', 'LLEG+',
    'EOG1', 'EOG2', 'ECG1', 'ECG2', 'SO1', 'SO2', 'Zy1', 'Zy2'
]

# Frequency bands
freq_bands = {
    'delta': (0.5, 4),
    'theta': (4, 8),
    'alpha': (8, 12),
    'beta': (12, 30),
    'gamma': (30, 40)
}

# Analysis parameters
method = "welch"
fmin = 0.5
fmax = 40
n_fft = 1024
n_per_seg = n_fft
n_overlap = int(n_per_seg / 2)
window = "hamming"

threshold = dict(eeg=400e-6)  # Threshold for bad epoch rejection

# Mind states (sleep stages)
mindstates = ['W', 'N1', 'N2', 'N3', 'R']

# Columns of interest
coi = [
    'sub_id', 'subtype', 'channel', 'mindstate',
    'abs_delta', 'abs_theta', 'abs_alpha', 'abs_beta', 'abs_gamma',
    'rel_delta', 'rel_theta', 'rel_alpha', 'rel_beta', 'rel_gamma'
]

cols_power = [
    'abs_delta', 'abs_theta', 'abs_alpha', 'abs_beta', 'abs_gamma',
    'rel_delta', 'rel_theta', 'rel_alpha', 'rel_beta', 'rel_gamma'
]

bigdic = {f: [] for f in coi}

# Define the montage before using it
montage_file = os.path.join(data_dir, 'Co-registered_average_positions.pos')
if os.path.exists(montage_file):
    montage_data = pd.read_csv(montage_file, sep='\t', header=None)
    electrode_names = montage_data[1].str.upper().tolist()
    coordinates = montage_data[[2, 3, 4]].values
    coordinates[:, [0, 1]] = coordinates[:, [1, 0]]  # Swap X and Y
    ch_pos = dict(zip(electrode_names, coordinates))
    montage = mne.channels.make_dig_montage(ch_pos=ch_pos, coord_frame='head')
else:
    print("Montage file not found. Proceeding without setting montage.")
    # Optionally, you can create a default montage or handle it differently
    montage = None  # Set to None if montage is not available

# %% Main Processing Loop
for subject_num in range(1, 30):
    subject_id = f'EPCTL{subject_num:02d}'
    subject_folder = os.path.join(raw_dir, subject_id)
    edf_file = os.path.join(subject_folder, f"{subject_id}.edf")
    sleep_stage_file = os.path.join(subject_folder, f"{subject_id}.txt")
    
    if not os.path.exists(edf_file) or not os.path.exists(sleep_stage_file):
        print(f"Missing data for subject {subject_id}, skipping.")
        continue
    
    print(f"Processing subject: {subject_id}")
    
    # Load the raw EEG data
    raw = mne.io.read_raw_edf(edf_file, preload=True)
    
    # Apply preprocessing steps
    raw.notch_filter(freqs=[50, 100], method='spectrum_fit')  # Notch filter at 50 and 100 Hz
    raw.filter(1, 45, fir_design='firwin')  # Band-pass filter between 1-45 Hz
    raw.resample(256, npad='auto')  # Resample to 256 Hz
    
    # Load the sleep stage annotations
    sleep_stages_df = pd.read_csv(sleep_stage_file, sep='\t', header=None)
    sleep_stages_df.columns = ['label', 'start_time', 'duration']
    sleep_stage_mapping = {'W': 0, 'N1': 1, 'N2': 2, 'N3': 3, 'R': 4}
    
    # Create annotations for sleep stages
    annotations = mne.Annotations(onset=[], duration=[], description=[])
    for idx, row in sleep_stages_df.iterrows():
        onset = float(row['start_time'])
        duration = float(row['duration'])
        label = row['label']
        annotations.append(onset=onset, duration=duration, description=label)
    
    raw.set_annotations(annotations)
    
    # Set montage using your electrode positions
    # Assuming you have 'electrode_df' and 'montage' from previous steps
    # You can include the code to create 'electrode_df' and 'montage' here
    # For brevity, assuming 'montage' is already created
    raw.rename_channels({ch_name: ch_name.replace('-Ref', '').upper() for ch_name in raw.ch_names if '-Ref' in ch_name})
    raw.set_montage(montage, on_missing='ignore')
    
    # Clean the data
    channels_to_exclude = [ch for ch in raw.ch_names if ch not in eeg_channels + exclude_channels]
    raw_eeg_cleaned = raw.copy().drop_channels(channels_to_exclude)
    
    # Apply average reference
    raw_eeg_cleaned.set_eeg_reference('average', projection=True)
    raw_eeg_cleaned.apply_proj()
    
    # Save the cleaned data for future use
    cleaned_file = os.path.join(preprocessed_dir, f"{subject_id}_cleaned_raw.fif")
    raw_eeg_cleaned.save(cleaned_file, overwrite=True)
    
    # Segment data into epochs based on sleep stages
    events, event_id = mne.events_from_annotations(raw_eeg_cleaned, event_id=sleep_stage_mapping)
    tmin = 0.0
    tmax = 30.0  # 30-second epochs
    
    epochs = mne.Epochs(raw_eeg_cleaned, events, event_id=event_id, tmin=tmin, tmax=tmax, baseline=None, preload=True)
    epochs.drop_bad(threshold)  # Drop bad epochs based on threshold
    
    # Loop over mind states
    for ms in mindstates:
        print(f'Processing mind state: {ms}')
        if ms not in epochs.event_id:
            continue  # Skip if mind state is not present
        
        epochs_ms = epochs[ms]
        if len(epochs_ms) == 0:
            continue
        
        # Compute PSD using compute_psd()
        psd = epochs_ms.compute_psd(
            method=method,
            fmin=fmin,
            fmax=fmax,
            n_fft=n_fft,
            n_overlap=n_overlap,
            n_per_seg=n_per_seg,
            window=window,
            picks=eeg_channels
        )
        
        # Extract PSD data
        psd_data = psd.get_data()
        freqs = psd.freqs
        
        # Compute band power
        for i_epoch in range(len(epochs_ms)):
            this_psd = psd_data[i_epoch]
            
            abs_bandpower_ch = {
                f"abs_{band}": np.nanmean(
                    this_psd[:, np.logical_and(freqs >= borders[0], freqs <= borders[1])],
                    axis=1
                )
                for band, borders in freq_bands.items()
            }
            
            total_power = np.sum(
                [abs_bandpower_ch[f"abs_{band}"] for band in freq_bands.keys()],
                axis=0
            )
            
            rel_bandpower_ch = {
                f"rel_{band}": abs_bandpower_ch[f"abs_{band}"] / total_power
                for band in freq_bands.keys()
            }
            
            for i_ch, channel in enumerate(eeg_channels):
                bigdic['sub_id'].append(subject_id)
                bigdic['subtype'].append('EPCTL')  # Update as needed
                bigdic['mindstate'].append(ms)
                bigdic['channel'].append(channel)
                for col in cols_power:
                    if col.startswith('abs'):
                        bigdic[col].append(abs_bandpower_ch[col][i_ch])
                    if col.startswith('rel'):
                        bigdic[col].append(rel_bandpower_ch[col][i_ch])
        
        # Clear variables to free memory
        del psd, psd_data
    # Clear variables to free memory
    del raw, raw_eeg_cleaned, epochs

# %% Create DataFrame and Save Results
df = pd.DataFrame.from_dict(bigdic)
bandpower_csv = os.path.join(results_dir, "bandpower.csv")
df.to_csv(bandpower_csv, index=False)
print(f"Band power data saved to {bandpower_csv}")


#%% Plotting Topographical Maps
#Assuming you have channel locations available
info = mne.create_info(ch_names=eeg_channels, sfreq=256, ch_types='eeg')
if montage is not None:
    info.set_montage(montage)

subtypes = ['EPCTL']  # Update as needed

for feature in cols_power:
    fig, ax = plt.subplots(
        nrows=1, ncols=len(subtypes), figsize=(16, 4), constrained_layout=True)
    columns = [feature, "channel"]
    
    vmin = df[columns][feature].min()
    vmax = df[columns][feature].max()
    
    for i_st, subtype in enumerate(subtypes):
        temp_values = df[columns].loc[df['subtype'] == subtype]
        temp_values = temp_values.groupby('channel', as_index=False).mean()
        temp_values['channel'] = pd.Categorical(
            temp_values['channel'],
            categories=eeg_channels,
            ordered=True
        )
        df_sorted = temp_values.sort_values('channel')
        values = df_sorted[feature].to_numpy()
        
        divider = make_axes_locatable(ax[i_st])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        im, cm = plot_topomap(
            data=values,
            pos=info,
            axes=ax[i_st],
            contours=2,
            cmap="Purples",
            vlim=(vmin, vmax),
            show=False
        )
        fig.colorbar(im, cax=cax, orientation='vertical')
        ax[i_st].set_title(f"{subtype}", fontweight="bold")
    
    fig.suptitle(f"{feature}\n", fontsize="xx-large", fontweight="bold")
    plt.show()

# %% Statistical Analysis with Linear Mixed Models
# Prepare data for LME
mean_df = df.groupby(
    ['sub_id', 'subtype', 'channel', 'mindstate'],
    as_index=False).mean()

interest = 'rel_alpha'  # Feature of interest
fdr_corrected = True  # Set to True if you want FDR correction

# Define the model formula
model_formula = f"{interest} ~ C(mindstate) + (1 | sub_id)"

# Initialize lists to store results
tvals = []
pvals = []
channels_list = []

for chan in eeg_channels:
    # Subset data for the current channel
    subdf = mean_df[mean_df['channel'] == chan]
    
    # Fit the linear mixed model
    md = smf.mixedlm(model_formula, subdf, groups=subdf['sub_id'])
    mdf = md.fit()
    tvals.append(mdf.tvalues[1])  # Index 1 corresponds to the mindstate effect
    pvals.append(mdf.pvalues[1])
    channels_list.append(chan)

# Apply FDR correction if desired
if fdr_corrected:
    _, corrected_pvals = fdrcorrection(pvals)
    display_pvals = corrected_pvals
else:
    display_pvals = pvals

# Plot the t-values on a topographical map
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
im, cm = plot_topomap(
    data=tvals,
    pos=info,
    axes=ax,
    contours=3,
    mask=np.array(display_pvals) <= 0.05,
    mask_params=dict(marker='o', markerfacecolor='w', markeredgecolor='k',
                     linewidth=0, markersize=6),
    cmap="viridis",
    show=False
)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title(f"T-values for {interest} (Effect of Mind State)", fontweight="bold")
plt.show()
