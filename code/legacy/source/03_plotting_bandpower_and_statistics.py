# plotting_and_statistics.py

import os
import glob
import mne
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Paths (adjust these paths as needed)
bandpower_dir = '/Users/borjan/code/healthy-sleep-project/results/bandpower'
preprocessed_dir = '/Users/borjan/code/healthy-sleep-project/data/preprocessed_sensor_space'
plots_dir = '/Users/borjan/code/healthy-sleep-project/results/bandpower/plots'

# Create the plots directory if it doesn't exist
os.makedirs(plots_dir, exist_ok=True)

# Load bandpower data
df = pd.read_csv(os.path.join(bandpower_dir, "bandpower.csv"))

# Load an example epochs file to get channel info
epochs_file = glob.glob(os.path.join(preprocessed_dir, "*_epochs.fif"))[0]
epochs = mne.read_epochs(epochs_file, preload=False)
info = epochs.info
channels = epochs.ch_names

# Columns of interest
cols_power = [
    'abs_delta', 'abs_theta', 'abs_alpha', 'abs_sigma', 'abs_beta', 'abs_gamma',
    'rel_delta', 'rel_theta', 'rel_alpha', 'rel_sigma', 'rel_beta', 'rel_gamma'
]

# Sleep states
sleep_states = ["W", "N1", "N2", "N3", "R"]

# Plotting
for feature in cols_power:
    fig, axes = plt.subplots(
        nrows=1, ncols=len(sleep_states), figsize=(20, 4), constrained_layout=True)

    vmin = df[feature].min()
    vmax = df[feature].max()

    for i_ss, ss in enumerate(sleep_states):
        temp_values = df[df['sleep_state'] == ss][['channel', feature]]

        temp_values = temp_values.groupby('channel', as_index=False).mean()
        temp_values['channel'] = pd.Categorical(
            temp_values['channel'],
            categories=channels,
            ordered=True
        )
        df_sorted = temp_values.sort_values('channel')

        values = df_sorted[feature].to_numpy()

        ax = axes[i_ss]
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        im, _ = mne.viz.plot_topomap(
            data=values,
            pos=info,
            axes=ax,
            contours=2,
            cmap="Purples",
            vmin=vmin,
            vmax=vmax,
            show=False
        )
        plt.colorbar(im, cax=cax, orientation='vertical')
        ax.set_title(f"{ss}", fontweight="bold")

    fig.suptitle(f"{feature}\n", fontsize="xx-large", fontweight="bold")
    plt.show()
    # Save figure
    fig.savefig(os.path.join(plots_dir, f"{feature}_topomap.png"), dpi=300)
    plt.close(fig)

# Statistical comparisons using Linear Mixed-Effects Models

# Group data by subject, sleep_state, and channel
mean_df = df.groupby(
    ['sub_id', 'sleep_state', 'channel'],
    as_index=False
).mean()

interest = 'rel_alpha'  # Change this to the band you are interested in
fdr_corrected = False  # Set to True if you want FDR correction

# We will compare each sleep state against W (Wake)
reference_state = 'W'
comparison_states = ['N1', 'N2', 'N3', 'R']

fig, axes = plt.subplots(
    nrows=1, ncols=len(comparison_states), figsize=(20, 4), constrained_layout=True)

for i, sleep_state in enumerate(comparison_states):
    temp_tval = []
    temp_pval = []
    chan_l = []

    cond_df = mean_df[mean_df['sleep_state'].isin([reference_state, sleep_state])]

    for chan in channels:
        subdf = cond_df[cond_df['channel'] == chan].dropna()

        if len(subdf['sleep_state'].unique()) < 2:
            # Not enough data to compare
            temp_tval.append(np.nan)
            temp_pval.append(1.0)
            chan_l.append(chan)
            continue

        model = f"{interest} ~ C(sleep_state, Treatment('{reference_state}'))"

        md = smf.mixedlm(model, subdf, groups=subdf['sub_id'], missing='drop')
        try:
            mdf = md.fit()
            t_value = mdf.tvalues[f"C(sleep_state, Treatment('{reference_state}'))[T.{sleep_state}]"]
            p_value = mdf.pvalues[f"C(sleep_state, Treatment('{reference_state}'))[T.{sleep_state}]"]
            temp_tval.append(t_value)
            temp_pval.append(p_value)
        except Exception as e:
            print(f"Could not fit model for channel {chan}: {e}")
            temp_tval.append(np.nan)
            temp_pval.append(1.0)

        chan_l.append(chan)

    temp_tval = np.array(temp_tval)
    temp_pval = np.array(temp_pval)

    # Handle NaN values in temp_tval
    if np.any(np.isnan(temp_tval)):
        nan_indices = np.where(np.isnan(temp_tval))[0]
        temp_tval[nan_indices] = np.nanmean(temp_tval)

    if fdr_corrected:
        _, corrected_pval = fdrcorrection(temp_pval)
        display_pval = corrected_pval
    else:
        display_pval = temp_pval

    ax = axes[i]
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    im, _ = mne.viz.plot_topomap(
        data=temp_tval,
        pos=info,
        axes=ax,
        contours=3,
        mask=np.asarray(display_pval) <= 0.05,
        mask_params=dict(marker='o', markerfacecolor='w', markeredgecolor='k',
                         linewidth=0, markersize=6),
        cmap="viridis",
        vlim=(np.percentile(temp_tval, 5), np.percentile(temp_tval, 95)),
        show=False
    )
    plt.colorbar(im, cax=cax, orientation='vertical')
    ax.set_title(f"{sleep_state} vs {reference_state}", fontweight="bold")

fig.suptitle(f"Statistical Comparison for {interest}", fontsize="xx-large", fontweight="bold")
plt.show()
# Save figure
fig.savefig(os.path.join(plots_dir, f"{interest}_statistical_comparison.png"), dpi=300)
plt.close(fig)
