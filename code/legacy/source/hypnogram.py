#%% -*- coding: utf-8 -*-
"""
Plot Hypnograms for Each Participant
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import yasa

# Directory paths
data_dir = '/Users/borjan/code/healthy-sleep-project/data/raw/'
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/sleep-results/hypno/'

# Create the results directory if it doesn't already exist
os.makedirs(results_dir, exist_ok=True)

# Define mapping from sleep stage labels to numerical codes
sleep_stage_mapping = {
    'W': 0,   # Wake
    'N1': 1,  # N1 sleep
    'N2': 2,  # N2 sleep
    'N3': 3,  # N3 sleep
    'R': 4,   # REM sleep
    # 'L': -1  # Decide how to handle 'L'; see explanation below
}

# Loop through participants 01 to 29
for subject_num in range(1, 30):
    # Format subject identifier with leading zeros (e.g., EPCTL01, EPCTL02, ..., EPCTL29)
    subject_id = f'EPCTL{subject_num:02d}'
    
    # Construct the file path
    file_path = os.path.join(data_dir, subject_id, f'{subject_id}.txt')
    
    # Check if the file exists
    if os.path.exists(file_path):
        print(f'Processing {subject_id}...')
        
        # Initialize list to store sleep stages
        sleep_stages = []
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if line == '':
                        continue  # Skip empty lines
                    # Split the line into columns
                    parts = line.split()
                    if len(parts) != 3:
                        print(f'Warning: Invalid line format in {subject_id}: {line}')
                        continue  # Skip invalid lines
                    stage_label, start_time, duration = parts
                    # Map the stage label to numerical code
                    if stage_label in sleep_stage_mapping:
                        stage_code = sleep_stage_mapping[stage_label]
                    else:
                        # Handle unknown stage labels
                        print(f'Warning: Unknown sleep stage label "{stage_label}" in {subject_id}')
                        continue  # Skip unknown stages
                    # Number of epochs for this duration
                    num_epochs = int(int(duration) / 30)  # Assuming duration is in seconds, and epochs are 30 seconds
                    # Append the stage code for each epoch
                    sleep_stages.extend([stage_code] * num_epochs)
            # Now sleep_stages is a list of stage codes for each 30-sec epoch
            # Convert to numpy array
            sleep_stages = np.array(sleep_stages)
            
            # Plot the hypnogram
            fig, ax = plt.subplots(1, 1, figsize=(10, 4), constrained_layout=True)
            ax = yasa.plot_hypnogram(sleep_stages, fill_color="lightblue", ax=ax)
            ax.set_title(f'Hypnogram for {subject_id}')
            
            # Save the figure
            fig_filename = os.path.join(results_dir, f'{subject_id}_hypnogram.png')
            plt.savefig(fig_filename)
            plt.close(fig)
            print(f'Hypnogram saved for {subject_id} at {fig_filename}\n')
        except Exception as e:
            print(f'Error processing {subject_id}: {e}')
    else:
        print(f'File not found for {subject_id}: {file_path}\n')

# %%
