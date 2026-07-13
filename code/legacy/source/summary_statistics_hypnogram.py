#%% -*- coding: utf-8 -*-
"""
Calculate and Plot Summary Statistics for Sleep Stages Including Wake Time
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
    # 'L': -1  # Decide how to handle 'L'; exclude or include as appropriate
}

# Initialize list to store percentages for all participants
all_percentages = []

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
            
            # Calculate total recording time (in minutes)
            total_time = len(sleep_stages) * 30 / 60  # in minutes
            
            # Calculate time spent in each sleep stage (in minutes)
            stage_times = {}
            for stage_code in [0, 1, 2, 3, 4]:  # Including Wake (0) and sleep stages
                time_in_stage = len(sleep_stages[sleep_stages == stage_code]) * 30 / 60  # in minutes
                stage_times[stage_code] = time_in_stage
            
            # Calculate percentage of total recording time for each stage
            percentages = {}
            for stage_code in stage_times:
                if total_time > 0:
                    percentages[stage_code] = (stage_times[stage_code] / total_time) * 100
                else:
                    percentages[stage_code] = 0.0  # Avoid division by zero
            
            # Append percentages to the list
            all_percentages.append(percentages)
            
        except Exception as e:
            print(f'Error processing {subject_id}: {e}')
    else:
        print(f'File not found for {subject_id}: {file_path}\n')

# After processing all participants, calculate group averages and standard deviations
# Initialize lists for each sleep stage
stages = [0, 1, 2, 3, 4]  # Including Wake (0)
stage_labels = ['Wake', 'N1', 'N2', 'N3', 'REM']
percentage_values = {stage: [] for stage in stages}

for participant_percentages in all_percentages:
    for stage in stages:
        percentage_values[stage].append(participant_percentages.get(stage, 0.0))

# Calculate mean and standard deviation for each stage
means = []
stds = []
for stage in stages:
    stage_percentages = percentage_values[stage]
    mean = np.mean(stage_percentages)
    std = np.std(stage_percentages)
    means.append(mean)
    stds.append(std)

# Plotting the average percentages with error bars
plt.figure(figsize=(10, 6))
x_pos = np.arange(len(stage_labels))
plt.bar(x_pos, means, yerr=stds, align='center', alpha=0.8, capsize=10, color='skyblue')
plt.xticks(x_pos, stage_labels)
plt.xlabel('Sleep Stages')
plt.ylabel('Percentage of Total Recording Time')
plt.title('Average Sleep Stage Distribution Including Wake Across Participants')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Save the figure
summary_fig_filename = os.path.join(results_dir, 'average_sleep_stage_distribution_including_wake.png')
plt.savefig(summary_fig_filename)
plt.show()

print(f'Average sleep stage distribution plot including wake saved at {summary_fig_filename}')

# %%
