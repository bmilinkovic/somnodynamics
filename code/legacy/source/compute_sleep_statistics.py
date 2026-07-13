#%% -*- coding: utf-8 -*-
"""
Compute Sleep Statistics for Each Participant and Create Summary Table and Plots
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yasa

# Directory paths
data_dir = '/Users/borjan/code/healthy-sleep-project/data/raw/'
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/sleep-results/hypno/'

# Create the results directory if it doesn't already exist
os.makedirs(results_dir, exist_ok=True)

# Define mapping from sleep stage labels to numerical codes compatible with YASA
sleep_stage_mapping = {
    'W': 0,    # Wake
    'N1': 1,   # N1 sleep
    'N2': 2,   # N2 sleep
    'N3': 3,   # N3 sleep
    'R': 4,    # REM sleep
    'L': -2    # Unscored or Lights Off (you can adjust this based on your data)
}

# Initialize a list to store the statistics for all participants
all_stats = []

# Loop through participants 01 to 29
for subject_num in range(1, 30):
    # Format subject identifier with leading zeros (e.g., EPCTL01)
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
                        # Handle unknown sleep stage labels
                        print(f'Warning: Unknown sleep stage label "{stage_label}" in {subject_id}')
                        continue  # Skip unknown stages
                    # Number of epochs for this duration
                    num_epochs = int(int(duration) / 30)  # Assuming duration is in seconds, and epochs are 30 seconds
                    # Append the stage code for each epoch
                    sleep_stages.extend([stage_code] * num_epochs)
            # Convert sleep stages to numpy array
            sleep_stages = np.array(sleep_stages)
            
            # Compute sleep statistics using YASA
            sf_hyp = 1 / 30  # Sampling frequency (1 value per 30 seconds)
            stats = yasa.sleep_statistics(sleep_stages, sf_hyp)
            
            # Add participant ID to the stats dictionary
            stats['Participant'] = subject_id
            
            # Append stats to the list
            all_stats.append(stats)
            
        except Exception as e:
            print(f'Error processing {subject_id}: {e}')
    else:
        print(f'File not found for {subject_id}: {file_path}\n')

# After processing all participants, create a DataFrame
stats_df = pd.DataFrame(all_stats)

# Reorder columns to place 'Participant' first
cols = ['Participant'] + [col for col in stats_df.columns if col != 'Participant']
stats_df = stats_df[cols]

# Save the DataFrame to a CSV file
summary_csv_filename = os.path.join(results_dir, 'sleep_statistics_summary.csv')
stats_df.to_csv(summary_csv_filename, index=False)
print(f'Sleep statistics summary saved at {summary_csv_filename}')

# Display the DataFrame
print(stats_df)

# Visualize some key sleep statistics across participants
# Example: Plotting Total Sleep Time (TST), Sleep Efficiency (SE), and Sleep Onset Latency (SOL)

# List of metrics to plot
metrics_to_plot = ['TST', 'SE', 'SOL', 'WASO', 'N1', 'N2', 'N3', 'REM']

# Create histograms for each metric
for metric in metrics_to_plot:
    plt.figure(figsize=(8, 6))
    plt.hist(stats_df[metric], bins=10, color='skyblue', edgecolor='black')
    plt.xlabel(metric)
    plt.ylabel('Number of Participants')
    plt.title(f'Distribution of {metric} Across Participants')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    # Save the figure
    fig_filename = os.path.join(results_dir, f'{metric}_distribution.png')
    plt.savefig(fig_filename)
    plt.show()
    print(f'{metric} distribution plot saved at {fig_filename}')

# Optional: Create boxplots to compare metrics
plt.figure(figsize=(12, 6))
stats_df.boxplot(column=metrics_to_plot[0:4])
plt.ylabel('Minutes')
plt.title('Boxplot of Sleep Statistics (Time in Minutes)')
plt.tight_layout()
# Save the figure
boxplot_filename = os.path.join(results_dir, 'sleep_statistics_boxplot.png')
plt.savefig(boxplot_filename)
plt.show()
print(f'Sleep statistics boxplot saved at {boxplot_filename}')

# %%
