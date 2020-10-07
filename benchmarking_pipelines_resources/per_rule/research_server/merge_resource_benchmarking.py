import pandas as pd
import os
from glob import glob
import re

# Function to remove filepath in rule names
def remove_filepath_rule(dirname):
    return dirname.rsplit('/',1)[1]

# Function to extract threads from the filepath
def extract_filepath_threads(dirname):
    return re.findall(r'\d+', dirname)[0]

# Function to remove filepath in pipeline names
def remove_filepath_pipeline(dirname):
    return dirname.rsplit('/',4)[1]

# Define the directories to search
benchmarking_directories = glob("./*/*/*/workflow/benchmarks/*/*")

appended_data = []

for benchmark in benchmarking_directories:
    # Read file
    df = pd.read_csv(benchmark, index_col = None, sep = '\t')

    # Add column with benchmarking repeat number
    repeat_col = ['1','2','3']
    df['repeat'] = repeat_col

    # Add column that defines 'sample' from filename
    df['sample'] = os.path.basename(benchmark)
    df['sample'] = df['sample'].str.replace('.tsv', '')

    # Add column that defines 'rule' from file directory
    df['rule'] = remove_filepath_rule(os.path.dirname(benchmark))

    # Add column that defines the thread number of each rule/sample from file directory
    df['threads_per_rule'] = extract_filepath_threads(os.path.dirname(benchmark))

    # Add column that defines pipeline from file directory
    df['pipeline'] = remove_filepath_pipeline(os.path.dirname(benchmark))

    # Add column that gives the relative filepath
    df['filepath'] = os.path.dirname(benchmark)

    # Append all dataframes
    appended_data.append(df)

# Merge into one dataframe
appended_data = pd.concat(appended_data)

# Write to csv
appended_data.to_csv('resource_benchmarking_merged.csv', index = False)
