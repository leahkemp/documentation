import pandas as pd
import os
from glob import glob
import re

# Function to remove filepath in rule names
def remove_filepath_rule(dirname):
    return dirname.rsplit('/',1)[1]

# Function to extract threads from the filepath
def extract_filepath_threads(dirname):
    return re.findall(r'\d+', dirname)[1]

# Function to extract pipeline run sample number from the filepath
def extract_filepath_no_samples(dirname):
    return re.findall(r'\d+', dirname)[0]

# Function to remove filepath in resource type names
def remove_filepath_resource_type(dirname):
    return dirname.rsplit('/',8)[1]

# Function to remove filepath in run type names
def remove_filepath_run_type(dirname):
    return dirname.rsplit('/',6)[1]

# Function to remove filepath in data type names
def remove_filepath_data_type(dirname):
    return dirname.rsplit('/',9)[1]

# Function to remove filepath in pipeline names
def remove_filepath_pipeline(dirname):
    return dirname.rsplit('/',4)[1]

# Function to get the node the rule/sample was run on from the contents of the log file
def get_node_from_log(dirname):
    log_file = open(dirname, "rt")
    contents = log_file.read()
    log_file.close()
    match =  re.search("kscprod-(\w*)\s", contents)

    if not match:
        return "NA"

    return match.group(0)


log_directories = glob("./exome/*/*/*/*/*/workflow/logs/*/*")

# Create empty dataframe to append data to
df = pd.DataFrame(columns = ['sample', 'rule', 'threads_per_rule', 'no_samples', 'resource_type', 'run_type', 'data_type', 'pipeline', 'filepath', 'node'])

appended_data = []

for log in log_directories:

    # Create variables that define contents to fill dataframe
    sample = os.path.basename(log)
    sample = sample.replace('.log', '')
    rule = remove_filepath_rule(os.path.dirname(log))
    threads_per_rule = extract_filepath_threads(os.path.dirname(log))
    no_samples = extract_filepath_no_samples(os.path.dirname(log))
    resource_type = remove_filepath_resource_type(os.path.dirname(log))
    run_type = remove_filepath_run_type(os.path.dirname(log))
    data_type = remove_filepath_data_type(os.path.dirname(log))
    pipeline = remove_filepath_pipeline(os.path.dirname(log))
    filepath = os.path.dirname(log)
    node = get_node_from_log(log)

    # Append variables to dataframe
    df = df.append({'sample': sample, 'rule': rule, 'threads_per_rule': threads_per_rule, 'no_samples': no_samples, 'resource_type': resource_type, 'run_type': run_type, 'data_type': data_type, 'pipeline': pipeline, 'filepath': filepath, 'node': node}, ignore_index=True)



# Write to csv
df.to_csv('nodes.csv', index = False)