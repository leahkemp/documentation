import pandas as pd

csv1 = pd.read_csv('resource_benchmarking_merged_production.csv')
csv2 = pd.read_csv('nodes.csv')
merged = pd.merge(csv1, csv2, how = 'left', on = ['sample', 'rule', 'threads_per_rule', 'no_samples', 'resource_type', 'run_type', 'data_type', 'pipeline'], suffixes=('_benchmarking_file', '_log_file'))
merged.to_csv("merged.csv", index = False)