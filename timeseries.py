import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/alpha_-30_-76_5/timeseries/alpha_-30_-76_5_abs_timeseries_concat.csv')



condition_identifier = ["" for x in range(df.shape[0])]
for i, row in df.iterrows():

    if 'Control' in row['FileNameConcat']:
        condition_identifier[i] = str('Control')

    elif 'Simon' in row['FileNameConcat']:
        condition_identifier[i] = 'Simon'

    elif 'Flanker' in row['FileNameConcat']:
        condition_identifier[i] = 'Flanker'

    elif 'Multi' in row['FileNameConcat']:
        condition_identifier[i] = 'MultiSource'

df['condition_identifier'] = condition_identifier
