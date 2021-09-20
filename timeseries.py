import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

filename = 'E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/gamma_-14_-88_5/timeseries/gamma_-14_-88_5_abs_timeseries_concat.csv'

df = pd.read_csv(filename)

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

df['condition'] = condition_identifier
df['mean'] = df.mean(axis=1)
# df['peak'] = df.iloc[:, 29:46].min(axis=1) # Alpha
df['peak'] = df.iloc[:, 27:34].max(axis=1) # Gamma

df.to_csv(filename)

time = np.linspace(-500, 1000, df.shape[1]-2)

plt.figure()
plt.plot(time, df[df['condition'] == 'Control'].mean(), label='Control')
plt.plot(time, df[df['condition'] == 'Simon'].mean(), label='Simon')
plt.plot(time, df[df['condition'] == 'Flanker'].mean(), label='Flanker')
plt.plot(time, df[df['condition'] == 'MultiSource'].mean(), label='MultiSource')
plt.axvline(x=0, color='black', alpha=.5)
plt.axvline(x=150, color='black', linestyle='--', alpha=.5)
plt.axvline(x=300, color='black', linestyle='--', alpha=.5)
plt.legend()
plt.savefig('E:/Data/MSIT_MIND/msit_visual/figures/gamma_-14_-88_5_abs.jpg', dpi=300)
plt.show()

