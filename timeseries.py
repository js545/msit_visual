import scipy
import numpy as np
import pandas as pd
import pingouin as pg
import matplotlib.pyplot as plt
from scipy.stats import f_oneway
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multicomp import pairwise_tukeyhsd

########################################################################################################################

filename = 'E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/gamma_30_-88_9/timeseries/gamma_30_-88_9_rel_timeseries_concat.csv'

scipy.stats.shapiro()


########################################################################################################################

# df = pd.read_csv(filename)
#
# condition_identifier = ["" for x in range(df.shape[0])]
# for i, row in df.iterrows():
#
#     if 'Control' in row['FileNameConcat']:
#         condition_identifier[i] = str('Control')
#
#     elif 'Simon' in row['FileNameConcat']:
#         condition_identifier[i] = 'Simon'
#
#     elif 'Flanker' in row['FileNameConcat']:
#         condition_identifier[i] = 'Flanker'
#
#     elif 'Multi' in row['FileNameConcat']:
#         condition_identifier[i] = 'MultiSource'
#
# df['condition'] = condition_identifier
# # df['mean'] = df.iloc[:, 29:46].mean(axis=1) # Alpha
# df['mean'] = df.iloc[:, 27:34].mean(axis=1) # Gamma
# # df['peak'] = df.iloc[:, 29:46].min(axis=1) # Alpha
# df['peak'] = df.iloc[:, 27:34].max(axis=1) # Gamma
#
# df.to_csv(filename, index=False)

# ANOVA

df = pd.read_csv(filename)

# aovrm = AnovaRM(df, 'mean', 'ID', within=['condition'])
# res = aovrm.fit()
# print(res)
aovrm = pg.rm_anova(dv='peak', within='condition', subject='ID', data=df, detailed=True)
print(aovrm)
# aovrm = pg.friedman(dv='mean', within='condition', subject='ID', data=df)
# print(aovrm)

# scipy.stats.ttest_rel(df[df['condition'] == 'Control']['mean'], df[df['condition'] == 'MultiSource']['mean'])

########################################################################################################################
# Two-way RM-ANOVA

aovrm = pg.rm_anova(dv='mean', within=['condition', 'group'], subject='ID', data=df, detailed=True)
print(aovrm.round(3))
pg.print_table(aovrm)

# tukey = pairwise_tukeyhsd(endog=df['peak'], groups=df['condition'], alpha=.05)
# print(tukey)


########################################################################################################################
# Plot figures

time = np.linspace(-500, 1000, df.shape[1]-4)

plt.figure()
plt.plot(time, df[df['condition'] == 'Control'].mean()[:61], label='Control')
plt.plot(time, df[df['condition'] == 'Simon'].mean()[:61], label='Simon')
plt.plot(time, df[df['condition'] == 'Flanker'].mean()[:61], label='Flanker')
plt.plot(time, df[df['condition'] == 'MultiSource'].mean()[:61], label='MultiSource')
plt.axvline(x=0, color='black', alpha=.5)
plt.axvline(x=150, color='black', linestyle='--', alpha=.5)
plt.axvline(x=300, color='black', linestyle='--', alpha=.5)
plt.legend()
plt.savefig('E:/Data/MSIT_MIND/msit_visual/figures/gamma_-14_-88_5_abs.jpg', dpi=300)
plt.show()

# Plot HIV vs Control Figures

df_hiv = df[df['group'] == 'HIV']
df_con = df[df['group'] == 'Control']

plt.figure()
plt.plot(time, df_hiv[df_hiv['condition'] == 'Control'].mean()[:61], label='HIV_Control', color='blue')
plt.plot(time, df_con[df_con['condition'] == 'Control'].mean()[:61], label='Control_Control', linestyle='--', color='blue')
plt.plot(time, df_hiv[df_hiv['condition'] == 'Simon'].mean()[:61], label='HIV_Simon', color='red')
plt.plot(time, df_con[df_con['condition'] == 'Simon'].mean()[:61], label='Control_Simon', linestyle='--', color='red')
plt.plot(time, df_hiv[df_hiv['condition'] == 'Flanker'].mean()[:61], label='HIV_Flanker', color='green')
plt.plot(time, df_con[df_con['condition'] == 'Flanker'].mean()[:61], label='Control_Flanker', linestyle='--', color='green')
plt.plot(time, df_hiv[df_hiv['condition'] == 'MultiSource'].mean()[:61], label='HIV_MS', color='purple')
plt.plot(time, df_con[df_con['condition'] == 'MultiSource'].mean()[:61], label='Control_MS', linestyle='--', color='purple')
plt.axvline(x=0, color='black', alpha=.5)
plt.axvline(x=150, color='black', linestyle='--', alpha=.5)
plt.axvline(x=300, color='black', linestyle='--', alpha=.5)
plt.legend()
plt.show()