import numpy as np
import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
import matplotlib.pyplot as plt
import scipy

########################################################################################################################
# Compare RT in a mixed model ANOVA
df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_rmanova.csv')

aovrm = pg.mixed_anova(dv='Response_Time', within='Condition', between='Group', subject='MIND_ID', data=df)
print(aovrm.round(3))
pg.print_table(aovrm)

df_con = df[df['Group'] == 'Control']
ms = df_con[df['Condition'] == 'MultiSource']['Response_Time'].values
flanker = df_con[df['Condition'] == 'Flanker']['Response_Time'].values
simon = df_con[df['Condition'] == 'Simon']['Response_Time'].values
control = df_con[df['Condition'] == 'Control']['Response_Time'].values

ms_con = scipy.stats.ttest_rel(ms, control)
ms_flanker = scipy.stats.ttest_rel(ms, flanker)
ms_simon = scipy.stats.ttest_rel(ms, simon)

flanker_con = scipy.stats.ttest_rel(flanker, control)
flanker_simon = scipy.stats.ttest_rel(flanker, simon)

simon_con = scipy.stats.ttest_rel(simon, control)

df_hiv = df[df['Group'] == 'HIV']
ms = df_hiv[df['Condition'] == 'MultiSource']['Response_Time'].values
flanker = df_hiv[df['Condition'] == 'Flanker']['Response_Time'].values
simon = df_hiv[df['Condition'] == 'Simon']['Response_Time'].values
control = df_hiv[df['Condition'] == 'Control']['Response_Time'].values

ms_con = scipy.stats.ttest_rel(ms, control)
ms_flanker = scipy.stats.ttest_rel(ms, flanker)
ms_simon = scipy.stats.ttest_rel(ms, simon)

flanker_con = scipy.stats.ttest_rel(flanker, control)
flanker_simon = scipy.stats.ttest_rel(flanker, simon)

simon_con = scipy.stats.ttest_rel(simon, control)

########################################################################################################################
# Compare Accuracy in a Friedman's test

df = pd.read_csv('E:/Data/MSIT_MIND/accuracy_129_friedman.csv')
df_con = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']
ftest = scipy.stats.friedmanchisquare(df_con[df_con['Condition'] == 'Control']['Accuracy'].values,
                              df_con[df_con['Condition'] == 'Simon']['Accuracy'].values,
                              df_con[df_con['Condition'] == 'Flanker']['Accuracy'].values,
                              df_con[df_con['Condition'] == 'MultiSource']['Accuracy'].values)


aovrm = pg.friedman(dv='mean', within='condition', subject='ID', data=df, method='chisq')
print(aovrm)
posthoc = sp.posthoc_conover(df_con, val_col='Accuracy', group_col='Condition', p_adjust='bonf')
print(posthoc)
